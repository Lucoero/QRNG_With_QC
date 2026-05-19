# -*- coding: utf-8 -*-
"""
Analyse data:
    Libreria para tomas las secuencias obtenidas por los QC para analizar las propiedades
    de aleatoriedad. Esto es:
        - Estimar Hmin y una entropía total (si es posible)
        
        - Probar los distintos tests estadísticos elegidos (ver si me basta con NIST o uso tambien Dieharder)
        - Realizar tests de orden dos (tests a los resultados de los tests). Usando el enfoque de Shen
        
    Los resultados se almacenaran en txts asociados a cada computador, concatenando con 
    lo que haya en esos documentos previamente. 
"""
#%% Librerias
import numpy as np
import os as os
import WriteFile as WF
import tqdm as tqdm
# Tests NIST (de stevenang: https://github.com/stevenang/randomness_testsuite)
from NIST.ApproximateEntropy import ApproximateEntropy as aet
from NIST.Complexity import ComplexityTest as ct
from NIST.CumulativeSum import CumulativeSums as cst
from NIST.FrequencyTest import FrequencyTest as ft
from NIST.RandomExcursions import RandomExcursions as ret
from NIST.RunTest import RunTest as rt
from NIST.Spectral import SpectralTest as st
from NIST.Universal import Universal as ut
import NIST.BinaryMatrix as bm
# Calculo de Complejidades
import gzip as gzip # Compresor de Huffmann https://docs.python.org/3/library/gzip.html
import lzip as lzip # Compresor de Cadenas de Markov https://www.nongnu.org/lzip/
import pyppmd as ppm # Compresor PPM https://pyppmd.readthedocs.io/en/latest/api_guide.html

from uncertainties import ufloat
#%% Variables globales
readStep = 10**6 # N of bytes for each read order
epsilon = 0.01 # Tolerance of tests
#%% Funciones

#%%% Ploteados
def Axe_Canvas(n,cost,ax, linesNames = [], multipleLines = False):
    ax.set_xlabel("n Bits")
    ax.set_ylabel("Cost ($)")
    ax.grid()
    if multipleLines:
        nLines = len(n)
        if len(linesNames) == 0:
            linesNames = [f"Line {i+1}" for i in range(nLines)]
        for i in range(nLines):
            ax.plot(n[i], cost[i], label = linesNames[i])
        ax.legend()
    else:
        ax.plot(n[i],cost[i])
    return ax

#%%% Analisis de Entropia
def Compute_Prob(path, showCount = False):
    """
    Estima las probabilidades de que se manifieste un 0 o un 1 en la secuencia
    de acuerdo a la frecuencia en las cadenas. Cuanto mas larga la secuencia, 
    mas confianza tendra esta estimacion.
    
    Esta funcion si se aplica showCount=True muestra tambien las diferencias entre 0s y 1s.
    """
    zCount = 0
    oCount = 0
    bitCount = 0
    
    try:
        with open(path,"r") as f: # Lo leemos a cachos para evitar overflow
            strArr = f.read(readStep)
            currLen = len(strArr)
            while currLen != 0:
                curr_oCount = strArr.count('1')
                oCount += curr_oCount
                zCount += currLen-curr_oCount
                
                bitCount += currLen
                currBit = f.read(readStep)
                currLen = len(currBit)
    except:
        print(f"ERROR: No se ha podido encontrar el archivo en el path \n {path}")
        return None,None
    if showCount:
        print(f"0s: {zCount}, 1s: {oCount}, bits: {bitCount}, diffPer: {(oCount-zCount)/bitCount}")
    # Calculamos las probabilidades
    zProb  = zCount/bitCount
    oProb = 1-zProb
    return zProb,oProb

def Estimate_Hmin(fileName):
    """
    Basta con tomar el suceso con mas probabilidades de la secuencia
    
    Usamos la def 
    
    Hmin = -log(Max(P_i))
    
    TODO: Ver como cuantizar la confianza de la estimacion
    """
    maxP = max(Compute_Prob(fileName))
    return -np.log2(maxP)

def Estimate_H(fileName):
    """
    Calculamos H de acuerdo a la definicion de entropia de Shannon
    
    TODO: cuantizar la confianza de la estimacion
    """
    probs = Compute_Prob(fileName)
    return -probs[0]*np.log2(probs[0]) - probs[1]*np.log2(probs[1])
#%%% Analisis de la complejidad (distintos compresores)
def Estimate_KComplexity(path):
    """
    Estimate_KComplexity:
        Dado un archivo con bits, comprobamos su compresion con diferentes algoritmos
        (Deflate// Huffman con gzip; Lempel-ziv-Cadenas de Markov con lzip y PPM con pyppmd)
        Tomamos como Complejidad de Kolmogorov como la longitud de la maxima compresion posible, y devolvemos
        tambien la funcion de deficiencia 
    """
    y = WF.ReadFile(path,readAsBit = True)
    #print("Estimating Kolmogorov_Complexity")
    initLength = len(y)
    deflate = len(gzip.compress(y, compresslevel = 9)) # Huffmann
    lempel = len(lzip.compress_to_buffer(y,level = 9)) # Lempel-Ziv
    PPM = len(ppm.compress(y, max_order = 64)) # PPM. maxorder = 64 es la compresion maxima
    maxComp = min(deflate,lempel,PPM)
    d = initLength - maxComp
    probD = 2**(-d)
    return maxComp,d,probD, probD > epsilon # Si superas epsilon, eres aleatoria
#%%% Calculo de la tasa de generado
def Compute_Creation_Ratio(filePath):
    """
    Compute_Creation_Ratio:
        Dados los tiempos almacenados en la carpeta Times, calculamos
        la media de generado con su error, y luego la tasa de generado por bit
        con su error
    """
    times = [] # Array de tiempos 
    bits = [] # Array de longitud de la cadena
    with open(filePath) as file:
        line = file.readline()
        while line != "":
            sides = line.split(":") # Separo la linea en lado derecho e izquierdo
            times.append(float(sides[1]))
            compData = sides[0].split("--")[1] # Elimino el nombre del computador
            aux = compData.split("-") # Separo en elementos
            bits.append(float(aux[0])*float(aux[1]))
            line = file.readline()
    # Con los tiempos y los bits calculamos las medias
    times = np.array(times,dtype = float)
    bits = np.array(bits,dtype = int)
    meanTime = np.mean(times)
    std = np.std(times)
    errMean = std/np.sqrt(len(times))
    
    totalBits = np.sum(bits)
    ratio = totalBits/meanTime # Velocidad bits/s
    errRatio = errMean/meanTime * ratio # Error abs
    return ratio,errRatio
#%%% Calculo de la funcion Costo 

def AWS_Cost_Function(nBits,cs,ct,sMax,sMin,nqbits):
    """
    Cost_Function:
        Calcula el coste 
    """
    # Tiene que funcionar entre sMin y sMax. En los bits intermedios debe ser infinito
    #if nBits > nqbits: nqbits = nBits
    nShots = nBits//nqbits
    nShots = np.where(nShots == 0, 1,nShots) # Si estamos en el caso de pedir menos bits que qbits
    T = nShots//sMax# Contamos que se debe empezar con 1 task
    R = nShots%sMax
    rLess_sMin = R < sMin
    # Donde R < sMin debo tener una task mas y los shots no cambian (solo que ahora es compatible por repartir smax + R en dos tasks)
    f = np.where(rLess_sMin,ct*(T+2) + cs*(sMax*T + R),ct*(T+1) + cs*(sMax*T + R)) # El +1 en las T es porque contamos la task 0
    f = np.where(nShots < 1, ct + cs*nShots ,f) # Caso concreto de pedir menos bits que qbits
    return f

def IBM_Cost_Function(nBits,ratio, cost_perSecond):
    """
    Dado que IBM cobra por minuto usado, calculamos el coste
    por nBits
    """
    # Calculo el tiempo necesario para generar esos bits
    timeNeeded = nBits/ratio 
    return cost_perSecond*timeNeeded

#%%% TEST SEGUNDO ORDEN
def Mix_Sequences(dirPr,dirRef,outDir = "Mix_Sequences/Mixed_IBM"):
    """
    Mix Sequences: Dados dos ficheros de bits, combinamos los ficheros con un xor
    """
    bitsPrArr = WF.Scan_Dir(dirPr)
    bitsRefArr = WF.Scan_Dir(dirRef)
    lenPr = len(bitsPrArr)
    lenRef = len(bitsRefArr)
    minLen = min(lenPr,lenRef)
    for j in tqdm.tqdm(range(minLen)):
        bitsPr = WF.ReadFile(bitsPrArr[j],readAsBit = True)
        bitsRef = WF.ReadFile(bitsRefArr[j],readAsBit = True)
        nPr = len(bitsPr)
        nRef = len(bitsRef)
        nOut = min(nPr,nRef)
        
        bitsOut = np.empty(nOut,dtype = str)
        for i in range(nOut):
            bitsOut[i] = bitsPr[i] ^ bitsRef[i] # ^ es la operacion XOR en bits 
        
        prBaseName = os.path.basename(bitsPrArr[j])
        refBaseName = os.path.basename(bitsRefArr[j])
        
        WF.Write_Bitfiles(outDir, prBaseName+"--"+refBaseName, bitsOut)
    return   

def Compute_Distribution(t,pValues):
    """
    Dada una distribucion de pvalues, calcula la distribucion numerica asociada
    """
    m = len(pValues)
    nPoints = len(t)
    out = np.zeros(nPoints)
    for i in range(nPoints): 
        aux = np.where(pValues <= t[i])
        out[i] = len(aux[0])
    return out/m

def KS_Shen_Test(dirPr,dirMix, outDir = "KS_Results",nPoints = 5000):
    """
    dirPr: Directorio donde se encuentran las secuencias del generador problema
    dirRef: Directorio donde se encuentran los p-values del generador de referencia
    nPoints: Numero de puntos de la discretizacion de las distribuciones G
    Pasamos por toda la bateria NIST
    """
    K = np.sqrt(-0.5 *np.log2(epsilon/2)) # Cota test segundo orden
    
    tests = [ft.monobit_test,Block_Frequency_Test,rt.run_test, rt.longest_one_block_test,st.spectral_test,aet.approximate_entropy_test,
             cst.cumulative_sums_test,ut.statistical_test,
             Binary_Matrix_Test,ct.linear_complexity_test]
    
    testsNames = ["Monobit Test","Frequency Within a Block Test", "Run Test", "Longest One Block Test",
                  "Spectral Test", "Approximate Entropy Test", "Cusum Test", "Maurer Test", "Binary Matrix Test",
                  "Linear Complexity Test"]
    
    #tests = [ft.monobit_test] # Debugging
    #testsNames = ["Monobit Test"]
    # run tests a parte
    nTests = len(tests)
    #%%%% BLOQUE I: Extraccion de informacion
    # Cargamos en memoria la lista de paths a cada archivo
    pathPr = WF.Scan_Dir(dirPr)
    pathMix = WF.Scan_Dir(dirMix)
    
    nPr = len(pathPr)
    nMix = len(pathMix)
    
    if nPr > nMix:
        minFile = nMix
    else:
        minFile = nPr
    
    Gpr = np.zeros((nTests,nPoints)) 
    GMix = np.zeros((nTests,nPoints))
    tMatrix = np.zeros(nTests,dtype = object)
    
    pr_pValues = np.zeros((nTests,minFile))
    mix_pValues = np.zeros((nTests,minFile))
    
    KSResults = np.zeros(nTests)
    passArr = np.zeros(nTests,dtype = bool)
    for i in tqdm.tqdm(range(nTests)):     # Vamos test a test calculando KS
        for j in tqdm.tqdm(range(minFile)):
            currMix = WF.ReadFile(pathMix[j])
            currPr = WF.ReadFile(pathPr[j])
            
            # Paso el test
            mix_pValues[i,j] = tests[i](currMix)[0] 
            pr_pValues[i,j] = tests[i](currPr)[0]
            
        # Compruebo cual es el maximo pvalue para hacer el dominio de las distribuciones
        maxMix = np.max([mix_pValues[i,:],pr_pValues[i,:]])
        # Calculo las distribuciones
        tMatrix[i] = np.linspace(0,maxMix,nPoints)
        Gpr[i,:] = Compute_Distribution(tMatrix[i],pr_pValues[i,:])
        GMix[i,:] = Compute_Distribution(tMatrix[i],mix_pValues[i,:])
        
        # Calculo KS
        KSResults[i] = np.max(abs(Gpr[i,:]-GMix[i,:]))
        passArr[i] = KSResults[i] < K
        
    # Trabajo a parte los casos de random excursions
    print("Running Random Excursions Test...")
    # Random excursion
    nRd = 8 # Los random excursions son 8 tests
    mixRD = np.zeros((nRd,minFile))
    prRD = np.zeros((nRd,minFile))
    tMatrixRd = np.zeros(nRd, dtype = object)
    KSRD = np.zeros(nRd)
    passRD = np.zeros(nRd,dtype = bool)
    
    GprRD = np.zeros((nRd,nPoints))
    GMixRD = np.zeros((nRd,nPoints))
    
    for j in tqdm.tqdm(range(minFile)):
        currMix = WF.ReadFile(pathMix[j])
        currPr = WF.ReadFile(pathPr[j])
        
        # Paso el test
        mixAux = ret.random_excursions_test(currMix) #[(... , pvalue,...), (...)]
        prAux = ret.random_excursions_test(currPr)
        for i in range(nRd):
            mixRD[i,j] = mixAux[i][3] 
            prRD[i,j] = prAux[i][3]

    for i in range(nRd):
        # Compruebo cual es el maximo pvalue para hacer el dominio de las distribuciones
        maxMix = np.max([mixRD[i,:],prRD[i,:]])
        # Calculo las distribuciones
        tMatrixRd[i] = np.linspace(0,maxMix,nPoints)
        GprRD[i,:] = Compute_Distribution(tMatrixRd[i],prRD[i,:])
        GMixRD[i,:] = Compute_Distribution(tMatrixRd[i],mixRD[i,:])
        
        # Calculo KS
        KSRD[i] = np.max(abs(GprRD[i,:]-GMixRD[i,:]))
        passRD[i] = KSRD[i] < K    
    
    
    # Variant 
    
    print("Running Random Excursions Variant Test...")
    
    nVar = 18 # Los variant corresponden a 18 tests 
    mixVar = np.ones((nVar,minFile))*2 # Se inicializan en 2 para no tener en cuenta los pasos que no existan
    prVar = np.ones((nVar,minFile))*2
    tMatrixVar = np.zeros(nVar,dtype = object)
    KSVar = np.zeros(nVar)
    passVar = np.zeros(nVar,dtype = bool)
    
    GprVar = np.zeros((nVar,nPoints))
    GMixVar = np.zeros((nVar,nPoints))
    
    for j in tqdm.tqdm(range(minFile)):
        currMix = WF.ReadFile(pathMix[j])
        currPr = WF.ReadFile(pathPr[j])
        
        # Paso el test
        mixAux = ret.variant_test(currMix) #[(... , pvalue,...), (...)]
        prAux = ret.variant_test(currPr)
        nMix = len(mixAux)
        nPr = len(prAux)
        for mx in range(nMix):
            pos = int(mixAux[mx][0][:-2]) # elimino el .0
            if pos > 0:
                index = pos + 8 # Se suma uno menos por quitar el 0 de la correspondencia
            else:
                index = pos + 9 # Pasamos de posiciones a indices
            # print(pos,index)
            mixVar[index,j] = mixAux[mx][3] 
        for pr in range(nPr):
            pos = int(prAux[pr][0][:-2])
            if pos > 0: 
                index = pos + 8 
            else:
                index = pos + 9
            prVar[index,j] = prAux[pr][3]
            
    for i in range(nVar):
        # Compruebo cual es el maximo pvalue para hacer el dominio de las distribuciones
        maxMix = np.max([mixVar[i,:],prVar[i,:]])
        # Calculo las distribuciones
        tMatrixVar[i] = np.linspace(0,maxMix,nPoints)
        GprVar[i,:] = Compute_Distribution(tMatrixVar[i],prVar[i,:])
        GMixVar[i,:] = Compute_Distribution(tMatrixVar[i],mixVar[i,:])
        
        # Calculo KS
        KSVar[i] = np.max(abs(GprVar[i,:]-GMixVar[i,:]))
        passVar[i] = KSVar[i] < K
    
    print("Done! Saving Results...")
    # Escribo en un fichero las distancas
    msg = []
    msg.append("SECOND ORDEN TESTS RESULTS (Based on KS distance and Shen's proposal)")
    msg.append("Files from directories:")
    msg.append(f"\t RNG: {dirPr}")
    msg.append(f"\t Mixed: {dirMix}")
    msg.append(f"epsilon = {epsilon} --> K = {K}")
    msg.append("")
    
    msg.append("")
    
    for i in range(nTests):
        msg.append(f"\t (o) {testsNames[i]}: {KSResults[i]} -- {passArr[i]}")
    
    # Run Test
    numbersRd = ["-4","-3","-2","-1","+1","+2","+3","+4"]
    msg.append("\t (o) Random Excursions Test:")
    for i in range(8):
        msg.append(f"\t \t {numbersRd[i]}: {KSRD[i]} -- {passRD[i]}")
    
    # Variant
    numbersVar = ["-9","-8","-7","-6","-5","-4","-3","-2","-1",
               "+1", "+2", "+3", "+4", "+5", "+6", "+7", "+8", "+9"]
    msg.append("\t (o) Random Excursions Variant Test:")
    for i in range(17):
        msg.append(f"\t \t {numbersVar[i]}: {KSVar[i]} -- {passVar[i]}")

    WF.Write_Analysis_Files("KS", msg, directory = outDir)
    
    
    # Escribo en una carpeta las distribuciones obtenidas
    WF.Write_Distribution_Files(tMatrix, Gpr, testsNames, directory = "KS_Results/Distributions/Pr")
    WF.Write_Distribution_Files(tMatrix, GMix, testsNames, directory = "KS_Results/Distributions/Mixed")
    # Escribo RunTest
    WF.Write_Distribution_Files(tMatrixRd, GprRD, numbersRd, directory = "KS_Results/Distributions/Pr/Run_Test")
    WF.Write_Distribution_Files(tMatrixRd, GMixRD, numbersRd, directory = "KS_Results/Distributions/Mixed/Run_Test")
    # Escribo test variant
    numbersVar = ["-9","-8","-7","-6","-5","-4","-3","-2","-1",
               "+1", "+2", "+3", "+4", "+5", "+6", "+7", "+8", "+9"]
    WF.Write_Distribution_Files(tMatrixVar,GprVar,numbersVar,directory = "KS_Results/Distributions/Pr/Variant_Test")
    WF.Write_Distribution_Files(tMatrixVar, GMixVar, numbersVar, directory = "KS_Results/Distributions/Mixed/Variant_Test")
    return 


#%%% Tests estadisticos NIST
def Construct_Sol_Message(testName,p,passMsg):
    return f"\t \t \t (o) {testName}: " + str(p) +" -- " +  str(passMsg)

def Binary_Matrix_Test(binary_data): # Implementado por Lucas Hernández Bellón
    # Como en NIST, usaremos M = Q = 32 (32 filas, 32 columnas)
    M = 32
    Q = 32
    nData = len(binary_data)
    # Dada la longitud de los datos, buscamos el numero de matrices que generaremos
    N = nData//(M**2)
    # Compruebo cuantos bits descarto
    # nLeft = nData%(M**2)
    # Construyo las matrices fila a fila
    #print(len(dataCrop))
    # Contador de cada rango
    rankM = 0
    rankM1 = 0
    rankLessM1 = 0
    for i in range(N):
        auxMatrix = np.zeros((M,Q),dtype = int)
        for j in range(M): # Relleno fila a fila
            for k in range(Q):    
                #print(i,j,k)
                #print(M*(i + j) + k)
                auxMatrix[j,k] = binary_data[M*Q*i + M*j + k]
                
        
        # Computo los rangos, clasificandolo en r = N, r = N-1 y r < N-1
        matrix = bm.BinaryMatrix(auxMatrix,M,Q)
        rank = matrix.compute_rank()
        if rank == M: rankM += 1
        elif rank == M-1: rankM1 += 1
        else: rankLessM1 += 1 
    # Computo el p-value del test con X^2. 
    """
    Notese que los valores de probabilidades teoricos se evaluan con:
        pM simeq prod_{j = 1}^{\infty} (1-2^{-j})
        p(M-1) simeq 2 pm ... 
        y sucesivo. Solo funciona si M = Q > 10
        Solo contamos los tres rangos primeros. Se podrian tener en cuenta mas para ser mas precisos, pero con 
        M y Q lo suficientemente grandes no hace falta (porque tendran valores < 0.005)
    """
    X2 = (rankM - 0.2888*N)**2/(0.2888*N) + (rankM1 - 0.5776*N)**2/(0.5776*N) + (rankLessM1 - 0.1336*N)**2/(0.13336*N)
    pval = np.exp(-X2/2)
    # Bool para ver si se ha superado el test
    #print(rankM,rankM1,rankLessM1)
    return pval, pval > epsilon

def Block_Frequency_Test(data):
    data_size = len(data) # Para elegir cuantos bloques
    """
    Eleccion de bloques:
        1. Deben ser de tamanno mayor a 20
        2. No debe haber mas de 100 bloques
    """
    block_size = max(data_size//99,21) 
    n_blocks = data_size//block_size
    return ft.block_frequency(data,block_size = block_size) # Segun Nist el numero de bloques no deberia ser mayor a 100. No se por que    

def NIST_Battery(QC,filePath):
    """
    Ejecutamos toda la bateria de tests programada por steve
    en la data de fileName, y guardamos los resultados en el txt asociado.
    """
    message = ["\t --> NIST Battery Test:"]
    # Data es un string directamente. Podria almacenar hasta un tamanno de 10^9. Si quiero mas, mejor tirar por bitarray.
    data = WF.ReadFile(filePath)
    #%%%% TESTS DE FRECUENCIAS EN CADENAS
    message.append("\t \t --> FREQUENCY IN SEQUENCES")
    # Monobit
    p,passMsg = ft.monobit_test(data)
    message.append(Construct_Sol_Message("Monobit_Test", p, passMsg))
    
    # Within a block
    data_size = len(data) # Para elegir cuantos bloques
    """
    Eleccion de bloques:
        1. Deben ser de tamanno mayor a 20
        2. No debe haber mas de 100 bloques
    """
    block_size = max(data_size//99,21) 
    n_blocks = data_size//block_size # Solo para el mensaje del txt 
    
    p,passMsg = Block_Frequency_Test(data)
    message.append(Construct_Sol_Message(f"Within a Block (used {block_size} bits for each block. Used {n_blocks} blocks)", p, passMsg))
    
    # runTest
    p,passMsg = rt.run_test(data)
    message.append(Construct_Sol_Message("runTest",p,passMsg))
    
    # Longest One Block
    p,passMsg = rt.longest_one_block_test(data)
    message.append(Construct_Sol_Message("Longest One Block",p,passMsg))
    
    # Discrete Fourier Transform (Spectral) Tests
    message.append("\t \t SPECTRAL TESTS")
    p,passMsg = st.spectral_test(data)
    message.append(Construct_Sol_Message("Spectral Test",p,passMsg))
    
    #%%%% TESTS DE MECANICA ESTADISTICA 
    message.append("\t \t STATISTICAL MECHANICS TESTS")
    
    # Approximate Entropy Test
    p,passMsg = aet.approximate_entropy_test(data)
    message.append(Construct_Sol_Message("Approximate Entropy Test",p,passMsg))
   
    # Cusum Test
    p,passMsg = cst.cumulative_sums_test(data)
    message.append(Construct_Sol_Message("Cusum Test",p,passMsg))
    
    # Random Excursions Test
    result = ret.random_excursions_test(data) # [pos,xVal,xObs, pos_pValue,passed]
    message.append("\t \t \t (o) Random Excursions Test:")
    for i in range(len(result)):
        message.append(f"\t \t \t \t {result[i][0]}: {result[i][3]} -- {result[i][4]}")
    
    """
    NOTA: POR ALGUNA RAZON SI NO HAY UNOS 10**7 BITS VARIANT TEST NO LLEGA DE -9 A 9
    """
    # Random Excursions Variant Test
    result = ret.variant_test(data) 
    # print(result)
    message.append("\t \t \t (o) Random Excursions Variant Test:")
    for i in range(len(result)):
        message.append(f"\t \t \t \t {result[i][0]}: {result[i][3]} -- {result[i][4]}")
    
    #%%%% TESTS DE COMPLEJIDAD DE KOLMOGOROV
    message.append("\t \t KOLMOGOROV COMPLEXITY TEST")
    # Maurer Universal Statistical Test
    p,passMsg = ut.statistical_test(data)
    message.append(Construct_Sol_Message("Universal Maurer Test",p,passMsg))
    
    # Binary Matrix Tests.
    """
    La parte de Kho esta incompleta. Falta construir las submatrices a partir de los bits
    y comprobar el rango de cada submatriz como criterio del test
    Decido hacerlo a mano con numpy
    """
    p,passMsg = Binary_Matrix_Test(data)
    message.append(Construct_Sol_Message("Binary Matrix Test",p,passMsg))

    # Linear Complexity Test
    p,passMsg = ct.linear_complexity_test(data)
    message.append(Construct_Sol_Message("Linear Complexity Test",p,passMsg))
    
    
    # Non-Overlapping Tests (Comprobar secuencias sospechosas) NO LO INCLUIMOS PORQUE NO TENEMOS SECUENCIAS SOPECHOSAS
    # Overlapping Test (IDEM)

    #%%%% END MESSAGE
    message.append("--> END NIST BATTERY TEST ----------")
    WF.Write_Analysis_Files(QC,message,showEndMessage=False)
    return

#%% MEAN RESULTS
def Get_Summary_Results(dataDir, outDir="Data_Summary"):
    """
    Get_Mean_Results:
        Dado un directorio con ficheros de los resultados de todos los tests,
        obtiene un valor medio de los resultados para el TFG
        Tomese IBM_Results como fichero de lectura de ejemplo para ver
        como funciona exactamente
    """
    dataFiles = WF.Scan_Dir(dataDir)
    for i in tqdm.tqdm(range(len(dataFiles))):
        currFile = dataFiles[i]
        currComp = os.path.basename(currFile).split("_")[0]
        
        #%%% DataLists
        minEntropyList = []
        entropyList = []
        
        
        maximumComplexityList = []
        DeficiencyFunctionList = []
        ComplexPValueList = []
        
        monobitList = []
        withinList = []
        runList =[]
        longestList = []
        
        spectralList = []
        approxEntList = []
        cusumList = []
        
        # Esto no tendria que haberlo hecho asi, la proxima vez un for
        rdExc0 = [] # -4
        rdExc1 = [] # -3
        rdExc2 = [] # -2
        rdExc3 = [] # -1
        rdExc4 = [] # +1
        rdExc5 = [] # +2
        rdExc6 = [] # +3
        rdExc7 = [] # +4
        
        varExc0 = [] # -9
        varExc1 = [] # -8
        varExc2 = [] # -7
        varExc3 = [] # -6
        varExc4 = [] # -5
        varExc5 = [] # -4
        varExc6 = [] # -3
        varExc7 = [] # -2
        varExc8 = [] # -1
        varExc9 = [] # +1
        varExc10 = [] # +2
        varExc11 = [] # +3
        varExc12 = [] # +4
        varExc13 = [] # +5
        varExc14 = [] # +6
        varExc15 = [] # +7
        varExc16 = [] # +8
        varExc17 = [] # +9
        
        maurerList = []
        matrixList = []
        linearList = []
        
        #%%% Seccion Lectura datos
        bitEndPos = -5
        with open(currFile,"r") as file:
            print(currFile)
            currLine = file.readline()
            while currLine != "":
                if "Entropy Estimations" in currLine: # Primera seccion
                    currLine = file.readline() 
                    aux = float(currLine.split("=")[1][:bitEndPos]) # Quito la palabra bits
                    minEntropyList.append(aux) 
                    
                    currLine = file.readline()
                    aux = float(currLine.split("=")[1][:bitEndPos]) # Quito la palabra bits
                    entropyList.append(aux)
                    
                    # Pasamos al bloque de complejidad
                    file.readline()
                    file.readline() # Estamos en linea "Complexity Estimations:"
                    
                    currLine = file.readline()
                    aux = float(currLine.split("=")[1][:bitEndPos]) # Quito la palabra bits
                    maximumComplexityList.append(aux)
                    
                    currLine = file.readline()
                    aux = float(currLine.split("=")[1]) # Quito la palabra bits

                    DeficiencyFunctionList.append(aux)
                    
                    currLine = file.readline()
                                        
                    aux = float(currLine.split("=")[1]) # Quito la palabra bits
                    ComplexPValueList.append(aux)
                    
                    # Pasamos a Tests NIST
                    file.readline()
                    file.readline()
                    file.readline()
                    file.readline() # Estamos en linea "FREQUENCY IN SEQUENCES"
                    
                    currLine = file.readline()
                    #print(currLine)
                    aux = float(currLine.split(":")[1].split("--")[0])
                    monobitList.append(aux)
                    
                    currLine = file.readline()
                    aux = float(currLine.split(":")[1].split("--")[0])
                    withinList.append(aux)
                    
                    currLine = file.readline()
                    aux = float(currLine.split(":")[1].split("--")[0])
                    runList.append(aux)
                    
                    currLine = file.readline()
                    aux = float(currLine.split(":")[1].split("--")[0])
                    longestList.append(aux)
                    
                    file.readline() # Linea "Spectral Tests"
                    currLine = file.readline()
                    aux = float(currLine.split(":")[1].split("--")[0])
                    spectralList.append(aux)
                    
                    file.readline() # Linea "Statistical Mechanics Tests"
                    currLine = file.readline()
                    aux = float(currLine.split(":")[1].split("--")[0])
                    approxEntList.append(aux)
                    
                    currLine = file.readline()
                    aux = float(currLine.split(":")[1].split("--")[0])
                    cusumList.append(aux)
                    
                    # Region Random Excursions
                    file.readline()
                    currLine = file.readline()
                    aux = float(currLine.split(":")[1].split("--")[0])
                    rdExc0.append(aux)
                    
                    currLine = file.readline()
                    aux = float(currLine.split(":")[1].split("--")[0])
                    rdExc1.append(aux)
                    
                    currLine = file.readline()
                    aux = float(currLine.split(":")[1].split("--")[0])
                    rdExc2.append(aux)
                    
                    currLine = file.readline()
                    aux = float(currLine.split(":")[1].split("--")[0])
                    rdExc3.append(aux)
                    
                    currLine = file.readline()
                    aux = float(currLine.split(":")[1].split("--")[0])
                    rdExc4.append(aux)
                    
                    currLine = file.readline()
                    aux = float(currLine.split(":")[1].split("--")[0])
                    rdExc5.append(aux)
                    
                    currLine = file.readline()
                    aux = float(currLine.split(":")[1].split("--")[0])
                    rdExc6.append(aux)
                    
                    currLine = file.readline()
                    aux = float(currLine.split(":")[1].split("--")[0])
                    rdExc7.append(aux)

                    # Seccion Variant:
                    file.readline()
                    # A veces el test variant no llega de -9 a 9, hay que ir con cuidado
                    currLine = file.readline() 
                    if "-9" in currLine: 
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc0.append(aux)
                    
                        currLine = file.readline()
                    if "-8" in currLine:
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc1.append(aux)
                    
                        currLine = file.readline()
                    if "-7" in currLine:
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc2.append(aux)
                    
                        currLine = file.readline()
                    if "-6" in currLine:
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc3.append(aux)
                        
                        currLine = file.readline()
                    if "-5" in currLine:
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc4.append(aux)
                    
                        currLine = file.readline()
                    if "-4" in currLine:
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc5.append(aux)
                    
                        currLine = file.readline()
                    if "-3" in currLine:
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc6.append(aux)
                    
                        currLine = file.readline()
                    if "-2" in currLine:
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc7.append(aux)
                    
                        currLine = file.readline()
                    if "-1" in currLine:
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc8.append(aux)
                    
                        currLine = file.readline()
                    if "+1" in currLine:
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc9.append(aux)
                    
                        currLine = file.readline()
                    if "+2" in currLine:
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc10.append(aux)
                    
                        currLine = file.readline()
                    if "+3" in currLine:
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc11.append(aux)
                        
                        currLine = file.readline()
                    if "+4" in currLine:
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc12.append(aux)
                        
                        currLine = file.readline()
                    if "+5" in currLine:
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc13.append(aux)
                        
                        currLine = file.readline()
                    if "+6" in currLine:
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc14.append(aux)
                        
                        currLine = file.readline()
                    if "+7" in currLine:
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc15.append(aux)
                        
                        currLine = file.readline()
                    if "+8" in currLine:    
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc16.append(aux)
                        
                        currLine = file.readline()
                    if "+9" in currLine:
                        aux = float(currLine.split(":")[1].split("--")[0])
                        varExc17.append(aux)
                        
                        file.readline() 
                    # Fin region larga
                    
                    # Estamos en linea "KOLMOGOROV COMPLEXITY TEST"
                    currLine = file.readline()
                    aux = float(currLine.split(":")[1].split("--")[0])
                    maurerList.append(aux)
                    
                    currLine = file.readline()
                    aux = float(currLine.split(":")[1].split("--")[0])
                    matrixList.append(aux)
                    
                    currLine = file.readline()
                    aux = float(currLine.split(":")[1].split("--")[0])
                    linearList.append(aux)
                    
                # Vamos al siguiente bloque
                currLine = file.readline()
        #%%% Seccion Calculos estadisticos
        nData = len(entropyList) # n de ficheros de este ordenador
        sqrData = np.sqrt(nData)
        # Ya hemos leido todos los datos. Procedemos a calcular medias
        meanMinH = np.mean(minEntropyList)
        errMinH = np.std(minEntropyList)/sqrData
        
        meanH = np.mean(entropyList)
        errH = np.std(minEntropyList)/sqrData
        
        meanMaxComplex = np.mean(maximumComplexityList)
        errMaxComplex = np.std(maximumComplexityList)/sqrData 
        
        meanDeficiency = np.mean(DeficiencyFunctionList)
        errDeficiency = np.std(DeficiencyFunctionList)/sqrData
        
        meanComplexPValue = np.mean(ComplexPValueList)
        errComplexPValue = np.std(ComplexPValueList)/sqrData
        
        meanMonobit = np.mean(monobitList)
        errMonobit = np.std(monobitList)/sqrData
        
        meanWithin = np.mean(withinList)
        errWithin = np.std(withinList)/sqrData
        
        meanRun = np.mean(runList)
        errRun = np.std(runList)/sqrData
        
        meanLongest = np.mean(longestList)
        errLongest = np.std(longestList)/sqrData
        
        meanSpectral = np.mean(spectralList)
        errSpectral = np.std(spectralList)/sqrData
        
        meanApprox = np.mean(approxEntList)
        errApprox = np.std(approxEntList)/sqrData
        
        meanCusum = np.mean(cusumList)
        errCusum = np.std(cusumList)/sqrData
        
        meanRD0 = np.mean(rdExc0)
        errRD0 = np.std(rdExc0)/sqrData
        
        meanRD1 = np.mean(rdExc1)
        errRD1 = np.std(rdExc1)/sqrData
        
        meanRD2 = np.mean(rdExc2)
        errRD2 = np.std(rdExc2)/sqrData
        
        meanRD3 = np.mean(rdExc3)
        errRD3 = np.std(rdExc3)/sqrData
        
        meanRD4 = np.mean(rdExc4)
        errRD4 = np.std(rdExc4)/sqrData
        
        meanRD5 = np.mean(rdExc5)
        errRD5 = np.std(rdExc5)/sqrData
        
        meanRD6 = np.mean(rdExc6)
        errRD6 = np.std(rdExc6)/sqrData
        
        meanRD7 = np.mean(rdExc7)
        errRD7 = np.std(rdExc7)/sqrData
        
        meanVar0 = np.mean(varExc0)
        errVar0 = np.std(varExc0)/sqrData
         
        meanVar1 = np.mean(varExc1)
        errVar1 = np.std(varExc1)/sqrData
        
        meanVar2 = np.mean(varExc2)
        errVar2 = np.std(varExc2)/sqrData
        
        meanVar3 = np.mean(varExc3)
        errVar3 = np.std(varExc3)/sqrData
         
        meanVar4 = np.mean(varExc4)
        errVar4 = np.std(varExc4)/sqrData
        
        meanVar5 = np.mean(varExc5)
        errVar5 = np.std(varExc5)/sqrData
        
        meanVar6 = np.mean(varExc6)
        errVar6 = np.std(varExc6)/sqrData
        
        meanVar7 = np.mean(varExc7)
        errVar7 = np.std(varExc7)/sqrData
        
        meanVar8 = np.mean(varExc8)
        errVar8 = np.std(varExc8)/sqrData
        
        meanVar9 = np.mean(varExc9)
        errVar9 = np.std(varExc9)/sqrData
        
        meanVar10 = np.mean(varExc10)
        errVar10 = np.std(varExc10)/sqrData
        
        meanVar11 = np.mean(varExc11)
        errVar11 = np.std(varExc11)/sqrData
        
        meanVar12 = np.mean(varExc12)
        errVar12 = np.std(varExc12)/sqrData
        
        meanVar13 = np.mean(varExc13)
        errVar13 = np.std(varExc13)/sqrData
        
        meanVar14 = np.mean(varExc14)
        errVar14 = np.std(varExc14)/sqrData
        
        meanVar15 = np.mean(varExc15)
        errVar15 = np.std(varExc15)/sqrData
        
        meanVar16 = np.mean(varExc16)
        errVar16 = np.std(varExc16)/sqrData
        
        meanVar17 = np.mean(varExc17)
        errVar17 = np.std(varExc17)/sqrData
        
        meanMaurer = np.mean(maurerList)
        errMaurer = np.std(maurerList)/sqrData
        
        meanMatrix = np.mean(matrixList)
        errMatrix = np.std(matrixList)/sqrData
        
        meanLinear = np.mean(linearList)
        errLinear = np.std(linearList)/sqrData
        
        #%%% Seccion de escritura
        msg = []
        msg.append(f"{currComp} Data Summary: Used {nData} files")
        
        msg.append("\t--> Entropy Estimations:")
        msg.append(f"\t\t --> H_min = {ufloat(meanMinH,errMinH)} bits")
        msg.append(f"\t\t\t --> Worst H_min = {np.min(minEntropyList)} bits")
        msg.append(f"\t\t --> H = {ufloat(meanH,errH)} bits")
        msg.append(f"\t\t\t --> Worst H = {np.min(entropyList)} bits")
        
        msg = msg + ["\t --> Complexity Estimations:",f"\t\t --> Maximum Complexity = {ufloat(meanMaxComplex,errMaxComplex)} bits", 
                     f"\t \t \t --> Worst Complexity =  {np.min(maximumComplexityList)} bits",
                    f"\t \t --> Deficiency Function = {ufloat(meanDeficiency,errDeficiency)}",
                    f"\t \t\t --> Worst Deficiency Function = {np.max(DeficiencyFunctionList)}",
                    f"\t \t --> p-value = {ufloat(meanComplexPValue,errComplexPValue)} -- {str(meanComplexPValue >= epsilon)}",
                    f"\t \t\t --> Worst p-value = {np.min(ComplexPValueList)} --  {np.min(ComplexPValueList) >= epsilon}"]
        
        msg.append("\t --> NIST Battery Test:")
        
        msg.append("\t \t --> FREQUENCY IN SEQUENCES")
        msg.append(Construct_Sol_Message("Monobit_Test", ufloat(meanMonobit,errMonobit), meanMonobit >= epsilon))
        msg.append(f"\t \t \t \t Worst p-value: {np.min(monobitList)} --  {np.min(monobitList) >= epsilon}")
        msg.append(Construct_Sol_Message("Within a Block", ufloat(meanWithin,errWithin), meanWithin >= epsilon))
        msg.append(f"\t \t \t \t Worst p-value: {np.min(withinList)} --  {np.min(withinList) >= epsilon}")
        msg.append(Construct_Sol_Message("Run Test", ufloat(meanRun,errRun), meanRun >= epsilon))
        msg.append(f"\t \t \t \t Worst p-value: {np.min(runList)} --  {np.min(runList) >= epsilon}")
        msg.append(Construct_Sol_Message("Longest One Block",ufloat(meanLongest,errLongest),meanLongest >= epsilon))
        msg.append(f"\t \t \t \t Worst p-value: {np.min(longestList)} --  {np.min(longestList) >= epsilon}")
        
        msg.append("\t \t SPECTRAL TESTS")
        msg.append(Construct_Sol_Message("Spectral Test",ufloat(meanSpectral,errSpectral),meanSpectral >= epsilon))
        msg.append(f"\t \t \t \t Worst p-value: {np.min(spectralList)} --  {np.min(spectralList) >= epsilon}")
        msg.append("\t \t STATISTICAL MECHANICS TESTS")
        msg.append(Construct_Sol_Message("Approximate Entropy Test",ufloat(meanApprox,errApprox),meanApprox >= epsilon))
        msg.append(f"\t \t \t \t Worst p-value: {np.min(approxEntList)} --  {np.min(approxEntList) >= epsilon}")
        msg.append(Construct_Sol_Message("Cusum Test",ufloat(meanCusum,errCusum),meanCusum >= epsilon))
        msg.append(f"\t \t \t \t Worst p-value: {np.min(cusumList)} --  {np.min(cusumList) >= epsilon}")
        
        msg.append("\t \t \t (o) Random Excursions Test:")
        msg.append(f"\t \t \t \t -4: {ufloat(meanRD0,errRD0)} -- {meanRD0 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {np.min(rdExc0)} --  {np.min(rdExc0) >= epsilon}")
        msg.append(f"\t \t \t \t -3: {ufloat(meanRD1,errRD1)} -- {meanRD1 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {np.min(rdExc1)} --  {np.min(rdExc1) >= epsilon}")
        msg.append(f"\t \t \t \t -2: {ufloat(meanRD2,errRD2)} -- {meanRD2 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {np.min(rdExc2)} --  {np.min(rdExc2) >= epsilon}")
        msg.append(f"\t \t \t \t  -1: {ufloat(meanRD3,errRD3)} -- {meanRD3 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {np.min(rdExc3)} --  {np.min(rdExc3) >= epsilon}")
        msg.append(f"\t \t \t \t +1: {ufloat(meanRD4,errRD4)} -- {meanRD4 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {np.min(rdExc4)} --  {np.min(rdExc4) >= epsilon}")
        msg.append(f"\t \t \t \t +2: {ufloat(meanRD5,errRD5)} -- {meanRD5 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {np.min(rdExc5)} --  {np.min(rdExc5) >= epsilon}")
        msg.append(f"\t \t \t \t +3: {ufloat(meanRD6,errRD6)} -- {meanRD6 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {np.min(rdExc6)} --  {np.min(rdExc6) >= epsilon}")
        msg.append(f"\t \t \t \t +4: {ufloat(meanRD7,errRD7)} -- {meanRD7 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {np.min(rdExc7)} --  {np.min(rdExc7) >= epsilon}")
        
        msg.append("\t \t \t (o) Random Excursions Variant Test:")
        msg.append(f"\t \t \t \t -9: {ufloat(meanVar0,errVar0)} -- {meanVar0 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc0,default = np.nan)} --  {min(varExc0,default = np.nan) >= epsilon}")
        msg.append(f"\t \t \t \t -8: {ufloat(meanVar1,errVar1)} -- {meanVar1 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc1,default = np.nan)} --  {min(varExc1,default = np.nan) >= epsilon}")
        msg.append(f"\t \t \t \t -7: {ufloat(meanVar2,errVar2)} -- {meanVar2 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc2,default = np.nan)} --  {min(varExc2,default = np.nan) >= epsilon}")
        msg.append(f"\t \t \t \t -6: {ufloat(meanVar3,errVar3)} -- {meanVar3 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc3,default = np.nan)} --  {min(varExc3,default = np.nan) >= epsilon}")
        msg.append(f"\t \t \t \t -5: {ufloat(meanVar4,errVar4)} -- {meanVar4 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc4,default = np.nan)} --  {min(varExc4,default = np.nan) >= epsilon}")
        msg.append(f"\t \t \t \t -4: {ufloat(meanVar5,errVar5)} -- {meanVar5 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc5,default = np.nan)} --  {min(varExc5,default = np.nan) >= epsilon}")
        msg.append(f"\t \t \t \t -3: {ufloat(meanVar6,errVar6)} -- {meanVar6 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc6,default = np.nan)} --  {min(varExc6,default = np.nan) >= epsilon}")
        msg.append(f"\t \t \t \t -2: {ufloat(meanVar7,errVar7)} -- {meanVar7 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc7,default = np.nan)} --  {min(varExc7,default = np.nan) >= epsilon}")
        msg.append(f"\t \t \t \t -1: {ufloat(meanVar8,errVar8)} -- {meanVar8 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc8,default = np.nan)} --  {min(varExc8,default = np.nan) >= epsilon}")
        msg.append(f"\t \t \t \t +1: {ufloat(meanVar9,errVar9)} -- {meanVar9 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc9,default = np.nan)} --  {min(varExc9,default = np.nan) >= epsilon}")
        msg.append(f"\t \t \t \t +2: {ufloat(meanVar10,errVar10)} -- {meanVar10 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc10,default = np.nan)} --  {min(varExc10,default = np.nan) >= epsilon}")
        msg.append(f"\t \t \t \t +3: {ufloat(meanVar11,errVar11)} -- {meanVar11 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc11,default = np.nan)} --  {min(varExc11,default = np.nan) >= epsilon}")
        msg.append(f"\t \t \t \t +4: {ufloat(meanVar12,errVar12)} -- {meanVar12 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc12,default = np.nan)} --  {min(varExc12,default = np.nan) >= epsilon}")
        msg.append(f"\t \t \t \t +5: {ufloat(meanVar13,errVar13)} -- {meanVar13 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc13,default = np.nan)} --  {min(varExc13,default = np.nan) >= epsilon}")
        msg.append(f"\t \t \t \t +6: {ufloat(meanVar14,errVar14)} -- {meanVar14 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc14,default = np.nan)} --  {min(varExc14,default = np.nan) >= epsilon}")
        msg.append(f"\t \t \t \t +7: {ufloat(meanVar15,errVar15)} -- {meanVar15 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc15,default = np.nan)} --  {min(varExc15,default = np.nan) >= epsilon}")
        msg.append(f"\t \t \t \t +8: {ufloat(meanVar16,errVar16)} -- {meanVar16 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc16,default = np.nan)} --  {min(varExc16,default = np.nan) >= epsilon}")
        msg.append(f"\t \t \t \t +9: {ufloat(meanVar17,errVar17)} -- {meanVar17 >= epsilon}")
        msg.append(f"\t \t \t \t \t Worst p-value: {min(varExc17,default = np.nan)} --  {min(varExc17,default = np.nan) >= epsilon}")
        
        msg.append("\t \t KOLMOGOROV COMPLEXITY TEST")
        msg.append(Construct_Sol_Message("Universal Maurer Test",ufloat(meanMaurer,errMaurer),meanMaurer >= epsilon))
        msg.append(f"\t \t \t \t Worst p-value: {np.min(maurerList)} --  {np.min(maurerList) >= epsilon}")
        msg.append(Construct_Sol_Message("Binary Matrix Test",ufloat(meanMatrix,errMatrix),meanMatrix >= epsilon))
        msg.append(f"\t \t \t \t Worst p-value: {np.min(matrixList)} --  {np.min(matrixList) >= epsilon}")
        msg.append(Construct_Sol_Message("Linear Complexity Test",ufloat(meanLinear,errLinear),meanLinear >= epsilon))
        msg.append(f"\t \t \t \t Worst p-value: {np.min(linearList)} --  {np.min(linearList) >= epsilon}")
        
        WF.Write_Analysis_Files(currComp,msg,directory=outDir,showEndMessage=False)
    return
    