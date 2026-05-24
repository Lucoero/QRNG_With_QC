# -*- coding: utf-8 -*-
"""
Get_OComplex:
    Calculamos la escala de la complejidad de los circuitos // tasa obtencion de bits
    HIPOTESIS: SERÁ LINEAL
"""
import numpy as np
import matplotlib.pyplot as plt
import WriteFile as WF
import MyCircuits as Circuits
import ExecuteCircuits as Execute
import tqdm as tqdm
import warnings as warnings
import os as os
import scipy.stats as stats
warnings.filterwarnings("ignore")
plt.close("all")

#%% Variables de entrada
startN = 10**4
endN = 9*10**4
nPoints = 10

timesOutput = "Times/Circuits_Complexity/QC"

OnlyPlots = True # Get Only plots
#%% Calculos
NArr = np.linspace(startN,endN, nPoints,dtype = int)
NArr = np.array([2*10**5], dtype = int)
# Qbits de cada computador
hQbitsIBM = 156 # Numero de qbits para el circuito de hadamard transform
hQbitsIonQ = 36
hQbitsQuEra = 256
hQbitsRigetti = 106
hQbitsSim = 20
hQbitsAQT = 12
# Cotas superiores de shots
lsIonQ = 5000
lsRigetti = 5*10**4
lsQuEra = 1000
lsAQT = 2000
# Cotas inferiores de shots
liIonQ = 100
liRigetti = 10
liQuEra = 1
liAQT = 1

nJobs = 1

IBMBackend = "ibm_kingston" 
IonQBackend = "Forte 1" # us-east-1
QuEraBackend = "Aquila" # us-east-1 # OJO, AQUILA USA OTRO PARADIGMA, se debe indicar en Run_On_AWS dejando el nombre del backend tal cual
RigettiBackend = "Cepheus-1-108Q" # 107 qbits. us-west-1
AQTBackend = "IBEX Q1" # Forte Averiado

#%%% Bucle de calculo
if not OnlyPlots:
    for i in tqdm.tqdm(range(len(NArr))):
        N = NArr[i]
        ibmShots = N//hQbitsIBM
        simShots = N//hQbitsSim
        rigShots = N//hQbitsRigetti
        QuEraShots = N//hQbitsQuEra
        IonQShots = N//hQbitsIonQ
        AQTShots = N//hQbitsAQT
        
        #SimpleCircuit = Circuits.Simple_Qbit(printCircuit=False), por eso restamos 1
        HIBM = Circuits.Hadamard_Transform(hQbitsIBM,printCircuit = False) # Esta el qbit 0
        HIonQ = Circuits.Hadamard_Transform(hQbitsIonQ,printCircuit = False)
        HRigetti = Circuits.Hadamard_Transform(hQbitsRigetti,printCircuit = False) 
        HSim = Circuits.Hadamard_Transform(hQbitsSim,printCircuit = False)
        HAQT = Circuits.Hadamard_Transform(hQbitsAQT)
    
        #IBMHadamard = Execute.Run_On_IBM(HIBM,IBMBackend,Nshots = ibmShots, extraName = f"_{IBMBackend}--{hQbitsIBM}", countReport=False,directory="Outputs/OComplexity",timeDirectory=timesOutput)
        #SimHadamard = Execute.Simulation(HSim,Nshots =  simShots, saveData = True, extraName = f"_Sim--{hQbitsSim}", countReport = False) # Overflow si usas mas de 20 qbits, ojo
        #RiggettiHadamard = Execute.Run_On_AWS(HRigetti, RigettiBackend,rigShots,liRigetti,lsRigetti,extraName = f"_{RigettiBackend}--{hQbitsRigetti}", subQC = "Rigetti",directory="Outputs/OComplexity",timeDirectory=timesOutput)
        #QuEraHadamard = Execute.Run_On_Aquila(hQbitsQuEra,QuEraShots,liQuEra,lsQuEra,extraName = f"_{QuEraBackend}--{hQbitsQuEra}", simulate=False,directory="Outputs/OComplexity",timeDirectory=timesOutput)
        #IonQHadamard = Execute.Run_On_AWS(HIonQ, IonQBackend,IonQShots,liIonQ,lsIonQ,extraName = f"_{IonQBackend}--{hQbitsIonQ}", subQC = "IonQ",directory="Outputs/OComplexity",timeDirectory=timesOutput)
        #UrandomHadamard = Execute.Run_Urandom(N,extraName = f"_CSPRNG--{N}",directory="Outputs/OComplexity",timeDirectory=timesOutput)
        #AQTHadamard = Execute.Run_On_AWS(HAQT, AQTBackend,AQTShots,liAQT,lsAQT,extraName = f"_{AQTBackend}--{hQbitsAQT}", subQC = "AQT",directory="Outputs/OComplexity",timeDirectory=timesOutput)
#%% Plots
plt.close("all")
fontsize = 10 #  text scale for plots
numbersize = [13,13] # scale for number axes
plt.rc('font', size=fontsize) # Para el escalado de texto en plots

def Axe_Canvas(Bits,Times,ax, linesNames = [],multipleLines = False):
    pSize = 10
    ax.grid()
    ax.set_ylabel("Time (s)")
    ax.set_xlabel("Bits")
    if multipleLines:
        nLines = len(Bits)
        if len(linesNames) == 0:
            linesNames = [f"Line {i+1}" for i in range(nLines)]
        for i in range(nLines):
            # Calculamos el ajuste
            pol = stats.linregress(Bits[i],Times[i])
            y = pol[1] + pol[0]*Bits[i] # Coeficientes de menor a mayor orden
            ec = f"y(x) = {round(pol[0],10)}x + {round(pol[1],4)}" + r"; $r^2$" + f"= {round(pol[2],4)}"  
            ax.scatter(Bits[i],Times[i],s = pSize)
            ax.plot(Bits[i],y,label = linesNames[i] + " (" + ec +")", linestyle = "dashed", linewidth = 1.5)
    else:
        ax.scatter(Bits,Times, s = pSize)
        pol = stats.linregress(Bits,Times)
        y = pol[1] + pol[0]*Bits # Coeficientes de menor a mayor orden
        ec = f"{round(pol[0],10)}*x + {round(pol[1],4)}" + r"; $r^2$" + f"= {round(pol[2],4)}" 
        ax.plot(Bits,y,linestyle = "dashed",label = ec)
    return
# Cargo en memoria los tiempos
pathArr = WF.Scan_Dir(timesOutput)
fig,ax = plt.subplots()
fig.suptitle("Comparison Between O(f(n)) times for each QRNG")
nFiles = len(pathArr)
bitList = []
timeList = []
namesList = []
for i in range(nFiles):
    currPath = pathArr[i]
    namesList.append(os.path.basename(currPath).split("_")[0])
    times,bits = WF.LoadTimes(currPath)
    timeList.append(times)
    bitList.append(bits)
    
Axe_Canvas(bitList,timeList,ax,linesNames = namesList, multipleLines = True)
ax.legend()
fig.show()


