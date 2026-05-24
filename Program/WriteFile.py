# -*- coding: utf-8 -*-
"""
WriteFile: Para escribir en ficheros los numeros
"""
import os as os
import numpy as np
import tqdm as tqdm
from decimal import Decimal
#%% Escritura de ficheros
def Write_Bitfiles(directory,origin ,results, overwrite_shots = False, shots = 0):
    if not os.path.exists(directory):
        os.makedirs(directory)
    N = np.size(results) # Results debe ser un array
    if type(results) == str:
        N = len(results)
    Nshots = N # Suponemos que los datos estan guardados a chunks de data
    if overwrite_shots: Nshots = shots
    filename = origin+f"-{Nshots:e}-0-.txt"
    path = directory+os.sep+filename
    # Vemos que no exista ya el fichero
    i = 1
    while (os.path.exists(path)) and i < 10**6:
        filename = origin+ f"-{Nshots:e}-{i}-.txt"
        path = directory+os.sep+filename
        i += 1
    try:
        with open(path,"a") as file:
            print("Saving Results...")
            if type(results) == str:
                file.write(results)
            else:
                for i in tqdm.tqdm(range(N)):
                    #print(str(bit))
                    file.write(results[i])
        print(f"\nRESULTS STORED IN: \n{path}\n")
    except:
        print(f"*****\nERROR: FILE {filename} HAS NOT BEEN CREATED IN:\n{path}\n*****")
        raise
        return "FILE_ERROR"
    return filename
def Combine_Bitfiles(inDir,outDir):
    """
    Combine_Bitfiles:
        Dado un directorio con varios ficheros de bits,
        creamos un fichero en el que concatenamos los bits 
        y lo guardamos en el outDir con el formato de nombre adecuado
        El formato de nombre se espera obtener del primer fichero de inDir.
    """
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    fileList = Scan_Dir(path = inDir) 
    nBits = 0 # inicializamos el contador de bits
    shots = 0
    name = "temp.txt"
    auxPath = outDir + os.sep + name
    with open(auxPath,"w") as outFile:
        for i in range(len(fileList)):
            with open(fileList[i],"r") as inFile:
                data = inFile.read()
                nData = len(data)
                nBits += nData
                
                auxShot = int(Decimal(os.path.basename(fileList[i]).split("--")[1].split("-")[1]))
                shots += auxShot
                outFile.write(data)
    
    # Renombramos el fichero al formato
    aux = os.path.basename(fileList[0])
    qc, rest = aux.split("--")
    qbits = rest.split('-')[0]
    origin = qc + "--" + qbits
    newName =  origin + f"-{shots:e}-0-.txt" # Suponemos que han usado el mismo numero de qbits
    path = outDir +os.sep+newName
    # Vemos que no exista ya el fichero
    i = 1
    while (os.path.exists(path)) and i < 10**6:
        newName = origin+ f"-{shots:e}-{i}-.txt"
        path = outDir+os.sep+newName
        i += 1
    os.rename(auxPath,path)
    print(f"FILES IN DIRECTORY {inDir} GOT MERGED IN FILE:\n {path}")
    print(f"Total Bits: {nBits}. Total shots: {shots}, qbits: {qbits}")
    return newName
def Delete_Directory(directory):
    """
    Para eliminar temp
    """
    fileList = Scan_Dir(path = directory) 
    for file in fileList:
        os.remove(file)
    os.rmdir(directory)
    return
def GetDirectory(path):
    words = path.split('/')[:-1]
    directory = '/'.join(words)
    return directory

def Write_Analysis_Files(QC, message, directory = "Data_Processing", showEndMessage = False):
    """
    Escribe el mensaje necesario en el path
    message corresponde a un array de strings que queremos escribir con saltos de linea
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    fileName = f"{QC}_Results.txt"
    path = directory + os.sep + fileName
    with open(path, "a") as f:
        if np.size(message) == 1: # Es solo un mensaje, no un array
            f.write(message)
        else:     
            for line in message:
                f.write(line)
                f.write("\n")
        f.write("\n")
    if showEndMessage:
        print(f"\nRESULTS STORED IN\n {path}")
    return 

def Write_Distribution_Files(tArr,GArr, nameArr, directory):
    """
    TArr debe ser una lista de arrays de la forma
    [[tArr1],[tArr2],...]
    Idem para GArr:
    [[GArr1],[GArr2],...]
    NameArr es el nombre de los ficheros:
    """
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    nFiles = len(nameArr)
    for i in tqdm.tqdm(range(nFiles)):
        filePath = directory + os.sep + nameArr[i] + ".txt"
        nCoords = len(tArr[i])
        with open(filePath,"w") as file:
            for j in range(nCoords):
                file.write(f"{tArr[i][j]} \t {GArr[i][j]}\n")
    return
    
#%% Lectura de ficheros
def Scan_Dir(path= "Outputs"):
    """
    Recolectamos todos los ficheros en carpetas y subcarpetas
    Consideramos que solo hay txts
    """
    fileList = []
    with os.scandir(path) as entries:
        for entry in entries:
            if entry.is_file():
                filename, file_extension = os.path.splitext(entry.path)
                if file_extension == ".txt": # Solo aceptamos txt
                    fileList.append(entry.path)
            elif entry.is_dir():
                fileList = fileList + Scan_Dir(path = entry.path)
    return fileList
def ReadFile(filePath, readAsBit = False):
    if readAsBit:
        readType = "rb"
    else:
        readType = "r"
    try:
        with open(filePath,readType) as f:
            data = f.read()
    except:
        print(f"ERROR: NO SE HA PODIDO leer el archivo en path \n {filePath}")
        return None
    return data

def Load_Distributions(directory):
    """
    Carga en memoria las distribuciones de KS almacenadas.
    En un directorio habra:
        -Seccion Mixed
        -Seccion Pr
    """
    
    
    return

def LoadTimes(filePath):
    times = []
    bits = []
    with open(filePath) as file:
        line = file.readline()
        while line != "":
            sides = line.split(":") # Separo la linea en lado derecho e izquierdo
            times.append(float(sides[1]))
            compData = sides[0].split("--")[1] # Elimino el nombre del computador
            aux = compData.split("-") # Separo en elementos
            bits.append(float(aux[0])*float(aux[1]))
            line = file.readline()
    times = np.array(times,dtype = float)
    bits = np.array(bits,dtype = int)
    return times,bits