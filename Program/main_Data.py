# -*- coding: utf-8 -*-
"""
main_Data

El fichero encargado de analizar todos los datos
"""
#%% Librerias
import numpy as np
import Analyse_Data as AD
import WriteFile as WF
import os as os
import tqdm as tqdm
import uncertainties as unc


#%% Variables generales
dataFolder = "Outputs/BCryptGenRandom"
#dataFolder = "Prueba"
timesFolder = "Times"

#%% Analisis de datos.
print(40*"*" + "\nSTARTING THE ANALYSIS OF THE BITFILES DATA\n" + "*"*40)
# Cogemos todos los ficheros txt
fileList = WF.Scan_Dir(dataFolder) # Ficheros de bitstrings
timeFileList = WF.Scan_Dir(timesFolder)

# Realizamos el analisis
"""
for i in tqdm.tqdm(range(len(fileList))):
    filepath = fileList[i]
    file = os.path.basename(filepath) 
    QC = file.split("_")[0]
    WF.Write_Analysis_Files(QC, f"ANALYSIS OF {file}:")
    
    Hmin = AD.Estimate_Hmin(filepath)
    
    H = AD.Estimate_H(filepath)
    
    WF.Write_Analysis_Files(QC, ["\t--> Entropy Estimations:", f"\t\t --> H_min = {Hmin} bits", f"\t\t --> H = {H} bits"], showEndMessage= False)
    
    maxC,d,probD,passed = AD.Estimate_KComplexity(filepath)
    WF.Write_Analysis_Files(QC,["\t --> Complexity Estimations:",f"\t\t --> Maximum Complexity = {maxC} bits", f"\t \t --> Deficiency Function = {d}",
                               f"\t \t --> Probability of that deficiency with Randomness Hypothesis = {probD}", f"\t \t --> Test Passed: {str(passed)}"])
    AD.NIST_Battery(QC,filepath)
"""
print(40*"*" + "\nGETTING THE MEAN OF THE COMPUTING TIME IN QCs\n" + "*"*40)
for i in tqdm.tqdm(range(len(timeFileList))):
    # Aqui vamos fichero a fichero calculando la media con su error
    #print(timeFileList[i])
    ratio,errRatio = AD.Compute_Creation_Ratio(timeFileList[i])
    bitSpeed = unc.ufloat(ratio,errRatio) # Maneja errores
    QC = os.path.basename(timeFileList[i]).split("_")[0]
    WF.Write_Analysis_Files("Times", f"{QC} Bits ratio (bits/s): {bitSpeed:2e}",directory = "Data_Processing")    
    pass
