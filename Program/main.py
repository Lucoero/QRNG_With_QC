# -*- coding: utf-8 -*-
"""
Main
"""
#%% Librerias
import numpy as np
import matplotlib.pyplot as plt
import WriteFile as WF
import MyCircuits as Circuits
import CallComputers as Computers
import ExecuteCircuits as Execute
import time as time
import tqdm as tqdm
import warnings as warnings
warnings.filterwarnings("ignore")
plt.close("all")
print("************")
print("MAIN STARTED")
print("************")
#%% Variables
N = 1*10**7 # Numero de bits que vamos a generar
simAquila = False # Elige si usar el simulador o no (Debugging) # Parece que no hay bug: Por Rydberg blockade mayor qbits los bits 0 y 1 dejan de ser equiprobables
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

nJobs = 50
#IBMBackend= "ibm_marrakesh"

IBMBackend = "ibm_kingston" 
IonQBackend = "Forte 1" # us-east-1
QuEraBackend = "Aquila" # us-east-1 # OJO, AQUILA USA OTRO PARADIGMA, se debe indicar en Run_On_AWS dejando el nombre del backend tal cual
RigettiBackend = "Cepheus-1-108Q" # 107 qbits. us-west-1
AQTBackend = "IBEX Q1" # Forte Averiado

t0 = time.time()
#%% Creacion de circuitos a usar
if N < hQbitsIBM: hQbitsIBM = N
if N < hQbitsIonQ: hQbitsIonQ = N
if N < hQbitsQuEra: hQbitsQuEra = N;
if N < hQbitsRigetti: hQbitsRigetti = N
if N < hQbitsSim: hQbitsSim = N

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
#%% Conseguimos los datos
#SimpleCount = Execute.Simulation(qc = SimpleCircuit, Nshots = N, saveData = True, countReport=False)

#IBMCount = Execute.Run_On_IBM(SimpleCircuit, IBMBackend, Nshots = 1, saveData=False)



print("RUNNING MAIN LOOP:\n")
for i in tqdm.tqdm(range(nJobs)):
    #IBMHadamard = Execute.Run_On_IBM(HIBM,IBMBackend,Nshots = ibmShots, extraName = f"_{IBMBackend}--{hQbitsIBM}", countReport=False)
    #SimHadamard = Execute.Simulation(HSim,Nshots =  simShots, saveData = True, extraName = f"_Sim--{hQbitsSim}", countReport = False) # Overflow si usas mas de 20 qbits, ojo
    #RiggettiHadamard = Execute.Run_On_AWS(HRigetti, RigettiBackend,rigShots,liRigetti,lsRigetti,extraName = f"_{RigettiBackend}--{hQbitsRigetti}", subQC = "Rigetti")
    #QuEraHadamard = Execute.Run_On_Aquila(hQbitsQuEra,QuEraShots,liQuEra,lsQuEra,extraName = f"_{QuEraBackend}--{hQbitsQuEra}", simulate=simAquila)
    #IonQHadamard = Execute.Run_On_AWS(HIonQ, IonQBackend,IonQShots,liIonQ,lsIonQ,extraName = f"_{IonQBackend}--{hQbitsIonQ}", subQC = "IonQ")
    UrandomHadamard = Execute.Run_Urandom(N,extraName = f"_CSPRNG--{N}")
    #AQTHadamard = Execute.Run_On_AWS(HAQT, AQTBackend,AQTShots,liAQT,lsAQT,extraName = f"_{AQTBackend}--{hQbitsAQT}", subQC = "AQT")
    pass

# Obtenemos bits por np Random

#%% Fin del programa
print("************")
print("PROGRAM ENDED")
print(f"TIME SPENT (s): {time.time()-t0}")
print("************")