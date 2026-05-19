# -*- coding: utf-8 -*-
"""
ExecuteCircuits
"""
import qiskit as qis
import numpy as np
import WriteFile as WF
import CallComputers as Computers

import time as time
import datetime as datetime
import braket.aws as aws # Para recuperar trabajos de aws

import tqdm as tqdm

from qiskit.primitives import StatevectorSampler
from qiskit.visualization import plot_histogram
from qiskit_ibm_runtime import QiskitRuntimeService

import os as os

import collections as collections # Extra para conteo sencillo

#%% Funciones

def Simulation(qc, Nshots, saveData = False, directory = "Outputs/Simulation", timeDirectory = "Times", extraName = "", countReport = True):
    """
    ReadSimulation:
        Given a Quanutm Circuit, simulates its outcome Nshots times and shows the results
        via terminal and histogram. Each shot result can be stored in path.
        
    Inputs:
        qc: Quantum Circuit to simulate
        Nshots: Number of simulations needed
        saveData: Choose to Save data in txt. Defaults False
        directory: directory Path
    Returns:
        counts: Dictionary that stores the counts of each output
        
    Proposal: Make the function run with multiple qc as inputs, so it can use Parallel
    provided by sampler.run.
    
    Usa PCG64
    """
    print("Running Simulation...") # Por la barra de progreso
    t0 = time.time()
    sampler = StatevectorSampler()
    job = sampler.run([qc], shots=Nshots)
    result = job.result()
    print("Job ended!")
    t1 = time.time()
    if saveData:
        shotResults = np.array(result[0].data.meas.get_bitstrings(),dtype = str) # Aqui tenemos el obtenido en cada shot
        # job metrics no se puede usar con simulaciones. Usamos time.time
        savedName = WF.Write_Bitfiles(directory,"Simulation" + extraName, shotResults)
        WF.Write_Analysis_Files("Simulation", f"{savedName}-{len(shotResults)}: {t1-t0}",timeDirectory)  
    counts = None
    print("----------------")
    if countReport:  
        counts = result[0].data.meas.get_counts()# Aqui un conteo
        print(f"Simulation result with {Nshots} shots: \n{counts}")
        plot_histogram(counts)
    print(f"Time on the task (s): {t1-t0}")
    print("----------------")
    return counts

def Run_Urandom(Nshots, saveData = True, directory = "Outputs/BCryptGenRandom", timeDirectory = "Times", extraName = ""):
    """
    Run_Urandom:
        Generamos una secuencia aleatoria de bits con PRNG criptografico del OS (en este caso, de Windows)
        
    """
    t0 = time.time()
    # Calculamos los bytes que necesitamos
    nBytes = (Nshots + 7) // 8 # Al pasarlo a bits obtendremos 8 bits por cada byte. El +7 mete un redondeo hacia arriba de la division
    # Obtenemos los bytes con BCryptGenRandom de windows (a paritr de os.urandom)
    bytesArr = os.urandom(nBytes)
    """
    # Pasamos a un array de strings de bits
    shotResults = np.zeros(nBytes,dtype = str)
    for i in range(nBytes):
        shotResults[i] = ''.join(f"{bytesArr[i]:08b}") # Puede ser que tengamos unos 8 bits mas, pero es despreciable
    """
    shotResults = ''.join(f'{b:08b}' for b in bytesArr)

    shotResults = shotResults[:Nshots]
    t1 = time.time()
    print(len(shotResults))
    #print(shotResults)
    if saveData:
        savedName = WF.Write_Bitfiles(directory,"BCryptGenRandom" + extraName , shotResults)
        WF.Write_Analysis_Files("BCryptGenRandom", f"{savedName}-{1}: {t1-t0}", timeDirectory)
    return

def Run_On_IBM(qc, backend, Nshots, saveData = True, directory = "Outputs/IBM", timeDirectory = "Times", extraName = "", countReport = False):
    t0 = time.time()
    if type(backend == str): print(f"Running on IBM backend {backend}... ")
    else: print(f"Running on IBM backend {backend.name}... ")

    job = Computers.CallIBM(qc,backend,Nshots)
    # Para ver la cola
    #jobs = QiskitRuntimeService().jobs(backend_name=backend,limit=50)

    #queued_jobs = [job for job in jobs if (job.status() == "QUEUED" or job.status() == "RUNNING")]

    #print(f"Current jobs in queue (with ours; sometimes wrong): {len(queued_jobs)}")
    print("Program will wait until job ends...")
    result = job.result()
    
    print("Job ended!")
    t1 = time.time()
    if saveData:
        shotResults = np.array(result[0].data.meas.get_bitstrings(),dtype = str) # Aqui tenemos el obtenido en cada shot
        timeSpent = job.metrics().get('usage', {})['quantum_seconds']

        savedName = WF.Write_Bitfiles(directory,"IBM" + extraName, shotResults)
        WF.Write_Analysis_Files("IBM",f"{savedName}-{len(shotResults)}-(s):{timeSpent}", directory = timeDirectory)        
    counts = None
    print("----------------")
    if countReport: 
        counts = result[0].data.meas.get_counts()# Aqui un conteo
        print(f"IBM result with {Nshots} shots: \n{counts}")
        plot_histogram(counts)
    print(f"Time on the QC (s): {timeSpent}")
    print(f"Time of the whole task (s): {t1-t0}")
    print("----------------")
    return counts    
#%% Funciones de AWS

def AWS_GetTimeSpent(jobId):
    """
    GetTimeSpent: 
        Dado un task de AWS obtenido por qiskit, recupera su objeto en AWS y calcula el tiempo total del task en segundos
    https://github.com/amazon-braket/qiskit-braket-provider/blob/main/qiskit_braket_provider/providers/braket_quantum_task.py
    """
    task = aws.AwsQuantumTask(jobId)
    result = task.result()
    metadata = result.task_metadata
    
    # No hay un valor exacto de QTime, asi que lo hacemos con el tiempo total del task quitandole el de la cola
    start = metadata.createdAt[:-1] # Elimino la "Z" del final porque el uso horario es cte
    end = metadata.endedAt[:-1]
    # Primero formateamos de tal forma que se pueda calcular rapido las diferencias
    startDate,startHours = start.split("T")
    startDate = np.array(startDate.split("-"),dtype = int)
    startHours = np.array(startHours.split(":"),dtype = float)
    
    endDate,endHours = end.split("T")
    endDate = np.array(endDate.split("-"),dtype = int)
    endHours = np.array(endHours.split(":"),dtype = float)
    
    # Obtenemos los microsegundos
    startMicroseconds = int((startHours[-1]-np.floor(startHours[-1]))*1000)
    endMicroseconds = int((endHours[-1]-np.floor(endHours[-1]))*1000)
    
    startHours[-1] = np.floor(startHours[-1])
    startHours = np.array(startHours,dtype = int)
    
    endHours[-1] = np.floor(endHours[-1])
    endHours = np.array(endHours,dtype = int)
    # Paso a libreria datetime para computar la diferencia
    # datetime usa year, month, day, hour, minute, second, microsecond , y tzinfo.
    startDateTime = datetime.datetime(startDate[0],startDate[1],startDate[2],startHours[0],startHours[1],startHours[2],startMicroseconds)
    endDateTime = datetime.datetime(endDate[0],endDate[1],endDate[2],endHours[0],endHours[1],endHours[2],endMicroseconds)
    
    timeSpent = (endDateTime-startDateTime).total_seconds()
    """
    except: 
        print("ERROR: TIME COULD NOT BE ESTIMATED")
        timeSpent = 0
    """
    return timeSpent

def Run_Task(qc,backend,Nshots, directory,extraName):
    job = Computers.CallAWS(qc,backend,Nshots) # Jobs en qiskit son las tasks de AWS
    t1 = time.time() # Tiempo para luego descartarlo del tiempo total
    print("Program will wait until job ends...")
    # Comprobamos cuando cambia de estado
    aux = 0
    while job.status() == 'QUEUED':
        if aux > 10000:
            print("#",end = "")
            aux = 0
        aux += 1
        pass
    t2 = time.time() # Se ejecutando el codigo
    result = job.result()
    print("Job ended!")
    shotResults = result.get_memory() # Aqui tenemos el obtenido en cada shot
    # Calculamos el tiempo estimado de la cola en el hardware
    queueTime = t2-t1
    timeSpent = AWS_GetTimeSpent(job.job_id())  # Ya en segundos
    # Quitamos el tiempo de la cola
    qTime = abs(timeSpent-queueTime)
    WF.Write_Bitfiles(directory,"AWS" + extraName, shotResults, overwrite_shots=True,shots = Nshots)
    return qTime

def Run_On_AWS(qc, backend, Nshots,li, ls, directory = "Outputs/AWS", timeDirectory = "Times", extraName = "", subQC = ""):
    """
    TODO: FIND BUG: Total bits stored doesnt quite correlate with total shots * qbits used.
    """
    
    if type(backend == str): print(f"Running on AWS backend {backend}... ")
    else: print(f"Running on AWS backend {backend.name}... ")
    T = Nshots//ls
    R = Nshots%ls
    
    if R < li and T == 0:
        print(f"Task was not done: We needed {R} shots, but minimun shots are {li}")        
        return None
    print(f"Full Tasks needed for this backend: {T}")
    print(f"Last task will have {R} shots")
    if subQC != "":
        directory += os.sep + subQC
        
    print("Executing task multiple times to fill the bitfile...")
    qTime = 0
    auxDir = "Temp" # Guardare en un directorio temporal
    for i in range(T):
        print(40*"-")
        print(f"Current Task: {i+1}")
        qTime_Task = Run_Task(qc, backend, ls,directory = auxDir,extraName = extraName)
        qTime += qTime_Task
        
    print(40*"-")
    if R < li:
        print(f"Shots removed: {R}")
        Nshots -= R
    else:
        # Una ultima ejecucion para el resto
        print("Last Execution for bits remaining...")
        qTime_Task = Run_Task(qc, backend, R,directory = auxDir,extraName = extraName)
        qTime += qTime_Task
    print("----------------")

    print("Merging in a single file...")
    savedName = WF.Combine_Bitfiles(auxDir, directory)
    # Eliminamos temp
    WF.Delete_Directory(auxDir)
    print("Saving Time Spent...")
    WF.Write_Analysis_Files(subQC,f"{savedName}-{Nshots}-(s):{qTime}", directory = timeDirectory)        
    print("Bitfile Created!")
    print("----------------")
    return None

#%% AQUILA (ES UN CASO A PARTE)

# TODO: PROGRMAR EL CODIGO EN CONCRETO PARA AQUILA (Dinamic Hamiltonian Simulator)
def Aquila_Bitstring(AquilaResult):
    meas = AquilaResult.measurements # Ojo, es una lista de arrays
    nShots = len(meas)
    nQbits = len(meas[0].post_sequence) # Suponemos que cuando se pierden e se sigue pillando qbits
    bitsArr = np.empty((nQbits*nShots), dtype = str)
    for i in range(nShots):
        bitsArr[i*nQbits:(i+1)*nQbits] = abs(meas[i].post_sequence.flatten()-1) # Invertimos el coteo
    return bitsArr

def Run_Task_Aquila(nQbits,Nshots, directory, extraName, simulate = False):
    # https://docs.aws.amazon.com/braket/latest/developerguide/braket-get-started-hello-ahs.html
    
    task = Computers.CallAquila(nQbits, Nshots, simulate = simulate)
    t1 = time.time()
    #print(f"Job created and returned as a variable. Job ID: {task.job_id()}")
    print("Program will wait until job ends...")
    result = task.result()   # Es un formato (nShots, nQbits)
    t2 = time.time()
    """
    Si tenemos un 1--> atomo detectado (ground)
    si tenemos un 0 --> atomo no detectado (o rydberg o se ha escapado)
    Nosotros invertiremos el valor de los bits para que cuadre con el resto de 
    computadores en la funcion Aquila_Bitstring
    """
    print("Job ended!")
    print("Creating bitArray...")
    bitstring = Aquila_Bitstring(result)
    print(collections.Counter(bitstring))
    # Calculamos el tiempo estimado de la cola en el hardware
    queueTime = t2-t1
    if simulate:
        timeSpent = 0
    else:
        timeSpent = AWS_GetTimeSpent(task.id)  # Ya en segundos
    # Quitamos el tiempo de la cola
    qTime = abs(timeSpent-queueTime)
    
    WF.Write_Bitfiles(directory,"Aquila" + extraName, bitstring, overwrite_shots = True, shots = Nshots) # bitstring no tiene los shots: lo metemos a mano
    return qTime

def Run_On_Aquila(nQbits,Nshots,li, ls, directory = "Outputs/AWS/Aquila", timeDirectory = "Times", extraName = "", simulate = False):
    print("Running on Aquila Backend from QuEra...")
    #print(Nshots)
    T = Nshots//ls # N de taks completos
    R = Nshots%ls # n de shots en el ultimo task

    if R < li and T == 0:
        print(f"Task was not done: We needed {R} shots, but minimun shots are {li}")        
        return None
    
    print(f"Full Tasks needed for this backend: {T}")
    print(f"Last task will have {R} shots")
    print("Executing task multiple times to fill the bitfile...")
    qTime = 0
    auxDir = "Temp" # Guardare en un directorio temporal
    for i in range(T):
        print(40*"-")
        print(f"Current Task: {i+1}")
        qTime_Task = Run_Task_Aquila(nQbits,ls,directory = auxDir,extraName = extraName, simulate = simulate)
        qTime += qTime_Task
        
    print(40*"-")
    if R < li:
        print(f"Shots removed: {R}")
        Nshots -= R
    else:
        # Una ultima ejecucion para el resto
        print("Last Execution for bits remaining...")
        print(R)
        qTime_Task = Run_Task_Aquila(nQbits,R,directory = auxDir,extraName = extraName, simulate = simulate)
        qTime += qTime_Task
    print("----------------")

    print("Merging in a single file...")
    savedName = WF.Combine_Bitfiles(auxDir, directory)
    # Eliminamos temp
    WF.Delete_Directory(auxDir)
    print("Saving Time Spent...")
    WF.Write_Analysis_Files("QuEra",f"{savedName}-{Nshots}-(s):{qTime}", directory = timeDirectory)        
    print("Bitfile Created!")
    print("----------------")
    return None
