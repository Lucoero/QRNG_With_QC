# -*- coding: utf-8 -*-
"""
Conectar con los distintos ordenadores cuanticos
"""
#%% Import libraries
import qiskit as qis
import numpy as np

import qiskit_braket_provider as qisbra # AWS
import braket.aws as aws # Para recuperar trabajos de aws con aws.AwsQuantumTask(arn)
import qiskit_ionq as qision #IonQ
import qiskit_ibm_runtime as qisibm #IBM
from qiskit_ibm_runtime import SamplerV2 as IBM_Sampler
import os as os
import dotenv as dotenv
import tqdm as tqdm
import matplotlib.pyplot as plt


# para aquila
from braket.ahs.atom_arrangement import AtomArrangement
from braket.timings.time_series import TimeSeries
from braket.ahs.driving_field import DrivingField
from braket.timings.time_series import TimeSeries
from braket.ahs.analog_hamiltonian_simulation import AnalogHamiltonianSimulation
from braket.devices import LocalSimulator # Debugging
from collections import Counter # Visualizar resultados

# To remove transpiler warning for IonQ (we already know qiskit uses lvl 2 instead of 1)
import warnings
from qiskit_ionq.exceptions import IonQTranspileLevelWarning
warnings.filterwarnings("ignore", category=IonQTranspileLevelWarning)

"""
Computadores escogidos

IBM: ibm-kingston (Usa Heron R2)
AWS:
    QuEra: Aquila (basado en atomos neutros atrapados con láseres, 256 qbits). us-east1. $0.30/tarea + $0.01/intento (bajo demanda)
        https://us-west-2.console.aws.amazon.com/braket/home?region=us-west-2#/devices/arn:aws:braket:us-east-1::device/qpu/quera/Aquila
    IonQ: Forte Enterprise 1 (iones atrapados)
        https://us-west-2.console.aws.amazon.com/braket/home?region=us-west-2#/devices/arn:aws:braket:us-east-1::device/qpu/ionq/Forte-Enterprise-1
    Rigetti:  Cepheus™-1-108Q $0.30/tarea + $0.000425/intento (bajo demanda). 
     <<las puertas CZ son más resistentes a los errores de fase, comunes en los sistemas superconductores>>
        https://us-west-2.console.aws.amazon.com/braket/home?region=us-west-2#/devices/arn:aws:braket:us-west-1::device/qpu/rigetti/Cepheus-1-108Q
    AWS Per se: 
    
Para settear
https://aws.amazon.com/es/blogs/quantum-computing/setting-up-your-local-development-environment-in-amazon-braket/
Paso 2:
https://aws.amazon.com/es/blogs/quantum-computing/introducing-the-qiskit-provider-for-amazon-braket/       
Braket:
    https://us-west-2.console.aws.amazon.com/braket/home?region=us-west-2#/dashboard
"""
#%% Load Providers
print("Loading and Connecting to dependencies of Computers...")

dotenv.load_dotenv("APIS.env")


#%%% This section is needed when this code is executed for the first time in your pc.
# Put your credentials in an env with those keywords in the root directory of the code
"""
qisibm.QiskitRuntimeService.save_account(
    token=os.getenv("IBM"),
    instance=os.getenv("IBM_Instance"),
    set_as_default = True,
    overwrite = True
    )    
"""
#%%% 
IBMProvider = qisibm.QiskitRuntimeService()
AWSProvider = qisbra.BraketProvider()

print("Dependencies Loaded!")
print("-"*40)
#%% Functions
def Status():
    """
    Check availability of all Computers, and return its backends
    """
    print("IBM COMPUTERS:\n")
    IBMBackends = IBMProvider.backends(simulator = False)
    for i in range(len(IBMBackends)):
        print("Name:", IBMBackends[i].name)
        if IBMBackends[i] == None:
            print("Cant access status")
        else:
            print('Status:')
            print('  Operational: ', IBMBackends[i].status().operational)
            print('  Pending jobs:', IBMBackends[i].status().pending_jobs)
            print('  Status message:',IBMBackends[i].status().status_msg)
        print("-" * 40)
    print("For a full description of each QC, check: \n https://quantum.cloud.ibm.com/computers \n")
    print("*"*40)
    
    print("AWS COMPUTERS:\n")
    AWSBackends = AWSProvider.backends()
    
    for i in range(len(AWSBackends)):
        print("Name:", AWSBackends[i].name)
        print("-" * 40)
    print("For a full description of each QC, check: \n https://quantum.cloud.ibm.com/computers \n")
    print("*"*40)
    
    print("Returned tuple of lists of each QC available in print order.")
    return IBMBackends,AWSBackends

def CallIBM(qc,backend,shots, lvl = 1,drawTranspile = False):
    """
    Call an IBM QC for a query. 
    https://quantum.cloud.ibm.com/docs/es/guides/hello-world
    """
    # Get Backend
    if type(backend) == str:
        backend = IBMProvider.backend(backend)
        
    print(f"Preparing IBM job at backend {backend.name}")
    # First  Optimize circuit for the backend choosed
    qcT = qis.compiler.transpile(qc,backend = backend, optimization_level= lvl) # Its exactly the same as generate_preset_pass_manager
    #https://quantum.cloud.ibm.com/docs/es/guides/visualize-circuits For saving circuits in pdf
    if drawTranspile:
        qcT.draw("mpl",idle_wires=False)
    
    # Then create and send job
    # Retrieve job details
    sampler = IBM_Sampler(mode = backend)
    sampler.options.default_shots = shots
    job = sampler.run([qcT])
    
    print(f"Job created and returned as a variable. Job ID: {job.job_id()}")
    print("Jobs can be checked at \n https://eu-de.quantum.cloud.ibm.com/workloads")
    return job


#%% AWS
def CallAWS(qc,backend,shots, lvl = 1, drawTranspile = False):
    """
    Call an AWS QC for a query
    region options:
       us-west-1 (rigetti)
       us-east-1 (QuEra, IonQ)
    """
    # Get Backend
    if type(backend) == str:
        backend = AWSProvider.get_backend(backend)
    print(f"Preparing AWS job at backend {backend.name}")
    # First  Optimize circuit for the backend choosed
    qcT = qis.compiler.transpile(qc,backend = backend, optimization_level=lvl) # Its exactly the same as generate_preset_pass_manager
    qcT.data = [inst for inst in qcT.data
    if inst[0].name != "barrier"] # Elimino las barreras (no compatibles con IonQ ni QuEra)
    #https://quantum.cloud.ibm.com/docs/es/guides/visualize-circuits For saving circuits in pdf
    if drawTranspile:
        qcT.draw("mpl",idle_wires=False)
    
    # Then create and send task
    # Retrieve task details
    task = backend.run(qcT,shots = shots)#,memory = True)
    
    print(f"Job created and returned as a variable. Job ID: {task.job_id()}")
    print("Jobs can be checked at \n https://us-west-1.console.aws.amazon.com/braket/home?region=us-west-1#/tasks  \n ")
    return task

#%% Aquila
def Visualize_Driving_Field(drive):
    # https://docs.aws.amazon.com/braket/latest/developerguide/braket-get-started-hello-ahs.html
    fig, axes = plt.subplots(3, 1, figsize=(12, 7), sharex=True)
    
    ax = axes[0]
    time_series = drive.amplitude.time_series
    ax.plot(time_series.times(), time_series.values(), '.-')
    ax.grid()
    ax.set_ylabel('Omega [rad/s]')
    
    ax = axes[1]
    time_series = drive.detuning.time_series
    ax.plot(time_series.times(), time_series.values(), '.-')
    ax.grid()
    ax.set_ylabel('Delta [rad/s]')
    
    ax = axes[2]
    time_series = drive.phase.time_series
    # Note: time series of phase is understood as a piecewise constant function
    ax.step(time_series.times(), time_series.values(), '.-', where='post')
    ax.set_ylabel('phi [rad]')
    ax.grid()
    ax.set_xlabel('time [s]')
    
    plt.show()
    return

def CallAquila(nQbits,shots, showDist = False, showField = False, simulate = False):
    # https://docs.aws.amazon.com/braket/latest/developerguide/braket-get-started-hello-ahs.html
    # TODO: Adaptar esto porque gran parte deberia estar en MyCircuits
    # Preparamos el grid de atomos
    grid = AtomArrangement()
    """
    Por los calculos del TFG, sabemos que podemos sacar la distancia entre atomos como 
    d = W/sqrt(nQbits) (um, si W en um)
    Donde W = 75 um la anchura del grid maximo (Altura simeq anchura)
    El numero de qbits debe ser par
    """
    W = 75 * 10**(-6)  # Anchura maxima grid (m)     H = 75 um # Altura maxima grid (realmente es 76)
    d = W/np.sqrt(nQbits) # en m
    nFilas = int(np.ceil(W/d)) # n qbits por fila
    nCol = nFilas # n qbits por columna
    # Vamos colocando uno a uno
    aux = 0 # Contador de cuantos qbits llevamos
    for i in range(nFilas):
        for j in range(nCol): # Voy rellenando columna a columna
            if aux < nQbits:
                grid.add([d*i, d*j])
                aux += 1
    # Si queremos vemos la distribucion
    if showDist:
        fig,ax = plt.subplots()
        xs, ys = [grid.coordinate_list(dim) for dim in (0, 1)]
        # print(len(xs))
        ax.plot(xs, ys, 'r.', ms=15)
        for idx, (x, y) in enumerate(zip(xs, ys)):
            ax.text(x, y, f" {idx}", fontsize=12)
        fig.show()

    # Grid Creado, inicializamos las variables
    omega = 1.58 * 10**6 # Float cte
    #omega = omega*10  # Frecuencia maxima, por si la queremos usar
    
    Omega = TimeSeries() # Time Series
    
    
    Delta = TimeSeries()
    Delta.put(0.,0.)
    
    Phi = TimeSeries()
    Phi.put(0.,0.)
    
    endTime = np.pi/(omega) # Segundos. Se obtiene de buscar t tal que Omega t = pi/2
    """
    dt = 5.0 * 10**(-8) # Minimo que acepta aquila
    times = np.arange(0,endTime + dt,dt) # Por alguna razon dice que la duracion entre pulsos no es igual
    """
    times = np.linspace(0,endTime,30) # Parecen pocos puntos, pero el intervalo es muy pequenno
    y = omega * np.sin(np.pi * times/endTime)**2
    y[-1] = 0.
    for i in range(len(times)):
        Omega.put(times[i],y[i])

    Delta.put(endTime,0.)
    Phi.put(endTime,0.)
    
    drive = DrivingField(amplitude = Omega, phase = Phi, detuning = Delta) 
    if showField: 
        Visualize_Driving_Field(drive)
    
    # Cerramos el AHS program
    program = AnalogHamiltonianSimulation(register=grid, hamiltonian=drive)
    
    # Creamos el task y lo enviamos
    if simulate:
        aquila = LocalSimulator("braket_ahs") # Local
        task = aquila.run(program,shots = shots)
    else:
        aquila = aws.AwsDevice("arn:aws:braket:us-east-1::device/qpu/quera/Aquila")
        program_discretized = program.discretize(aquila) # Redondea para que aquila pueda procesarlo
        task = aquila.run(program_discretized,shots = shots)
        print(f"Job created and returned as a variable. Job ID: {task.id}")
        print("Jobs can be checked at \n https://us-west-1.console.aws.amazon.com/braket/home?region=us-west-1#/tasks  \n ")
    return task


    
