# -*- coding: utf-8 -*-
"""
Prueba de simulacion
"""
import numpy as np
import matplotlib.pyplot as plt
import qiskit as qis
import WriteFile as WF
from qiskit.primitives import StatevectorSampler
from qiskit.visualization import plot_histogram

plt.close("all")
#%% Variables
N = 10**6 # Ojo, luego se va a hacer un array de este size, cuidado. Con 10**7 ya tarda bastante
if N > 10**7:
    print("Reducimos N a 10**7 para evitar overflow")
    N = 10**7
path = r"Outputs"
secType = "Simulado"
#%% Creacion del circuito
qc = qis.QuantumCircuit(1) # Creo un circuito con un qbit (inicializado por defecto en 0)
qc.h(0)
qc.measure_all()

qc.draw(output = "mpl") #mpl = matplotlib
#%% Lectura del circuito (aqui es donde entra la naturaleza de simulacion)
sampler = StatevectorSampler()
job = sampler.run([qc], shots=N)
result = job.result()

shotResults = np.array(result[0].data.meas.get_bitstrings(),dtype = int) # Aqui tenemos el obtenido en cada shot
WF.WriteResults(directory = path,sec_origin = secType, results = shotResults)

counts = result[0].data.meas.get_counts()# Aqui un conteo
print(counts)
#%% Histogramas
plot_histogram(counts)

