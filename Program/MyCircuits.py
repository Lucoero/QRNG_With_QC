# -*- coding: utf-8 -*-
"""
MyCircuits
"""
import qiskit as qis
import matplotlib.pyplot as plt

def Simple_Qbit(printCircuit = False):
    """
    Simple_Qbit:
        Create a 1 Qbit circuit with 50/50 prob of 1 or 0
    """
    
    qc = qis.QuantumCircuit(1) # Creo un circuito con un qbit (inicializado por defecto en 0)
    qc.h(0)
    qc.measure_all()

    if printCircuit:
        qc.draw(output = "mpl") #mpl = matplotlib
    return qc

def Hadamard_Transform(n, printCircuit = False):
    """
    Hadamard_Transform:
        Create a circuit of n qbits with a Hadamard Transform 
        (n Hadamard gates in parallel)
    """
    qc = qis.QuantumCircuit(n)
    for i in range(n):
        qc.h(i)
    qc.measure_all()
    if printCircuit:
        qc.draw(output = "mpl")
    return qc