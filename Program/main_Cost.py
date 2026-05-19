# -*- coding: utf-8 -*-
"""
Main_Cost

Script para analizar los costes de las operaciones.
"""
import Analyse_Data as AD
import numpy as np
import matplotlib.pyplot as plt
#%% Variables
Cepheus = [0.000425,0.3,50000,10,107] # coste por shot, coste por task, lim sup de shot, lim inf de shot, qbits
Forte = [0.08,0.3, 5000 ,100,36]
Aquila = [0.01,0.3,1000,1,256] 
AQT = [0.0235,0.3,2000,1,12]

vs = 5.309 * 10**7 # Ratio IBM
Ibm = [vs,96/60] # Ibm solo cobra por tiempo. coste dolares/s

minBit, maxBit = 10**4,10**10
#%% Obtengamos las funciones
nBits = np.linspace(minBit,maxBit, 10**6)

FCep = AD.AWS_Cost_Function(nBits, *Cepheus)
FForte = AD.AWS_Cost_Function(nBits,*Forte)
FAquila = AD.AWS_Cost_Function(nBits,*Aquila)
FAQT = AD.AWS_Cost_Function(nBits,*AQT)

FIbm = AD.IBM_Cost_Function(nBits, *Ibm)

# Presupuesto
maxCost = 66
#%% Ploteamos
plt.close("all")
fontsize = 13 #  text scale for plots
numbersize = [13,13] # scale for number axes
plt.rc('font', size=fontsize) # Para el escalado de texto en plots
#%%% Entre los de AWS
fig,ax = plt.subplots()
fig.suptitle("Cost Between AWS Computers")
ax.grid()
ax.set_xlabel("bits")
ax.set_ylabel(r"Cost ($ \$ $)")
#ax.semilogx(nBits,FCep, label = "Cepheus™-1-108Q")
#ax.semilogx(nBits, FForte,label = "Forte-1")
#ax.semilogx(nBits,FAquila, label = "Aquila")
ax.plot([minBit,maxBit],[maxCost,maxCost], linestyle = "dashed", label = f"Budget ({maxCost} $\$ $)", color = "black")
ax.semilogx(nBits, FAQT,label = "IBEX Q1")
ax.legend()
fig.show()
#%%% Para IBM
fig3,ax3 = plt.subplots()
fig3.suptitle("Approximate Cost of IBM")
ax3.grid()
ax3.set_xlabel("bits")
ax3.set_ylabel(r"Cost ($ \$ $)")
ax3.semilogx(nBits,FIbm)
fig3.show()
#%%% Para todos
fig2,ax2 = plt.subplots()
fig2.suptitle("Cost Between Computers")
ax2.set_ylim(0,np.max(FIbm))
ax2.grid()
ax2.set_xlabel("bits")
ax2.set_ylabel(r"Cost ($ \$ $)")
ax2.semilogx(nBits,FCep, label = "Cepheus™-1-108Q")
ax2.semilogx(nBits, FAQT,label = "Ibex")
ax2.semilogx(nBits,FAquila, label = "Aquila")
ax2.plot(nBits, FIbm, label = "ibm_kingston")
ax2.plot([minBit,maxBit],[maxCost,maxCost], linestyle = "dashed", label = f"AWS Budget ({maxCost} $\$ $)", color = "black")
ax2.legend()
fig2.show()