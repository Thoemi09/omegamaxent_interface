#! /usr/bin/env python

import OmegaMaxEnt_TRIQS as OT
from pytriqs.archive import HDFArchive as HA
from matplotlib import pyplot as plt
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

A=HA("G.h5",'r')
G=A['G']

GR=OT.compute_GfReFreq(G,interactive_mode=False)

file=open("spectr_0_0.dat",'r')
A00_data=np.loadtxt(file)
file.close()

file=open("spectr_0_1.dat",'r')
A01_data=np.loadtxt(file)
file.close()

file=open("spectr_1_1.dat",'r')
A11_data=np.loadtxt(file)
file.close()

w=A00_data[:,0]
A00=A00_data[:,1]
A01=A01_data[:,1]
A11=A11_data[:,1]

with HA("G_Re_freq.h5", 'r') as A:
    G=A['G']

G00=G['0','0']
G01=G['0','1']
G11=G['1','1']

wG=np.array([om.value for om in G00.mesh])

G00_i=G00.data.imag
G01_i=G01.data.imag
G11_i=G11.data.imag

lgd=["$A_{00}$","-2Im$[G_{00}]$"]
plt.figure(1)
plt.plot(w,A00,'b-')
plt.plot(wG,-2*G00_i,'r-')
plt.legend(lgd)

lgd=["$A_{01}$","-2Im$[G_{01}]$"]
plt.figure(2)
plt.plot(w,A01,'b-')
plt.plot(wG,-2*G01_i,'r-')
plt.legend(lgd)

lgd=["$A_{11}$","-2Im$[G_{11}]$"]
plt.figure(3)
plt.plot(w,A11,'b-')
plt.plot(wG,-2*G11_i,'r-')
plt.legend(lgd)

plt.show()