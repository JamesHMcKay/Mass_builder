#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special

plt.figure()

A=np.genfromtxt('mass_splittings/deltam_Q.txt',usecols=[0,1])
x=A[:,0]
ylower=A[:,1]


plt.plot(x,ylower,'-')

#A=np.genfromtxt('../Figures/data/mass_splittings_explicit.txt',usecols=[0,1,2,3,4,5])
#x=A[:,0]
#ylower=A[:,1]
#ymid=A[:,2]
#yupper1=A[:,3]
#yupper2=A[:,4]
#yupper3=A[:,5]


xlabel(r"Q",fontsize=14)
ylabel(r"$\Delta m$ (Gev)",fontsize=14)

plt.xscale('log')
plt.ylim([0,0.5])
plt.xlim([0,13170.3])

plt.savefig("mass_splittings/test.eps")
