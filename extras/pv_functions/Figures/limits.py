#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special

plt.figure()



A=np.genfromtxt('../Figures/data/limits.txt',usecols=[0,1,2,3,4,5])
x=A[:,0]
y=A[:,1]
y2=A[:,2]
y3=A[:,3]
y4=A[:,4]
y5=A[:,5]


plt.plot(x,y4,'-.',color='red',label='$r=1.001$',linewidth=1.5) #
plt.plot(x,y5,'--',color='red',label='$r=1.0001$',linewidth=1.5) #
plt.plot(x,y,'-',color='black',label='$r=1$',linewidth=1.5) #
plt.plot(x,y2,'--',color='blue',label='$r=0.9999$',linewidth=1.5) #
plt.plot(x,y3,'-.',color='blue',label='$r=0.99$',linewidth=1.5) #


xlabel(r"$M$ ",fontsize=16)
ylabel(r"$\frac{M}{\pi}[B_0(rM,M,m_1)-B_0(rM,M,m_2)]$ ",fontsize=16)

plt.tick_params(labelsize=14)

plt.xscale('log')
plt.legend(loc=2)
plt.ylim([-10,120])
plt.xlim([min(x),max(x)])

plt.savefig("../Figures/Figures/limits.pdf")
