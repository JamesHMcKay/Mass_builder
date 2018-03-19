#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special

plt.figure()



A=np.genfromtxt('../Figures/data/derivatives.txt',usecols=[0,1,2])
x=A[:,0]
y=A[:,1]
y2=A[:,2]

#y4=A[:,4]
#y5=A[:,5]



plt.plot(x,y,'-',color='black',label='numerical') #
plt.plot(x,y2,'--',color='black',label='integral') #
#plt.plot(x,y4,'-',color='black',label='$r=0.9$') #
#plt.plot(x,y5,'--',color='black',label='$r=0.6$') #


xlabel(r"$M$ ",fontsize=14)
ylabel(r"$\Sigma(p)$ ",fontsize=14)


plt.xscale('log')

plt.xscale('log')
plt.legend(loc=2)
#plt.xlim([min(x),1e8])
#plt.ylim([0,105])

plt.savefig("../Figures/Figures/derivatives.eps")