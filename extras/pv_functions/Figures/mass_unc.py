#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special

plt.figure()



A=np.genfromtxt('../Figures/data/mass_unc.txt',usecols=[0,1,2])
x=A[:,0]
y=A[:,1]
y2=A[:,2]

plt.plot(x,y,'--',color='black',label='Implicit pole mass') #
plt.plot(x,y2,'-',color='black',label='Explicit pole mass') #

xlabel(r"$Q$ (GeV)",fontsize=14)
ylabel(r"$M_{pole}$ (GeV)",fontsize=14)
#plt.xlim([min(x),10e3])
#plt.ylim([min(sigma_v(s_plot)),max(sigma_v(s_plot))])

#plt.yscale('log')
plt.legend(loc=2)


plt.savefig("../Figures/Figures/mass_unc.eps")
