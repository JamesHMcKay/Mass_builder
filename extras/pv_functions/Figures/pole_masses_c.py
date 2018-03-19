#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special

plt.figure()



A=np.genfromtxt('../Figures/data/pole_masses_c.txt',usecols=[0,1,2,3])
x=A[:,0]
y=A[:,1]
y_simple=A[:,3]

plt.plot(x,x,'-',color='black',label='tree level mass') #
plt.plot(x,y,'-',color='green',label='pole mass') #
plt.plot(x,y_simple,'-',color='red',label='pole mass simple') #

xlabel(r"Degenerate mass (GeV)",fontsize=14)
ylabel(r"$M_{pole}$ (Gev)",fontsize=14)
#plt.xlim([min(x),1e5])
#plt.ylim([min(sigma_v(s_plot)),max(sigma_v(s_plot))])

plt.xscale('log')
plt.yscale('log')
plt.legend(loc=2)


plt.savefig("../Figures/Figures/pole_masses_c.eps")