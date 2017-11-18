#!/usr/bin/python

# Mass Builder
 
# James McKay
# Sep 2016
 
#--- plot_MSSM.py ---
 
# Simple routine to plot output from mass splittings example


import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special

from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
from matplotlib import rc

rc('text', usetex=True)
plt.rc('font', family='Computer Modern Roman',weight='normal')

fig=plt.figure()

ax = fig.add_subplot(1,1,1)


ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

A=np.genfromtxt('mass_splittings/data2/deltam_M_test.txt',usecols=[0,1,2,3,4])
x=A[:,0]
y1=A[:,1]*1000
y2=A[:,2]*1000
y3=A[:,3]*1000
y4=A[:,4]*1000


plt.plot(x,y1,'-',color='red',label='1-loop non-iterative') #

plt.plot(x,y2,'-',color='green',label='1-loop iterative') #

plt.plot(x,y3,'--',color='red') #

plt.plot(x,y4,'--',color='green') #

plt.xlim([10,1e5])
#plt.ylim([0,200])

leg = plt.legend(loc='lower left')
ax.set_xscale('log')

xlabel(r"Degenerate mass $M$ (GeV)",fontsize=18)
ylabel(r"$\Delta M$ (Mev)",fontsize=18)

plt.savefig("mass_splittings/figures2/plot_test.eps")


