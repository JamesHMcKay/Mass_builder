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


fig=plt.figure()

ax = fig.add_subplot(1,1,1)


ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

A=np.genfromtxt('examples/QED_test.txt',usecols=[0,1,2])
x=A[:,0]
y1=A[:,1]
y2=A[:,2]


print x

plt.plot(x,y1,'-',color='green',label='off-shell') #
plt.plot(x,y2,'--',color='red',label='on-shell') #


#ax.set_yscale('log')
#ax.set_xscale('log')

#plt.plot(x,c0*1000.*ones(size(y1)),'--',color='black',label='Limit') #

xlabel(r"momentum",fontsize=16)
ylabel(r"amplitude",fontsize=16)

plt.savefig("QED_test.eps")
