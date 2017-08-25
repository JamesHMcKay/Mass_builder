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

A=np.genfromtxt('models/MSSM/output/mass_splittings.txt',usecols=[0,1,2])
x=A[:,0]
y1=A[:,1]*1000
y2=A[:,2]*1000



plt.plot(x,y2,'--',color='red',label='1-loop') #

plt.plot(x,y1,'-',color='green',label='2-loop') #

plt.xlim([90,4000])
plt.ylim([145,172])

#plt.scatter(x,y4,color='green',label='1-loop iterative',s=2) #
#plt.scatter(x,y3,color='red',label='1-loop non-iterative',s=2) #

#plt.plot(x,y3,'-',color='green',label='1-loop non-iterative') #

#plt.plot(x,y4,'--',color='green',label='1-loop iterative') #

#ax.set_yscale('log')
ax.set_xscale('log')

#plt.plot(x,c0*1000.*ones(size(y1)),'--',color='black',label='Limit') #

xlabel(r"Degenerate mass $M$ (GeV)",fontsize=16)
ylabel(r"$\Delta M$ (Mev)",fontsize=16)



plt.savefig("mass_splittings_MSSM.eps")

# determine a polynomial fit through the 2-loop mass splitting

fig=plt.figure()

ax = fig.add_subplot(1,1,1)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

fx=np.polyfit(log(x),y1,4)
f=np.poly1d(fx)

print(np.poly1d(f))

r=np.linspace(90,4000)


plt.plot(r,f(log(r)),'-',color='black')

plt.plot(x,y1,'x',color='red')
plt.plot(x,y1,'-',color='green',label='2-loop') #


plt.xlim([90,4000])
plt.ylim([145,172])
ax.set_xscale('log')

plt.savefig("interpolation.eps")


