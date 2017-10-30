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

A=np.genfromtxt('mass_splittings/data2/deltam_M.txt',usecols=[0,1,2,3])
x=A[:,0]
y1=A[:,1]*1000
y2=A[:,2]*1000
y3=A[:,3]*1000


plt.plot(x,y1,'-',color='red',label='1-loop') #

plt.plot(x,y2,'-',color='green',label='2-loop') #

plt.plot(x,y3,'--',color='green',label='2-loop FS') #

plt.xlim([90,5000])
plt.ylim([150,180])
leg = plt.legend(loc='lower left')
ax.set_xscale('log')

xlabel(r"Degenerate mass $M$ (GeV)",fontsize=18)
ylabel(r"$\Delta M$ (Mev)",fontsize=18)

plt.savefig("mass_splittings/figures2/mass_splittings.eps")

# determine a polynomial fit through the 2-loop mass splitting

fig=plt.figure()

ax = fig.add_subplot(1,1,1)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

fx=np.polyfit(log(x),y2,4)
f=np.poly1d(fx)

print(np.poly1d(f))

r=np.linspace(log(90),log(4000))


plt.plot(exp(r),f(r),'-',color='black')

plt.plot(x,y2,'x',color='red',label="1-loop")
plt.plot(x,y2,'-',color='green',label='2-loop') #

gx = [-0.181509, 5.41949, -60.8831 , 305.383, -413.315]
g=np.poly1d(gx)
plt.plot(exp(r),g(r),'-',color='yellow')


plt.xlim([90,4000])
plt.ylim([145,172])
ax.set_xscale('log')

plt.savefig("mass_splittings/figures2/interpolation.eps")


