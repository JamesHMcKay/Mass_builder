#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
import pylab
from scipy import interpolate
from scipy import integrate
from scipy import special

from matplotlib import ticker

from matplotlib import rc
rc('text', usetex=True)
plt.rc('font', family='Computer Modern Roman',weight='normal')

fig=plt.figure()

ax = fig.add_subplot(1,1,1)
ax.set_xscale('log')
ax.set_yscale('log')

#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
ax.set_yticks([1e2, 1e4, 1e6,1e8,1e10,1e12])


#ax.yaxis. #set_major_formatter(ticker.FormatStrFormatter("%d"))

A=np.genfromtxt('mass_splittings/data/decays_iterative.txt',usecols=[0,1,2])

x=A[:,0]
ylower=A[:,1]*2.9979e11
yupper=A[:,2]*2.9979e11

plt.fill_between(x, ylower, yupper ,color='grey',label="iterative")

A=np.genfromtxt('mass_splittings/data/decays_explicit.txt',usecols=[0,1,2])

x=A[:,0]
ylower=A[:,1]*2.9979e11
yupper=A[:,2]*2.9979e11

plt.fill_between(x, ylower, yupper ,color='lightgrey',label="non-iterative")


xlabel(r"Degenerate mass $\hat{M}$ $($GeV$)$",fontsize=16)
ylabel(r"$\tau$ $($mm/$c$$)$",fontsize=16)

#plt.xscale('log')
#plt.yscale('log')
#plt.ylim([10,1e12])
plt.xlim([min(x),max(x)])



# get the value of Xi
Mass=np.genfromtxt('models/MSSM/input.txt',usecols=[1])
Xi = Mass[12]

plt.legend(loc='upper left',fontsize=16)
plt.savefig("mass_splittings/figures/decays.pdf")




