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

# This routine plots Figure 6 in arXiv:1710.01511

fig=plt.figure()

ax = fig.add_subplot(1,1,1)
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_yticks([1e2, 1e4, 1e6,1e8,1e10,1e12])

A=np.genfromtxt('mass_splittings/data/decays_iterative_1.txt',usecols=[0,1,2])

x=A[:,0]
ylower=A[:,1]*2.9979e11
yupper=A[:,2]*2.9979e11

plt.fill_between(x, ylower, yupper ,color='#7fbf7f',label="iterative")


A=np.genfromtxt('mass_splittings/data/decays_explicit_1.txt',usecols=[0,1,2])

x=A[:,0]
ylower=A[:,1]*2.9979e11
yupper=A[:,2]*2.9979e11

plt.fill_between(x, ylower, yupper ,color='grey',label="non-iterative")


xlabel(r"Degenerate mass $\hat{M}$ $($GeV$)$",fontsize=18)
ylabel(r"$\tau$ $($mm/$c$$)$",fontsize=18)

plt.xlim([min(x),max(x)])

plt.tick_params(labelsize=14)

plt.legend(loc='upper left',fontsize=16)
plt.savefig("mass_splittings/figures/decays.pdf")
