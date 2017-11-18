#!/usr/bin/python

import numpy as np
import pylab
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special
from scipy.interpolate import spline

from matplotlib import rc
rc('text', usetex=True)
plt.rc('font', family='Computer Modern Roman',weight='normal')

fig = plt.figure()

ax = fig.add_subplot(111)

# 2-loop curves
A=np.genfromtxt('mass_splittings/data2/decays_BR_2loop_MDM.txt',usecols=[0,3,4,5])

x=A[:,0]
BRpi=A[:,1]
BRmu=A[:,2]
BRe=A[:,3]

BRpi2 = BRpi[BRpi>0]
xpi = x[BRpi>0]

BRmu2 = BRmu[BRmu>0]
xmu = x[BRmu>0]

BRe2 = BRe[BRe>0]
xe = x[BRe>0]


ax.plot(xpi,BRpi2,'-',color='blue',label=r'$\chi^+\rightarrow\chi^0e^+\nu_e$',linewidth=1,zorder=500) #
ax.plot(xmu,BRmu2,'-',color='red',label=r'$\chi^+\rightarrow\chi^0\mu^+\nu_{\mu}$',linewidth=1,zorder=500) #
ax.plot(xe,BRe2,'-',color='black',label=r'$\chi^+\rightarrow\chi^0\pi^+$',linewidth=1,zorder=500) #

# 1-loop curves

A=np.genfromtxt('mass_splittings/data2/decays_BR_MDM.txt',usecols=[0,3,4,5])

x=A[:,0]
BRpi=A[:,1]
BRmu=A[:,2]
BRe=A[:,3]

BRpi2 = BRpi[BRpi>0]
xpi = x[BRpi>0]

BRmu2 = BRmu[BRmu>0]
xmu = x[BRmu>0]

BRe2 = BRe[BRe>0]
xe = x[BRe>0]


ax.plot(xpi,BRpi2,'--',color='blue',linewidth=1,zorder=500) #
ax.plot(xmu,BRmu2,'--',color='red',linewidth=1,zorder=500) #
ax.plot(xe,BRe2,'--',color='black',linewidth=1,zorder=500) #




ax.set_xlabel(r"$\mathrm{Degenerate\ mass}$ $\hat{M}$ $(\mathrm{GeV})$",fontsize=18)
ax.set_ylabel(r"$\mathrm{Branching\ ratio}$",fontsize=18)

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_yticks([1e-3,1e-2,1e-1,1])
ax.set_yticklabels(['$\mathrm{10}^{\mathrm{-3}}$','$\mathrm{0.01}$','$\mathrm{0.1}$','$\mathrm{1}$'],fontsize=18)



ax.set_xticks([1,100,1e3,1e4])
ax.set_xticklabels(['$\mathrm{1}$','$\mathrm{10}$','$\mathrm{100}$','$\mathrm{10}^{\mathrm{3}}$','$\mathrm{10}^{\mathrm{4}}$'],fontsize=18)

ax.set_xlim([1,1e4])
ax.set_ylim([1e-3,1.2])

leg = ax.legend(loc='lower left',frameon=False)

plt.annotate(r"$\mathrm{MDM}$", xy=(2,0.3), xytext=(2,0.3),fontsize=16)
plt.annotate(r"$Q = 173.34\ (\mathrm{GeV})$", xy=(2,0.2), xytext=(2,0.2),fontsize=16)



plt.savefig("mass_splittings/figures2/BR.pdf",bbox_inches='tight')



