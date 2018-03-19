#!/usr/bin/python

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

A=np.genfromtxt('../Figures/data/mass_splittings.txt',usecols=[0,1,2,3,4])
x=A[:,0]
y=A[:,1]*1000
y2=A[:,2]*1000

plt.plot(x,y,'-',color='black',label='Explicit pole mass') #

plt.plot(x,y2,'--',color='black',label='Implicit pole mass') #

#ax.set_yscale('symlog')
ax.set_xscale('log')

xlabel(r"Degenerate mass $M$ (GeV)",fontsize=16)
ylabel(r"$\Delta M$ (Mev)",fontsize=16)
#plt.ylim([0,710])
#plt.xlim([1,3e4])
#plt.ylim([min(sigma_v(s_plot)),max(sigma_v(s_plot))])


plt.legend(loc=2)

#plt.annotate('$M_{\chi^{++}}-M_{\chi^0}$', xy=(300, 600), xytext=(300, 600),fontsize=16)

#plt.annotate('$M_{\chi^+}-M_{\chi^0}$', xy=(70, 200), xytext=(70, 200),fontsize=16)

#plt.plot(x,y2,'-',color='black',label='Explicit pole mass') #
#plt.plot(x,y,'--',color='black',label='Implicit pole mass') #



#plt.axvspan((9.4e3-0.47e3), (9.4e3+0.47e3), alpha=0.8, color='green')

#plt.plot((9.4e3,9.4e3),(0,700),color='red')

plt.savefig("../Figures/Figures/mass_splittings.eps")

