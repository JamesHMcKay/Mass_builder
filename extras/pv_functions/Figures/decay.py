#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special

from matplotlib import ticker

fig=plt.figure()

ax = fig.add_subplot(1,1,1)
ax.set_xscale('log')
ax.set_yscale('log')

ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
ax.set_yticks([1e2, 1e4, 1e6,1e8,1e10,1e12])


#ax.yaxis. #set_major_formatter(ticker.FormatStrFormatter("%d"))


A=np.genfromtxt('../Figures/data/decay.txt',usecols=[0,1,2])
x=A[:,0]
y=A[:,1]*2.9979e11
y2=A[:,2]*2.9979e11

B=np.genfromtxt('../Figures/data/decay2.txt',usecols=[0,1,2])
z=B[:,1]*2.9979e11
z2=B[:,2]*2.9979e11


plt.plot(x,y2,'-',color='black',label='Explicit pole mass') #
plt.plot(x,y,'--',color='black',label='Implicit pole mass') #



plt.axvspan((9.4e3-0.47e3), (9.4e3+0.47e3), alpha=0.8, color='green')

plt.plot((9.4e3,9.4e3),(0,1e12),color='red')










#plt.plot(x,z2,'-',color='black') #
#plt.plot(x,z,'--',color='black') #
#plt.plot(x,theory,'--',color='blue',label='theory') #

xlabel(r"Degenerate mass $M$ (GeV)",fontsize=16)
ylabel(r"$\tau$ (mm/$c$)",fontsize=16)

#plt.ylim([min(sigma_v(s_plot)),max(sigma_v(s_plot))])

#plt.xscale('log')
#plt.yscale('log')
plt.legend(loc=2)
plt.ylim([10,1e12])
plt.xlim([1,2e4])
#plt.annotate('$M_{\chi^{++}}$', xy=(10,8 ), xytext=(10,8),fontsize=16)

#plt.annotate('$M_{\chi^+}$', xy=(40, 0.5e4), xytext=(40, 0.5e4),fontsize=16)




plt.savefig("../Figures/Figures/decay.eps")