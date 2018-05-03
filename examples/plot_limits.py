#!/usr/bin/python

import numpy as np
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special


from matplotlib import rc
rc('text', usetex=True)
plt.rc('font', family='Computer Modern Roman',weight='normal')


plt.figure()



A=np.genfromtxt('examples/limits.txt',usecols=[0,1,2,3,4,5])
x=A[:,0]
y=A[:,1]*1000
y2=A[:,2]*1000
y3=A[:,3]*1000
y4=A[:,4]*1000
y5=A[:,5]*1000


plt.plot(x,y4,'-.',color='red',label='$r=1.001$',linewidth=1.5) #
plt.plot(x,y5,'--',color='red',label='$r=1.0001$',linewidth=1.5) #
plt.plot(x,y,'-',color='black',label='$r=1$',linewidth=1.5) #
plt.plot(x,y2,'--',color='blue',label='$r=0.9999$',linewidth=1.5) #
plt.plot(x,y3,'-.',color='blue',label='$r=0.99$',linewidth=1.5) #


xlabel(r"$\hat{M}$ $(\mathrm{GeV})$ ",fontsize=20)
#ylabel(r"$\frac{M}{\pi}[\mathbf{B}(rM,M,m_1)-\mathbf{B}(rM,M,m_2)]$ ",fontsize=16)
#ylabel(r"$\frac{M}{\pi}[s_W^2B_0(rM,M,0)+c_W^2B_0(rM,M,m_Z)-B_0(rM,M,m_W)]$ $(\mathrm{MeV})$ ",fontsize=14)
ylabel(r"$\Delta B$ $(\mathrm{MeV})$ ",fontsize=20)

plt.tick_params(labelsize=20)

plt.xscale('log')
plt.legend(loc=2,fancybox=True, framealpha=0.5)
#plt.ylim([-10,120])
plt.ylim([-10,220])
plt.xlim([min(x),max(x)])

plt.savefig("examples/limits.pdf",bbox_inches='tight')
