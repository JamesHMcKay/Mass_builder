#!/usr/bin/python

# Mass Builder
 
# James McKay
# Sep 2016
 
#--- plot_VDM.py ---
 
# plot mass splitting for VDM model


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
ax.set_xscale('log')

ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

A=np.genfromtxt('models/VDM/output/mass_splittings.txt',usecols=[0,1,2,3])
x=A[:,0]
y1=A[:,1]*1000
y2=A[:,2]*1000
y3=A[:,3]*1000

i = complex(0.0, 1.0)
e = 0.3134

pi = 4.0*arctan(1.0)

mw = 80.385
mz = 91.1876

CW = mw/mz
CW2 = power(CW,2)
sw2 = 1.0-CW2
sw = power(sw2,0.5)

print "sw2 = ", sw2


delta = -0.577216 + log(4*pi);

F = (mw+mz*(sw2-1.0))


c0 = real((1./2.)*(5./16.)*( e**2 / (pi*sw2 ) ) * F )

c1 = real(i* (-1./8.) * (5./576.)*(i* e**4 / ((pi**3)*(sw2*2) ) ) * ( - ( 93 + 16 * i ) )  * F )

c2 = real( i*(1./16.) * (5./27648.)*(  e**6 / ((pi**5)*(sw2**3) ) ) * (( 93 + 16 * i )**2) * F )

c3 = real( (-5./128.) * (-5./1492992.)*(i* e**8 / ((pi**7)*(sw2**4) ) ) * (( - ( 87 + 16 * i ) )**3) * F )

c4 = real( (7./256.) * (-25./429981696.)*( e**10 / ((pi**9)*(sw2**5) ) ) * (( 87 + 16 * i )**4) * F )

c5 = real( (-21./1024.) * (5./5159780352.)*( e**12 / ((pi**11)*(sw2**6) ) ) * (( - ( 87 + 16 * i ) )**5) * F )


print c0
print c1
print c2
print c3
print c4
print c5


plt.plot(x,y2,'-',color='grey',label='Explicit pole mass (fixed $\mu$)') #
plt.plot(x,y3,'-',color='black',label='Explicit pole mass ($\mu = \hat{M}$)') #
plt.plot(x,y1,'--',color='black',label='Iterative pole mass') #


plt.plot(x,(0.210901 )*1000.*ones(size(x)),'-.',color='red') #

#plt.plot(x,(c0)*1000.*ones(size(x)),'-.',color='blue')



xlabel(r"$\hat{M}_V$ (GeV)",fontsize=16)
ylabel(r"$\Delta M=M_V^+-M_V^0$ (MeV)",fontsize=16)
plt.ylim([0,310])
plt.xlim([30,100000])
plt.legend(loc=2,fontsize=12)


plt.savefig("mass_splittings.eps")

# plot pole masses


fig=plt.figure()

ax = fig.add_subplot(1,1,1)
ax.set_xscale('log')
#ax.set_yscale('log')

ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

A=np.genfromtxt('models/VDM/output/masses.txt',usecols=[0,1,2])
x=A[:,0]
y1=A[:,1]
y2=A[:,2]


plt.plot(x,y1-x,'-',color='black',label='MVp') #

plt.plot(x,y2-x,'--',color='black',label='MV0') #

#plt.plot(x,x,'--',color='black',label='Tree-level M') #



xlabel(r"MS-bar mass $MV_p$, $MV_0$ (GeV)",fontsize=16)
ylabel(r"Pole masses MV_p, MV_0$ (Gev)",fontsize=16)

plt.legend(loc=2)


plt.savefig("pole_masses.eps")
