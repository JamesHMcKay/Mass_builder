#!/usr/bin/python

# Mass Builder
 
# James McKay
# May 2017
 
#--- plot_VDM.py ---
 
# plot mass splitting for VDM model


import numpy as np
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

from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter

def smooth(A,pts,myorder):	
	x=A[:,0]
	B = zeros([pts,6]);
	xnew = np.linspace(x.min(),x.max(),pts)
	B[:,1:6] = spline(log10(x),A[:,1:6],log10(xnew),order=myorder)
	B[:,0] = xnew
	return B

fig=plt.figure()

ax = fig.add_subplot(1,1,1)
ax.set_xscale('log')

ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

plt.figure()

A=np.genfromtxt('models/VDM/output/mass_splittings.txt',usecols=[0,1,2,3,4,5])

x=A[:,0]
y1=A[:,1]*1000
y2=A[:,2]*1000
y3=A[:,3]*1000
y4=A[:,4]*1000
y5=A[:,5]*1000

xnew1 = np.linspace(x.min(),x.max(),5000)
y1_smooth = spline(log10(x),y1,log10(xnew1))
y2_smooth = spline(log10(x),y2,log10(xnew1))
y3_smooth = spline(log10(x),y3,log10(xnew1))
y4_smooth = spline(log10(x),y4,log10(xnew1))
y5_smooth = spline(log10(x),y5,log10(xnew1))


plt.plot(x,y1,'--',color='red',label="$Q=2 m_t$",linewidth=1.2) #
#plt.plot(x,y2,'--',color='black',label="$Q=2 m_t$",linewidth=1.2) #
#plt.plot(x,y3,'--',color='yellow',label="$Q=2 m_t$",linewidth=1.2) #
#plt.plot(x,y4,'--',color='blue',label="$Q=2 m_t$",linewidth=1.2) #
#plt.plot(x,y5,'--',color='#ff7f00',label="$Q=2 m_t$",linewidth=1.2) #

plt.plot(xnew1,y2_smooth,'--',color='black',label="$Q=m_t/2$",linewidth=1.2)
plt.plot(xnew1,y3_smooth,'--',color='yellow',label="$Q=\hat{M}/2$",linewidth=1.2)
plt.plot(xnew1,y4_smooth,'--',color='blue',label="$Q= \hat{M}$",linewidth=1.2)
plt.plot(xnew1,y5_smooth,'--',color='#ff7f00',label="$Q=2 \hat{M}$",linewidth=1.2)

# create vector with maximum and minimum values

#A[:,4] = A[:,5]
#A = smooth(A,5000,3)

x=A[:,0]

# create vector with maximum and minimum values

maximum = zeros(size(x))
minimum = zeros(size(x))

for i in range(0,len(x)):
	maximum[i] = max(A[i,1:6])*1000
	minimum[i] = min(A[i,1:6])*1000
	
plt.fill_between(x, minimum,maximum,color='green',alpha=0.4)


# plot the one-loop result

A=np.genfromtxt('models/VDM/output/mass_splittings.txt',usecols=[0,6])
x=A[:,0]
y1=A[:,1]*1000
plt.plot(x,y1,'-',color='black',linewidth=1.3)




plt.xscale('log')

xlabel(r"$\hat{M}_V$ (GeV)",fontsize=18)
ylabel(r"$\Delta M=M_V^+-M_V^0$ (MeV)",fontsize=18)

#plt.ylim([-250,-100])
plt.ylim([100,250])
plt.xlim([30,1e4])

#leg = plt.legend(loc='lower right')

#for legobj in leg.legendHandles:
#    legobj.set_linewidth(2.0)


plt.savefig("mass_splittings.pdf")
