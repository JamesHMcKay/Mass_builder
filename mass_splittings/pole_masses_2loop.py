#!/usr/bin/python

import numpy as np
import pylab
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

A=np.genfromtxt('mass_splittings/data2/deltam_2loop_it.txt',usecols=[0,1,2,3,4,5])

x=A[:,0]
y1=A[:,1]
y2=A[:,2]
y3=A[:,3]
y4=A[:,4]
y5=A[:,4]

k1 = 0
k2 = 0
k3 = 0
k4 = 0
k5 = 0
z1 = zeros(size(x))
z2 = zeros(size(x))
z3 = zeros(size(x))
z4 = zeros(size(x))
z5 = zeros(size(x))
x1 = zeros(size(x))
x2 = zeros(size(x))
x3 = zeros(size(x))
x4 = zeros(size(x))
x5 = zeros(size(x))

for i in range(0,len(x)):
	
	if (y1[i] != 0):
		z1[k1] = y1[i]
		x1[k1] = x[i]
		k1 = k1 + 1
#else:
		#z1 = zeros(size(x))
		#x1 = zeros(size(x))
	if (y2[i] != 0):
		z2[k2] = y2[i]
		x2[k2] = x[i]
		k2 = k2 + 1
	#else:
		#z2 = zeros(size(x))
		#x2 = zeros(size(x))		
	if (y3[i] != 0):
		z3[k3] = y3[i]
		x3[k3] = x[i]	
		k3 = k3 + 1	
	#else:
		#z3 = zeros(size(x))
		#x3 = zeros(size(x))	
	if (y4[i] != 0):
		z4[k4] = y4[i]
		x4[k4] = x[i]	
		k4 = k4 + 1
	#else:
		#z4 = zeros(size(x))
		#x4 = zeros(size(x))	
	if (y5[i] != 0):
		z5[k5] = y5[i]
		x5[k5] = x[i]	
		k5 = k5 + 1
	#else:
		#z5 = zeros(size(x))
		#x5 = zeros(size(x))			
		

z1 = np.trim_zeros(z1)
z2 = np.trim_zeros(z2)
z3 = np.trim_zeros(z3)
z4 = np.trim_zeros(z4)
z5 = np.trim_zeros(z5)
x1 = np.trim_zeros(x1)
x2 = np.trim_zeros(x2)
x3 = np.trim_zeros(x3)
x4 = np.trim_zeros(x4)
x5 = np.trim_zeros(x5)

plt.plot(x1,z1,'-',color='yellow',linewidth=1.2) #
plt.plot(x2,z2,'-',color='red',linewidth=1.2) #
plt.plot(x3,z3,'-',color='black',linewidth=1.2) #
plt.plot(x4,z4,'-',color='blue',linewidth=1.2) #
plt.plot(x5,z5,'-',color='#ff7f00',linewidth=1.2) #


xlabel(r"Degenerate mass, $\hat{M}$ (GeV)",fontsize=18)
ylabel(r"$\Delta M$",fontsize=18)
plt.xscale('log')
plt.yscale('log')

#plt.xlim([10,1e4])
plt.ylim([0.0,100])

plt.tick_params(labelsize=14)

#plt.annotate(r'Fried-Yennie gauge $(\xi = 3)$', xy=(20,0.87), xytext=(20,0.87),fontsize=16)
plt.savefig("mass_splittings/figures2/deltam_it.pdf")





########################################################

########################################################

########################################################

########################################################
#### Compare iterative and explicit two-loop pole masses
########################################################

########################################################

########################################################

########################################################

plt.figure()


A=np.genfromtxt('mass_splittings/data2/pole_mass_unc_2loop_c_it.txt',usecols=[0,1,2,3,4,5])

x=A[:,0]
y1=A[:,1]
y2=A[:,2]
y3=A[:,3]
y4=A[:,4]
y5=A[:,4]

k1 = 0
k2 = 0
k3 = 0
k4 = 0
k5 = 0
z1 = zeros(size(x))
z2 = zeros(size(x))
z3 = zeros(size(x))
z4 = zeros(size(x))
z5 = zeros(size(x))
x1 = zeros(size(x))
x2 = zeros(size(x))
x3 = zeros(size(x))
x4 = zeros(size(x))
x5 = zeros(size(x))

for i in range(0,len(x)):
	
	if (y1[i] != 0):
		z1[k1] = y1[i]
		x1[k1] = x[i]
		k1 = k1 + 1
#else:
		#z1 = zeros(size(x))
		#x1 = zeros(size(x))
	if (y2[i] != 0):
		z2[k2] = y2[i]
		x2[k2] = x[i]
		k2 = k2 + 1
	#else:
		#z2 = zeros(size(x))
		#x2 = zeros(size(x))		
	if (y3[i] != 0):
		z3[k3] = y3[i]
		x3[k3] = x[i]	
		k3 = k3 + 1	
	#else:
		#z3 = zeros(size(x))
		#x3 = zeros(size(x))	
	if (y4[i] != 0):
		z4[k4] = y4[i]
		x4[k4] = x[i]	
		k4 = k4 + 1
	#else:
		#z4 = zeros(size(x))
		#x4 = zeros(size(x))	
	if (y5[i] != 0):
		z5[k5] = y5[i]
		x5[k5] = x[i]	
		k5 = k5 + 1
	#else:
		#z5 = zeros(size(x))
		#x5 = zeros(size(x))			
		

z1 = np.trim_zeros(z1)
z2 = np.trim_zeros(z2)
z3 = np.trim_zeros(z3)
z4 = np.trim_zeros(z4)
z5 = np.trim_zeros(z5)
x1 = np.trim_zeros(x1)
x2 = np.trim_zeros(x2)
x3 = np.trim_zeros(x3)
x4 = np.trim_zeros(x4)
x5 = np.trim_zeros(x5)

plt.plot(x1,z1/x1,'-',color='yellow',linewidth=1.2) #
plt.plot(x2,z2/x2,'-',color='red',linewidth=1.2) #
plt.plot(x3,z3/x3,'-',color='black',linewidth=1.2) #
plt.plot(x4,z4/x4,'-',color='blue',linewidth=1.2) #
plt.plot(x5,z5/x5,'-',color='#ff7f00',linewidth=1.2) #

A=np.genfromtxt('mass_splittings/data2/pole_mass_unc_2loop_c_AB.txt',usecols=[0,1,2,3,4,5])

x=A[:,0]
ylower=A[:,1]
ymid=A[:,2]
yupper1=A[:,3]
yupper2=A[:,4]
yupper3=A[:,5]

plt.plot(x,ylower/x,'--',color='yellow',linewidth=1.2) #
plt.plot(x,ymid/x,'--',color='red',linewidth=1.2) #
plt.plot(x,yupper1/x,'--',color='black',linewidth=1.2) #
plt.plot(x,yupper2/x,'--',color='blue',linewidth=1.2) #
plt.plot(x,yupper3/x,'--',color='#ff7f00',linewidth=1.2) #


# create vector with maximum and minimum values

maximum = zeros(size(x))
minimum = zeros(size(x))

for i in range(0,len(x)):
	maximum[i] = max(A[i,1:6])
	minimum[i] = min(A[i,1:6])


plt.fill_between(x, minimum/x,maximum/x,color='grey',alpha=0.5,label='non-iterative')


xlabel(r"Degenerate mass, $\hat{M}$ (GeV)",fontsize=18)
ylabel(r"$M^0_p /\hat{M}$",fontsize=18)
plt.xscale('log')

plt.xlim([10,1e4])
plt.ylim([0.85,1.10])

plt.tick_params(labelsize=14)
plt.legend(loc=2)

#plt.annotate(r'Fried-Yennie gauge $(\xi = 3)$', xy=(20,0.87), xytext=(20,0.87),fontsize=16)
plt.savefig("mass_splittings/figures2/pole_masses_c_it_compare.pdf")




#######################
### make the legend ###


# create a figure for the data
figData = pylab.figure()
ax = pylab.gca()

# plot something, it doesn't matter what
pylab.plot(x,x,'-',color='#7fbf7f',label="iterative uncertainty") #
pylab.plot(x,x,'-',color='grey',alpha=0.5,label="non-iterative uncertainty") #

# create a second figure for the legend
figLegend = pylab.figure(figsize = (8.5,0.6))

# produce a legend for the objects in the other figure
leg = pylab.figlegend(*ax.get_legend_handles_labels(), loc = 'upper center',ncol=5,frameon=False,fontsize=16)

for legobj in leg.legendHandles:
    legobj.set_linewidth(8.0)


# save the two figures to files
figLegend.savefig("mass_splittings/figures//legend_3.pdf")

