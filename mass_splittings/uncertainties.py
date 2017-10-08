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


# for xi = 3 we turn threshold corrections off, it doesn't change
# the plot except for some strange bumps at low M
#
# all plots made using
# ./mass_builder -g -m xiMSSM -i models/xiMSSM/lists/ChiChi.txt
#
# input values are
# mw 80.385
# mz 91.1876
#
# different gauges are plotted by chaning Xi in models/MSSM/input.txt


rc('text', usetex=True)
plt.rc('font', family='Computer Modern Roman',weight='normal')

plt.figure()

# get the value of Xi
Mass=np.genfromtxt('models/MSSM/input.txt',usecols=[1])
Xi = Mass[11]

print("Xi = " +  str(Xi))

if (Xi == 0):
	A=np.genfromtxt('mass_splittings/data/uncertainties_0.txt',usecols=[0,1,2,3,4])
if (Xi == 1):
	A=np.genfromtxt('mass_splittings/data/uncertainties_1.txt',usecols=[0,1,2,3,4])
if (Xi == 3):
	A=np.genfromtxt('mass_splittings/data/uncertainties_3.txt',usecols=[0,1,2,3,4])	


x=A[:,0]
iterative_lower=A[:,1]*1000
iterative_upper=A[:,2]*1000
non_iterative_lower=A[:,3]*1000
non_iterative_upper=A[:,4]*1000


plt.fill_between(x, iterative_lower, iterative_upper,color='#7fbf7f')


plt.fill_between(x, non_iterative_lower, non_iterative_upper,color='grey')#color='#cc4c33')

# plot curves over top of uncertainty bands

if (Xi == 0):
	A=np.genfromtxt('mass_splittings/data/mass_splittings_iterative_0.txt',usecols=[0,1,2,3,4,5])
if (Xi == 1):
	A=np.genfromtxt('mass_splittings/data/mass_splittings_iterative_1.txt',usecols=[0,1,2,3,4,5])
if (Xi == 3):
	A=np.genfromtxt('mass_splittings/data/mass_splittings_iterative_3.txt',usecols=[0,1,2,3,4,5])


x22=A[:,0]
ylower=A[:,1]*1000
ymid=A[:,2]*1000
yupper1=A[:,3]*1000
yupper2=A[:,4]*1000
yupper3=A[:,5]*1000


y1=A[:,1]*1000
y2=A[:,2]*1000
y3=A[:,3]*1000
y4=A[:,4]*1000
y5=A[:,5]*1000

k1 = 0
k2 = 0
k3 = 0
k4 = 0
k5 = 0
z1 = zeros(size(x22))
z2 = zeros(size(x22))
z3 = zeros(size(x22))
z4 = zeros(size(x22))
z5 = zeros(size(x22))
x1 = zeros(size(x22))
x2 = zeros(size(x22))
x3 = zeros(size(x22))
x4 = zeros(size(x22))
x5 = zeros(size(x22))

for i in range(0,len(x22)):
	if (y1[i] != 0):
		z1[k1] = y1[i]
		x1[k1] = x22[i]
		k1 = k1 + 1
	if (y2[i] != 0):
		z2[k2] = y2[i]
		x2[k2] = x22[i]
		k2 = k2 + 1		
	if (y3[i] != 0):
		z3[k3] = y3[i]
		x3[k3] = x22[i]	
		k3 = k3 + 1	
	if (y4[i] != 0):
		z4[k4] = y4[i]
		x4[k4] = x22[i]	
		k4 = k4 + 1
	if (y5[i] != 0):
		z5[k5] = y5[i]
		x5[k5] = x22[i]	
		k5 = k5 + 1		
		

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


plt.plot(x1,z1,'-',color='red',linewidth=1.2) #
plt.plot(x2,z2,'-',color='black',linewidth=1.2) #
plt.plot(x3,z3,'-',color='yellow',linewidth=1.2) #
plt.plot(x4,z4,'-',color='blue',linewidth=1.2) #
plt.plot(x5,z5,'-',color='#ff7f00',linewidth=1.2) #



#plt.plot(x2,ylower,'-',color='red',label="$Q=2 m_t$",linewidth=1.2) #
#plt.plot(x2,ymid,'-',color='black',label="$Q=0.5 m_t$",linewidth=1.2) #
#plt.plot(x2,yupper1,'-',color='yellow',label="$Q=0.5 M$",linewidth=1.2) #
#plt.plot(x2,yupper2,'-',color='blue',label="$Q=1 M$",linewidth=1.2) #
#plt.plot(x2,yupper3,'-',color='#ff7f00',label="$Q=2 M$",linewidth=1.2) #

if (Xi == 0):
	A=np.genfromtxt('mass_splittings/data/mass_splittings_explicit_0.txt',usecols=[0,1,2,3,4,5])
if (Xi == 1):
	A=np.genfromtxt('mass_splittings/data/mass_splittings_explicit_1.txt',usecols=[0,1,2,3,4,5])
if (Xi == 3):
	A=np.genfromtxt('mass_splittings/data/mass_splittings_explicit_3.txt',usecols=[0,1,2,3,4,5])


x2=A[:,0]
ylower=A[:,1]*1000
ymid=A[:,2]*1000
yupper1=A[:,3]*1000
yupper2=A[:,4]*1000
yupper3=A[:,5]*1000


xlabel(r"Degenerate mass, $\hat{M}$ $($GeV$)$",fontsize=18)
ylabel(r"$\Delta M = M_{\mathrm{pole}}^+ - M_{\mathrm{pole}}^0$  $($MeV$)$",fontsize=18)


plt.tick_params(labelsize=14)

plt.xscale('log')
plt.xlim([min(x),max(x)])
plt.ylim([0,250])





if (Xi == 3):
	plt.annotate(r'Fried-Yennie gauge $(\xi = 3)$', xy=(20,20), xytext=(20,20),fontsize=16)
	plt.savefig("mass_splittings/figures/deltam_uncertainties_xi_3.pdf")

if (Xi == 1):
	plt.annotate(r"Feynman-'t Hooft gauge $(\xi = 1)$", xy=(20,20), xytext=(20,20),fontsize=16)
	plt.savefig("mass_splittings/figures/deltam_uncertainties_xi_1.pdf",transparent=True)

if (Xi == 0):
	plt.annotate(r'Landau gauge $(\xi = 0)$', xy=(20,20), xytext=(20,20),fontsize=16)
	plt.savefig("mass_splittings/figures/deltam_uncertainties_xi_0.pdf")


#######################
### make the first legend ###


# create a figure for the data
figData = pylab.figure()
ax = pylab.gca()

plt.plot(x2,ylower,'-',color='red',label="$Q=2 m_t$") #
plt.plot(x2,ymid,'-',color='black',label="$Q=m_t/2$")
plt.plot(x2,yupper1,'-',color='yellow',label="$Q=\hat{M}/2$")#
plt.plot(x2,yupper2,'-',color='blue',label="$Q= \hat{M}$")
plt.plot(x2,yupper3,'-',color='#ff7f00',label="$Q=2 \hat{M}$")



# create a second figure for the legend
figLegend = pylab.figure(figsize = (8.5,0.6))

# produce a legend for the objects in the other figure
leg = pylab.figlegend(*ax.get_legend_handles_labels(), loc = 'upper left',ncol=5,frameon=False)

for legobj in leg.legendHandles:
    legobj.set_linewidth(3.0)


# save the two figures to files
figLegend.savefig("mass_splittings/figures//legend_1.pdf",transparent=True)



#######################
### make the second legend ###


# create a figure for the data
figData = pylab.figure()
ax = pylab.gca()

# plot something, it doesn't matter what
pylab.plot(x2,ylower,'-',color='#7fbf7f',label="iterative uncertainty") #
pylab.plot(x2,ymid,'-',color='grey',label="non-iterative uncertainty") #

# create a second figure for the legend
figLegend = pylab.figure(figsize = (8.5,0.6))

# produce a legend for the objects in the other figure
leg = pylab.figlegend(*ax.get_legend_handles_labels(), loc = 'upper center',ncol=5,frameon=False,fontsize=16)

for legobj in leg.legendHandles:
    legobj.set_linewidth(8.0)


# save the two figures to files
figLegend.savefig("mass_splittings/figures//legend_2.pdf",transparent=True)


