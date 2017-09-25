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


plt.figure()

A=np.genfromtxt('mass_splittings/data2/deltam_2loop_it.txt',usecols=[0,1,2,3])

x=A[:,0]
y1=A[:,1]*1000
y2=A[:,2]*1000
y3=A[:,3]*1000

xnew1 = np.linspace(x.min(),x.max(),5000)
xnew2 = np.linspace(2000,x.max(),5000)
xnew3 = np.linspace(6000,x.max(),5000)
y1_smooth = spline(log10(x),y1,log10(xnew1))
y2_smooth = spline(log10(x),y2,log10(xnew2))
y3_smooth = spline(log10(x),y3,log10(xnew3))


#plt.plot(x,y1,'--',color='black',linewidth=1.2)
#plt.plot(x,y2,'--',color='blue',linewidth=1.2)
#plt.plot(x,y3,'--',color='red',linewidth=1.2)

plt.plot(xnew1,y1_smooth,'-',color='blue',linewidth=1.2,label="$Q=10$ GeV")
plt.plot(xnew2,y2_smooth,'-',color='red',linewidth=1.2,label="$Q=100$ GeV")
plt.plot(xnew3,y3_smooth,'-',color='black',linewidth=1.2,label="$Q=1000$ GeV")



xlabel(r"Degenerate mass, $\hat{M}$ $($GeV$)$",fontsize=18)
ylabel(r"$\Delta M = M_p^+ - M_p^0$  $($MeV$)$",fontsize=18)

plt.xscale('log')
#plt.yscale('log')

#plt.xlim([10,1e4])
plt.ylim([0.0,16])

plt.tick_params(labelsize=14)


leg = plt.legend(loc='upper right')

for legobj in leg.legendHandles:
    legobj.set_linewidth(2.0)

plt.annotate(r"Feynman-'t Hooft gauge $(\xi = 1)$", xy=(1100,15), xytext=(1100,15),fontsize=16)
plt.annotate(r"partial two-loop", xy=(1100,14), xytext=(1100,14),fontsize=16)
plt.savefig("mass_splittings/figures2/deltam_2loop_iterative.pdf")




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

