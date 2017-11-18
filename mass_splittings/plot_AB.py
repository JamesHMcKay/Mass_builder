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
import matplotlib.gridspec as gridspec


from matplotlib import rc
rc('text', usetex=True)
plt.rc('font', family='Computer Modern Roman',weight='normal')



# plot the mass splitting uncertainties

fig = plt.figure()

ax = fig.add_subplot(111)


A=np.genfromtxt('mass_splittings/data2/deltam_1loop_MSSM_AB.txt',usecols=[0,1,2])

x=A[:,0]
lower_1loop=A[:,1]*1000
upper_1loop=A[:,2]*1000

xnew1 = np.linspace(x.min(),x.max(),5000)
lower_1loop_smooth = spline(log10(x),lower_1loop,log10(xnew1))
upper_1loop_smooth = spline(log10(x),upper_1loop,log10(xnew1))

ax.fill_between(xnew1, lower_1loop_smooth,upper_1loop_smooth,color='green',zorder=500,linewidth=1,label=r'$\mathrm{One\!-\!loop}$')	


A=np.genfromtxt('mass_splittings/data2/deltam_2loop_MSSM_AB.txt',usecols=[0,1,2])

x=A[:,0]
lower_2loop=A[:,1]*1000
upper_2loop=A[:,2]*1000

xnew1 = np.linspace(x.min(),x.max(),5000)
lower_2loop_smooth = spline(log10(x),lower_2loop,log10(xnew1))
upper_2loop_smooth = spline(log10(x),upper_2loop,log10(xnew1))

ax.fill_between(xnew1, lower_2loop_smooth,upper_2loop_smooth,color='red',alpha=0.3,zorder=100,linewidth=0,label=r'$\mathrm{Two\!-\!loop}$')	



ax.set_xlabel(r"$\mathrm{Degenerate\ mass}$ $\hat{M}$ $(\mathrm{GeV})$",fontsize=18)
ylabel(r"$\Delta M = M_{\mathrm{pole}}^+ - M_{\mathrm{pole}}^0$  $(\mathrm{MeV})$",fontsize=18)

ax.set_xscale('log')

ax.set_xlim([100,1e4])
ax.set_ylim([150,175])


ax.set_xticks([100,1e3,1e4])
ax.set_xticklabels(['$\mathrm{100}$','$\mathrm{10}^{\mathrm{3}}$','$\mathrm{10}^{\mathrm{4}}$'],fontsize=18)
plt.tick_params(labelsize=18)

plt.annotate(r"$\mathrm{Wino\ limit}$", xy=(2250,157), xytext=(2250,157),fontsize=18)
plt.annotate(r"$Q = 173.34\ (\mathrm{GeV})$", xy=(2250,155.5), xytext=(2250,155.5),fontsize=18)
leg = ax.legend(loc='lower right',frameon=False,fontsize=18)


plt.savefig("mass_splittings/figures2/deltam_AB_mssm.pdf",bbox_inches='tight')


fig = plt.figure()

ax = fig.add_subplot(111)


A=np.genfromtxt('mass_splittings/data2/deltam_1loop_MDM_AB.txt',usecols=[0,1,2])

x=A[:,0]
lower_1loop=A[:,1]*1000
upper_1loop=A[:,2]*1000

xnew1 = np.linspace(x.min(),x.max(),5000)
lower_1loop_smooth = spline(log10(x),lower_1loop,log10(xnew1))
upper_1loop_smooth = spline(log10(x),upper_1loop,log10(xnew1))

ax.fill_between(xnew1, lower_1loop_smooth,upper_1loop_smooth,color='green',zorder=500,linewidth=1,label=r'$\mathrm{One\!-\!loop}$')	


A=np.genfromtxt('mass_splittings/data2/deltam_2loop_MDM_AB.txt',usecols=[0,1,2])

x=A[:,0]
lower_2loop=A[:,1]*1000
upper_2loop=A[:,2]*1000

xnew1 = np.linspace(x.min(),x.max(),5000)
lower_2loop_smooth = spline(log10(x),lower_2loop,log10(xnew1))
upper_2loop_smooth = spline(log10(x),upper_2loop,log10(xnew1))

ax.fill_between(xnew1, lower_2loop_smooth,upper_2loop_smooth,color='red',alpha=0.3,zorder=100,linewidth=0,label=r'$\mathrm{Two\!-\!loop}$')	



ax.set_xlabel(r"$\mathrm{Degenerate\ mass}$ $\hat{M}$ $(\mathrm{GeV})$",fontsize=18)
ylabel(r"$\Delta M = M_{\mathrm{pole}}^+ - M_{\mathrm{pole}}^0$  $(\mathrm{MeV})$",fontsize=18)

ax.set_xscale('log')

ax.set_xlim([100,1e4])
ax.set_ylim([150,175])


ax.set_xticks([100,1e3,1e4])
ax.set_xticklabels(['$\mathrm{100}$','$\mathrm{10}^{\mathrm{3}}$','$\mathrm{10}^{\mathrm{4}}$'],fontsize=18)
plt.tick_params(labelsize=18)

plt.annotate(r"$\mathrm{Minimal\ dark\ matter}$", xy=(1500,157), xytext=(1500,157),fontsize=18)
plt.annotate(r"$Q = 173.34\ (\mathrm{GeV})$", xy=(1500,155.5), xytext=(1500,155.5),fontsize=18)
leg = ax.legend(bbox_to_anchor=(0.91,0.22),frameon=False,fontsize=18)

#for legobj in leg.legendHandles:
#    legobj.set_linewidth(1.0)
    
    
plt.savefig("mass_splittings/figures2/deltam_AB_mdm.pdf",bbox_inches='tight')
