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



# plot the BRs

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

A=np.genfromtxt('mass_splittings/data2/decays_BR_1loop_MDM.txt',usecols=[0,3,4,5])

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



ax.set_xticks([1,10,100,1000,10000])
ax.set_xticklabels(['$\mathrm{1}$','$\mathrm{10}$','$\mathrm{100}$','$\mathrm{10}^{\mathrm{3}}$','$\mathrm{10}^{\mathrm{4}}$'],fontsize=18)

ax.set_xlim([1,10000])
ax.set_ylim([1e-3,1.2])

leg = ax.legend(loc='lower left',frameon=False,fontsize=18)

plt.annotate(r"$\mathrm{Minimal\ dark\ matter}$", xy=(1.5,0.016), xytext=(1.5,0.016),fontsize=18)
plt.annotate(r"$Q = 173.34\ (\mathrm{GeV})$", xy=(1.5,0.01), xytext=(1.5,0.01),fontsize=18)



plt.savefig("mass_splittings/figures2/BR_MDM.pdf",bbox_inches='tight')




# plot the BRs for the second mass splitting

fig = plt.figure()

ax = fig.add_subplot(111)

# 2-loop curves
A=np.genfromtxt('mass_splittings/data2/decays2_BR_2loop_MDM.txt',usecols=[0,3,4,5,6])

x=A[:,0]
BRpi=A[:,1]
BRmu=A[:,2]
BRe=A[:,3]
BRka=A[:,4]

BRpi2 = BRpi[BRpi>0]
xpi = x[BRpi>0]

BRmu2 = BRmu[BRmu>0]
xmu = x[BRmu>0]

BRe2 = BRe[BRe>0]
xe = x[BRe>0]

BRka2 = BRe[BRka>0]
xka = x[BRka>0]

ax.plot(xpi,BRpi2,'-',color='blue',label=r'$\chi^{++}\rightarrow\chi^+e^+\nu_e$',linewidth=1,zorder=500) #
ax.plot(xmu,BRmu2,'-',color='red',label=r'$\chi^{++}\rightarrow\chi^+\mu^+\nu_{\mu}$',linewidth=1,zorder=500) #
ax.plot(xe,BRe2,'-',color='black',label=r'$\chi^{++}\rightarrow\chi^+\pi^+$',linewidth=1,zorder=500) #
ax.plot(xka,BRka2,'-',color='orange',label=r'$\chi^{++}\rightarrow\chi^+k^+$',linewidth=1,zorder=500) #

# 1-loop curves

A=np.genfromtxt('mass_splittings/data2/decays2_BR_1loop_MDM.txt',usecols=[0,3,4,5,6])

x=A[:,0]
BRpi=A[:,1]
BRmu=A[:,2]
BRe=A[:,3]
BRka=A[:,4]

BRpi2 = BRpi[BRpi>0]
xpi = x[BRpi>0]

BRmu2 = BRmu[BRmu>0]
xmu = x[BRmu>0]

BRe2 = BRe[BRe>0]
xe = x[BRe>0]

BRka2 = BRe[BRka>0]
xka = x[BRka>0]


ax.plot(xpi,BRpi2,'--',color='blue',linewidth=1,zorder=500) #
ax.plot(xmu,BRmu2,'--',color='red',linewidth=1,zorder=500) #
ax.plot(xe,BRe2,'--',color='black',linewidth=1,zorder=500) #
ax.plot(xka,BRka2,'--',color='orange',linewidth=1,zorder=500) #

ax.set_xlabel(r"$\mathrm{Degenerate\ mass}$ $\hat{M}$ $(\mathrm{GeV})$",fontsize=18)
ax.set_ylabel(r"$\mathrm{Branching\ ratio}$",fontsize=18)

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_yticks([1e-3,1e-2,1e-1,1])
ax.set_yticklabels(['$\mathrm{10}^{\mathrm{-3}}$','$\mathrm{0.01}$','$\mathrm{0.1}$','$\mathrm{1}$'],fontsize=18)



ax.set_xticks([1,10,100,1000,10000])
ax.set_xticklabels(['$\mathrm{1}$','$\mathrm{10}$','$\mathrm{100}$','$\mathrm{10}^{\mathrm{3}}$','$\mathrm{10}^{\mathrm{4}}$'],fontsize=18)

ax.set_xlim([1,10000])
ax.set_ylim([1e-3,1.2])

leg = ax.legend(loc='lower right',frameon=False,fontsize=18)

plt.annotate(r"$\mathrm{Minimal\ dark\ matter}$", xy=(180,0.045), xytext=(180,0.045),fontsize=18)
plt.annotate(r"$Q = 173.34\ (\mathrm{GeV})$", xy=(180,0.03), xytext=(180,0.03),fontsize=18)



plt.savefig("mass_splittings/figures2/BR2_MDM.pdf",bbox_inches='tight')






# plot the mass splitting uncertainties

fig = plt.figure()

ax = fig.add_subplot(111)


A=np.genfromtxt('mass_splittings/data2/uncertainties_MDM.txt',usecols=[0,1,2,3,4,5,6,7,8])

x=A[:,0]
lower_1loop=A[:,1]*1000
upper_1loop=A[:,2]*1000
lower2_1loop=A[:,3]*1000
upper2_1loop=A[:,4]*1000

lower_2loop=A[:,5]*1000
upper_2loop=A[:,6]*1000
lower2_2loop=A[:,7]*1000
upper2_2loop=A[:,8]*1000

xnew1 = np.linspace(x.min(),x.max(),5000)
lower_1loop_smooth = spline(log10(x),lower_1loop,log10(xnew1))
upper_1loop_smooth = spline(log10(x),upper_1loop,log10(xnew1))
lower_2loop_smooth = spline(log10(x),lower_2loop,log10(xnew1))
upper_2loop_smooth = spline(log10(x),upper_2loop,log10(xnew1))


lower2_1loop_smooth = spline(log10(x),lower2_1loop,log10(xnew1))
upper2_1loop_smooth = spline(log10(x),upper2_1loop,log10(xnew1))
lower2_2loop_smooth = spline(log10(x),lower2_2loop,log10(xnew1))
upper2_2loop_smooth = spline(log10(x),upper2_2loop,log10(xnew1))


ax.fill_between(xnew1, lower_1loop_smooth,upper_1loop_smooth,color='green',alpha=0.5,zorder=100,linewidth=0,label=r'$\mathrm{One\!-\!loop}$')	
ax.fill_between(xnew1, lower_2loop_smooth,upper_2loop_smooth,color='red',alpha=0.5,zorder=100,linewidth=0,label=r'$\mathrm{Two\!-\!loop}$')	


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


plt.savefig("mass_splittings/figures2/deltam_mdm.pdf",bbox_inches='tight')

## plot the second mass splitting

fig = plt.figure()

ax = fig.add_subplot(111)

ax.fill_between(xnew1, lower2_1loop_smooth,upper2_1loop_smooth,color='green',alpha=0.5,zorder=100,linewidth=0,label=r'$\mathrm{One\!-\!loop}$')	
ax.fill_between(xnew1, lower2_2loop_smooth,upper2_2loop_smooth,color='red',alpha=0.5,zorder=100,linewidth=0,label=r'$\mathrm{Two\!-\!loop}$')	


ax.set_xlabel(r"$\mathrm{Degenerate\ mass}$ $\hat{M}$ $(\mathrm{GeV})$",fontsize=18)
ylabel(r"$\Delta M = M_{\mathrm{pole}}^{++} - M_{\mathrm{pole}}^0$  $(\mathrm{MeV})$",fontsize=18)

ax.set_xscale('log')

ax.set_xlim([100,1e4])
ax.set_ylim([600,700])


ax.set_xticks([100,1e3,1e4])
ax.set_xticklabels(['$\mathrm{100}$','$\mathrm{10}^{\mathrm{3}}$','$\mathrm{10}^{\mathrm{4}}$'],fontsize=18)
plt.tick_params(labelsize=18)

plt.annotate(r"$\mathrm{Minimal\ dark\ matter}$", xy=(1500,628), xytext=(1500,628),fontsize=18)
plt.annotate(r"$Q = 173.34\ (\mathrm{GeV})$", xy=(1500,622), xytext=(1500,622),fontsize=18)
leg = ax.legend(bbox_to_anchor=(0.91,0.22),frameon=False,fontsize=18)


plt.savefig("mass_splittings/figures2/deltam2_mdm.pdf",bbox_inches='tight')



###### plot decay lifetimes


gs = gridspec.GridSpec(20, 20)


fig=plt.figure()
ax = fig.add_subplot(111)


ax.set_xscale('log')
ax.set_yscale('log')

#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

A=np.genfromtxt('mass_splittings/data2/decays_MDM.txt',usecols=[0,1,2,3,4,5,6,7,8])


x=A[:,0]
lower_1loop=A[:,1]*2.9979e11
upper_1loop=A[:,2]*2.9979e11
lower2_1loop=A[:,3]*2.9979e11
upper2_1loop=A[:,4]*2.9979e11


lower_2loop=A[:,5]*2.9979e11
upper_2loop=A[:,6]*2.9979e11
lower2_2loop=A[:,7]*2.9979e11
upper2_2loop=A[:,8]*2.9979e11


plt.fill_between(x, lower_1loop,upper_1loop,color='green',alpha=0.5,zorder=100,linewidth=0,label=r'$\mathrm{One\!-\!loop}$')	
plt.fill_between(x, lower_2loop,upper_2loop,color='red',alpha=0.5,zorder=100,linewidth=0,label=r'$\mathrm{Two\!-\!loop}$')	

plt.fill_between(x, lower2_1loop,upper2_1loop,color='green',alpha=0.5,zorder=100,linewidth=0)	
plt.fill_between(x, lower2_2loop,upper2_2loop,color='red',alpha=0.5,zorder=100,linewidth=0)	


plt.annotate(r"$\chi^{++}$", xy=(5,1e3), xytext=(5,1e3),fontsize=18)

plt.annotate(r"$\chi^{+}$", xy=(20,10e4), xytext=(20,10e4),fontsize=18)


plt.annotate(r"$\mathrm{Minimal\ dark\ matter}$", xy=(300,3), xytext=(300,3),fontsize=18)
plt.annotate(r"$Q=173.15$ $(\mathrm{GeV})$", xy=(300,1), xytext=(300,1),fontsize=18)

xlabel(r"$\mathrm{Degenerate\ mass}$ $\hat{M}$ $(\mathrm{GeV})$",fontsize=18)
ylabel(r"$\tau$ $(\mathrm{mm}/\mathrm{c})$",fontsize=18)

plt.xlim([3,1e4])
plt.ylim([0.1,1e6])


ax.set_xticks([1,10,100,1e3,1e4])
ax.set_xticklabels(['$\mathrm{1}$','$\mathrm{10}$','$\mathrm{100}$','$\mathrm{10}^{\mathrm{3}}$','$\mathrm{10}^{\mathrm{4}}$'],fontsize=18)


ax.set_yticks([10,100,1000,1e4,1e5,1e6])
ax.set_yticklabels(['$\mathrm{10}$','$\mathrm{100}$','$\mathrm{1000}$','$\mathrm{10}^{\mathrm{3}}$','$\mathrm{10}^{\mathrm{4}}$','$\mathrm{10}^{\mathrm{5}}$','$\mathrm{10}^{\mathrm{6}}$'],fontsize=18)

leg = ax.legend(loc='lower left',frameon=False,fontsize=18)


plt.tick_params(labelsize=18)



### do a zoomed in plot of large M limit

ax2 = fig.add_subplot(gs[1:5, 13:19])
ax2.set_xscale('log')
#ax.set_yscale('log')

ax2.fill_between(x, lower_1loop,upper_1loop,color='green',alpha=0.5,zorder=100,linewidth=0)	
ax2.fill_between(x, lower_2loop,upper_2loop,color='red',alpha=0.5,zorder=100,linewidth=0)	

#ax2.set_xlabel(r"$\hat{M}$ $(\mathrm{GeV})$",fontsize=16)
ax2.set_ylabel(r"$\tau$ $(\mathrm{mm}/\mathrm{c})$",fontsize=14)

ax2.set_xlim([200,1e4])
ax2.set_ylim([45,70])

ax2.set_xticks([])

ax2.set_yticks([50,60,70])
ax2.set_yticklabels(['$\mathrm{50}$','$\mathrm{60}$','$\mathrm{70}$'],fontsize=16)
ax2.set_title(r'$\chi^{+}$',fontsize=14)


### do a zoomed in plot of large M limit

ax3 = fig.add_subplot(gs[6:10, 13:19])
ax3.set_xscale('log')
#ax.set_yscale('log')

ax3.fill_between(x, lower2_1loop,upper2_1loop,color='green',alpha=0.5,zorder=100,linewidth=0)	
ax3.fill_between(x, lower2_2loop,upper2_2loop,color='red',alpha=0.5,zorder=100,linewidth=0)	

ax3.set_xlabel(r"$\hat{M}$ $(\mathrm{GeV})$",fontsize=14,labelpad=0.1)
ax3.set_ylabel(r"$\tau$ $(\mathrm{mm}/\mathrm{c})$",fontsize=14)

ax3.set_xlim([200,1e4])
ax3.set_ylim([0.35,0.5])

ax3.set_yticks([0.3,0.4,0.5])
ax3.set_yticklabels(['$\mathrm{0.3}$','$\mathrm{0.4}$','$\mathrm{0.5}$'],fontsize=16)

ax3.set_title(r'$\chi^{++}$',fontsize=14)


plt.savefig("mass_splittings/figures2/decays_MDM.pdf",bbox_inches='tight')
