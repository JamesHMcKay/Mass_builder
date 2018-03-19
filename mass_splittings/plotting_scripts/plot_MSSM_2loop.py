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
A=np.genfromtxt('mass_splittings/data2/decays_BR_2loop_MSSM.txt',usecols=[0,3,4,5])

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


ax.plot(xpi,BRpi2,'-',color='blue',label=r'$\chi^+\rightarrow\chi^0\pi^+$',linewidth=1,zorder=500) #
ax.plot(xmu,BRmu2,'-',color='red',label=r'$\chi^+\rightarrow\chi^0\mu^+\nu_{\mu}$',linewidth=1,zorder=500) #
ax.plot(xe,BRe2,'-',color='black',label=r'$\chi^+\rightarrow\chi^0e^+\nu_e$',linewidth=1,zorder=500) #

# 1-loop curves

A=np.genfromtxt('mass_splittings/data2/decays_BR_1loop_MSSM.txt',usecols=[0,3,4,5])

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


ax.plot(xpi,BRpi2,':',color='blue',linewidth=1.5,zorder=500) #
ax.plot(xmu,BRmu2,':',color='red',linewidth=1.5,zorder=500) #
ax.plot(xe,BRe2,':',color='black',linewidth=1.5,zorder=500) #




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

plt.annotate(r"$\mathrm{Wino\ limit}$", xy=(1.5,0.016), xytext=(1.5,0.016),fontsize=18)
plt.annotate(r"$Q = 173.34\ (\mathrm{GeV})$", xy=(1.5,0.01), xytext=(1.5,0.01),fontsize=18)



plt.savefig("mass_splittings/figures/BR_MSSM.pdf",bbox_inches='tight')

# plot the mass splitting uncertainties

fig = plt.figure()

ax = fig.add_subplot(111)


A=np.genfromtxt('mass_splittings/data2/uncertainties_MSSM.txt',usecols=[0,1,2,3,4])

x=A[:,0]
lower_1loop=A[:,1]*1000
upper_1loop=A[:,2]*1000
lower_2loop=A[:,3]*1000
upper_2loop=A[:,4]*1000

xnew1 = np.linspace(x.min(),x.max(),5000)
lower_1loop_smooth = spline(log10(x),lower_1loop,log10(xnew1))
upper_1loop_smooth = spline(log10(x),upper_1loop,log10(xnew1))
lower_2loop_smooth = spline(log10(x),lower_2loop,log10(xnew1))
upper_2loop_smooth = spline(log10(x),upper_2loop,log10(xnew1))

ax.fill_between(xnew1, lower_1loop_smooth,upper_1loop_smooth,color='green',alpha=0.5,zorder=100,linewidth=0,label=r'$\mathrm{One\!-\!loop}$')	
ax.fill_between(xnew1, lower_2loop_smooth,upper_2loop_smooth,color='red',alpha=0.5,zorder=100,linewidth=0,label=r'$\mathrm{Two\!-\!loop}$')	

# plot estimated uncertainties

A=np.genfromtxt('mass_splittings/data2/estimated_uncertainties_MSSM.txt',usecols=[0,1,2])

x=A[:,0]
lower_1loop=A[:,1]*1000
upper_1loop=A[:,2]*1000

xnew1 = np.linspace(x.min(),x.max(),5000)
lower_1loop_smooth = spline(log10(x),lower_1loop,log10(xnew1))
upper_1loop_smooth = spline(log10(x),upper_1loop,log10(xnew1))

ax.fill_between(xnew1, lower_1loop_smooth,upper_1loop_smooth,color='green',alpha=0.2,zorder=10,linewidth=0)	





ax.set_xlabel(r"$\mathrm{Degenerate\ mass}$ $\hat{M}$ $(\mathrm{GeV})$",fontsize=18)
ylabel(r"$\Delta M = M_{\mathrm{pole}}^+ - M_{\mathrm{pole}}^0$  $(\mathrm{MeV})$",fontsize=18)

ax.set_xscale('log')

ax.set_xlim([100,1e4])
ax.set_ylim([150,175])


ax.set_xticks([100,1e3,1e4])
ax.set_xticklabels(['$\mathrm{100}$','$\mathrm{10}^{\mathrm{3}}$','$\mathrm{10}^{\mathrm{4}}$'],fontsize=18)
plt.tick_params(labelsize=18)

plt.annotate(r"$\mathrm{Wino\ limit}$", xy=(2250,157), xytext=(2250,157),fontsize=18)
plt.annotate(r"$m_t/2 \leq Q \leq 2m_t$", xy=(2250,155.5), xytext=(2250,155.5),fontsize=18)
leg = ax.legend(loc='lower right',frameon=False,fontsize=18)


plt.savefig("mass_splittings/figures/deltam_mssm.pdf",bbox_inches='tight')



###### plot decay lifetimes


gs = gridspec.GridSpec(10, 10)


fig=plt.figure()
ax = fig.add_subplot(111)


ax.set_xscale('log')
ax.set_yscale('log')

#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

A=np.genfromtxt('mass_splittings/data2/decays_MSSM.txt',usecols=[0,1,2,3,4])


x=A[:,0]
lower_1loop=A[:,1]*2.9979e11
upper_1loop=A[:,2]*2.9979e11
lower_2loop=A[:,3]*2.9979e11
upper_2loop=A[:,4]*2.9979e11


A=np.genfromtxt('mass_splittings/data2/estimated_decays_MSSM.txt',usecols=[0,1,2])


x2=A[:,0]
unc_lower_1loop=A[:,1]*2.9979e11
unc_upper_1loop=A[:,2]*2.9979e11


plt.fill_between(x, lower_1loop,upper_1loop,color='green',alpha=0.5,zorder=100,linewidth=0,label=r'$\mathrm{One\!-\!loop}$')	
plt.fill_between(x2, unc_lower_1loop,unc_upper_1loop,color='green',alpha=0.2,zorder=10,linewidth=0)
plt.fill_between(x, lower_2loop,upper_2loop,color='red',alpha=0.5,zorder=100,linewidth=0,label=r'$\mathrm{Two\!-\!loop}$')	


plt.annotate(r"$\mathrm{Wino\ limit}$", xy=(4,300), xytext=(4,300),fontsize=18)
plt.annotate(r"$m_t/2 \leq Q \leq 2m_t$", xy=(4,150), xytext=(4,150),fontsize=18)

xlabel(r"$\mathrm{Degenerate\ mass}$ $\hat{M}$ $(\mathrm{GeV})$",fontsize=18)
ylabel(r"$\tau$ $(\mathrm{mm}/\mathrm{c})$",fontsize=18)

plt.xlim([3,1e4])
plt.ylim([10,1e6])


ax.set_xticks([10,100,1e3,1e4])
ax.set_xticklabels(['$\mathrm{10}$','$\mathrm{100}$','$\mathrm{10}^{\mathrm{3}}$','$\mathrm{10}^{\mathrm{4}}$'],fontsize=18)


ax.set_yticks([10,100,1000,1e4,1e5,1e6])
ax.set_yticklabels(['$\mathrm{10}$','$\mathrm{100}$','$\mathrm{1000}$','$\mathrm{10}^{\mathrm{3}}$','$\mathrm{10}^{\mathrm{4}}$','$\mathrm{10}^{\mathrm{5}}$','$\mathrm{10}^{\mathrm{6}}$'],fontsize=18)

leg = ax.legend(loc='lower left',frameon=False,fontsize=18)


plt.tick_params(labelsize=18)



### do a zoomed in plot of large M limit

ax2 = fig.add_subplot(gs[1:5, 5:9])
ax2.set_xscale('log')
#ax.set_yscale('log')

ax2.fill_between(x, lower_1loop,upper_1loop,color='green',alpha=0.5,zorder=100,linewidth=0)	
ax2.fill_between(x2, unc_lower_1loop,unc_upper_1loop,color='green',alpha=0.2,zorder=10,linewidth=0)	
ax2.fill_between(x, lower_2loop,upper_2loop,color='red',alpha=0.5,zorder=100,linewidth=0)	

ax2.set_xlabel(r"$\mathrm{Degenerate\ mass}$ $\hat{M}$ $(\mathrm{GeV})$",fontsize=16)
ax2.set_ylabel(r"$\tau$ $(\mathrm{mm}/\mathrm{c})$",fontsize=16)

ax2.set_xlim([200,1e4])
ax2.set_ylim([50,60])


plt.savefig("mass_splittings/figures/decays_MSSM.pdf",bbox_inches='tight')


#######   determine polynomial fits through the data


print 'Polynomial fit to two-loop result for deltam1'
A=np.genfromtxt('mass_splittings/data2/deltam_2loop_MSSM.txt',usecols=[0,1])
x=A[:,0]
deltam=A[:,1]*1000
x_fit = x[x>100]
y_fit = deltam[x>100]
x_fit2  = x_fit[x_fit<4000]
y_fit2 = y_fit[x_fit<4000]
fx=np.polyfit(log(x_fit2),y_fit2,4)
f=np.poly1d(fx)
print(np.poly1d(f))




