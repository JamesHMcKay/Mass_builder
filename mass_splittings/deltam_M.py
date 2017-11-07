#!/usr/bin/python

# Mass Builder
 
# James McKay
# Sep 2016
 
#--- plot_MSSM.py ---
 
# Simple routine to plot output from mass splittings example


import numpy as np
import pylab
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy import special

from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
from matplotlib import rc

#rc('text', usetex=True)
#plt.rc('font', family='Computer Modern Roman',weight='normal')

fig=plt.figure()

ax = fig.add_subplot(1,1,1)


#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

A=np.genfromtxt('mass_splittings/data2/deltam_M.txt',usecols=[0,1,2,3])
x=A[:,0]
y1=A[:,1]*1000
y2=A[:,2]*1000
y3=A[:,3]*1000


plt.plot(x,y1,'-',color='red',label='1-loop') #

plt.plot(x,y2,'-',color='green',label='2-loop') #

plt.plot(x,y3,'--',color='green',label='2-loop FS') #

#plt.xlim([100,4000])
plt.ylim([150,180])
#plt.ylim([660,680])
leg = plt.legend(loc='upper left')
ax.set_xscale('log')

xlabel(r"Degenerate mass $M$ (GeV)",fontsize=18)
ylabel(r"$\Delta M$ (Mev)",fontsize=18)

plt.savefig("mass_splittings/figures2/mass_splittings.pdf")



# determine a polynomial fit through the 2-loop mass splitting

fig=plt.figure()

subplots_adjust(hspace=0.3,wspace=0.)

ax1 = fig.add_subplot(2,1,1)
#ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

fx=np.polyfit(log(x),y2,4)
f=np.poly1d(fx)

print(np.poly1d(f))

r=np.linspace(log(100),log(4000))

ax1.plot(x,y2,'-',color='black',linewidth=2,label="$\mathrm{This\ work}$")

gx = [-0.181509, 5.41948, -60.8831 , 305.383, -413.315]
g=np.poly1d(gx)
ax1.plot(exp(r),g(r),'--',color='blue',linewidth=2,label="$\mathrm{Ibe\ et\ al.\ Eq.\ (17)}$")


ax1.set_xlim([200,4000])
ax1.set_ylim([158,165])
ax1.set_xscale('log')

ax2 = fig.add_subplot(4,1,3)

gx = [-0.181509, 5.41948, -60.8831 , 305.383, -413.315]
g=np.poly1d(gx)

leg = ax1.legend(loc='lower right')
ax1.set_xscale('log')



ax2.fill_between(x, -0.02*ones(size(y2)), 0.02*ones(size(y2)),color='grey',alpha=0.4)


ax2.plot(x,((g(log(x))-y2)/y2)*100,'-',color='black',linewidth=2)

#ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

ax1.set_ylabel('$\Delta M = M_{\mathrm{pole}}^+ - M_{\mathrm{pole}}^0$ $(\mathrm{MeV})$',fontsize=12)
ax2.set_xlabel('$M^0_{\mathrm{pole}}$ $(\mathrm{GeV})$',fontsize=12)
ax2.set_ylabel('$\mathrm{Difference}$ $(\mathrm{\%}$)$',fontsize=12)



ax2.set_xlim([200,4000])
#ax2.set_ylim([-0.05,0.05])
ax2.set_xscale('log')

ax1.set_xticks([200,1000,4000])
ax1.set_xticklabels([])

ax2.set_xticks([200,1000,4000])
ax2.set_xticklabels(['$\mathrm{200}$','$\mathrm{1000}$','$\mathrm{4000}$'])

ax2.set_yticks([-0.05,-0.02,0.02,0.05])
#ax2.set_yticklabels([-0.3,-0.2,-0.1,0,0.1])


plt.savefig("mass_splittings/figures2/difference.pdf",bbox_inches='tight')




