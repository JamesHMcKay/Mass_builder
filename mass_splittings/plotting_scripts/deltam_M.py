#!/usr/bin/python

# Mass Builder
 
# James McKay
# Sep 2016

# This routine makes Figure 4 in arXiv:1712.00968
# it compares the two-loop MSSM mass splitting
# with and without light quarks and fits
# a polynomial to each

# Required routines to run are
#  void Figures<T>::plot_M(Data data)
#  and an equivalent one without light quark masses (TODO)


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

rc('text', usetex=True)
plt.rc('font', family='Computer Modern Roman',weight='normal')

fig=plt.figure()

ax = fig.add_subplot(1,1,1)


#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

A=np.genfromtxt('mass_splittings/data2/deltam_M.txt',usecols=[0,1,2,3,4])
x=A[:,0]
y1=A[:,1]*1000
y2=A[:,2]*1000
x_2=A[:,3]
y2_2=A[:,4]*1000


A=np.genfromtxt('mass_splittings/data2/deltam_M_lightquarks.txt',usecols=[0,1,2,3,4])
xl=A[:,0]
y1l=A[:,1]*1000
y2l=A[:,2]*1000
x_2l=A[:,3]
y2_2l=A[:,4]*1000


# determine a polynomial fit through the 2-loop mass splitting

fig=plt.figure()

subplots_adjust(hspace=0.3,wspace=0.)

ax1 = fig.add_subplot(2,1,1)


r=np.linspace(log(100),log(4000))

ax1.plot(x,y2,'-',color='blue',linewidth=1.1,label="$\mathrm{1}\!-\!\mathrm{loop\ RGEs}$",zorder=10)
ax1.plot(x_2,y2_2,'-',color='red',linewidth=1.1,label="$\mathrm{2}\!-\!\mathrm{loop\ RGEs}$",zorder=10)

gx = [-0.181509, 5.41948, -60.8831 , 305.383, -413.315]
g=np.poly1d(gx)
ax1.plot(exp(r),g(r),'--',color='black',linewidth=1.1,label="$\mathrm{Ibe\ et\ al.\ Eq.\ (17)}$",zorder=100)


ax1.set_xlim([100,4000])
ax1.set_ylim([158,165])
ax1.set_xscale('log')

ax2 = fig.add_subplot(4,1,3)

gx = [-0.181509, 5.41948, -60.8831 , 305.383, -413.315]
g=np.poly1d(gx)

ax1.set_xscale('log')

leg = ax1.legend(loc='lower right',frameon=False,fontsize=16)



ax2.fill_between(x, -0.02*ones(size(y2)), 0.02*ones(size(y2)),color='grey',alpha=0.3)


ax2.plot(x,((y2-g(log(x)))/y2)*100,'-',color='blue',linewidth=1.1)
ax2.plot(x_2,((y2_2-g(log(x_2)))/y2_2)*100,'-',color='red',linewidth=1.1)


# plot the light quark lines
ax2.plot(xl,((y2l-g(log(xl)))/y2l)*100,'--',color='blue',linewidth=1.1)
ax2.plot(x_2l,((y2_2l-g(log(x_2l)))/y2_2l)*100,'--',color='red',linewidth=1.1)

ax2.plot(x_2,zeros(size(x_2)),'--',color='black',linewidth=1.1)


ax1.set_ylabel('$\Delta M = M_{\mathrm{pole}}^+ - M_{\mathrm{pole}}^0$ $(\mathrm{MeV})$',fontsize=12)
ax2.set_xlabel('$M^0_{\mathrm{pole}}$ $(\mathrm{GeV})$',fontsize=12)
ax2.set_ylabel('$\mathrm{Difference}$ $(\mathrm{\%}$)$',fontsize=12)

ax1.annotate(r'$\mathrm{Wino\ model}$', xy=(110,164), xytext=(110,164),fontsize=16)
ax1.annotate(r'$Q=163.3$ $\mathrm{GeV}$', xy=(110,163), xytext=(110,163),fontsize=16)


ax2.set_xlim([100,4000])
ax2.set_ylim([-0.05,0.05])
ax2.set_xscale('log')

ax1.set_xticks([100,200,1000,4000])
ax1.set_xticklabels([])

ax2.set_xticks([100,200,1000,4000])
ax2.set_xticklabels(['$\mathrm{100}$','$\mathrm{200}$','$\mathrm{1000}$','$\mathrm{4000}$'])

ax2.set_yticks([-0.05,-0.02,0.02,0.05])
#ax2.set_yticklabels([-0.3,-0.2,-0.1,0,0.1])


plt.savefig("mass_splittings/figures/difference.pdf",bbox_inches='tight')



x_fit = x[x>100]
y2_fit = y2[x>100]


x_fit2  = x_fit[x_fit<4000]
y2_fit2 = y2_fit[x_fit<4000]



fx=np.polyfit(log(x_fit2),y2_fit2,4)
f=np.poly1d(fx)

print '1-loop RGEs and negligible light quark masses'
print(np.poly1d(f))


x_2_fitl = x_2l[x_2l>100]
y2_2_fitl = y2_2l[x_2l>100]


x_2_fit2l  = x_2_fitl[x_2_fitl<4000]
y2_2_fit2l = y2_2_fitl[x_2_fitl<4000]



fx=np.polyfit(log(x_2_fit2l),y2_2_fit2l,4)
f=np.poly1d(fx)

print '2-loop RGEs and non-zero light quark masses'
print(np.poly1d(f))



