from __future__ import division,print_function

from scipy.optimize import curve_fit
ipy = False
try:
    ipy = __IPYTHON__
except NameError:
    pass
import numpy as np
import matplotlib

if not ipy:
    matplotlib.use('Qt4Agg')
    matplotlib.rcParams.update({"figure.autolayout": "true"})
    try:
        matplotlib.rcParams.update({"usetex": "true"})
    except KeyError:
        matplotlib.rcParams.update({"text.usetex": "true"})
    from matplotlib import font_manager
    font_manager.USE_FONTCONFIG = True
import matplotlib.pyplot as plt
from sys import path
path.append("/home/micha/phd")
from pyNoN import logbin,logfit
from scipy import odr
import figureManager2 as fm2
"""
This plots the critical scaling of backbone and conductivity for r=0,1,2,3
"""
def get_color(**kwargs):
    colors = ['r', 'g' , 'b', 'c']
    r = kwargs['r']
    return colors[r]
def direct_fit(x,y,exp0=0.5):
    odr_dat = odr.RealData(x,y)
    odr_mod = odr.models.Model(lambda B, x: B[1]* (x ** B[0]))
    odr_inst = odr.ODR(data=odr_dat,model=odr_mod,beta0=[exp0,1],maxit=1000)
    odr_out = odr_inst.run()
    odr_out.pprint()
    return odr_out

def powerlaw(beta,x):
    return beta[1] * np.array(x)**beta[0]

def get_crit(fitS=True,showFit=True,fit_min_above_pc=0,fit_max_above_pc=1,**kwargs):
    p,S,N = fm2.calculate_critical_exponent(**kwargs)
    if fitS:
        N=S
    #lb = logbin.LogBin(zip(p,N),resolution=0.04)
    #plt.loglog(lb.xavg,lb.yavg,'-s',color=get_color(**kwargs),label="$r = %i$"%kwargs['r'])
    fitp = [i for i in p if fit_min_above_pc <  i < fit_max_above_pc]
    fitN = [N[i] for i,pval in enumerate(p) if fit_min_above_pc<pval<fit_max_above_pc]
    odr_out = direct_fit(fitp,fitN)
    if showFit:
        plt.plot(fitp,powerlaw(odr_out.beta,fitp),lw=2,color='k')
    # fit = lambda p,A,beta : A* p**beta
    #[A,beta],res = curve_fit(fit,p,N,p0=(0.5,5.0/36))
    #print (A,beta,res)
    #plt.loglog(p,fit(p,A,beta),lw=3)
    return odr_out

def plot_points(fitS=True, **kwargs):
    p,S,N = fm2.calculate_critical_exponent(**kwargs)
    if fitS:
        plt.loglog(p,S,'o',alpha=0.8,ms=3.5,mew=0,color=get_color(**kwargs),label="$r = %i$"%kwargs["r"])

    else:
        plt.loglog(p,N,'o',alpha=0.8,ms=3.5,mew=0,color=get_color(**kwargs),label="$r = %i$"%kwargs["r"])

    return p,S,N


L=1500

plt.figure()
plt.ion()
for i in range(4):
    p,S,N = plot_points(fitS=True,avg_pc=True,q=1,r=i,L=L,depType=3,min_above_pc=2e-4,max_above_pc=0.05)
critS = []
critS.append(get_crit(q=1,r=3,L=L,depType=3,fitS=True,avg_pc=True,min_above_pc=2e-4,fit_min_above_pc=0.0004,fit_max_above_pc=0.0029))
critS.append(get_crit(q=1,r=2,L=L,depType=3,fitS=True,avg_pc=True,min_above_pc=2e-4,fit_min_above_pc=0.0005,fit_max_above_pc=0.00655))
critS.append(get_crit(q=1,r=1,L=L,depType=3,fitS=True,avg_pc=True,min_above_pc=2e-4,fit_min_above_pc=0.0005,fit_max_above_pc=0.0148))
critS.append(get_crit(q=1,r=0,L=L,depType=3,fitS=True,avg_pc=True,min_above_pc=2e-4,fit_min_above_pc=0.0005,fit_max_above_pc=0.0385))

#plt.axis('scaled')
plt.axis('tight')
plt.xlabel(r"$p - p_c$")
plt.ylabel(r"$\sigma$")
plt.legend(loc=4)
plt.tight_layout()
plt.show()
plt.savefig("/tmp/fig5a.pdf")

plt.figure()
plt.ion()
for i in range(4):
    plot_points(fitS=False,avg_pc=True,q=1,r=i,L=L,depType=3,min_above_pc=2e-4,max_above_pc=0.05)
critN = []
critN.append(get_crit(q=1,r=3,L=L,depType=3,fitS=False,avg_pc=True,min_above_pc=2e-4,fit_min_above_pc=3e-4,fit_max_above_pc=0.0029))
critN.append(get_crit(q=1,r=2,L=L,depType=3,fitS=False,avg_pc=True,min_above_pc=2e-4,fit_min_above_pc=0.0003,fit_max_above_pc=0.00467))
critN.append(get_crit(q=1,r=1,L=L,depType=3,fitS=False,avg_pc=True,min_above_pc=2e-4,fit_min_above_pc=0.0003,fit_max_above_pc=0.00851))
critN.append(get_crit(q=1,r=0,L=L,depType=3,fitS=False,avg_pc=True,min_above_pc=2e-4,fit_min_above_pc=0.0003,fit_max_above_pc=0.027))

#plt.axis('scaled')
plt.xlabel(r"$p - p_c$")
plt.ylabel(r"$B_\infty$")
plt.axis('tight')
plt.legend(loc=4)
plt.tight_layout()
plt.show()
plt.savefig("/tmp/fig5b.pdf")
for i,(S_,N_) in enumerate(reversed(zip(critS,critN))):
    print("%i : %.6f +/- %.3f (%i samples) \t %.6f +/- %.3f (%i samples)"%
        (i,S_.beta[0], S_.sd_beta[0], len(S_.y),
            N_.beta[0], N_.sd_beta[0], len(N_.y)))

