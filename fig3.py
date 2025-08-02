from __future__ import division

ipy = False
try:
    ipy = __IPYTHON__
except NameError:
    pass
import numpy as np
import matplotlib

if not ipy:
    matplotlib.use('Qt4Agg')
    matplotlib.rcParams.update({"usetex": "true"})
    from matplotlib import font_manager

    font_manager.USE_FONTCONFIG = True
import matplotlib.pyplot as plt
import figureManager as fm1
import figureManager2 as fm2

L = 1500
d = fm2.get_NS(q=1, r=0, L=L, type=3, topology="grid")
S, N = zip(*sorted(zip(fm1.get_vecs(d, ['N', 'S']))))[0]  #careful, * operator transposes order
N = np.array(N) / (L * L)
SofN = fm2.fit_S_to_N(N, S)
plt.figure()
plt.ion()
#plt.plot(N, S, ',',alpha=0.3)
n = np.linspace(0, 1, 1500)
plt.plot(n, SofN(n),lw=2)
plt.xlabel(r"$B_\infty$")
plt.ylabel(r"$\sigma$")
plt.axis([0, 1, 0, 1])
plt.show()
if not ipy:
    plt.savefig("BackboneSigma.pdf")
    plt.close()