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
    #matplotlib.rcParams.update({"usetex": "true"})
    from matplotlib import font_manager

    font_manager.USE_FONTCONFIG = True
import matplotlib.pyplot as plt
import figureManager as fm1
import figureManager2 as fm2

L = 1500
rl = [0, 1, 2, 3, 4, 5, 6]  # ,1000]
plt.figure()
plt.ion()
fm1.plot_all_runs(L=L, rlist=rl)
p = np.linspace(0.5, 1, 3000)

S, N = list(zip(*sorted(zip(fm1.get_vecs(fm2.get_NS(q=1, r=0, L=L, type=3, topology="grid"), ['N', 'S'])))))[0]  # careful, * operator transposes order
N = np.array(N) / (L * L)
SofN = fm2.fit_S_to_N(N, S)
d3_q0 = fm2.get_pinf(q=0, r=1000, L=L, depType=3, topology="grid")
if d3_q0:
    fm2.plot_theory(p, d3_q0, 'red', SofN)
plt.axis([0.58, 1, 0, 1])

