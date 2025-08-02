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
    matplotlib.rcParams.update({"figure.autolayout": "true"})
    try:
        matplotlib.rcParams.update({"usetex": "true"})
    except KeyError:
        matplotlib.rcParams.update({"text.usetex": "true"})
    matplotlib.rcParams.update({"font.size": "24"})
    matplotlib.rcParams.update({"legend.fontsize": "19"})
    from matplotlib import font_manager

    font_manager.USE_FONTCONFIG = True
import matplotlib.pyplot as plt
import figureManager as fm1
import figureManager2 as fm2
import theory

r_list = [0, 1, 2, 3, 4, 5]#, 6, 7], 1000]
# r_list = [5]
q = 1
depType = 3
L = 1000
plotBackbone = False
p = np.linspace(0.5, 1, 3000)
plt.figure()
plt.ion()
y_key = "N" if plotBackbone else "S"
y_label = r"$B_\infty$" if plotBackbone else r"$\sigma$"
for r in r_list:
    fp = fm2.get_interpolated_curve(p,r=r,L=L,q=q,y_key=y_key,depType=depType)
    # systems = fm2.get_systems(q=q, r=r, L=L, depType=depType)
    # pinf = fm2.get_pinf(q=q, r=r, L=L, depType=depType)
    # pinf_filtered = fm2.filter_runs(pinf, systems)
    # value_matrix = fm2.get_interpolated_value_matrix(p, pinf_filtered, pSquared=True, y_key=y_key, extend=False)
    # f = fm2.reduce_value_matrix_to_function(p, value_matrix, ignoreZeros=True)
    # fp = np.nan_to_num(f(p))
    if r == 0:
        label_text = "single net"
    elif r < L/2 :
        label_text = "$r=%i$" % r
    else:
        label_text = r"$r=\infty$"

    plt.plot(p, fp, lw=2, label=label_text)
if max(r_list) > L/2:

    d3_q0 = fm2.get_pinf(q=0, r=1000, L=1000, depType=3, topology="grid")
    value_matrix = fm2.get_interpolated_value_matrix(p, d3_q0)
    p_inf_measured = fm2.reduce_value_matrix_to_function(p, value_matrix)
    p_inf_interdep = theory.get_interdep_pinf(p, p_inf_measured)
    if plotBackbone:
        s_of_n = lambda x : x
    else:
        S, N = zip(*sorted(zip(fm1.get_vecs(fm2.get_NS(q=1, r=0, L=L, type=3, topology="grid"), ['N', 'S']))))[0]  #careful, * operator transposes order
        N = np.array(N) / (L * L)
        s_of_n = fm2.fit_S_to_N(N, S)
    plt.plot(p, s_of_n(p_inf_interdep(p)), 'r-.', label="theory")

plt.axis([0.58,1,0,1])
plt.xlabel("$p$")
plt.ylabel(y_label)
plt.legend(loc=2, ncol=1,labelspacing=0.2)
plt.show() if ipy else plt.savefig("/tmp/fig4_%s.pdf"%y_key)
plt.close()