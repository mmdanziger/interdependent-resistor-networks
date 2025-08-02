from __future__ import division
import figureManager2 as fm2
import figureManager as fm1
import theory
from scipy.optimize import curve_fit
from scipy.stats import linregress
import matplotlib.pyplot as plt
import numpy as np

'''This figure plots the scaling above p_c for the r=\infty case, for backbone and sigma, using measured or fit data.'''
L = 1000
q = 1
r = 1000
bb_from = 0
bb_to = 50
s_from = 0
s_to = 50

theoretical = False
def get_pNS(res_list,pc,Nc,Sc):
    return map(lambda x:(x['p']-pc,x['N']-Nc,x['S']-Sc), res_list)

def logfit(x,y):
    from scipy.stats import linregress
    x = np.log10(x)
    y = np.log10(y)
    beta, A, r, p_test, err = linregress(x,y)
    print("Beta = %.5f"%beta)
    return lambda x: 10**A*x**beta

p = np.linspace(0.7, 0.85, 30000)

'''experimental data'''
if not theoretical:
    systems = fm2.get_systems(q=q, r=r, L=L, depType=3)
    pinf = fm2.get_pinf(q=q, r=r, L=L, depType=3)
    pinf_filtered = fm2.filter_runs(pinf, systems)
    value_matrix = fm2.get_interpolated_value_matrix(p, pinf_filtered, pSquared=True, y_key="N", extend=False)
    f = fm2.reduce_value_matrix_to_function(p, value_matrix, ignoreZeros=True)
    fp = np.nan_to_num(f(p))
    pc, Nc = [(i, j) for i, j in zip(p, fp) if j > 0][0]


    value_matrix = fm2.get_interpolated_value_matrix(p, pinf_filtered, pSquared=True, y_key="S", extend=False)
    f = fm2.reduce_value_matrix_to_function(p, value_matrix, ignoreZeros=True)
    fp = np.nan_to_num(f(p))

    pc, Sc = [(i, j) for i, j in zip(p, fp) if j > 0][0]
    xs0 = p[p > pc] - pc
    ys0 = fp[fp > Sc] - Sc

    pNs=[]
    for res in pinf_filtered.itervalues():
        pNs.extend(get_pNS(res,pc,Nc,Sc))
    p,N,S = np.array(zip(*sorted(pNs)))
    x0 = p
    y0 = N
    ys0 = S
'''theoretical curve'''
if theoretical:
    d3_q0 = fm2.get_pinf(q=0, r=1000, L=1000, depType=3, topology="grid")
    value_matrix = fm2.get_interpolated_value_matrix(p, d3_q0)
    p_inf_measured = fm2.reduce_value_matrix_to_function(p, value_matrix)
    p_inf_interdep = theory.get_interdep_pinf(p, p_inf_measured)
    pi = p_inf_interdep(p)
    pc, pic = [(i, j) for i, j in zip(p, pi) if j > 0][0]
    x0 = p[p > pc] - pc
    y0 = pi[pi > pic] - pic
    S, N = zip(*sorted(zip(fm1.get_vecs(fm2.get_NS(q=1, r=0, L=L, type=3, topology="grid"), ['N', 'S']))))[
    0]  # careful, * operator transposes order
    N = np.array(N) / (L * L)
    s_of_n = fm2.fit_S_to_N(N, S)
    s = s_of_n(p_inf_interdep(p))
    pc, sc = [(i, j) for i, j in zip(p, s) if j > 0][0]
    xs0 = p[p > pc] - pc
    ys0 = s[s > sc] - sc

plt.figure()
plt.ion()
plt.loglog(x0, y0, '.')
f = lambda x, A, beta: A * x ** beta
x = x0[bb_from:bb_to]
y = y0[bb_from:bb_to]
[A, beta], cov = curve_fit(f, x, y)
beta, A, r, p_test, err = linregress(np.log10(x), np.log10(y))
f2 = lambda x, A, beta: (10 ** A) * x ** beta
plt.loglog(x, f2(x, A, beta), label="$\\tilde{P}_\infty$ : $\\beta = %.4f$" % beta)


xs = x0[s_from:s_to]
ys = ys0[s_from:s_to]
[As, betas], cov = curve_fit(f, x, ys)
betas, As, r, p_test, err = linregress(np.log10(x), np.log10(ys))
plt.loglog(xs0, ys0, '.')
plt.loglog(xs, f2(xs, As, betas), label="$\sigma$ : $\\beta = %.4f$" % betas)
plt.show()