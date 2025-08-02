from __future__ import division
import figureManager2 as fm2
import theory
import numpy as np
from itertools import product as iproduct
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rcParams

rcParams.update({'figure.autolayout': True})

N1 = 30000
N2 = 1000


def get_val(f, pA, pB):
    val = 0
    if (f(pA, pB) > 0.0):
        val += 2
    if pA > theory.pc:
        val += 1
    return val


'''
This figure shows the difference in co-existence space between structural and functional dependency
'''
p = np.linspace(0.5, 1, N1)
d2_q0 = fm2.get_pinf(q=0, r=1000, L=1000, depType=2, topology="grid")
value_matrix2 = fm2.get_interpolated_value_matrix(p, d2_q0)
p_inf_measured2 = fm2.reduce_value_matrix_to_function(p, value_matrix2)


sharp_pinf = theory.sharp_pinf


d3_q0 = fm2.get_pinf(q=0, r=1000, L=1000, depType=3, topology="grid")
value_matrix = fm2.get_interpolated_value_matrix(p, d3_q0)
p_inf_measured = fm2.reduce_value_matrix_to_function(p, value_matrix)

sharp_p1 = sharp_pinf(p_inf_measured, 0.48)
sharp_p2 = sharp_pinf(p_inf_measured2, 5 / 36)

f = lambda pA, pB: 1 - theory.pc / pB - sharp_p1(pA)
# f2 = lambda pA, pB: 1 - theory.pc / pB - p_inf_measured2(pA)
f2 = lambda pA, pB: 1 - theory.pc / pB - sharp_p2(pA)

p = np.linspace(0.5, 1, N2)
Z = np.zeros((len(p), len(p)))
for i, j in iproduct(range(len(p)), repeat=2):
    pA = p[i]
    pB = p[j]
    Z[i, j] = get_val(f, pA, pB)
p = np.linspace(0.5, 1, N2)
Z2 = np.zeros((len(p), len(p)))
for i, j in iproduct(range(len(p)), repeat=2):
    pA = p[i]
    pB = p[j]
    Z2[i, j] = get_val(f2, pA, pB)
plt.figure()
msk = lambda Z: np.ma.masked_array(data=Z, mask=Z == 0)
plt.subplot(121)
# plt.pcolormesh(p, p, msk(Z), cmap=cm.jet_r)
plt.imshow(msk(Z), extent=[p.min(), p.max(), p.min(), p.max()], interpolation="none", origin='lower')
# plt.axis([0.55, 1, 0.55, 0.7])
plt.axis("tight")
plt.xlabel("$p_A$")
plt.ylabel("$p_B$")
plt.title("Functional")

plt.subplot(122)
# plt.pcolormesh(p, p, msk(Z2), cmap=cm.jet_r)
plt.imshow(msk(Z2), extent=[p.min(), p.max(), p.min(), p.max()], interpolation="none", origin='lower')
# plt.axis([0.55, 1, 0.55, 0.7])
plt.axis("tight")
plt.xlabel("$p_A$")
plt.ylabel("$p_B$")
plt.title("Structural")

plt.show()
