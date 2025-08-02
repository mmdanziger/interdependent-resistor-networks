from __future__ import division
from pylab import *
from sys import argv
import simplejson as json
import os
from matplotlib import rcParams
from glob import glob
rcParams.update({'figure.autolayout': True})

fname = glob(argv[1])[0]
d = json.load(open(fname))

eps = 0.01
correct_type = lambda x: x if x!=1 else 2
msk = lambda Z: np.ma.masked_array(data=Z, mask=abs(Z) < eps).T
plot_theory = True

A = dict(((i['pA'], i['pB']), i['finalNa'] / (d['L'] * d['L'])) for i in d['runs'])
B = dict(((i['pA'], i['pB']), i['finalNb'] / (d['L'] * d['L'])) for i in d['runs'])

figure(figsize=(12, 9), dpi=400)
subplot(2, 2, 1) if plot_theory else subplot(1, 2, 1)
pA = np.array(sorted(list(set([i['pA'] for i in d['runs']]))))
pB = np.array(sorted(list(set([i['pB'] for i in d['runs']]))))

pB_ind = dict(((i, ind) for ind, i in enumerate(pB)))
pA_ind = dict(((i, ind) for ind, i in enumerate(pA)))

A_mat = np.zeros((len(pA), len(pB)))
B_mat = np.zeros((len(pA), len(pB)))

for [x, y], z in A.iteritems():
    A_mat[pA_ind[x], pB_ind[y]] = z
imshow(msk(A_mat), interpolation='none', extent=[pA.min(), pA.max(), pB.min(), pB.max()], origin='lower')
# if z>0:
# scatter(x,y,c=cm.jet(z),marker='s',lw=0,cmap=cm.jet)
# else:
# scatter(x,y,marker='s',lw=0,c='k')
xlabel("$p_A$")
ylabel("$p_B$")
# axis([0.5,1,0.5,1])
axis("tight")
title("A")
subplot(2, 2, 2) if plot_theory else subplot(1, 2, 2)
for [x, y], z in B.iteritems():
    B_mat[pA_ind[x], pB_ind[y]] = z
    # if z > 0:
    # scatter(x, y, c=cm.jet(z), marker='s', lw=0, cmap=cm.jet)
    # else:
    #     scatter(x, y, marker='s', lw=0, c='k')
imshow(msk(B_mat), interpolation='none', extent=[pA.min(), pA.max(), pB.min(), pB.max()], origin='lower')
xlabel("$p_A$")
ylabel("$p_B$")
# axis([0.5,1,0.5,1])
axis("tight")
title("B")
if plot_theory:
    import theory

    if d['topology'] == 'ER':
        pinf = theory.P_inf_factory(d['k'])
    elif d['topology'] == 'grid':
        pinf = theory.grid_pinf(linspace(0.5, 1, 2000), q=0, r=1000, L=1000, topology="grid", depType=correct_type(d['type']))
    g = lambda x: pinf(x) / x if x > 0 else 0
    if d['resolver'] == 'equiprobable':
        pwinA, pwinB = [0.5, 0.5]
    elif d['resolver'] == 'firstpriority':
        pwinA, pwinB = [0, 1]
    else:
        pwinA, pwinB = [0.5, 0.5]
    fA = lambda a, b: theory.conflict_cascade(a, b, g, pwinA=pwinA, pwinB=pwinB, qA=d['q'])[0]
    fB = lambda a, b: theory.conflict_cascade(a, b, g, pwinA=pwinA, pwinB=pwinB, qA=d['q'])[1]

    ZA = theory.fill_matrix(pA, pB, fA)
    ZB = theory.fill_matrix(pA, pB, fB)

    subplot(223)
    imshow(msk(ZA), interpolation='none', extent=[pA.min(), pA.max(), pB.min(), pB.max()], origin='lower')
    xlabel("$p_A$")
    ylabel("$p_B$")
    axis("tight")

    subplot(224)
    imshow(msk(ZB), interpolation='none', extent=[pA.min(), pA.max(), pB.min(), pB.max()], origin='lower')
    xlabel("$p_A$")
    ylabel("$p_B$")
    axis("tight")

show()
savefig("/tmp/%s"%os.path.basename(fname.strip(".json"))+".pdf")
