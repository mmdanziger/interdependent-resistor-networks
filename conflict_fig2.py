from __future__ import division
import theory
import json
from sys import argv
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

fname = glob(argv[1])[0]
try:
    i = int(argv[2])
except:
    i = 0

d = json.load(open(fname))
correct_type =  lambda x : x if x!=1 else 2
L = d["L"]
if d['topology'] == 'ER':
    pinf = theory.P_inf_factory(d['k'])
elif d['topology'] == 'grid':
    pinf = theory.grid_pinf(np.linspace(0.5, 1, 2000), q=0, r=1000, L=1000, topology="grid", depType=correct_type(d['type']))
g = lambda x: pinf(x) / x if x > 0 else 0
sim_cascadeA = lambda run : list(map(lambda x: x["NA"]/(L*L), run["history"]))
sim_cascadeB = lambda run : list(map(lambda x: x["NB"]/(L*L), run["history"]))
theory_cascade = lambda run: theory.conflict_cascade(run["pA"],run["pB"],g,returnAll=True)

run = d["runs"][i]
theoryA,theoryB = theory_cascade(run)
plt.figure()
plt.subplot(121)
plt.plot(theoryA,label="theory")

plt.plot(sim_cascadeA(run),label="sim")
plt.ylabel("$N_A$")
plt.legend()
plt.subplot(122)
plt.plot(theoryB,label="theory")
plt.plot(sim_cascadeB(run),label="sim")
plt.ylabel("$N_B$")
plt.legend()
plt.suptitle("$p_A = %.2f , p_B = %.2f$"%(run["pA"],run["pB"]))