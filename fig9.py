from __future__ import division
import numpy as np
import sys
sys.path.append("/home/micha/phd")
ipy = False
try:
    ipy = __IPYTHON__
except NameError:
    pass
import numpy as np
import matplotlib

if not ipy:
#    matplotlib.use('Qt4Agg')
    matplotlib.use('PDF')
    matplotlib.rcParams.update({"figure.autolayout": "true"})
    try:
        matplotlib.rcParams.update({"usetex": "true"})
    except KeyError:
        matplotlib.rcParams.update({"text.usetex": "true"})
    from matplotlib import font_manager
    font_manager.USE_FONTCONFIG = True
import matplotlib.pyplot as plt
import figureManager2 as fm2
import figureManager as fm1
import theory
import json
'''
    This figure plots the cascade backbone and conductivity for r = 5 and r = inf.
'''
if False:
    rinfsystems = [529,595,1926,593,1927,571,593,284]
    r5systems = [1551,1576,1673,1394]
    rinflist,r5list=[[],[]]
    for systemid in rinfsystems:
        rinflist.append(fm2.get_system_by_id(systemid))
    json.dump(rinflist,"/tmp/rinf.json")
    for systemid in r5systems:
        r5list.append(fm2.get_system_by_id(systemid))
    json.dump(r5list,"/tmp/r5.json")


#d2_q0 = fm2.get_pinf(q=0, r=1000, L=1000, depType=2, topology="grid")
#d2_q1 = fm2.get_pinf(q=1, r=1000, L=1000, depType=2, topology="grid")
result_list = fm2.get_pinf(q=0, r=1000, L=1000, depType=3, topology="grid")
#d3_q1 = fm2.get_pinf(q=1, r=1000, L=1000, depType=3, topology="grid")
p_vec = np.linspace(0.5, 1, 3000)
L = 1000
S, N = zip(*sorted(zip(fm1.get_vecs(fm2.get_NS(q=1, r=0, L=L, type=3, topology="grid"), ['N', 'S']))))[0]  #careful, * operator transposes order
N = np.array(N) / (L * L)
SofN = fm2.fit_S_to_N(N, S)
value_matrix = fm2.get_interpolated_value_matrix(p_vec, result_list)
p_inf_measured = fm2.reduce_value_matrix_to_function(p_vec, value_matrix)
p_inf_interdep = theory.get_interdep_pinf(p_vec, p_inf_measured)

g = lambda p: p_inf_measured(p)/p
pc = 0.8045668668546102
pc = 0.804648786321
bcascade = theory.interdependent_cascade_simple(pc,g,True)
scascade = list(map(SofN,bcascade))
t = list(range(len(bcascade)))

### rinf
rinffname =  "/home/micha/rnet_results/RNetDynamics-type=3-L=1000-r-1000-k=4.0000-l=-1.0000-q=1.0000-126377105.json"
rinfplots = None# [-7,-5,-3,0,2]
plt.figure()
fm1.plotDynamicsFromJson(rinffname,plot_key="S",window=3,add_title=False,explicit_plots=rinfplots)
plt.plot(t[::5],scascade[::5],'Dk',fillstyle="none",mew=1.3)
plt.savefig("/tmp/CascadeSigmaL1000r1000q1.pdf")
plt.close()

plt.figure()
fm1.plotDynamicsFromJson(rinffname,plot_key="Ng",window=3,add_title=False,explicit_plots=rinfplots)
plt.plot(t[::5],bcascade[::5],'sk',fillstyle="none",mew=1.3)
plt.savefig("/tmp/CascadeBackboneL1000r1000q1.pdf")
plt.close()

sys.exit()
#r5fname =  "/storage/home/micha/rnet_results/RNetDynamics-type=3-L=1000-r-5-k=4.0000-l=-1.0000-q=1.0000-827944040.json"
#r5fname =  "/storage/home/micha/rnet_results/RNetDynamics-type=3-L=1000-r-5-k=4.0000-l=-1.0000-q=1.0000-248229026.json"
r5fname =  "/storage/home/micha/rnet_results/RNetDynamics-type=3-L=1000-r-5-k=4.0000-l=-1.0000-q=1.0000-317116022.json"
r5plots = [-7,-5,-3,0,2]
plt.figure()
fm1.plotDynamicsFromJson(r5fname,plot_key="S",window=7,add_title=False,explicit_plots=r5plots)
plt.savefig("/tmp/CascadeSigmaL1000r5q1.pdf")
plt.close()

plt.figure()
fm1.plotDynamicsFromJson(r5fname,plot_key="Ng",window=7,add_title=False,explicit_plots=r5plots)
plt.savefig("/tmp/CascadeBackboneL1000r5q1.pdf")
plt.close()
