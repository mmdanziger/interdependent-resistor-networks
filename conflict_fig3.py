from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import matplotlib
from sys import argv
from glob import glob
import json
import itertools
matplotlib.rcParams.update({'font.size': 16})


#color constants
sat=1
val=0.9
hmin=0
hmax=0.9
eps = 1e-9

msk3 = lambda Z: np.ma.masked_array(data=Z, mask=abs(Z) < eps).T

fname = glob(argv[1])[0]
pA = float(argv[2])
pB = float(argv[3])


to_check = filter(lambda x: x["pA"]==pA and x["pB"]==pB,
                  json.load(open(fname)))[0]
arrlist = to_check["history"]

#dig through to get unique vals to make universal discrete colormap
merged = list(itertools.chain(*arrlist))
merged = list(itertools.chain(*merged))
merged = list(itertools.chain(*merged))
total_unique_vals = sorted(list(set( merged )))
del merged

#create hsv optimized discrete colormap
color_list = ['black']+[clrs.hsv_to_rgb([h,sat,val]) for h in np.linspace(hmin,hmax,num=(len(total_unique_vals)))]
cmap = clrs.ListedColormap(color_list)
bounds = total_unique_vals
norm = clrs.BoundaryNorm([i+eps for i in bounds], cmap.N)

for idx,[A,B] in enumerate(to_check['history']):
    fig = plt.figure(figsize=(19,10),dpi=400)
    ax1=fig.add_axes([0.1,0.1,0.3,0.9])
    ax2=fig.add_axes([0.5,0.1,0.3,0.9])
    cax=fig.add_axes([0.85,0.33,0.05,0.42])
    A=np.array(A); B=np.array(B);
    im=ax1.imshow(msk3(A),interpolation="none", origin="lower", cmap=cmap, norm=norm)
    ax2.imshow(msk3(B),interpolation="none", origin="lower", cmap=cmap, norm=norm)
    plt.suptitle("%i"%idx)
    cb = fig.colorbar(im, cax=cax,boundaries=bounds,cmap=cmap,norm=norm)
    cb.set_ticks(np.array(bounds)[:-1] + 0.5*np.diff(bounds))
    cb.set_ticklabels(total_unique_vals[1:])
    plt.show()
    plt.savefig("/tmp/competing_graph_%.4f_%.4f_%i.png"%(pA,pB,idx))
    plt.close()
