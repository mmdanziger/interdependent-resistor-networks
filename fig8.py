from __future__ import division,print_function
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
        matplotlib.rcParams.update({"font.size": "26"})
        matplotlib.rcParams.update({"legend.fontsize": "22"})
    from matplotlib import font_manager
    font_manager.USE_FONTCONFIG = True
import figureManager2 as fm2
from matplotlib import pyplot as plt

"""
This figure plots the difference between pc for structural and process based dependency,
over a range of r values
"""
r_range = range(50)

pc2_avg=[]
pc2_err=[]
pc3_avg=[]
pc3_err=[]
Nj2_avg=[]
Nj2_err=[]
Nj3_avg=[]
Nj3_err=[]
L=1000

for r in r_range:
    pc_2 = fm2.get_system_averages_by_key(fm2.get_systems(q=1,r=r,depType=2,L=L),['pc','Njump'],apply_rescale_p=True)
    pc2_avg.append(pc_2['pc']['mean'])
    pc2_err.append(pc_2['pc']['std'])
    Nj2_avg.append(pc_2['Njump']['mean']/(L*L))
    Nj2_err.append(pc_2['Njump']['std']/(L*L))


    pc_3 = fm2.get_system_averages_by_key(fm2.get_systems(q=1,r=r,depType=3,L=L),['pc','Njump'],apply_rescale_p=True)
    pc3_avg.append(pc_3['pc']['mean'])
    pc3_err.append(pc_3['pc']['std'])
    Nj3_avg.append(pc_3['Njump']['mean']/(L*L))
    Nj3_err.append(pc_3['Njump']['std']/(L*L))

#plot everything
plt.errorbar(r_range,pc3_avg,pc3_err, elinewidth=2, lw=3, color="green" , label=r"process-based")
plt.errorbar(r_range,pc2_avg,pc2_err, elinewidth=2, lw=3, color = "darkslateblue", label=r"structural")
plt.xlabel(r"$r$")
plt.ylabel(r"$p_c$")
plt.gcf().set_size_inches(6,6)
plt.legend(loc=4)
plt.show() if ipy else plt.savefig("/tmp/func_vs_process_pc.pdf")

plt.figure()
plt.errorbar(r_range,Nj3_avg,Nj3_err, elinewidth=2, lw=3, color="green" , label=r"process-based")
plt.errorbar(r_range,Nj2_avg,Nj2_err, elinewidth=2, lw=3, color = "darkslateblue", label=r"structural")
plt.xlabel(r"$r$")
plt.ylabel(r"Order parameter at $p_c$")
plt.gcf().set_size_inches(6,6)
plt.legend(loc=4)
plt.show() if ipy else plt.savefig("/tmp/func_vs_process_Njump.pdf")


