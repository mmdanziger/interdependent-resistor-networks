from __future__ import division
ipy=False
try:
    ipy=__IPYTHON__
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

import matplotlib.pyplot as plt
import figureManager2 as fm2
import figureManager as fm1

def get_idx_in_sorted_list(mylist,val,offset=0):
    '''

    :param mylist: sorted! list - only looks at first idx this version will break on nonnested lists
    :param val: value abovewhich to return value
    :return: idx to list
    '''
    for idx,listval in enumerate(mylist[offset:]):
        if listval[0] > val:
            return idx+offset

if __name__ == "__main__":

    # p = np.linspace(0.5, 1, 3000)
    # L = 1000
    # St2q0 = fm2.get_interpolated_curve(p,r=1000,L=L,q=0,y_key="S",depType=2)
    # St2q1 = fm2.get_interpolated_curve(p,r=1000,L=L,q=1,y_key="S",depType=2)
    # St3q0 = fm2.get_interpolated_curve(p,r=1000,L=L,q=0,y_key="S",depType=3)
    # St3q1 = fm2.get_interpolated_curve(p,r=1000,L=L,q=1,y_key="S",depType=3)
    # Nt2q0 = fm2.get_interpolated_curve(p,r=1000,L=L,q=0,y_key="N",depType=2)
    # Nt2q1 = fm2.get_interpolated_curve(p,r=1000,L=L,q=1,y_key="N",depType=2)
    # Nt3q0 = fm2.get_interpolated_curve(p,r=1000,L=L,q=0,y_key="N",depType=3)
    # Nt3q1 = fm2.get_interpolated_curve(p,r=1000,L=L,q=1,y_key="N",depType=3)
    #
    # plt.figure()
    # plt.plot(p,St2q0)
    # plt.plot(p,Nt2q0)
    # plt.plot(p,Nt2q0)



    d2_q0 = fm2.get_pinf(q=0, r=1000, L=1000, depType=2, topology="grid")
    d2_q1 = fm2.get_pinf(q=1, r=1000, L=1000, depType=2, topology="grid")
    d3_q0 = fm2.get_pinf(q=0, r=1000, L=1000, depType=3, topology="grid")
    d3_q1 = fm2.get_pinf(q=1, r=1000, L=1000, depType=3, topology="grid")
    p = np.linspace(0.5, 1, 3000)
    L = 1000
    S, N = list(zip(*sorted(zip(fm1.get_vecs(fm2.get_NS(q=1, r=0, L=L, type=3, topology="grid"), ['N', 'S'])))))[0]  #careful, * operator transposes order
    N = np.array(N) / (L * L)
    SofN = fm2.fit_S_to_N(N, S)


    plt.figure()
    plt.ion()
    #Smeasured, Sinterdep = fm2.plot_theory(p, d3_q0, 'red', SofN)
    Nmeasured, Ninterdep = fm2.plot_theory(p, d3_q0, 'blue')
    N2measured, N2interdep = fm2.plot_theory(p, d2_q0, 'green')

    lw=3
    dblue,mgreen,lred=["#000033","#0f800f","red"]
    plt.figure()
    plt.plot(p,N2measured(p),"--",lw=lw,label=r"$P_\infty$",color=dblue,)
    plt.plot(p,Nmeasured(p),lw=lw,label=r"$B_\infty$", color=mgreen)
    plt.xlabel("$p$")
    plt.ylabel("Node fraction")
    plt.plot(p,SofN(Nmeasured(p)),"-.",lw=lw,label=r"$\sigma$", color=lred)
    plt.legend(loc=2)
    plt.axis([0.5,1,0,1])
    plt.twinx()
    plt.ylabel("$\sigma$")
    plt.axis([0.5,1,0,1])
    plt.show() if ipy else plt.savefig("/tmp/PBSigmaq0.pdf")

    plt.figure()
    plt.plot(p,N2interdep(p),"--",lw=lw,label=r"$P_\infty$",color=dblue)
    plt.plot(p,Ninterdep(p),lw=lw,label=r"$B_\infty$", color=mgreen)
    plt.xlabel("$p$")
    plt.ylabel("Node fraction")
    plt.plot(p,SofN(Ninterdep(p)),"-.",lw=lw,label=r"$\sigma$",color=lred)
    allp,allN,allS = [[],[],[]]
    for res in d3_q1.values():
        p_, N_, S_ = fm1.get_vecs(res, ['p', 'N', 'S'])
        allp.extend(p_)
        allN.extend(N_)
        allS.extend(S_)
    p_desired = np.arange(0.5,1,0.001)
    all_res = zip(allp,allN,allS)
    last_ind=0
    idx_to_plot=[]
    for pd in p_desired:
        idx=get_idx_in_sorted_list(all_res,pd,last_ind)
        last_ind=idx
        idx_to_plot.append(idx)
    idx_to_plot = sorted(list(set(idx_to_plot)))
    p_plot,N_plot,S_plot = zip(*[all_res[i] for i in idx_to_plot])
    plt.plot(np.array(p_plot) ** 2, N_plot, 'kD', fillstyle="none", ms=7, mew=2)
    plt.plot(np.array(p_plot) ** 2, S_plot, 'ks',fillstyle="none", ms=7, mew=2)
    plt.legend(loc=2)
    plt.axis([0.5,1,0,1])
    plt.twinx()
    plt.ylabel("$\sigma$")
    plt.axis([0.5,1,0,1])
    plt.show() if ipy else plt.savefig("/tmp/PBSigmaq1.pdf")





    #plt.figure()
    #plt.plot(p, p, 'yellow')
    # for res in d2_q1.values():
    #     p_, N_ = fm1.get_vecs(res, ['p', 'N'])
    #     plt.plot(np.array(p_) ** 2, N_, ',')
    # for res in d3_q1.values():
    #     p_, N_, S_ = fm1.get_vecs(res, ['p', 'N', 'S'])
    #     plt.plot(np.array(p_) ** 2, N_, ',')
    #     plt.plot(np.array(p_) ** 2, S_, ',')
    #
    # if ipy:
    #     plt.show()
    # else:
    #     plt.savefig("/tmp/StructFunctTheoryAndResults.pdf")
    #     plt.close()
