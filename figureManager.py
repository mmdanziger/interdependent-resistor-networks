from __future__ import division,print_function
import requests
try:
    import simplejson as json
except ImportError:
    import json
from collections import defaultdict
from numpy import mean, std, array, zeros
import matplotlib.pyplot as plt
import numpy as np

fo = {}
fo[1] = {}
fo[1][2] = 990  # 8
fo[1][3] = 1000000


def getAllRuns(q, r, L, deptype):
    url = "http://physionet.ph.biu.ac.il/~micha/rnets/api/get_run.php"
    qstring = "?q=%.5f&r=%i&L=%i&type=%i" % (q, r, L, deptype)
    r = requests.request("GET", url + qstring)
    if len(r.text) > 0:
        return json.loads(r.text)


def getPcByR(q, L, deptype):
    url = "http://physionet.ph.biu.ac.il/~micha/rnets/api/get_pc_by_r.php"
    qstring = "?q=%.5f&L=%i&type=%i" % (q, L, deptype)
    r = requests.request("GET", url + qstring)
    if len(r.text) > 0:
        return json.loads(r.text)


def getDynamicsVector(q=1, L=1000, deptype=3):
    url = "http://physionet.ph.biu.ac.il/~micha/rnets/api/dynamics_vector.php?q=%.2f&L=%i" % (q, L)
    d = json.loads(requests.request("GET", url).text)
    return d


def rescale_p(p, q):
    return p * (1 - q + q * p)


def get_vecs(res_list, keylist=["p", "S"]):
    vecs = defaultdict(list)
    for row in res_list:
        for k in keylist:
            vecs[k].append(row[k])
    return [vecs[k] for k in keylist]


def get_pc(run):
    for i, row in enumerate(run[:-1]):
        if run[i + 1]['S'] > 0:
            break
    return (i, dict(p=row["p"], S=run[i + 1]["S"]))


def average_runs(q=1, L=1000, deptype=3, rlist=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]):
    meanpc = {}
    meandeltaS = {}
    results = {}
    for r in rlist:
        run_stats = defaultdict(list)
        d = getAllRuns(q, r, L, deptype)
        for run in d.values():
            i, rowCrit = get_pc(run)
            pc = rowCrit["p"]
            deltaS = rowCrit["S"]
            # print deltaS
            if r > fo[q][deptype]:
                p = [i["p"] - pc for i in run if i["p"] >= pc]
                S = [i["S"] - deltaS for i in run if i["p"] >= pc]
            else:
                p = [i["p"] for i in run if i["p"] >= pc]
                S = [i["S"] for i in run if i["p"] >= pc]
            run_stats["pc"].append(pc)
            run_stats["p"].append(p)
            run_stats["S"].append(S)
            run_stats["deltaS"].append(deltaS)
        meanpc[r] = mean(run_stats["pc"])
        meandeltaS[r] = mean(run_stats["deltaS"])
        x = meanpc[r] if r > fo[q][deptype] else 0
        y = meandeltaS[r] if r > fo[q][deptype] else 0
        plist = []
        Slist = []
        for p in run_stats["p"]:
            plist.append(array([p_ + x for p_ in p]))
        for S in run_stats["S"]:
            Slist.append(array([S_ + y for S_ in S]))
        results[r] = zip(plist, Slist)
    return results


def average_and_shift_runs(q=1, L=1000, deptype=3, rlist=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]):
    meanpc = {}
    meandeltaS = {}
    results = {}
    for r in rlist:
        run_stats = defaultdict(list)
        d = getAllRuns(q, r, L, deptype)
        for run in d.values():
            i, rowCrit = get_pc(run)
            pc = rowCrit["p"]
            deltaS = rowCrit["S"]
            # print deltaS
            if r > fo[q][deptype]:
                p = [i["p"] - pc for i in run if i["p"] >= pc]
                S = [i["S"] - deltaS for i in run if i["p"] >= pc]
            else:
                p = [i["p"] for i in run if i["p"] >= pc]
                S = [i["S"] for i in run if i["p"] >= pc]
            run_stats["pc"].append(pc)
            run_stats["p"].append(p)
            run_stats["S"].append(S)
            run_stats["deltaS"].append(deltaS)
        meanpc[r] = mean(run_stats["pc"])
        meandeltaS[r] = mean(run_stats["deltaS"])
        x = meanpc[r] if r > fo[q][deptype] else 0
        y = meandeltaS[r] if r > fo[q][deptype] else 0
        plist = []
        Slist = []
        for p in run_stats["p"]:
            plist.append(array([p_ + x for p_ in p]))
        for S in run_stats["S"]:
            Slist.append(array([S_ + y for S_ in S]))
        results[r] = zip(plist, Slist)
    return results


def plot_all_runs(q=1, L=1000, deptype=3, rlist=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]):
    d = average_runs(q, L, deptype, rlist)
    f = plt.figure()
    ax = f.add_subplot(111)
    allcols = ["blue", "green", "red", "cyan", "magenta", "orchid", "mediumpurple", "darkgoldenrod", "purple", "teal", "slateblue",
           "sienna"]
    np.random.shuffle(allcols)
    for r in d:
        try:
            color = next(ax._get_lines.color_cycle)
        except AttributeError:
            color = allcols[0]
            del allcols[0]
        labelled = False
        for p, S in d[r]:
            ax.plot(rescale_p(p, q), S, '.', color=color) if labelled else ax.plot(rescale_p(p, q), S, '.', color=color,
                                                                                   label="r=%i" % r)
            ax.plot(rescale_p(p, q), S, '-', alpha=0.5, color=color)
            labelled = True
    plt.legend(fontsize="xx-small", loc=2)
    plt.show()


def plot_pc_of_r(q=1, L=1000, deptypes=[2, 3], rescale=True):
    rvec = {}
    pc = {}
    pcstd = {}
    for deptype in deptypes:
        d = getPcByR(q, L, deptype)
        try:
            d.keys()
        except:
            print("No results obtained for type %i" % deptype)
            continue
        if rescale:
            rvec[deptype] = [r for r in sorted(d.keys(), key=lambda x: int(x))]
            pc[deptype] = [rescale_p(mean(d[r]["pc"]), q) for r in sorted(d.keys(), key=lambda x: int(x))]
            pcstd[deptype] = [std(rescale_p(array(d[r]["pc"]), q)) for r in sorted(d.keys(), key=lambda x: int(x))]
        else:
            rvec[deptype] = [r for r in sorted(d.keys(), key=lambda x: int(x))]
            pc[deptype] = [mean(d[r]["pc"]) for r in sorted(d.keys(), key=lambda x: int(x))]
            pcstd[deptype] = [std(array(d[r]["pc"])) for r in sorted(d.keys(), key=lambda x: int(x))]

    import matplotlib.pyplot as plt

    plt.figure()
    for t in deptypes:
        plt.errorbar(rvec[t], pc[t], pcstd[t], label="pc type %i" % t)
    plt.legend(fontsize="x-small")
    plt.show()


def plot_NOI_of_r(q=1, L=1000, deptypes=[2, 3], scatter=True):
    rvec = {}
    pc = {}
    pcstd = {}
    NOI = {}
    NOIstd = {}
    for deptype in deptypes:
        d = getPcByR(q, L, deptype)
        try:
            d.keys()
        except:
            print("No results obtained for type %i" % deptype)
            continue
        skeys = sorted(d.keys(), key=lambda x: int(x))
        rvec[deptype] = [r for r in skeys]
        pc[deptype] = [rescale_p(mean(d[r]["pc"]), q) for r in skeys]
        NOI[deptype] = [mean(d[r]["NOI"]) for r in skeys]
        NOIstd[deptype] = [std(d[r]["NOI"]) for r in skeys]
        pcstd[deptype] = [std(rescale_p(array(d[r]["pc"]), q)) for r in skeys]
    import matplotlib.pyplot as plt

    plt.figure()
    for t in deptypes:
        plt.errorbar(rvec[t], pc[t], pcstd[t], label="pc type %i" % t)
    plt.legend(fontsize="x-small")
    plt.twinx()

    for t in deptypes:
        if not scatter:
            plt.errorbar(rvec[t], NOI[t], NOIstd[t], label="NOI type %i" % t)
        else:
            for x in rvec[t]:
                for y in d[x]["NOI"]:
                    plt.scatter(x, y)
    plt.legend(fontsize="x-small", loc=0)
    plt.show()


def plotDynamicsFromJson(fname, window=3, plot_key="S", add_title=True, explicit_plots=None):
    x0 = json.load(open(fname))
    x = sorted(x0['runs'], key=lambda i: float(i['p']))
    for pc_ind, d in enumerate(x):
        if d['finalS'] > 0:
            break
    pc_ind -= 1
    if pc_ind <= 0:
        print("No pc data")
        return
    start_ind = pc_ind - window if pc_ind >= window else 0
    stop_ind = pc_ind + window if pc_ind + window <= len(x) else len(x)
    indices_to_plot = [pc_ind + i for i in explicit_plots] if explicit_plots else range(start_ind, stop_ind)
    norm = 1 if plot_key == "S" else x0['L'] * x0['L']
    for ind in indices_to_plot:
        plt.plot([i[plot_key] / norm for i in x[ind]["history"]], label="$p = %.5f$" % rescale_p(x[ind]["p"], x0['q']))
    plt.legend(ncol=2, fontsize="x-small")
    import os

    title_string = " ".join(os.path.basename(fname).strip(".json").split('-')[1:-1])
    ylabel_string = r"$\sigma$" if plot_key == "S" else r"$B_\infty$"
    xlabel_string = r"$t$"
    if add_title:
        plt.title(title_string)
    plt.xlabel(xlabel_string)
    plt.ylabel(ylabel_string)
    plt.show()
    plt.ion()
    # ofname = fname.replace(".json", "-window=%i.png" % window)
    #plt.savefig(ofname)
    #plt.close()


def get_damage_to_step(d, pvec, L, step, zeroNorm=True):
    Z = zeros((len(d), len(d[0])))
    zeronormvec = [i.values()[-1][step] if step < len(i.values()[-1]) else i.values()[-1][-1] for i in
                   sorted(d[0].values(), reverse=True)]
    for i, pvec in enumerate(d):  # iterate over r
        for j, this_p in enumerate(sorted(pvec, reverse=True)):
            N = pvec[this_p].values()[0]
            normval = zeronormvec[j] if zeroNorm else L * L
            Z[i, j] = N[step] if step < len(N) else N[-1]
            Z[i, j] /= normval
    return Z


def plotDynamicDamageArrays(q=1, L=1000, deptype=3, includeMax=True, stepsList=range(50), zeroNorm=True):
    d = getDynamicsVector(q, L, deptype)
    rvec = range(len(d))
    pvec = sorted(d[0].keys(), reverse=True)
    pvec = sorted(pvec, reverse=True)
    maxiter = max([max([len(i.values()[-1]) for i in j.values()]) for j in d])
    if includeMax:
        stepsList += [maxiter]
    spacing = 5
    for step in stepsList:
        plt.imshow(get_damage_to_step(d, pvec, L, step, zeroNorm).T, aspect=len(rvec) / len(pvec))
        plt.colorbar()
        ticks = plt.yticks(range(0, len(pvec), spacing) + [len(pvec) - 1],
                           ["%.3f" % rescale_p(float(i), q) for en, i in enumerate(sorted(pvec, reverse=True)) if
                            en % spacing == 0] + ["%.3f" % rescale_p(float(min(pvec)), q)])
        plt.xlabel("r")
        plt.ylabel("p0")
        plt.title("L=%i q=%.2f toStep %i" % (L, q, step))
        plt.savefig("/tmp/toStep_%i.png" % step)
        plt.close()


