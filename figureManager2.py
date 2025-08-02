from __future__ import division,print_function
import figureManager as fm1
import requests
try:
    import simplejson as json
except ImportError:
    import json
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from  scipy import interp
import theory


def get_pinf(q=1, r=0, L=1000, k=4.0, lam=0, depType=3, topology="grid"):
    url = "http://physionet.ph.biu.ac.il/~micha/rnets/api/v2/get_pinf.php"
    query_string = {"grid": "?q=%.5f&r=%i&L=%i&type=%i&topology=grid" % (q, r, L, depType,),
                    "ER": "?q=%.5f&L=%i&k=%.8f&type=%i&topology=ER" % (q, L, k, depType),
                    "lambda": "?q=%.5f&L=%i&k=%.8f&lambda=%.8f&type=%i&topology=lambda" % (q, L, k, lam, depType)
    }
    r = requests.request("GET", url + query_string[topology])
    if len(r.text) > 0:
        return json.loads(r.text)


def get_systems(q=1, r=0, L=1000, k=4.0, lam=0, depType=3, topology="grid"):
    url = "http://physionet.ph.biu.ac.il/~micha/rnets/api/v2/get_systems.php"
    query_string = {"grid": "?q=%.5f&r=%i&L=%i&type=%i&topology=grid" % (q, r, L, depType,),
                    "ER": "?q=%.5f&L=%i&k=%.8f&type=%i&topology=ER" % (q, L, k, depType),
                    "lambda": "?q=%.5f&L=%i&k=%.8f&lambda=%.8f&type=%i&topology=lambda" % (q, L, k, lam, depType)
    }
    r = requests.request("GET", url + query_string[topology])
    if len(r.text) > 0:
        return json.loads(r.text)


def get_NS(**kwargs):  # q=1, r=None, L=1000, k=None, lam=None, depType=3, topology="grid"):
    url = "http://physionet.ph.biu.ac.il/~micha/rnets/api/v2/get_NS.php"
    # params = {"q" : q , "r" : r , "L" : L, "k": k, "lambda": lam , "type": depType, "topology":topology}
    # query_string = {"grid": "?q=%.5f&r=%i&L=%i&type=%i&topology=grid" % (q, r, L, depType,),
    # "ER": "?q=%.5f&L=%i&k=%.8f&type=%i&topology=ER" % (q, L, k, depType),
    # "lambda": "?q=%.5f&L=%i&k=%.8f&lambda=%.8f&type=%i&topology=lambda" % (q, L, k, lam, depType)
    # }
    r = requests.request("GET", url, params=kwargs)
    if len(r.text) > 0:
        return json.loads(r.text)


'''
Returns a function which linearly interpolates p and y_key based on the data in result_list.
Values for 0 and 1 are added automatically to make all functions correctly cover the entire domain.
'''


def interpolate_pinf(result_list, y_key="N", extend=True, pSquared=False):
    p, y = fm1.get_vecs(result_list, ["p", y_key])
    if pSquared:
        p = map(lambda x: x * x, p)
    if extend:
        return theory.interp_func([0] + p + [1], [0] + y + [1])
    else:
        return theory.interp_func(p + [1], y + [1], left=0, right=0)


def get_interpolated_value_matrix(x_values, results_dict, y_key="N", extend=True, pSquared=False):
    function_list = [interpolate_pinf(res, y_key=y_key, pSquared=pSquared, extend=extend) for res in
                     sorted(results_dict.values())]
    return np.array([f(x_values) for f in function_list])


def reduce_value_matrix_to_function(x_values, value_matrix, ignoreZeros=False):
    if ignoreZeros:
        return theory.interp_func(x_values, np.sum(value_matrix, axis=0) / np.sum(value_matrix > 0, axis=0), left=0,
                                  right=1)
    else:
        return theory.interp_func(x_values, np.mean(value_matrix, axis=0))


def plot_theory(p_vec, result_list, color='blue', s_of_n=None, interdepOnly=False):
    value_matrix = get_interpolated_value_matrix(p_vec, result_list)
    p_inf_measured = reduce_value_matrix_to_function(p_vec, value_matrix)
    p_inf_interdep = theory.get_interdep_pinf(p_vec, p_inf_measured)
    if not s_of_n:
        s_of_n = lambda x: x
    if not interdepOnly:
        plt.plot(p_vec, s_of_n(p_inf_measured(p_vec)), '-.', color=color, label="single")
    plt.plot(p_vec, s_of_n(p_inf_interdep(p_vec)), color=color, label="interdependent")
    plt.xlim(xmax=1)
    plt.show()
    return p_inf_measured,p_inf_interdep


def fit_S_to_N(N, S):
    from scipy.optimize import curve_fit

    # p5 = lambda x, p0, p1, p2, p3: p0 * x ** 5 + p1 * x ** 4 + p2 * x ** 3 + p3 * x ** 2 + (1.001 - p0 - p1 - p2 - p3) * x
    p6 = lambda x, p0, p1, p2, p3, p4: p0 * x ** 6 + p1 * x ** 5 + p2 * x ** 4 + p3 * x ** 3 + p4 * x ** 2 + ( 1.001 - p0 - p1 - p2 - p3 - p4) * x
    # p7 = lambda x, p0, p1, p2, p3, p4, p5: p0 * x ** 7 + p1 * x ** 6 + p2 * x ** 5 + p3 * x ** 4 + p4 * x ** 3 + p5 * x ** 2 + (1.001 - p0 - p1 - p2 - p3 - p4 - p5) * x
    # p8 = lambda x, p0, p1, p2, p3, p4, p5, p6: p0 * x ** 8 + p1 * x ** 7 + p2 * x ** 6 + p3 * x ** 5 + p4 * x ** 4 + p5 * x ** 3 + p6 * x ** 2 + (1.001 - p0 - p1 - p2 - p3 - p4 - p5 - p6) * x
    p = p6
    popt, pcov = curve_fit(p, N, S)
    print("Coefficients are",end='\t')
    print(popt)
    print("Norm : %.4f" % np.linalg.norm(popt))
    print("With variance",end='\t')
    print(np.diag(pcov))
    print("Norm : %.4f" % np.linalg.norm(np.diag(pcov)))
    return lambda x: p(x, *popt)


def get_system_averages_by_key(systems_dict, key_list=['pc', 'Njump', 'Sjump'], apply_rescale_p=False):
    '''
    @systems_dict : output from get_systems in this file
    '''
    stat_dict = defaultdict(list)
    for system in systems_dict.values():
        for key in key_list:
            if key == 'pc' and apply_rescale_p:
                stat_dict[key].append(fm1.rescale_p(system[key], system['q']))
            else:
                stat_dict[key].append(system[key])
    try:
        return dict((key, {'mean': np.mean(val), 'std': np.std(val), 'median': np.median(val), 'max': max(
            val), 'min': min(val)} ) for key, val in stat_dict.iteritems())
    except TypeError:
        print("Hmmm. Did you forget to run insert_derived_critical_values on physionet?")


def filter_runs(pinf_dict, systems_dict):
    stat_dict = get_system_averages_by_key(systems_dict)
    return dict((id,
                 filter(lambda sim: sim['p'] > stat_dict['pc']['mean'] and sim['N'] > 0, pinf)) for id, pinf in
                pinf_dict.iteritems())


def calculate_critical_exponent(q=1, r=0, L=1000, depType=3, min_above_pc=1e-4, max_above_pc=0.05, force_pc=None,
                                avg_pc=False, avg_Sc=False, avg_Nc=False):
    systems = get_systems(q=q, r=r, L=L, depType=depType)
    sys_stats = get_system_averages_by_key(systems)
    if not systems:
        print("No systems found. Try different parameters.")
        return None, None
    pinf = get_pinf(q=q, r=r, L=L, depType=depType)
    p, S, N = [[], [], []]
    if avg_pc:
        avg_pc = np.mean([i['pc'] for i in systems.values()])
    # print avg_pc

    Sc = sys_stats["Sjump"]["mean"] if avg_Sc else 0
    Nc = sys_stats["Njump"]["mean"] / (L * L) if avg_Nc else 0
    for sid, sim in pinf.iteritems():
        p_, S_, N_ = fm1.get_vecs(sim, ['p', 'S', 'N'])
        if force_pc:
            pc = force_pc
        elif avg_pc:
            pc = avg_pc
        else:
            pc = systems[sid]['pc']
        p.extend([val - pc for i, val in enumerate(p_) if (pc + min_above_pc) < val < (pc + max_above_pc)])
        S.extend([val - Sc for i, val in enumerate(S_) if (pc + min_above_pc) < p_[i] < (pc + max_above_pc)])
        N.extend([val - Nc for i, val in enumerate(N_) if (pc + min_above_pc) < p_[i] < (pc + max_above_pc)])
    return zip(*sorted(filter(lambda x: all(map(lambda y: y > 0, x)), zip(p, S, N))))


def get_pinf_from_json(fname):
    d = json.load(open(fname))
    L = d["L"]
    N = L * L
    p, Nga, Ngb, Sa, Sb = [[] for i in range(5)]
    for run in sorted(d["runs"], key=lambda x: x['p']):
        p.append(run['p'])
        even = len(run['history']) % 2 == 0
        Nga.append(run['history'][-1]['Ng'] / N if not even else run['history'][-2]['Ng'] / N)
        Ngb.append(run['history'][-1]['Ng'] / N if even else run['history'][-2]['Ng'] / N)
        if d['type'] == 3:
            Sa.append(run['history'][-1]['S'] if not even else run['history'][-2]['S'])
            Sb.append(run['history'][-1]['S'] if even else run['history'][-2]['S'])
    return p, Nga, Ngb, Sa, Sb


def get_interpolated_curve(p, r=1, L=1500, q=1, y_key="S", depType=3):
    systems = get_systems(q=q, r=r, L=L, depType=depType)
    pinf = get_pinf(q=q, r=r, L=L, depType=depType)
    pinf_filtered = filter_runs(pinf, systems)
    value_matrix = get_interpolated_value_matrix(p, pinf_filtered, pSquared=True, y_key=y_key, extend=False)
    f = reduce_value_matrix_to_function(p, value_matrix, ignoreZeros=True)
    fp = np.nan_to_num(f(p))
    return fp


def bin_data(x, y, spacing='log', bin_count=None):
    from collections import Counter

    if not bin_count:
        bin_count = len(x) // 10
    if spacing == 'log':
        epsilon = 0.01 * np.log10(max(x)) / np.log10(min(x))
        bin_x = np.logspace(np.log10(min(x)) - epsilon, np.log10(max(x)) + epsilon, bin_count)
    elif spacing == 'linear':
        epsilon = 0.01 * (max(x) - min(x))
        bin_x = np.linspace(min(x) - epsilon, max(x) + epsilon, bin_count)
    else:
        raise ValueError("Please select 'linear' or 'log' bin spacing, not '" + str(spacing) + "'")
    bin_centers = ((np.pad(bin_x, (0, 1), mode="constant") + np.pad(bin_x, (1, 0), mode="constant") ) / 2)[1:-1]
    bin_centers = np.concatenate((np.asarray([bin_x[0]]), bin_centers[:-1]))
    bin_y = bin_x * np.nan
    bin_inds = np.digitize(x, bin_x)
    norm = dict(Counter(bin_inds))
    for b_ind, y_point in zip(bin_inds, y):
        if np.isnan(bin_y[b_ind]):
            bin_y[b_ind] = y_point
        else:
            bin_y[b_ind] += y_point
    bin_y = map(lambda x: x[1] / norm.get(x[0], 1), enumerate(bin_y))

    return np.asarray(zip(*filter(lambda x: not np.any(np.isnan(x)), zip(bin_centers, bin_y))))


def get_system_by_id(systemid):
    url = "http://physionet.ph.biu.ac.il/~micha/rnets/api/v2/get_dynamics.php"
    r = requests.request("GET",url,params={"systemID" :  systemid})
    try:
        return json.loads(r.text)
    except json.JSONDecodeError:
        print(r.text)
        raise ValueError("Something went wrong with the json.  probably a server error.")

