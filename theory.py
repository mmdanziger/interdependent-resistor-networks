from __future__ import division, print_function
import numpy as np
from scipy.optimize import fsolve
from scipy import interp
import scipy

verbose = False
pc = 0.5927
num_eps = 1e-5
'''
Dummy percolation function, looks kind of like real function to test methods
'''
theta = lambda x: 0.5 * ((np.abs(x) / x) + 1)
dummy = lambda p: np.sqrt(theta(p - pc) * (p - pc)) / np.sqrt(1 - pc) / p  # if p>pc else 0


def P_inf_factory(k):
    '''
    ER network percolation function returns the function, not the value
    :param k: average degree
    :return: percolation profile function
    '''
    from scipy.optimize import fsolve

    def f(p):
        return fsolve(lambda x: x - p * (1 - np.exp(-k * x)), 0.5, full_output=False)

    return f


def grid_pinf(p, sharpen=True, **kwargs):
    '''
    Downloads empirical percolation profile for lattices
    :param p: np array upon which interpolation is carried out
    :param sharpen: should the function be sharpened to correct critical exponent at pc?
    :param kwargs: arguments passed to fm2.get_pinf
    :return: p_inf function
    '''
    import figureManager2 as fm2

    d = fm2.get_pinf(**kwargs)
    value_matrix = fm2.get_interpolated_value_matrix(p, d)
    p_inf_measured = fm2.reduce_value_matrix_to_function(p, value_matrix)
    if sharpen:
        beta = 5 / 36 if kwargs['depType'] == 2 else 0.48
        return sharp_pinf(p_inf_measured, beta)
    else:
        return p_inf_measured


def interp_func(xp, yp, **kwargs):
    # return scipy.interpolate.interp1d(xp,yp)
    return lambda x: interp(x, xp, yp, **kwargs)


def conflict_cascade(pA, pB, gA, gB=None, pwinB=0.5, pwinA=0.5, qA=1, qB=None, returnAll=False):
    '''

    :param pA:
    :param pB:
    :param gA:
    :param gB:
    :param pwinB:
    :param pwinA:
    :param qA:
    :param qB:
    :param returnAll:
    :return:
    '''
    if not qB:
        qB = qA
    if not gB:
        gB = gA
    psi_tag = [pA]
    phi_tag = [pB]
    while phi_tag[-1] > num_eps and psi_tag[-1] > num_eps:
        phi_last = phi_tag[-1]
        psi_last = psi_tag[-1]
        if verbose:
            print("(%.4f,%.4f)" % (phi_last, psi_last))
            raw_input()
        psi_tag.append(pA * (1 - qA * pwinB * pB * gB(phi_last)))
        phi_tag.append(pB * (1 - qB * pwinA * pA * gA(psi_last)))
        if abs(phi_tag[-2] - phi_tag[-1]) < num_eps and abs(psi_tag[-2] - psi_tag[-1]) < num_eps:
            break
        if phi_tag[-1] * gB(phi_tag[-1]) < num_eps and psi_tag[-1] * gA(psi_tag[-1]) < num_eps:
            break
    if returnAll:
        return zip(*map(lambda x, y: (x * gA(x), y * gB(y)), psi_tag, phi_tag))
    else:
        return [psi_tag[-1] * gA(psi_tag[-1]), phi_tag[-1] * gB(phi_tag[-1])]


def interdependent_cascade_simple(p, g,  returnAll=False):
    '''

    :param pA:
    :param pB:
    :param gA:
    :param gB:
    :param pwinB:
    :param pwinA:
    :param qA:
    :param qB:
    :param returnAll:
    :return:
    '''

    psi_tag = [p]
    phi_tag = [p]
    while phi_tag[-1] > num_eps and psi_tag[-1] > num_eps:
        phi_last = phi_tag[-1]
        psi_last = psi_tag[-1]
        if verbose:
            print("(%.4f,%.4f)" % (phi_last, psi_last))
            raw_input()
        psi_tag.append(p * g(phi_last))
        phi_tag.append(p * g(psi_last))
        if abs(phi_tag[-2] - phi_tag[-1]) < num_eps and abs(psi_tag[-2] - psi_tag[-1]) < num_eps:
            break
        if phi_tag[-1] * g(phi_tag[-1]) < num_eps and psi_tag[-1] * g(psi_tag[-1]) < num_eps:
            break
    if returnAll:
        return list(map(lambda x: x*g(x), psi_tag))
        #return zip(*map(lambda x, y: (x * g(x), y * g(y)), psi_tag, phi_tag))
    else:
        #return [psi_tag[-1] * g(psi_tag[-1]), phi_tag[-1] * g(phi_tag[-1])]
        return psi_tag[-1] * g(psi_tag[-1])


def new_conflict_cascade(pA, pB, gA, gB=None, pwinB=0.5, pwinA=0.5, qA=1, qB=None, returnAll=False):
    if not qB:
        qB = qA
    if not gB:
        gB = gA
    psi_tag = [pA*gA(pA)]
    phi_tag = [pB*gB(pB)]
    C = [1]
    while phi_tag[-1] > num_eps and psi_tag[-1] > num_eps:
        phi_last = phi_tag[-1]
        psi_last = psi_tag[-1]
        C.append(qA * psi_last * gA(psi_last) * phi_last * gB(phi_last))
        psi_tag.append(psi_last - C[-1] * pwinB)
        phi_tag.append(phi_last - C[-1] * pwinA)
        if abs(phi_tag[-2] - phi_tag[-1]) < num_eps and abs(psi_tag[-2] - psi_tag[-1]) < num_eps:
            break
    if returnAll:
        return zip(*map(lambda x, y, z: (x * gA(x), y * gB(y), z), psi_tag, phi_tag, C))
    else:
        return [psi_tag[-1] * gB(psi_tag[-1]), phi_tag[-1] * gA(phi_tag[-1])]


def conflict_self_consistent(pA, pB, gA, gB=None, pwinB=0.5, pwinA=0.5, returnAll=False):
    if gB is None:
        gB=gA
    conditionA = lambda x,y: x - pA*(1 - pwinB*pB*gB(y))
    conditionB = lambda x,y: y - pB*(1 - pwinA*pA*gA(x))
    def to_solve(xy):
        x,y=xy
        return np.c_[conditionA(x,y),conditionB(x,y)].flatten()
    sol, info, ier, msg = fsolve(to_solve,[0.5,0.5],full_output=1)
    x,y=sol
    return np.c_[x*gA(x),y*gB(y)].flatten()

def get_interdep_pinf_at_p(p, p_inf, **kwargs):
    if p == 1:
        return 1
    if p <= 0:
        return 0
    to_solve = lambda x: x - np.sqrt(p * p_inf(x))
    sol, info, ier, msg = fsolve(to_solve, 0.9, full_output=1, **kwargs)
    if ier == 1:
        return p_inf(sol.max())
    else:
        return 0


def get_interdep_pinf(p_vec, p_inf, **kwargs):
    p_inf_vec = 0 * p_vec
    for i, p in enumerate(p_vec):
        p_inf_vec[i] = get_interdep_pinf_at_p(p, p_inf, **kwargs)
    return interp_func(p_vec, p_inf_vec)


def graphical_solver(p, p_inf):
    import matplotlib.pyplot as plt

    x = np.linspace(0, 1, 1000)
    plt.plot(x, np.sqrt(p * p_inf(x)))
    plt.plot(x, x)
    plt.show()
    plt.savefig("/tmp/fig.pdf")


def fill_matrix(vec1, vec2, f):
    '''

    :param vec1: M x 1 vector of values
    :param vec2: N x 1 vector of values
    :param f:  function of two variables
    :return: matrix valued at grid of vec1,vec2 values
    '''
    Z = np.zeros((len(vec1), len(vec2)))
    for i, val1 in enumerate(vec1):
        for j, val2 in enumerate(vec2):
            Z[i, j] = f(val1, val2)
    return Z


def sharp_pinf(p_inf, beta, pc=0.5927, match_at=0.6):
    '''

    :param p_inf: percolation profile function
    :param beta: critical exponent to sharpen to
    :param pc: pc
    :param match_at: the point where you cross back to the pinf function from the forced scaling
    :return: new percolation function, matched at 0.6
    '''

    def f(p):

        if p < pc:
            return 0
        elif p < match_at:
            return p_inf(match_at) * ((p - pc) / (match_at - pc)) ** beta
        else:
            return p_inf(p)

    return f