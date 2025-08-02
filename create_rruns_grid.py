from __future__ import division
from rrun import rrun, process_failed_list
from itertools import product as iproduct


def do_params(r):
    global q, L, deptype
    print [q, r, L]
    rr = rrun(q=q, r=r, L=L, deptype=deptype, submit=False)
    #print pcr.get_pc_guess(1)
    rr.get_timestamp()
    rr.kwargs["mode"] = 1
    rr.kwargs["p"] = 0.77
    rr.kwargs["dp"] = 0.0002
    rr.kwargs["dpc"] = 2
    rr.kwargs["steps"] = -1
    rr.random = False
    rr.dynamics = True
    rr.submitDynamics = True
    rr.run_program()

    return 0


if __name__ == '__main__':
    from sys import argv

    q = 1
    L = 1500
    deptype = argv[1]
    rmin = 0
    rmax = 2
    repititions = 8
    r_vals = [i for i in range(rmin,rmax+1)]*repititions
    try:
        from multiprocessing import Process, Pool, cpu_count
        pool = Pool(processes=int(argv[2]))
        res = pool.map_async(do_params, r_vals)
        print res.get()
    except ImportError:
        print "Not utilizing multiprocessing"
        map(do_params, r_vals)
    process_failed_list()
