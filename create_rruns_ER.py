from __future__ import division
from rrun import rrun,process_failed_list
from itertools import product as iproduct


def do_params(kavg):
    global q,L,deptype
    print [q,kavg,L]
    rr=rrun(q=q,k=kavg,L=L,deptype=deptype,submit=False)
    #print pcr.get_pc_guess(1)
    rr.get_timestamp()
    rr.kwargs["mode"] = 1
    rr.kwargs["p"] = 0.2
    rr.kwargs["dp"] = 0.01
    rr.kwargs["dpc"] = 50
    rr.kwargs["steps"] = -1
    rr.random = True
    rr.dynamics = True
    rr.submitDynamics = True
    rr.run_program()
    
    return 0



if __name__ == '__main__':
    from sys import argv
    q=0
    L=500
    deptype=argv[1]
    kmin=4
    kmax=6
    dk=0.1
    ksteps = int((kmax-kmin)/dk)
    k_vals = [kmin + i*dk for i in range(ksteps)]
    try:
        from multiprocessing import Process,Pool,cpu_count 
        pool = Pool(processes=int(argv[2]))
        res = pool.map_async(do_params,k_vals)
        print res.get()
    except ImportError:
        print "Not utilizing multiprocessing"    
        map(do_params,param_list)
    process_failed_list()
