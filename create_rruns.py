from __future__ import division
from rrun import rrun
from itertools import product as iproduct

def do_params(x):
    global q,L,deptype
    print [q,x,L]
    rr=rrun(q=q,k=x[0],lam=x[1],L=L,deptype=deptype,submit=False)
    #print pcr.get_pc_guess(1)
    rr.get_timestamp()
    rr.kwargs["mode"] = 1
    rr.kwargs["p"] = 1
    rr.kwargs["dp"] = 0.01
    rr.kwargs["dpc"] = 100
    rr.kwargs["steps"] = -1
    rr.dynamics = True
    #rr.submitDynamics = True
    rr.run_program()
    return 0



if __name__ == '__main__':
    from sys import argv
    q=1
    L=400
    deptype=argv[1]
    kmin=4
    kmax=5
    dk=0.1
    ksteps = int((kmax-kmin)/dk)
    k_vals = [kmin + i*dk for i in range(ksteps)]
    lambda_vals = [1/i for i in range(4,15)]
    param_list = iproduct(k_vals,lambda_vals)
    try:
        from multiprocessing import Process,Pool,cpu_count 
        pool = Pool(processes=int(argv[2]))
        res = pool.map_async(do_params,param_list)
        print res.get()
    except ImportError:
        print "Not utilizing multiprocessing"    
        map(do_params,param_list)
    
