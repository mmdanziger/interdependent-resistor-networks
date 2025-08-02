from __future__ import division,print_function
from rrun import rrun,process_failed_list
from itertools import product as iproduct
import numpy as np
from numpy.random import shuffle

def do_params(lam):
    global q,L,kavg,deptype
    print([q,kavg,lam,L])
    rr=rrun(q=q,k=kavg,L=L,deptype=deptype,submit=False)
    #print pcr.get_pc_guess(1)
    rr.get_timestamp()
    rr.kwargs["mode"] = 1
    rr.kwargs["lambda"] = lam
    rr.kwargs["p"] = 1 / kavg
    rr.kwargs["dp"] = 0.005
    rr.kwargs["dpc"] = 5
    rr.kwargs["steps"] = -1
    rr.random = False
    rr.dynamics = True
    rr.history = False
    rr.submitDynamics = True
    rr.run_program()
    
    return 0



if __name__ == '__main__':
    from sys import argv
    q=0.5
    L=1000
    deptype=argv[1]
    kavg=4
    #zeta_list =  [36,37,38,39] + [40,41,42,43,44,45,46,47,48,49,50]*2
    zeta_list = [0.2,0.5,1,2,5,10,11,12,13,14,15,16,17,18,19,20]
    lambda_list = [1 / i for i in zeta_list]*2
    shuffle(lambda_list) #so that the results come in in different orders, to get a better picture
    try:
        from multiprocessing import Process,Pool,cpu_count 
        pool = Pool(processes=int(argv[2]))
        res = pool.map_async(do_params,lambda_list)
        print(res.get())
    except ImportError:
        print( "Not utilizing multiprocessing")
        list(map(do_params,param_list))
    process_failed_list()
