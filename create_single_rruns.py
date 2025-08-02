from __future__ import division
from rrun import rrun


def do_r(x):
    global q,L,deptype
    print [q,x,L]
    rr=rrun(q=q,r=x,L=L,deptype=deptype,submit=False)
    #print pcr.get_pc_guess(1)
    rr.get_timestamp()
    rr.kwargs["mode"] = 3
    rr.kwargs["p"] = 0.84
    rr.kwargs["dp"] = 0.004
    rr.kwargs["steps"] = 30
    rr.dynamics = True
    rr.submitDynamics=True
    rr.run_program()
    return 0



if __name__ == '__main__':
    from sys import argv
    q=1
    r_list=range(25)
    L=2000
    deptype=argv[1]
    try:
        from multiprocessing import Process,Pool,cpu_count 
        pool = Pool(processes=int(argv[2]))
        res = pool.map_async(do_r,r_list)
        print res.get()
    except ImportError:
        print "Not utilizing multiprocessing"    
        map(do_r,r_list)
    
