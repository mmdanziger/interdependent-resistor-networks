import execute.find_pc as fp


def do_r(x):
    global q,L
    print [q,x,L]
    pcr=fp.pc_run(q=q,r=x,L=L,m=3,submit=False)
    #print pcr.get_pc_guess(1)
    pcr.get_timestamp()
    pcr.pc_guess["pc"]=0.74
    pcr.default_step=0.001
    pcr.run_program()
#    pcr.load_results()#also submits
    return 0



if __name__ == '__main__':
    q=0.7
    r_list=range(1)
    L=4000
    todo = {0.7:range(50)}
    for q in todo:
        r_list = todo[q]
        try:
            from multiprocessing import Process,Pool,cpu_count 
            pool = Pool(processes=6)
            res = pool.map_async(do_r,r_list)
            print res.get()
        except ImportError:
            print "Not utilizing multiprocessing"    
            map(do_r,r_list)
        
import execute.load as l
from glob import glob

for fname in glob("*L=4000*.output"):
    r=l.submit_results(fname)#asynchronous makes the db choke