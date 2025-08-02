from __future__ import division
import pyNoN.execute.load as ld
import sh
from collections import defaultdict
from itertools import product as iproduct
import numpy as np

def spontaneous_collapse(dtype=3,globname="/home/micha/rnet_results/multiplex/DependentLattices-type=%i-L=200-lambda=*"):
    zero=[]
    nonzero=[]
    data=[]
    if type(sh.glob(globname%dtype)) == list:
        for fname in sh.glob(globname%dtype):
            wcout = sh.wc(fname,"-l").stdout.strip().split()
            md = ld.get_metadata(wcout[1],keys=["k","lambda","L"])
            if wcout[0] == '3':
                zero.append((md["k"],md["lambda"]))
            else:
                nonzero.append((md["k"],md["lambda"]))
                s=ld.load_smart_summary(fname)
                aNg=np.array([i['Ng'] for i in s])
                dNg=np.diff(aNg)
                md["jump"]=dNg.max() 
                md["flag"]="OK"
                try:
                    md["pc"]=s[np.where(dNg==md["jump"])[0]]["p"]
                except TypeError:
                    md["pc"]=0
                    md["flag"]="CHECK_PC"
                try:
                    for i in s:
                        S=i["S"]
                        if S>0:
                            break
                except KeyError:
                    S=0
                md["sJump"] = S
                md["NOImax"]=max([i["NOI"] for i in s])
                md["run"]=s
                md["fname"]=fname
                md["jump"]/=(md["L"]*md["L"])
                data.append(md)

    return [sorted(zero),sorted(nonzero),data]
def both_spontaneous_collapse(globname="/home/micha/rnet_results/multiplex/DependentLattices-type=%i-L=200-lambda=*"):
    [z2,nz2,d2]=spontaneous_collapse(dtype=2,globname=globname)
    [z3,nz3,d3]=spontaneous_collapse(dtype=3,globname=globname)
    z={2:z2,3:z3}
    nz={2:nz2,3:nz3}
    d={2:d2,3:d3}
    return [z,nz,d]
def plot_spontaneous_collapse(zero,nonzero,dtype=3):
    from matplotlib import pyplot as plt
    z = zero[dtype]
    nz = nonzero[dtype]
    plt.plot( [i[0] for i in z], [1/i[1] for i in z], 'rs', ms=20,label="==0")
    plt.plot( [i[0] for i in nz], [1/i[1] for i in nz], 'bs', ms=20,label="!=0")
    plt.axis([2.4,5.1,0,10])
    plt.xlabel(r"$\langle k \rangle$")
    plt.ylabel(r"$\frac{1}{\lambda}$")
    plt.title("Type %i"%dtype)
    plt.legend()
'''
Returns value rescaled to interval [0,1] based on allvals
'''
def rescale_val(val,allvals):
    return (val - min(allvals))/( max(allvals) - min(allvals) )

def getMaxConsec(data,checkVal=True):
    from itertools import groupby
    groups=[]
    for k,g in groupby(data):
        if k == checkVal:
            groups.append(len(list(g)))
    return max(groups) if groups else 0

def plot_param(d_bothtypes,dtype=3,param="jump"):
    from matplotlib import pyplot as plt
    from matplotlib import cm
    d=d_bothtypes[dtype]
    allvals = [i[param] for i in d]
    plt.figure()
    colors=plt.imshow([allvals]*2)
    plt.clf()
    for i in d:
        plt.plot(i["k"],1/i["lambda"], 's', ms=20, color=cm.jet(rescale_val(i[param],allvals)))
    plt.colorbar(colors)
    plt.xlabel(r"$\langle k \rangle$")
    plt.ylabel(r"$\frac{1}{\lambda}$")
    xmin,xmax,ymin,ymax=plt.axis("tight")
    intx=xmax-xmin
    inty=ymax-ymin
    padding=0.05
    plt.axis([xmin-padding*intx,xmax+padding*intx,ymin-padding*inty,ymax+padding*inty])
    plt.title("Type %i : %s"%(dtype,param))


class MTXReader:
    def __init__(self,fname):
        if fname.find("_b.mtx") >0:
            self.vector_fname = fname
            self.matrix_fname = fname[:-5]+"SPD.mtx"
        else:
            self.matrix_fname=fname
            self.vector_fname=fname[:-7]+"b.mtx"
        self.load_matrix()
        self.load_vector()

    def load_matrix(self):
        from scipy import sparse
        f=open(self.matrix_fname)
        f.readline()
        md= f.readline().strip().split()
        M=int(md[0])
        N=int(md[1])
        nnz=int(md[2])
        mdata=[[int(j) for j in i.strip().split()] for i in f.readlines()]
        csrij = np.array([[i[0]-1,i[1]-1] for i in mdata]).T
        csrdata = [i[2] for i in mdata]
        self.L = sparse.csr_matrix((csrdata,csrij))
    
    def load_vector(self):
        self.v = np.array([float(i.strip()) for i in open(self.vector_fname).readlines()[2:]]) 
            