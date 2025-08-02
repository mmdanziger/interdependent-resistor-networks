import pyNoN.execute.load as l
from glob import glob
import matplotlib.pyplot as plt

def keyv(d,k):
    return [i[k] for i in d]

flist = glob("/home/micha/phd/resistornets/protoresults/DependentLattices-type=2-L=1000-q=1*.summary")
flist = glob("/tmp/protoresults/DependentLattices-type=2-L=1000-q=1*.summary")

def plot_flist(flist):
    results=[]
    for fname in flist:
        d = l.get_metadata(fname, ["q","r","L","type"])
        results_ =  l.load_smart_summary(fname)
        d["results"] = results_
        results.append(d)
    results.sort(key= lambda x :x["r"])
    #i want it sorted right...
    plt.figure()
    for r_dict in results:
        results_ = r_dict["results"]    
        plt.plot(keyv(results_,"p"),keyv(results_,"S"),'.-',label="r = %i"%r_dict["r"])
    plt.legend()
    plt.title("q=%.4f L=%i type=%i"%(d["q"],d["L"],d["type"]))
    return results
