from __future__ import division
from scipy.interpolate import interp1d
from numpy import linspace,mean,std
import requests
try: import simplejson as json
except ImportError: import json

class run_averager:
    def __init__(self,runList=None):
        self.runList = runList
        self.interp_min = 0
        self.interp_max = 1
        self.interp_steps = 1000
        self.interp_funclist = []
    def get_interp_vec(self):
        return linspace(self.interp_min,self.interp_max,self.interp_steps)
        
    def download_qrl(self,q=0.7,r=3,L=2900):
        url = "http://physionet.ph.biu.ac.il/~micha/api/qrl_all_runs.php?q=%.2f&r=%i&L=%i"%(q,r,L)
        req = requests.request("GET",url)
        d = json.loads(req.text)
        self.runList = [ [v["p"],v["N"]] for v in d.values()]
    def get_interp_funclist(self,force=False):
        if len(self.interp_funclist)>0 and not force:
            return
        for x,y in self.runList:
            self.interp_funclist.append(interp1d(x,y))
            #self.min_list.append(min(x))
            #self.max_list.append(max(x))
    def scatter_all_runs(self):
        import matplotlib.pyplot as plt
        plt.figure()
        for x,y in self.runList:
            plt.plot(x,y,'.')
        plt.show()
    def average_runs(self):
        ivec = self.get_interp_vec()
        self.interp_vals=[]
        self.get_interp_funclist()
        for i in ivec:
            self.interp_vals.append(0)
            N=0
            for f in self.interp_funclist:
                try:
                    self.interp_vals[-1]+=f(i)
                    N+=1
                except ValueError:
                    pass
            if N>0:
                self.interp_vals[-1]/=N
    def interp_evaluate(self,value):
        self.get_interp_funclist()
        evaluations=[]
        for f in self.interp_funclist:
            try:
                evaluations.append(f(value))
            except ValueError:
                pass
        return [mean(evaluations),std(evaluations),len(evaluations)]
    def get_average_curve(self,interp_min=None,interp_max=None,interp_steps=None):
        if interp_min != None:
            self.interp_min = interp_min
        if interp_max != None:
            self.interp_max = interp_max
        if interp_steps != None:
            self.interp_steps = interp_steps
        self.average_runs()
        return [self.get_interp_vec(),self.interp_vals]