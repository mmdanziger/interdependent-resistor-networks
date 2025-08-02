from __future__ import division
from networkx import *
import numpy as np
try:
    from networkx.linalg import laplacian_matrix
except:
    pass
from itertools import product as iproduct
eps=1e-9

def continue_attack(p,ao=None,N=50):
    """
    Use like this:
    ao=continue_attack(0.7,N=70)
    for p in [0.67,0.65,0.63,0.6]:
        continue_attack(p,ao,N=70)

    """
    rn = resnet.resnet(N=N)
    if ao is None:
        rn.initializeAttackOrder()
    else:
        rn.loadAttackOrder(ao)
    rn.doAttack(p)
    rn.createColoredDeathMap()
    from matplotlib import pyplot as plt
    plt.savefig("/tmp/RNetIllustration_N%i_p%i.pdf"%(N,int(p*100)))
    return rn.attackOrder



class resnet:
    def __init__(self,N=50, backend='amg'):
        self.N = N
        self.backend = backend
        self.generateNNLattice()
        self.generateNNNLattice()
        self.aliveA=dict((i,1) for i in self.nnLattice.nodes())
        #self.alibeB={}
        self.nnClusters=[]
        self.nnnClusters=[]
    def generateNNLattice(self):
        self.nnLattice = grid_graph([self.N,self.N])
        for node in self.nnLattice.nodes_iter():
            self.nnLattice.node[node]["label"]=str(node)
    def generateNNNLattice(self):
        N=self.N
        self.nnnLattice = grid_graph([self.N,self.N])
        x=[( (i,j),(i+1,j+1) ) for i,j in iproduct(range(N-1),range(N-1)) ]
        y=[( (i,j),(i+1,j-1) ) for i,j in iproduct(range(N-1),range(1,N)) ]
        self.nnnLattice.add_edges_from(x+y)
    def getNNClusters(self):
        self.nnClusters=[]
        seen=[v for v in self.nnLattice if v not in self.aliveA]
        for source in self.nnLattice:
            if source not in seen:
                cluster=[]
                nextLevel = {source}
                while nextLevel:
                    thisLevel=nextLevel
                    nextLevel = {}
                    for v in thisLevel:
                        if v not in seen:
                            seen.append(v)
                            cluster.append(v)
                            nextLevel.update(self.nnLattice[v])
                self.nnClusters.append(cluster)
    def trimToBackbone(self):
        self.nnClusters=[]
        seen=[v for v in self.nnLattice if v not in self.aliveA]
        for source in self.nnLattice:
            if source not in seen:
                cluster=[]
                nextLevel = {source}
                while nextLevel:
                    thisLevel=nextLevel
                    nextLevel = {}
                    for v in thisLevel:
                        if v not in seen:
                            seen.append(v)
                            cluster.append(v)
                            filtered_d = dict( (j,k) for j,k in self.nnLattice[v].iteritems() if 'current' in k and abs(k['current'])>eps)
                            nextLevel.update(filtered_d)
                self.nnClusters.append(cluster)
        giant = self.getNNGiant()
        self.aliveA = dict( (k,1) for k in giant)
        for v in self.nnLattice:
            if v not in giant and self.nnLattice.node[v]['dead'] == '0':
                self.nnLattice.node[v]['dead']='current'
    def restoreAlive(self):
        self.generateNNLattice()
        self.generateNNNLattice()
        self.aliveA = dict((i,1) for i in self.nnLattice.nodes())
        
    def getNNGiant(self):
        if self.nnClusters == []:
            self.getNNClusters()
        self.nnClusters = sorted(self.nnClusters, key = lambda x: len(x), reverse=True)
        self.nnGiant = self.nnClusters[0]
        return self.nnGiant
    def trimToGiant(self):
        giant = self.getNNGiant()
        self.aliveA = dict( (k,1) for k in giant)
        for v in self.nnLattice:
            if v not in giant and self.nnLattice.node[v]['dead'] == '0':
                self.nnLattice.node[v]['dead']='perc'
                   
    def getNNNClusters(self):
        self.nnnClusters=[]
        seen=self.aliveA.keys()
        for source in self.nnnLattice:
            if source not in seen:
                cluster=[]
                nextLevel = {source}
                while nextLevel:
                    thisLevel=nextLevel
                    nextLevel = {}
                    for v in thisLevel:
                        if v not in seen:
                            seen.append(v)
                            cluster.append(v)
                            nextLevel.update(self.nnnLattice[v])
                self.nnnClusters.append(cluster)
    def getNNNGiant(self):
        if self.nnnClusters == []:
            self.getNNNClusters()
        self.nnnClusters = sorted(self.nnnClusters, key = lambda x: len(x), reverse=True)
        self.nnnGiant = self.nnnClusters[0]
        return self.nnnGiant
    def initializeAttackOrder(self):
        from numpy.random import rand,shuffle
        self.attackOrder =  list(self.nnLattice.nodes())
        np.random.shuffle(self.attackOrder)
    def loadAttackOrder(self,attackOrder):
        if attackOrder and len(attackOrder) == self.N*self.N:
            self.attackOrder = attackOrder
            return
        else: 
            raise ValueError("Bad attack order passed")
        
    def randAttack(self,p):
        try:
            len(self.attackOrder)
        except AttributeError:
            self.initializeAttackOrder()
        self.restoreAlive()
        until_index = self.N*self.N - int(p*self.N*self.N)
        for i in range(until_index):
            v = self.attackOrder[i]
            if v in self.aliveA:
                del self.aliveA[v] 
            self.nnLattice.node[v]['dead']='attack'
        for i in range(until_index,self.N*self.N):
            v = self.attackOrder[i]
            self.nnLattice.node[v]['dead']='0'
        
    def generateAliveMatrix(self):
        from numpy import zeros
        self.A = zeros((self.N,self.N)) 
        for v in self.aliveA:
            self.A[v]=1
    def generateClusterNumberMatrix(self):
        from numpy import zeros
        self.A = zeros((self.N,self.N))
        if self.nnnClusters == []:
            self.getNNNClusters()
        for cIndex, cluster in enumerate(self.nnnClusters):
            for node in cluster:
                self.A[node]=cIndex+1
    def findVoltageCurrentSource(self):
        from numpy import array,matrix,zeros
        v0 = zeros(len(self.aliveA))
        ground=[]
        for node in self.nnLattice:
            self.nnLattice.node[node]["voltage"]=-1
        foundsource=0
        foundsink=0
        for i,node in enumerate(sorted(self.aliveA.keys())):
            if node[0]==0:
                v0[i]=1
                foundsource=1
            if node[0]==self.N-1:
                ground.append(i)
                foundsink=1
        if foundsource*foundsink == 0:
            print "No spanning cluster"
            return
        L = laplacian_matrix(self.nnLattice,nodelist=sorted(self.aliveA.keys()))
        L=L[:-len(ground),:-len(ground)]
        v0 = v0[:-len(ground)]
        from scipy import sparse
        Lsparse = sparse.csr_matrix(L)
        if self.backend == 'amg':
            import pyamg
            solver = pyamg.smoothed_aggregation_solver(Lsparse)
            self.v = solver.solve(v0, tol=1e-19, accel='cg')
        else:
            from scipy.sparse.linalg import spsolve
            self.v = spsolve(Lsparse,v0)
        
        for i,node in enumerate(sorted(self.aliveA.keys())):
            try:
                self.nnLattice.node[node]['voltage'] = self.v[i]
            except IndexError:
                self.nnLattice.node[node]['voltage'] = 0
        self.v0 = v0
        self.ground = ground
        self.Laplacian = L
    def findVoltageVoltageSource(self):
        from numpy import array,matrix,zeros,dot
        v0 = zeros(len(self.aliveA))
        ground=[]
        for node in self.nnLattice:
            self.nnLattice.node[node]["voltage"]=-1
            self.nnLattice.node[node]["label"]=str(node)
        sources=0
        sinks=0
        for i,node in enumerate(sorted(self.aliveA.keys())):
            if node[0]==0:
                v0[i]=1
                sources+=1
            if node[0]==self.N-1:
                ground.append(i)
                sinks+=1
        if sources*sinks == 0:
            print "No spanning cluster"
            return False
        L = laplacian_matrix(self.nnLattice,nodelist=sorted(self.aliveA.keys()))
        
        b = (-dot(L,v0)).A[0]
        L=L[sources:-sinks,sources:-sinks]
        b = b[sources:-sinks]
        from scipy import sparse
        Lsparse = sparse.csr_matrix(L)
        if self.backend == 'amg':
            import pyamg
            solver = pyamg.smoothed_aggregation_solver(Lsparse)
            self.v = solver.solve(b, tol=1e-19, accel='cg')
        else:
            from scipy.sparse.linalg import spsolve
            self.v = spsolve(Lsparse,b)
        v0[sources:-sinks]=self.v
        self.v=v0
        for i,node in enumerate(sorted(self.aliveA.keys())):
            try:
                self.nnLattice.node[node]['voltage'] = self.v[i]
            except IndexError:
                self.nnLattice.node[node]['voltage'] = 0
        self.v0 = v0
        self.ground = ground
        self.Laplacian = L
        return True

    def findCurrent(self):
        try:
            self.v
        except:
            if not self.findVoltageVoltageSource():
                return False
        from numpy import dot,resize,array
        nodelist = [node for i,node in enumerate(sorted(self.aliveA.keys()))]
        edgelist = [(i,j) for i,j in self.nnLattice.edges() if i in nodelist and j in nodelist]
        A = incidence_matrix(self.nnLattice, nodelist=nodelist, edgelist=edgelist, oriented=True)
        
        v=resize(self.v,len(self.aliveA))
        v[-len(self.ground):]=0
        self.current = array(dot(v,A))[0]
        if not (len(self.current) == len(edgelist)):
            print "Edgelist %i " % len(edgelist)
            print "current vec",; print self.current.shape
        for i,(s,t) in enumerate(edgelist):
            self.nnLattice.edge[s][t]['current'] = self.current[i]
        self.Ainc = A
        self.v = v
    def doAttack(self,p):
        self.randAttack(p)
        self.trimToGiant()
        if self.findVoltageVoltageSource():
            self.findCurrent()
            self.trimToBackbone()
            return True
        else:
            return False
    def getPartitionSizes(self):
        from collections import defaultdict
        d=defaultdict(int)
        for n in self.nnLattice:
            d[self.nnLattice.node[n]['dead']]+=1
        return d
    def findTotalOutCurrent(self):
        try: 
            self.current
        except:
            self.findCurrent()
        out_current = 0
        nodelist = [node for i,node in enumerate(sorted(self.aliveA.keys()))]
        edgelist = [(i,j) for i,j in self.nnLattice.edges() if i in nodelist and j in nodelist]
        for i,(s,t) in enumerate(edgelist):
            if s[0] == self.N-1 or t[0] == self.N-1:
                print (s,t),
                print self.nnLattice.edge[s][t]['current']
                out_current+=abs(self.nnLattice.edge[s][t]['current'])
        self.out_current = out_current
        return out_current
    def findTotalInCurrent(self):
        try: 
            self.current
        except:
            self.findCurrent()
        in_current=0
        nodelist = [node for i,node in enumerate(sorted(self.aliveA.keys()))]
        edgelist = [(i,j) for i,j in self.nnLattice.edges() if i in nodelist and j in nodelist]
        for i,(s,t) in enumerate(edgelist):
            if s[0] == 0 or t[0] == 0:
                print (s,t),
                print self.nnLattice.edge[s][t]['current']
                in_current+=abs(self.nnLattice.edge[s][t]['current'])
        self.in_current = in_current
        return in_current
    def generateVoltageMatrix(self):
        try:
            self.v
        except:
            self.findVoltage()
        from numpy import ones
        self.Av = ones((self.N,self.N))*-1
        livenodes=sorted(self.aliveA.keys())
        for i,voltage in enumerate(self.v):
            self.Av[livenodes[i]] = voltage
    def writeForGephi(self):
        '''we need to make the ids ints...'''
        outG=Graph()
        nodelist = [i for i in sorted(self.nnLattice.nodes())]
        nodeinfo = [(ind+1, self.nnLattice.node[i]) for ind,i in enumerate(sorted(self.nnLattice.nodes()))]
        outG.add_nodes_from(nodeinfo)
        edgelist = [(nodelist.index(i)+1,nodelist.index(j)+1,self.nnLattice.edge[i][j]) for i,j in self.nnLattice.edges_iter()]
        outG.add_edges_from(edgelist)
        for i in outG.nodes_iter():
            try:
                outG.node[i]['voltage'] = float(outG.node[i]['voltage'])
            except:
                pass
        for i,j in outG.edges_iter():
            try:
                outG.edge[i][j]['current'] = float(outG.edge[i][j]['current'])
            except:
                pass
        return outG
    def getDeathMap(self):
        from numpy import zeros
        A = zeros((self.N,self.N),dtype=str)
        for i,j in  iproduct(range(self.N),range(self.N)):
            A[i,j] = self.nnLattice.node[(i,j)]['dead']
        return A
    def colorDeathMap(self,A,cdict={'p':[0,1,0],'a':[0,0,0],'0':[1,1,1],'c':[0.3,0.3,1]}):
        from numpy import zeros
        from numpy.random import rand
        types = list(set(A.flatten()))
        numdict = {'0':0,'c':1,'p':2,'a':3}
        for t in types:
            if t not in cdict:
                cdict[t]=[rand(),rand(),rand()] 
        import matplotlib
        cmap=matplotlib.colors.ListedColormap([(1,1,1,1),(0,1,0,1),(0.3,0.3,1,1),(0,0,0,1)])
        bounds=[0,1,2,3,4]
        norm=matplotlib.colors.BoundaryNorm(bounds,cmap.N)
        B = zeros((self.N,self.N))
        for i,j in iproduct(range(self.N),range(self.N)):
            B[i,j]=numdict[A[i,j]]
        #B = zeros((self.N,self.N,3))
        #for i,j in iproduct(range(self.N),range(self.N)):
            #B[i,j]=cdict[A[i,j]]
        return (B,cmap,norm)
    def createColoredDeathMap(self):
        '''Expects that you have already initialized and run self.doAttack(p)'''
        A=self.getDeathMap()
        (B,cmap,norm)=self.colorDeathMap(A)
        import matplotlib.pyplot as plt
        from numpy import arange
        plt.figure()
        plt.pcolormesh(arange(self.N),arange(self.N),B,cmap=cmap,norm=norm)
        plt.axis("image")
        plt.xticks([])
        plt.yticks([])
        '''save if you want'''
        
        
        