from __future__ import division 
from glob import glob
from sys import argv
from os import path

flist = glob(argv[1]+"/*.matrix")

def entries(fname):
    with open(fname) as f:
        for i,line in enumerate(f.xreadlines()):
            pass
    return i-1

def newfname(fname):
    return fname[:-7]+"_SPD.mtx"

for fname in flist:
    entry_count = entries(fname)
    f = open(fname)
    g = open(newfname(fname),'w')
    g.write(f.readline())
    g.write(f.readline().strip()[:-1]+str(entry_count)+"\n")
    for line in f.readlines():
        g.write(line)
    f.close()
    g.close()