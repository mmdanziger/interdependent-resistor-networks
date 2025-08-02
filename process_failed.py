from __future__ import division
import rrun
from sys import argv
from os import path

if __name__ == "__main__":
    basepath = argv[1]
    listfname = path.join(basepath,rrun.RetryList)
    rrun.process_failed_list(listfname)
