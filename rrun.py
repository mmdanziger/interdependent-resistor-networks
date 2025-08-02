from __future__ import division, print_function
from random import random, randrange
from subprocess import call, STDOUT
from time import strftime, sleep, time
from itertools import product as iproduct
from os import remove
import sys, pickle, os

sys.path += ['/home/micha/phd', '/home/micha/phd/resistornets']
from pyNoN.execute.load import *

try:
    import simplejson as json
except ImportError:
    import json

RetryList = "FailedSubmissions.txt"


def add_to_retry_list(fname):
    global RetryList
    with open(RetryList, 'a') as f:
        f.write(fname + '\n')


def process_failed_list(listfname=None):
    if listfname is None:
        listfname = RetryList
    try:
        flist = [i.strip() for i in open(listfname)]
    except IOError:
        print("No failed submissions found.")
        return
    still_failed = []
    for fname in flist:
        try:
            open(fname).close()
            if rrun().submit_dynamics(fname) == 0:
                still_failed.append(fname)
        except IOError:
            print("File not found: %s " % fname)
            still_failed.append(fname)
    if len(still_failed) != 0:
        with open(listfname, 'w') as f:
            f.write("\n".join(still_failed))
    else:
        remove(listfname)


class rrun:
    def __init__(self, q=1, k=4, r=-1, lam=-1, L=1000, deptype=3, submit=True):
        self.kwargs = {}
        self.kwargs["q"] = q
        self.kwargs["L"] = L
        self.kwargs["k"] = k
        self.kwargs["lambda"] = lam
        self.kwargs["type"] = deptype
        self.kwargs["dp"] = 0.001
        self.kwargs["dpc"] = 100
        self.kwargs["steps"] = 100
        self.kwargs["p"] = 0.75
        self.kwargs["r"] = r
        self.history = True
        self.dynamics = True
        self.profile = False
        self.random = False
        self.error = False
        self.pc_guess = {"pc": 0.75, "std": 0}
        self.submit = submit
        self.submitDynamics = False

    def get_pc_guess(self, getmean=1):
        url = "http://physionet.ph.biu.ac.il/~micha/api/pc.php"
        params = {"q": self.kwargs["q"], "r": self.kwargs["r"], "L": self.kwargs["L"], "mean": getmean}
        try:
            import requests

            r = requests.request("GET", url, params=params)
            if len(r.text) == 0:
                url = "http://physionet.ph.biu.ac.il/~micha/api/pc_bounds.php"
                r = requests.request("GET", url, params=params)
            if len(r.text) == 0:
                print("Failed to get response from server")
                return None
            try:
                d = json.loads(r.text)
            except:
                print("Unable to parse response " + r.text)
                return None
        except ImportError:
            import urllib2, urllib

            query = "?" + urllib.urlencode(params)
            r = urllib2.urlopen(urllib2.Request(url + query))
            if len(r.read()) == 0:
                url = "http://physionet.ph.biu.ac.il/~micha/api/pc_bounds.php"
                r = urllib2.urlopen(urllib2.Request(url + query))
            if len(r.read()) == 0:
                print("Failed to get response from server")
                return None
            try:
                d = json.loads(r.read())
                self.pc_guess = d
            except:
                print("Unable to parse response " + r.read())
                return None

        self.pc_guess = d
        self.kwargs["p"] = d["pc"]
        if self.pc_guess["std"] > 0:
            self.kwargs["dp"] = self.pc_guess["std"]
        return d

    def get_timestamp(self):
        if "timestamp" in self.kwargs:
            return self.kwargs["timestamp"]
        t = time()
        if (t - int(t)) > 0:
            t = str(t - int(t)).split('.')[-1]
        else:
            t = str(int(t))
        t = t[:9]  # for the 32 bit longs...
        if t[0] == '0':
            t = t[1:]
        self.kwargs["timestamp"] = t
        return t

    def get_kw_array(self):
        kwarray = ["idresnet.bin"]
        for k, v in self.kwargs.items():
            kwarray.append("--" + k)
            kwarray.append(str(v))
        if self.history:
            kwarray.append("--history")
        if self.dynamics:
            kwarray.append("--dynamics")
        if self.profile:
            kwarray.append("--profile")
        if self.error:
            kwarray.append("--error")
        if self.random:
            kwarray.append("--random")
        return kwarray

    def run_program(self):
        call(self.get_kw_array())
        if self.submit:
            self.submit_smart_summary()
        if self.submitDynamics:
            self.submit_dynamics()

    def submit_dynamics(self, fname=None):
        from glob import glob

        if fname is None:
            file_list = glob("RNetDynamics*-%s.json" % self.kwargs["timestamp"])
            if len(file_list) != 1:
                print("Cannot find file. Found:"),
                print(file_list)
                return 0
            fname = file_list[0]
            metadict = self.kwargs.copy()
        else:
            metadict = dict()
        metadict["filename"] = fname
        metadata = json.dumps(metadict)

        with open(fname) as f:
            try:
                results = json.load(f)
            except:
                print("Unable to load %s" % fname)
                return 0
        import requests

        url = "http://physionet.ph.biu.ac.il/~micha/rnets/api/v2/add_run.php"
        headers = {'Content-type': 'application/json', 'Accept': 'text/plain'}
        r = requests.request("POST", url, data=json.dumps({"metadata": metadata, "results": results}), headers=headers)
        try:
            server_report = json.loads(r.text)
            if server_report["status"] != "success":
                print("Server failed with message status %s" % server_report["status"])
                add_to_retry_list(fname)
                return 0

        except:
            print("Unable to parse server response: %s" % r.text)
            add_to_retry_list(fname)
            return 0

        return 1

    def submit_smart_summary(self):
        from glob import glob

        globme = "DependentLattices*" + self.kwargs["timestamp"] + ".summary"
        print(globme)
        try:
            fname = glob(globme)[0]
        except:
            print("No file found matching: " + globme)
            return
        d = get_metadata(fname, ["q", "r", "L", "type"])
        if self.already_submitted():
            return
        results = load_smart_summary(fname)
        d['double'] = 1
        d['precision'] = self.kwargs['dpc']
        d['file'] = os.path.basename(fname)
        d['results'] = results
        url = "http://physionet.ph.biu.ac.il/~micha/rnets/add_run.php"
        headers = {'Content-type': 'application/json', 'Accept': 'text/plain'}
        try:
            import requests

            r = requests.request("POST", url, data=json.dumps([d]), headers=headers)
            print(r.text)
        except ImportError:
            import urllib, urllib2

            r = urllib2.urlopen(urllib2.Request(url, json.dumps([d]), headers))
            print(r.read())

    def already_submitted(self):
        '''
        This will eventually be a program that queries the server to see if
        a given file has been submitted to avoid the overhead of loading and transferring the information when it may be redundant
        '''
        return False

    def load_summary(self):
        from glob import glob

        globme = "DependentLattices-L=*" + self.kwargs["time_stamp"] + ".summary"
        print(globme)
        for fname in glob(globme):
            print("Found " + fname)
            submit_summary(fname)

    def load_results(self):
        from glob import glob

        for fname in glob("DependentLattices*" + self.kwargs["time_stamp"] + ".output"):
            submit_results(fname)
