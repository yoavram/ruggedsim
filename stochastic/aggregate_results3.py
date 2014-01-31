import sys
import glob
import os
import io
import time
import gzip
import ujson as json
import pandas as pd
import re
import numpy as np

if len(sys.argv) == 0:
    print "Please give a jobname!"
    sys.exit()
jobname = sys.argv[1] 
print "job name:", jobname

files = glob.glob("output/"+jobname+"/*.data")
print 'files:', len(files), files[0]

data = []
for fname in files:
    with io.BufferedReader(gzip.open(fname)) as f:
        d = json.load(f)
        data.append(d)
print 'records', len(data)

df = []
for d in data:
    c = {}
    c['U'] = d['U']
    c['beta'] = d['beta']
    c['s'] = d['s']
    c['pop_size'] = d['pop_size']
    c['G'] = d['G']
    c['H'] = d['H']
    c['tau1'] = d['tau1']
    c['pi1'] = d['pi1']
    c['tau2'] = d['tau2']
    c['pi2'] = d['pi2']
    c['invasion_rate'] = d.get('invasion_rate', 0.5)
    c['app_time'] = d['app']['t'] - d['msb']['t']    
    c['fix_time'] = d['fix']['t'] - d['app']['t']
    c['fix'] = d['fix']['success']
    c['invader'] = d['fix']['invader']
    df.append(c)
df = pd.DataFrame(df)
print 'dataframe shape:', df.shape

with gzip.open("df_%s.csv.gz" % jobname, 'wb') as f:
	df.to_csv(f)
print "DataFrame written to df_%s.csv.gz" % jobname
