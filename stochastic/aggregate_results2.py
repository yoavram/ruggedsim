import sys
import glob
import os
import io
import time
import cPickle
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
    c['tau'] = d['tau']
    c['pi'] = d['pi']
    c['app_time'] = d['app']['t'] - d['msb']['t']
    c['app_who'] = np.argmax(d['app']['p'][2])
    c['fix_time'] = d['fix']['t'] - d['app']['t']
    c['fix'] = d['fix']['success']
    if c['fix']:
        c['fix_who'] = np.argmax(d['fix']['p'][2])
    else:
        c['fix_who'] = None
    app_p_ab = d['app']['p'][0]
    c['app_click'] = app_p_ab[0] == 0   
    c['app_almost_click'] = app_p_ab[0] < app_p_ab[1]
    fix_p_AB = d['fix']['p'][2]
    c['fix_click'] = np.argmax(d['fix']['p'][2]) > np.argmax(d['app']['p'][2])
    df.append(c)
df = pd.DataFrame(df)
print 'dataframe shape:', df.shape

with gzip.open("df_%s.csv.gz" % jobname, 'wb') as f:
	df.to_csv(f)
print "DataFrame written to df_%s.csv.gz" % jobname
