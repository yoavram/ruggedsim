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
    c['tau'] = d['tau']
    c['pi'] = d['pi']    

    c['msb_W'] = d['msb']['W']

    c['app_time'] = d['app']['t'] - d['msb']['t']    
    c['app_W'] = d['app']['W']
    c['app_AB0'] = d['app']['p'][2][0]
    c['app_AB1'] = d['app']['p'][2][1]
    c['app_AB2'] = d['app']['p'][2][2]
    c['app_AB3'] = d['app']['p'][2][3]
    c['app_AB4'] = d['app']['p'][2][4]
    c['app_best'] = np.array(d['app']['p'][2]).nonzero()[0].min()

    c['fix_time'] = d['fix']['t'] - d['app']['t']
    c['fix'] = d['fix']['success']
    c['fix_W'] = d['fix']['W']
    c['fix_AB0'] = d['fix']['p'][2][0]
    c['fix_AB1'] = d['fix']['p'][2][1]
    c['fix_AB2'] = d['fix']['p'][2][2]
    c['fix_AB3'] = d['fix']['p'][2][3]
    c['fix_AB4'] = d['fix']['p'][2][4]
    if c['fix']:
        c['fix_best'] = np.array(d['fix']['p'][2]).nonzero()[0].min()
    else:
        c['fix_best'] = None

    df.append(c)
df = pd.DataFrame(df)
print 'dataframe shape:', df.shape

fname = 'df_fix_%s.csv.gz' % jobname
with gzip.open(fname, 'wb') as f:
	df.to_csv(f)
print "DataFrame written to", fname
