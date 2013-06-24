import sys
import glob
import os
import time
import cPickle
import gzip
import ujson as json
import pandas as pd
import re
import numpy as np

filename_pattern = re.compile('^pop_(?P<pop>\d+)_G_(?P<G>\d+)_s_(?P<s>\d\.?\d*)_H_(?P<H>\d\.?\d*)_U_(?P<U>\d\.?\d*)_beta_(?P<beta>\d\.?\d*)_pi_(?P<pi>\d\.?\d*)_tau_(?P<tau>\d\.?\d*)_(?P<date>\d{4}-\w{3}-\d{1,2})_(?P<time>\d{2}-\d{2}-\d{2}-\d{6}).(?P<extension>\w+)$')

def parse_filename(fname):
	m = filename_pattern.match(fname)
	if m:
		return m.groupdict()
	else:
		return dict()

def process_data_file(fname):
	fpath = folder + fname
	params = parse_filename(fname)
	if not params:
		print "Failed parsing file name", fpath
		return {},[],[]
	with gzip.open(fpath) as f:
		data = json.load(f,precise_float=True)
	if not data:
		print "Failed reading data", fpath
		return {},[],[]
	data.update(params)
	W = data.pop('W')
	p = data.pop('p')
	apps = None
	try:
		apps =  data.pop('apps')
	except KeyError:
		pass
	data['fname'] = fname
	for k in ['tau',  'G',  'H',  'pop',  'beta',  'U',  'T',  'pop_size',  's',  'pi']:
		if str == type(data[k]):
			data[k] = eval(data[k])  
	return data, W, p, apps

if __name__ == '__main__':
	if len(sys.argv) == 0:
		print "Please give a jobname!"
		sys.exit()
	jobname = sys.argv[1] 
	print "job name:", jobname
	folder = 'output/%s/' % jobname
	print "Working folder:", jobname
	tic = time.clock()
	file_list = glob.glob1(folder, '*.data')
	all_data = [None] * len(file_list)
	apps_data = {}
	all_p = {}
	all_W = {}
	print "Processing", len(file_list), "data files"
	if len(file_list) == 0:
		sys.exit()

	for i,fname in enumerate(file_list):
	    data,W,p,apps = process_data_file(fname)
	    all_data[i] = data
	    if apps != None:
	        apps_data[data['fname']] = apps
	    all_p[data['fname']] = p
	    all_W[data['fname']] = W
	toc = time.clock()
	print "Processed all files in", (toc-tic), "seconds"
	print "all_data length:", len(all_data)
	if len(all_data) == 0:
		sys.exit()

	with gzip.open( "all_data_%s.pickle.gz" % jobname, "w") as f:
		cPickle.dump(all_data,f)
	with gzip.open( "apps_data_%s.pickle.gz" % jobname, "w") as f:
		cPickle.dump(apps_data,f)
	with gzip.open( "all_p_%s.pickle.gz" % jobname, "w") as f:
		cPickle.dump(all_p,f)
	with gzip.open( "all_W_%s.pickle.gz" % jobname, "w") as f:
		cPickle.dump(all_W,f)

	print "apps_data length:", len(apps_data)
	if len(apps_data) == 0:
		sys.exit()

    	difs_data = {}
	difs = []
	for k,v in apps_data.items():
		if len(v):
			difs_data[k] = [v[0]] + [v[i] - v[i-1] for i in range(1,len(v))]
			difs.extend(difs_data[k])
	print np.mean(difs), np.min(difs), np.max(difs)

	df_list = []
	for res in all_data:
		fname = res['fname']
		for i,a in enumerate(apps_data[fname]):
			d = difs_data[fname][i]
			rec = res.copy()
			rec['app'] = a
			rec['dif'] = d
			df_list.append(rec)
	df = pd.DataFrame(data=df_list)
	fname = 'appearances_%s.csv.gz' % jobname
	with gzip.open(fname, 'w') as f:
    		df.to_csv(f)
    	print "Saved appearances to", fname