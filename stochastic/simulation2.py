from time import clock
from os import makedirs, rename
from os.path import sep, exists, dirname
from datetime import datetime
import json
import gzip

import numpy as np
from numpy import e

from model import *

VERSION = '3'
# local peaks with distance 2
# int   bin                     fitness H(s=0.05)
#17 (1, 0, 0, 0, 1, 0, 0, 0) 0.973 0.554984583762
#129 (1, 0, 0, 0, 0, 0, 0, 1) 0.896 2.32142857143
#9 (1, 0, 0, 1, 0, 0, 0, 0) 0.83 4.09638554217
#3 (1, 1, 0, 0, 0, 0, 0, 0) 0.83 4.09638554217
#192 (0, 0, 0, 0, 0, 0, 1, 1) 0.945 1.16402116402
wildtype = 9

# utility functions

def make_path(filename):
	path = dirname(filename)
	if not exists(path):
		print("Creating path: %s" % path)
		makedirs(path)
	return exists(path)


def cat_file_path(extension):
	return output_dir + sep + job_name + sep + simulation_id + extension

## Setting up the simulation infrastructure

# load parameters to global namespace
import args, params
args_and_params = args.args_and_params()
globals().update(args_and_params)
if not 'simulation_id' in args_and_params:
	simulation_id = "pop_%d_G_%d_s_%g_U_%g_beta_%g_pi_%g_tau_%g_" % (pop_size,G,s,U,beta,pi,tau)
	simulation_id += datetime.now().strftime('%Y-%b-%d_%H-%M-%S-%f')
	args_and_params['simulation_id'] = simulation_id
params_filename = cat_file_path(params_ext)
make_path(params_filename)
params.save(params_filename, args_and_params)

# load logging
import log
log_filename = cat_file_path(log_ext)
make_path(log_filename)
log.init(log_filename, console, debug)
logger = log.get_logger('simulation')

# log initial stuff
logger.info("Simulation ID: %s", simulation_id)
logger.info("Logging to %s", log_filename)
logger.info("Parametes from file and command line: %s", params.to_string(args_and_params, short=True))
logger.info("Parameters saved to file %s", params_filename)

def run():
	tic = clock()

	# init population
	w = rugged_fitness(s, 256, G, wildtype=wildtype)
	w[0,:] = 0
	print "*** W[%d]: %.3f, W[0]:, %.3f" % (wildtype, w[wildtype,0], w[0,0])
	if SIMe:
		logger.info("SIMe mode - using pi=0 and tau=1 until MSB is reached")
		mutation_rates = mutation_rates_matrix(U, 0, 1, w)
	else:
		mutation_rates = mutation_rates_matrix(U, pi, tau, w)
	Mm = big_mutation_matrix(mutation_rates, 256, small_background_mutation_matrix)
	mutation_rates2 = mutation_rates.copy()
	if unloaded:
		mutation_rates2[:,1:] = 0
	Mu = big_mutation_matrix((mutation_rates2 * beta).transpose(), G, small_strain_mutation_matrix)

	p = mutation_free_population(256, G, start_strain=wildtype)
	shape = p.shape
	W = mean_fitness(p,w)

	logger.info("Starting simulation V.%s", VERSION)
	
	tick = 0

	## MSB

	while tick < 100:
		# selection
		p = w * p
		
		# strain mutations
		p = Mu.dot( p.flatten(order="F") )
		p = p.reshape(shape, order="F")

		# background mutations 
		p = Mm.dot( p.flatten(order="C") )
		p = p.reshape(shape, order="C")

		p /= p.sum()

		# drift
		if pop_size > 0:
			p = np.random.multinomial(pop_size, p.flatten()) / np.float64(pop_size)
			p = p.reshape(shape)

		# mean fitness
		W = mean_fitness(p,w)
	
		# monitoring and logging

		if tick_interval != 0 and tick % tick_interval == 0:
			logger.debug("Tick %d", tick)
		tick += 1

	msb_dict = {'p': p.tolist(), 'W': W, 't': tick}
	logger.info("MSB reached at tick %d with mean fitness %.4g", tick, W)
	w = rugged_fitness(s, 256, G, wildtype=wildtype)
	print "*** W[%d]: %.3f, W[0]:, %.3f" % (wildtype, w[wildtype,0], w[0,0])
	print "*** H =", (w[0,0] - 1)/s
	mutation_rates = mutation_rates_matrix(U, pi, tau, w)
	Mm = big_mutation_matrix(mutation_rates, 256, small_background_mutation_matrix)
	mutation_rates2 = mutation_rates.copy()
	if unloaded:
		raise Warning("Not implemented!")
		mutation_rates2[:,1:] = 0
	Mu = big_mutation_matrix((mutation_rates2 * beta).transpose(), G, small_strain_mutation_matrix)

	## Double mutant appearance

	while p[0,:].sum() == 0:
		# selection
		p = w * p
		
		# strain mutations
		p = Mu.dot( p.flatten(order="F") )
		p = p.reshape(shape, order="F")

		# background mutations 
		p = Mm.dot( p.flatten(order="C") )
		p = p.reshape(shape, order="C")

		p /= p.sum()

		# drift
		if pop_size > 0:
			p = np.random.multinomial(pop_size, p.flatten()) / np.float64(pop_size)
			p = p.reshape(shape)

		# mean fitness
		W = mean_fitness(p,w)
	
		# monitoring and logging

		if tick_interval != 0 and tick % tick_interval == 0:
			logger.debug("Tick %d", tick)
		tick += 1

	app_dict = {'p': p.tolist(), 'W': W, 't': tick}

	logger.info("Double mutant appeared at tick %d with mean fitness %.4g", tick, W)
	AB0,AB1,AB2,AB3 = p[0,0],p[0,1],p[0,2],p[0,3]
	logger.info("AB/0 %.4g, AB/1 %.4g, AB/2 %.4g, AB/3 %.4g", AB0, AB1, AB2, AB3)


	## Double mutant fixation

	while p[0,:].sum() > 0 and p[0,:].sum() < 1:
		# selection
		p = w * p
		
		# NO strain mutations
		# p = Mu.dot( p.flatten(order="F") )
		# p = p.reshape(shape, order="F")

		# background mutations 
		p = Mm.dot( p.flatten(order="C") )
		p = p.reshape(shape, order="C")

		p /= p.sum()

		# drift
		if pop_size > 0:
			p = np.random.multinomial(pop_size, p.flatten()) / np.float64(pop_size)
			p = p.reshape(shape)

		# mean fitness
		W = mean_fitness(p,w)
	
		# monitoring and logging

		if tick_interval != 0 and tick % tick_interval == 0:
			logger.debug("Tick %d", tick)
		tick += 1

	fix_dict = {'p': p.tolist(), 'W': W, 't': tick, 'success': bool(p[0,:].sum() > 0)}
	
	if fix_dict['success']:
		logger.info("Fixation at tick %d with mean fitness %.4g and AB frequency %.4g", tick, W, p[0,:].sum())
		AB0,AB1,AB2,AB3 = p[0,0],p[0,1],p[0,2],p[0,3]
		logger.info("AB/0 %.4g, AB/1 %.4g, AB/2 %.4g, AB/3 %.4g", AB0, AB1, AB2, AB3)
	else:
		logger.info("Extinction at tick %d with mean fitness %.4g", tick, W)

	# wrap up 
	toc = clock()
	logger.info("Simulation V.%s finished, %d ticks, time elapsed %.3f seconds", VERSION, tick, (toc - tic))

	# output file
	output_filename = cat_file_path(data_ext)
	make_path(output_filename)
	data = { 'simulation_id':simulation_id,
		'G':G, 
		'pop_size':pop_size, 
		's':s, 
		'U':U, 
		'beta':beta, 
		'pi':pi, 
		'tau':tau, 
		'msb':msb_dict,
		'app':app_dict,
		'fix':fix_dict 
	}
	with gzip.open(output_filename, 'w') as f:
		json.dump(data, f, indent=4, separators=(',', ': '))
	logger.info("Saved output to %s", output_filename)
	
	return output_filename


if __name__=="__main__":
	fname = run()
