from time import clock
from os import makedirs, rename
from os.path import sep, exists, dirname
from datetime import datetime
import json
import gzip

import numpy as np
from numpy import e

from model import *

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
	simulation_id = "pop_%d_G_%d_s_%g_H_%g_U_%g_beta_%g_pi_%g_tau_%g_" % (pop_size,G,s,H,U,beta,pi,tau)
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
	w = smooth_fitness(s, H, 3, G)
	mutation_rates = mutation_rates_matrix(U, pi, tau, w)
	Mm = big_mutation_matrix(mutation_rates, 3, small_background_mutation_matrix)
	Mu = big_mutation_matrix((mutation_rates * beta).transpose(), G, small_strain_mutation_matrix)

	p = mutation_free_population(3, G)
	shape = p.shape
	W = [mean_fitness(p,w)]

	logger.info("Starting simulation with %d ticks", ticks)
	
	fixation_count = 0
	tick = 0
	change_tick = -1
	appearances = []

	while tick < ticks and fixation_count < 2000:		
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
		W += [mean_fitness(p,w)]
	
		# monitoring and logging

		if tick_interval != 0 and tick % tick_interval == 0:
			logger.debug("Tick %d", tick)

		if W[-1] < e ** (-(1 + beta) * U) and change_tick == -1:
			logger.debug("Changing fitness landscape at tick %d with mean fitness %f" % (tick, W[-1]))
			w = rugged_fitness(s, H, 3, G)
			change_tick = tick

		if p[2,0] > 0 and change_tick >= 0:
			if fixation_count == 0:
				logger.debug("Started counting fixation at tick %d" % tick)
				appearances.append(tick)
			fixation_count += 1
		else:
			if fixation_count > 0:
				logger.debug("Stopped counting fixation at tick %d" % tick)
				fixation_count = 0
		tick += 1
	toc = clock()
	logger.info("Simulation finished, %d ticks, time elapsed %.3f seconds",tick, (toc - tic))
	
	# output file
	output_filename = cat_file_path(data_ext)
	make_path(output_filename)
	W = [x for i,x in enumerate(W) if i % stats_interval == 0]
	data = {'G':G, 
		'pop_size':pop_size, 
		's':s, 
		'H':H, 
		'U':U, 
		'beta':beta, 
		'pi':pi, 
		'tau':tau, 
		'p':p.tolist(), 
		'W':W, 
		'ticks':tick, 
		'apps':appearances
	}
	if len(appearances) > 0 and fixation_count > 0:
		data['T'] = appearances[-1]  
	else:
		data['T'] = None
	with gzip.open(output_filename, 'w') as f:
		json.dump(data, f, indent=4, separators=(',', ': '))
	logger.info("Saved output to %s", output_filename)
	
	return p, W, output_filename


if __name__=="__main__":
	p,W,f = run()
