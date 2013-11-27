from time import clock
from os import makedirs, rename
from os.path import sep, exists, dirname
from datetime import datetime
import json
import gzip

import numpy as np
from numpy import e

from model import *

VERSION = '3.0'

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
import args3, params
args_and_params = args3.args_and_params()
globals().update(args_and_params)
if not 'simulation_id' in args_and_params:
	simulation_id = "pop_%d_G_%d_s_%g_H_%g_U_%g_beta_%g_pi1_%g_tau1_%g_pi2_%g_tau2_%g_" % (pop_size,G,s,H,U,beta,pi1,tau1,pi2,tau2)
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
	
	# resident
	mutation_rates = mutation_rates_matrix(U, pi1, tau1, w)
	Mm1 = big_mutation_matrix(mutation_rates, 3, small_background_mutation_matrix)
	mutation_rates2 = mutation_rates.copy()
	if unloaded:
		mutation_rates2[:,1:] = 0
	Mu1 = big_mutation_matrix((mutation_rates2 * beta).transpose(), G, small_strain_mutation_matrix)

	p1 = mutation_free_population(3, G) 

	# invader
	mutation_rates = mutation_rates_matrix(U, pi2, tau2, w)
	Mm2 = big_mutation_matrix(mutation_rates, 3, small_background_mutation_matrix)
	mutation_rates2 = mutation_rates.copy()
	if unloaded:
		mutation_rates2[:,1:] = 0
	Mu2 = big_mutation_matrix((mutation_rates2 * beta).transpose(), G, small_strain_mutation_matrix)
	
	p2 = mutation_free_population(3, G) 

	# go on...
	shape = p1.shape
	W1 = mean_fitness(p1,w)
	W2 = mean_fitness(p2,w)
	W = W1 + W2

	logger.info("Starting simulation V.%s", VERSION)
	
	tick = 0

	## MSB

	while tick < 5000:
		# selection
		p1 = w * p1
		p2 = w * p2
		
		# strain mutations
		p1 = Mu1.dot( p1.flatten(order="F") )
		p1 = p1.reshape(shape, order="F")
		p2 = Mu2.dot( p2.flatten(order="F") )
		p2 = p2.reshape(shape, order="F")

		# background mutations 
		p1 = Mm1.dot( p1.flatten(order="C") )
		p1 = p1.reshape(shape, order="C")
		p2 = Mm2.dot( p2.flatten(order="C") )
		p2 = p2.reshape(shape, order="C")
		
		p1 /= p1.sum()
		p2 /= p2.sum()

		# drift
		if pop_size > 0:
			if p1.sum() > 0:
				p1 = np.random.multinomial(int(p1.sum()*pop_size), p1.flatten()/p1.sum()) / np.float64(pop_size)
				p1 = p1.reshape(shape)
			if p2.sum() > 0:
				p2 = np.random.multinomial(int(p2.sum()*pop_size), p2.flatten()/p2.sum()) / np.float64(pop_size)
				p2 = p2.reshape(shape)

		# mean fitness
		W1 = mean_fitness(p1,w)
		W2 = mean_fitness(p2,w)
		W = W1 + W2
	
		# monitoring and logging

		if tick_interval != 0 and tick % tick_interval == 0:
			logger.debug("Tick %d", tick)
		tick += 1


	msb_dict = {'p1': p1.tolist(), 'p2': p2.tolist(), 'W': W, 't': tick}
	logger.info("MSB reached at tick %d with mean fitness %.4g", tick, W)
	logger.info("Mixing resident and invader (%.2f) and changing the fitness landscape", invasion_rate)
	w = rugged_fitness(s, H, 3, G)

	# resident
	mutation_rates = mutation_rates_matrix(U, pi1, tau1, w)
	Mm1 = big_mutation_matrix(mutation_rates, 3, small_background_mutation_matrix)
	mutation_rates2 = mutation_rates.copy()
	if unloaded:
		mutation_rates2[:,1:] = 0
	Mu1 = big_mutation_matrix((mutation_rates2 * beta).transpose(), G, small_strain_mutation_matrix)

	p1 *= (1 - invasion_rate)

	# invader
	mutation_rates = mutation_rates_matrix(U, pi2, tau2, w)
	Mm2 = big_mutation_matrix(mutation_rates, 3, small_background_mutation_matrix)
	mutation_rates2 = mutation_rates.copy()
	if unloaded:
		mutation_rates2[:,1:] = 0
	Mu2 = big_mutation_matrix((mutation_rates2 * beta).transpose(), G, small_strain_mutation_matrix)
	
	p2 *= invasion_rate

	## Double mutant appearance

	while p1[2,:].sum() == 0 and p2[2,:].sum() == 0:
		# selection
		p1 = w * p1
		p2 = w * p2
		
		# strain mutations
		p1 = Mu1.dot( p1.flatten(order="F") )
		p1 = p1.reshape(shape, order="F")
		p2 = Mu2.dot( p2.flatten(order="F") )
		p2 = p2.reshape(shape, order="F")

		# background mutations 
		p1 = Mm1.dot( p1.flatten(order="C") )
		p1 = p1.reshape(shape, order="C")
		p2 = Mm2.dot( p2.flatten(order="C") )
		p2 = p2.reshape(shape, order="C")

		psum = p1.sum() + p2.sum()
		p1 /= psum
		p2 /= psum

		# drift
		if pop_size > 0:
			if p1.sum() > 0:
				p1 = np.random.multinomial(int(p1.sum()*pop_size), p1.flatten()/p1.sum()) / np.float64(pop_size)
				p1 = p1.reshape(shape)
			if p2.sum() > 0:
				p2 = np.random.multinomial(int(p2.sum()*pop_size), p2.flatten()/p2.sum()) / np.float64(pop_size)
				p2 = p2.reshape(shape)

		# mean fitness
		W1 = mean_fitness(p1,w)
		W2 = mean_fitness(p2,w)
		W = W1 + W2
	
		# monitoring and logging

		if tick_interval != 0 and tick % tick_interval == 0:
			logger.debug("Tick %d", tick)
		tick += 1

	app_dict = {'p1': p1.tolist(), 'p2': p2.tolist(), 'W': W, 't': tick}

	logger.info("Double mutant appeared at tick %d with mean fitness %.4g", tick, W)
	AB0,AB1,AB2,AB3 = p1[2,0]+p2[2,0],p1[2,1]+p2[2,1],p1[2,2]+p2[2,2],p1[2,3]+p2[2,3]
	logger.info("AB/0 %.4g, AB/1 %.4g, AB/2 %.4g, AB/3 %.4g", AB0, AB1, AB2, AB3)

	logger.info("Resident: %.2f, Invader: %.2f", p1.sum(), p2.sum())

	## Double mutant fixation

	while (p1[2,:].sum() > 0 and p1[2,:].sum() < 1) or (p2[2,:].sum() > 0 and p2[2,:].sum() < 1):
		# selection
		p1 = w * p1
		p2 = w * p2
		
		# NO strain mutations
		# p = Mu.dot( p.flatten(order="F") )
		# p = p.reshape(shape, order="F")

		# background mutations 
		p1 = Mm1.dot( p1.flatten(order="C") )
		p1 = p1.reshape(shape, order="C")
		p2 = Mm2.dot( p2.flatten(order="C") )
		p2 = p2.reshape(shape, order="C")

		psum = p1.sum() + p2.sum()
		p1 /= psum
		p2 /= psum

		# drift
		if pop_size > 0:
			if p1.sum() > 0:
				p1 = np.random.multinomial(int(p1.sum()*pop_size), p1.flatten()/p1.sum()) / np.float64(pop_size)
				p1 = p1.reshape(shape)
			if p2.sum() > 0:
				p2 = np.random.multinomial(int(p2.sum()*pop_size), p2.flatten()/p2.sum()) / np.float64(pop_size)
				p2 = p2.reshape(shape)

		# mean fitness
		W1 = mean_fitness(p1,w)
		W2 = mean_fitness(p2,w)
		W = W1 + W2
	
		# monitoring and logging

		if tick_interval != 0 and tick % tick_interval == 0:
			logger.debug("Tick %d", tick)
		tick += 1

	fix_dict = {'p1': p1.tolist(), 
				'p2': p2.tolist(), 
				'W': W, 
				't': tick, 
				'success': bool(p1[2,:].sum()+p2[2,:].sum() > 0), 
				'invader': p2.sum()
				}
	
	if fix_dict['success']:
		logger.info("Fixation at tick %d with mean fitness %.4g and AB frequency %.4g", tick, W, p1[2,:].sum()+p2[2,:].sum())
		AB0,AB1,AB2,AB3 = p1[2,0]+p2[2,0],p1[2,1]+p2[2,1],p1[2,2]+p2[2,2],p1[2,3]+p2[2,3]
		logger.info("AB/0 %.4g, AB/1 %.4g, AB/2 %.4g, AB/3 %.4g", AB0, AB1, AB2, AB3)
	else:
		logger.info("Extinction at tick %d with mean fitness %.4g", tick, W)

	logger.info("Invader frequency %.4g", p2.sum())
	
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
		'H':H, 
		'U':U, 
		'beta':beta, 
		'pi1':pi1, 
		'tau1':tau1, 
		'pi2':pi2, 
		'tau2':tau2, 
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
