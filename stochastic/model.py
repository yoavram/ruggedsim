import numpy as np
import random
from math import floor
from scipy.spatial.distance import cdist, hamming
import cython_load
import model_c

import log
logger = log.get_logger('model')


def create_uniform_mutation_load_population(pop_size, num_classes):
	return np.random.multinomial(pop_size, [1.0 / num_classes] * num_classes)


def create_mutation_free_population(pop_size, num_classes):
	population = np.zeros(num_classes, dtype=np.int)
	population[0] = pop_size
	return population


def create_rates(basic_rate, num_classes):
	return np.array([basic_rate] * num_classes)


def create_muation_rates(mu, num_classes):
	return create_rates(mu, num_classes)


def create_recombination_rates(r, genomes, fitness, s, num_loci):
	return create_rates(r, genomes.shape[0])


def create_mutation_rates_with_modifiers(mu, genomes, fitness, s, num_loci, pi, tau):
	rates = create_rates(mu, genomes.shape[0])
	hypers = fitness < pi
	rates[hypers] *= tau
	return rates


def create_target_genome(num_loci):
	return np.array(np.zeros(num_loci), dtype=np.int)


def hamming_fitness_genome(genome, target_genome, s, H,num_loci):
	load = hamming(genome[:num_loci], target_genome) * target_genome.shape[0]
	fitness = (1 - s) ** load
	if target_genome[:2].sum() == 2:
		genotype = genome[:2].sum()
		if  genotype == 2:
			fitness *= (1 + s * H)
		elif genotype == 0:
			fitness /= (1 - s) ** 2
	return fitness


def hamming_fitness_genomes(genomes, target_genome, s, H, num_loci):
	fitness = np.apply_along_axis(hamming_fitness_genome, 1, genomes, target_genome, s, H, num_loci)
	return fitness


def genome_to_num(genome, num_loci):
	return genome[:num_loci].nonzero()[0]


def genomes_to_nums(genomes, num_loci):
	nums = [genome_to_num(g, num_loci) for g in genomes]
	return nums


def find_row_nums(nums, target):
	# cython version is slightly faster, last time I checked
	for i,n in enumerate(nums):
		if np.array_equal(n, target):
			return i
	return -1


def drift(population):
	pop_size = population.sum()
	p = population / float(pop_size)
	population[:] = np.random.multinomial(pop_size, p)
	return population


def selection(population, fitness):
	pop_size = population.sum()
	p = population * fitness.reshape(population.shape)
	p[:] = p / p.sum()
	population[:] = np.random.multinomial(pop_size, p)
	return population


def mutation(population, mutation_rates, beta):
	### background mutations
	mutations = np.random.poisson(population * mutation_rates)
	total_mutations = mutations.sum()	
	mutations  = np.array((mutations, population)).min(axis=0) # no more than one mutation per individual
	# DEBUG STUFF
	if total_mutations > mutations.sum():
		logger.debug("Reduced %.4f of mutations from %d to %d" % ((1 - mutations.sum() / float(total_mutations)), total_mutations, mutations.sum()))
	# DEBUG END
	for strain in arange(population.shape[0]): 
		strain_mutations = mutations[strain]
		strain_mutations[-1] = 0
		population[strain] -= strain_mutations
		strain_mutations = np.roll(strain_mutations, 1)
		population[strain] += strain_mutations

	### adaptive mutations
	mutations = np.random.poisson(mutation_rates[:3,:] * beta * population[:3,:])
	for strain,load in zip(*mutations.nonzero()):
		#print "strain",strain,"load",load,population[strain,load], mutations[strain,load]
		if strain > 0:
			population[strain,load] -= mutations[strain,load]
			population[3,load] += mutations[strain,load]
			#print population[3,load]
		elif strain == 0:
			# double mutations
			strain_size = population[strain,load]
			individual_mutations = np.random.multinomial(mutations[strain,load], [1./strain_size] * strain_size)
			for individual, num_mutations in enumerate(individual_mutations):
			    	#print "individual",individual,"num_mutations",num_mutations
				if num_mutations == 1:
					new_strain = 1 + (individual % 2)
					population[strain, load] -= 1
					population[new_strain, load] += 1
				elif num_mutations > 1:
					population[strain, load] -= 1
					population[3, load] += 1
			#print  population[3,load]
	return population
