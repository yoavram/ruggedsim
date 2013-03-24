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


def mutation(population, genomes, mutation_rates, num_loci, target_genome, nums, beta):
	mutations = np.random.binomial(n=population, p=mutation_rates, size=population.shape)
	total_mutations = mutations.sum()	
	mutations  = np.array((mutations, population)).min(axis=0) # no more than one mutation per individual
	# DEBUG STUFF
	if total_mutations > mutations.sum():
		logger.debug("Reduced %.4f of mutations from %d to %d" % ((1 - mutations.sum() / float(total_mutations)), total_mutations, mutations.sum()))
	# DEBUG END
	mutations_cumsum = mutations.cumsum()
	total_mutations = mutations_cumsum[-1]
	loci = np.random.randint(0, num_loci, total_mutations)
	loci_split = np.split(loci, mutations.cumsum())[:-1] # split by strain
	new_counts = {}
	new_genomes = {}

	for strain in range(population.shape[0]):
		population[strain] = population[strain] - mutations[strain]
		assert population[strain] >= 0  # ASSERT

		_loci = loci_split[strain]	
		if len(_loci) == 0:
			continue
		allele_change = np.random.binomial(1, 1 - beta, len(_loci))
		new_alleles = (target_genome[_loci] + allele_change) % 2 

		for i, locus in enumerate(_loci):
			new_allele = new_alleles[i]
			key = (strain, locus, new_allele)
			if key in new_counts:
				new_counts[key] += 1
			else:
				new_genome = genomes[strain, :].copy()				
				new_genome[locus] = new_allele
				new_counts[key] = 1
				new_genomes[key] = new_genome

	if len(new_genomes) > 0:
		for key, new_genome in new_genomes.items():
			n = genome_to_num(new_genome, num_loci)
			index = model_c.find_row_nums(nums, n)
			if index != -1:
				new_genomes.pop(key)
				population[index] += new_counts.pop(key)

	if len(new_genomes) > 0:
		population = np.append(population, new_counts.values())
		genomes = np.vstack((genomes, new_genomes.values()))

	return population, genomes
