import log
logger = log.get_logger('model')

import numpy as np
from scipy.stats import poisson
from scipy.linalg import block_diag


def mprint(M, precision=3):
    '''mprint(matrix, precision=int) -> None
    pretty prints a matrix with specified precision.
    '''
    for row in M:
        for cell in row:
            print ("%." +str(precision) + "f") % cell,
        print


def mutation_free_population(strains, G):
	p = array([ [0] * G for _ in range(strains) ], dtype=np.float64)
	p[0,0] = 1
	assert np.allclose(p.sum(), 1)
	return p


def uniform_population(strains, G):
	x = 1.0 / (strains * G)
	p =  ones((strains, G), dtype=np.float64) * x
	assert np.allclose(p.sum(), 1)
	return p


def smooth_fitness(s, H, strains, G):
	w = array([ [ (1 - s ) ** (k + i) for k in range(G)] for i in range(strains) ])
	return w


def rugged_fitness(s, H, strains, G):
	w = smooth_fitness(s, H, strains, G)
	w[-1,:] *= (1 + s * H)/((1 - s) ** 2) # fix fitness of the double mutant strain
	return w


def mean_fitness(p, w):
	return (p*w).sum()


def mutation_rates_matrix(U, pi, tau, w):
	mutation_rates = np.ones(w.shape) * U
	mutation_rates[w < pi] *= tau
	return mutation_rates


def big_mutation_matrix(mutation_rates, repeats, small_mutation_matrix_function):
	M = np.zeros((0,0))
	for i in range(repeats):
		m = small_mutation_matrix_function(mutation_rates[i,:])
		M = block_diag(M, m)
	assert allclose(M.sum(axis=0),1)
	return M


def small_background_mutation_matrix(mutation_rates):
	assert mutation_rates.shape[0] == len(mutation_rates)
	mutation_rvs = poisson(mutation_rates)
	m = diag(mutation_rvs.pmf(0))
	for k in range(1,mutation_rates.shape[0]):
		m += np.diag(mutation_rvs.pmf(k)[:-k],-k)
	# absorb further mutations in the last class
	for j in range(mutation_rates.shape[0]):
		m[-1,j] = 1 - mutation_rvs.cdf(mutation_rates.shape[0] - 2 - j)[j]
	return m


def small_strain_mutation_matrix(mutation_rates):
	assert mutation_rates.shape[0] == len(mutation_rates)
	assert mutation_rates.shape[1] == 3
	mu = mutation_rates
	u = array([ [ (1 - mu[0]) ** 2, 0, 0 ], [2 * mu[0] * (1 - mu[0]) , 1 - mu[1], 0], [ mu[0] ** 2, mu[1], 1 ] ])
	return u


def find_fixation_time(p, W, w, U):
	'''a heuristic to calculate the fixation time from the results of a simulation
	'''
	w0 = w[p==p.max()] * e**(-U)
	if w0 < 1:
		return np.infty,0
	W0 = W >= w0
	t0 = W0.argmax()
	return t0, W[t0]