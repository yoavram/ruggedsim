import log
logger = log.get_logger('model')

import numpy as np
from numpy import e
from scipy.stats import poisson
from scipy.linalg import block_diag
import pickle

FITNESS_LANDSCAPE_FILENAME = "Franke2011_landscape.dict"

def mprint(M, precision=3):
    '''mprint(matrix, precision=int) -> None
    pretty prints a matrix with specified precision.
    '''
    for row in M:
        for cell in row:
            print ("%." +str(precision) + "f") % cell,
        print


def strain2int(strain):    
    return sum([g*2**i for i,g in enumerate(strain)])


def mutation_free_population(strains, G, start_strain):
	p = np.array([ [0] * G for _ in range(strains) ], dtype=np.float64)
	p[start_strain,0] = 1
	assert np.allclose(p.sum(), 1)
	return p


def msb_population(strains, G, start_strain, U, s):
	p = np.array([ [0] * G for _ in range(strains) ], dtype=np.float64)
	p[start_strain,:] = poisson(U/s).pmf(range(G))
	p /= p.sum()
	assert np.allclose(p.sum(), 1)
	return p


def uniform_population(strains, G):
	x = 1.0 / (strains * G)
	p =  ones((strains, G), dtype=np.float64) * x
	assert np.allclose(p.sum(), 1)
	return p


def smooth_fitness(s, strains, G):
	w = np.array([ [ (1 - s ) ** k for k in range(G)] for i in range(strains) ])	
	return w


def load_fitness(filename=FITNESS_LANDSCAPE_FILENAME, intkey=True, cache=None):
	if cache == None:
		with open(filename) as f:
			cache = pickle.load(f)
	if intkey:
		cache2 = {}
		for strain,w in cache.items():
			strain = strain2int(strain)
			cache2[strain] = w
		return cache2
	else:
		return cache


def rugged_fitness(s, strains, G, wildtype):
	fitness = load_fitness()	
	w = smooth_fitness(s, strains, G)
	for strain in range(strains):
		w[strain,:] *= fitness.get(strain, 0)
	w /= w[wildtype,0]
	return w


def mean_fitness(p, w):
	return (p*w).sum()


def hamming(x,y):
    assert len(x) == len(y)
    return sum([1 for i in range(len(x)) if x[i] != y[i]])


def mutation_rates_matrix(U, pi, tau, w):
	if pi < 0:
		return mutation_rates_matrix_simk(U, abs(pi), tau, w)
	mutation_rates = np.ones(w.shape) * U
	mutation_rates[w < pi] *= tau
	return mutation_rates


def mutation_rates_matrix_simk(U, k, tau, w):
	wk = w**k
	wk[wk>1] = 1
	return tau * U - U * (tau - 1) * wk


def big_mutation_matrix(mutation_rates, repeats, small_mutation_matrix_function):
	assert mutation_rates.shape[0] == 256 or mutation_rates.shape[1] == 256
	M = np.zeros((0,0))
	for i in range(repeats):
		m = small_mutation_matrix_function(mutation_rates[i,:])
		M = block_diag(M, m)
	assert np.allclose(M.sum(axis=0),1)
	return M


def small_background_mutation_matrix(mutation_rates):
	assert mutation_rates.shape[0] == len(mutation_rates)	
	mutation_rvs = poisson(mutation_rates)
	m = np.diag(mutation_rvs.pmf(0))
	for k in range(1,mutation_rates.shape[0]):
		m += np.diag(mutation_rvs.pmf(k)[:-k],-k)
	# absorb further mutations in the last class
	for j in range(mutation_rates.shape[0]):
		m[-1,j] = 1 - mutation_rvs.cdf(mutation_rates.shape[0] - 2 - j)[j]
	return m


def small_strain_mutation_matrix(mutation_rates):
	assert mutation_rates.shape[0] == len(mutation_rates)
	mu = mutation_rates
	fitness = load_fitness(intkey=False)
	strains = fitness.keys()
	m = len(strains)
	M = np.diag([1.0] * m)
	for i in range(m):
	    for j in range(m):
	        gi,gj = strains[i],strains[j]        
	        if hamming(gi,gj) == 1:
	            M[j,i] = mu[i]
	            M[i,i] -= mu[i]
	assert np.allclose(M.sum(1),1)	
	return M
	#u = np.array([ [ (1 - mu[0]) ** 2, 0, 0 ], [2 * mu[0] * (1 - mu[0]) , 1 - mu[1], 0], [ mu[0] ** 2, mu[1], 1 ] ])

def find_fixation_time(p, W, w, U):
	'''a heuristic to calculate the fixation time from the results of a simulation
	'''
	w0 = w[p==p.max()] * e**(-U)
	if w0 < 1:
		return np.infty,0
	W0 = W >= w0
	t0 = W0.argmax()
	return t0, W[t0]