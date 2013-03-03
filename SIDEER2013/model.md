# Model

Following [@Hadany2002].

## The model

There are two haploid loci with alleles *a/A* and *b/B*.
The individual fitness landscape has two peaks and a valley between them, and the population occupies the lower peak [@Crow1990]:

1. The frequent type at the begining is *ab*.
1. The relative fitness of *ab* is 1
1. The single mutants *Ab* and *aB* have a lower fitness *1-s* 
1. The double mutant *AB* has a higher fitness *1+sH*

Forward and backward mutation rates in both loci are equal to $\mu$ and there is no recombination.

### Waiting time for double mutant
Denote the population size as *N* and assume that $\frac{s}{\mu} < N << (\frac{s}{\mu})^2$ so that drift due to small population size is unlikely and production of *AB* is not too common. For example, if $\mu=10^{-6}$ and $s=10^{-2}$ then $10^4 < N << 10^8$.

In the absence of *AB*, the frequencies of *Ab* and *aB* should approach the mutation-selection balance frequency $p=\frac{\mu}{s}$ and the probability that a random new individual is a *AB* individuals is approximated by $q_n$. The condition on *N* guarantees that $Nq_n$ is very small. Therefore the probability that an *AB* individual would appear in a population without *AB* , $1-(1-q_n)^N$, can be approximated by $Nq_n$, using the Binomial series expansion:
$$
1-(1-q_n)^N = 1 - (1 + Nq_n + O(Nq_n^2)) = 
Nq_n + O(Nq_n^2) \Rightarrow
1-(1-q_n)^N \approx Nq_n
$$
The waiting time for the appearance of *AB* follows a geometric distribution with expectation:
$$
E[T] \approx 1/Nq_n
$$
We therefore need to find $q_n$.

The frequency of *ab* at equilibrium is $1-2p$ and that of *Ab* and *aB* is *p* each. The transition from each of these genotypes to *ab* due to mutation is $(1-2p)\mu^2$ and $2p\mu(1-\mu)(1-s)$. Together these transition sum up to $q_n$:
$$
q_n = (1-2p)\mu^2 + 2p\mu(1-\mu)(1-s) = 
\frac{\mu^2}{s}(2+2s\mu-s-4\mu) = 
\frac{\mu^2}{s}(2-s) + O(\mu^3) = 
$$
Therefore 
$$
E[T] \approx \frac{s}{N\mu^2(2-s)}
$$
This is an increasing function of *s* ($\frac{dE[T]}{ds}=\frac{2}{N\mu^2 (s-2)^2}>0$) and a decreasing function of $\mu$ and *N* ($\frac{d E[T] }{d \mu} = \frac{2s}{N(s-2)u^3} < 0, \frac{dE[T]}{dN}=\frac{s}{N^2 \mu^2 (s-2)} < 0$).

Rewriting $q_n$ with two mutation rates, where $\mu_0$ and $\mu_1$ are the mutation rates of the wild type and the single mutants ($p=\mu_0/s$) and assuming $1 >> \mu_1 \ge \mu_0 $:
$$
q_n = (1-2p)\mu_0^2 + 2p\mu_1(1-\mu_1)(1-s) = 
(1-2\mu_0/s)\mu_0^2 + 2\mu_0 \mu_1/s(1-\mu_1)(1-s) = 
\mu_0^2 + 2\mu_0 \mu_1(1-s)/s + O(\mu_1^3) =
\mu_0(\mu_0 + 2 \mu_1 \frac{1-s}{s}) + O(\mu_1^3) 
$$
Nowincreasing the mutation rates decreases the waiting time but the effect of the two rates is different in size: