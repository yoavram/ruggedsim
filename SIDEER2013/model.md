# Model

Following [@Hadany2002].

## The model

There are two haploid loci with alleles *a/A* and *b/B* as well as infinitely other loci with infinitely many alleles. The population evolved in a fitness landscape in which *ab* was the prefered genotype, and we start the dynamics immediately after a new fitness peak has risen for the genotype *AB*, thus the individual fitness landscape has two peaks and a valley between them, and the population occupies the lower peak [@Crow1990]:

1. The frequent type at the begining is *ab*.
1. The relative fitness of *ab* is 1
1. The single mutants *Ab* and *aB* have a lower fitness *1-s* 
1. The double mutant *AB* has a higher fitness *1+sH*

Forward mutation rates in both loci are equal to $\mu_k$ and the genomic mutation rate is $U_k$, where *k* denotes the genotype. Denote by *U* and $\mu$ the basal rates without stress-induction which are also the ones used before the environmental change. There is no recombination. There are also no back mutations, so we assume that individuals that acquired deleterious mutations in loci other that *A/a* and *B/b* are evolutionary dead ends and will not contribute to the adaptation process. Therefore there are three ways to shift from the lower peak to the higher peak - directly, by two mutations in a single generations, and in two steps via *Ab* or *aB*.

### Waiting time for double mutant
In the absence of *AB*, the fractions $f_0$ and $f_1$ of the wild type *ab* and the single mutants (*Ab* and *aB*) should approach the mutation-selection balance frequencies $f_1=\mu/s$ and $f_0=1-U/s$, assuming *s* is the selection coefficient of deleterious mutations. Transition from these genotypes to *ab* due to mutation is $f_0 e^{-U_0} \mu_0^2$ and $2f_1 (1-s)e^{-U_1}\mu_1$, where the $e^{-U_k}$ terms express the load imposed by deleterious mutations. Together these transitions sum up to the probability that a new born individual is a double mutant $q_n$:
$$
q_n = (1-\frac{U}{s})e^{-U_0}\mu_0^2 + 2 \frac{\mu}{s} e^{-U_1}\mu_1 \approx \\\\ 
\frac{1}{s}((s-U)(1-U_0)\mu_0^2+2\mu \mu_1(1-U_1)(1-s))
$$
Denote $\mu_k = \tau_k \mu$ and $U_k = \tau_k U$,
$$
q_n \approx \frac{\mu^2}{s}((s-U)(1-\tau_0 U)\tau_0^2 + 2\tau_1 (1-\tau_1 U)(1-s))
$$
Setting $\tau=\tau_0=\tau_1$:
$$
q_n \approx \frac{\tau\mu^2}{s}(1-\tau U)(\tau(s-U) + 2-2s)
$$
Differentiating by $\tau$:
$$
\frac{d q_n}{d \tau} \approx \frac{\mu^2}{s}( 2+3 \tau U(\tau U -2)+s(\tau U (4-3 \tau)-2))
$$
Now, setting $\tau = 1$ for a non-mutator strain,
$$
\frac{d q_n(\tau=1)}{d \tau} \approx \frac{\mu^2}{s}( 2+3 U(U -2)+s(U-2)) = \\\\
\frac{\mu^2}{s}(2+(3U+s)(U-2)) \approx \\\\
\frac{2 \mu^2}{s}(1-3U+s)
$$
Therefore, if $1-3U+s>0$ then increasing the mutation rates of wild-types and single mutants would increase the appearance rate of the adaptive double mutants.This condition can be formulated as $1+s>3U$ or as a simpler, narrower condition for low values of *s*, $U<\frac{1}{3}$ which is a very reasonable condition for bacteria and yeast.

Denote the population size as *N* and assume that $\frac{s}{\mu} < N << (\frac{s}{\mu})^2$ so that drift due to small population size is unlikely and production of *AB* is not too common. For example, if $\mu=10^{-6}$ and $s=10^{-2}$ then $10^4 < N << 10^8$.

The condition on *N* guarantees that $Nq_n$ is very small. Therefore the probability that an *AB* individual would appear in a population without *AB* , $1-(1-q_n)^N$, can be approximated by $Nq_n$, using the Binomial series expansion:
$$
1-(1-q_n)^N = 1 - (1 + Nq_n + O(Nq_n^2)) = 
Nq_n + O(Nq_n^2) \Rightarrow
1-(1-q_n)^N \approx Nq_n
$$
The waiting time for the appearance of *AB* follows a geometric distribution with expectation:
$$
E[T] \approx 1/Nq_n
$$
Therefore 
$$
E[T] \approx \frac{s}{N\mu^2((s-U)(1-\tau_0 U)\tau_0^2 + 2\tau_1 (1-\tau_1 U)(1-s))}
$$
Or if we set $\tau=\tau_0=\tau_1$:
$$
E[T] \approx \frac{s}{N(1-\tau U)\tau\mu^2((s-U)\tau + 2 (1-s))} = \\\\
\frac{s}{N(1-\tau U)\tau\mu^2(\tau s-\tau U + 2-2s)} = \\\\
\frac{s}{N(1-\tau U)\tau\mu^2(2 + s(\tau-2) -\tau U)}
$$
