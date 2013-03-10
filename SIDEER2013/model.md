# Adaptation with stress-induced mutation

Following [@Hadany2003b].

There are two haploid loci with alleles *a/A* and *b/B* as well as infinitely other loci with infinitely many alleles. The population evolved in a fitness landscape in which *ab* was the prefered genotype, and we start the dynamics immediately after a new fitness peak has risen for the genotype *AB*, thus the individual fitness landscape has two peaks and a valley between them, and the population occupies the lower peak [@Crow1990]:

1. The frequent type at the begining is *ab*.
1. The relative fitness of *ab* is 1
1. The single mutants *Ab* and *aB* have a lower fitness *1-s* 
1. The double mutant *AB* has a higher fitness *1+sH*
1. Mutations in other loci reduce the fitness in a multiplicative manner by *1-s*

Mutation rates in the focus loci are equal to $\mu_k$ and the genomic mutation rate is $U_k$, where *k* denotes the genotype. *U* and $\mu$ are the basal mutation rates used before the environmental change. There is no recombination. There are also no back mutations, so we assume that individuals that acquired deleterious mutations in loci other that *A/a* and *B/b* are evolutionary dead ends and will not contribute to the adaptation process. Therefore there are three ways to shift from the lower peak to the higher peak - directly, by two mutations in a single generations, and in two steps via *Ab* or *aB*.

## Waiting time for double mutant

### Probability of appearance

In the absence of *AB*, the fractions $f_{ab}$, $f_{Ab}$ and $f_{Ab}$ of the wild type *ab* and the single mutants (*Ab* and *aB*) should approach the mutation-selection balance frequencies $f_{aB}=f_{Ab}=\mu/s$ and $f_ab=1-U/s$, assuming *s* is the selection coefficient of deleterious mutations. Transition from these genotypes to the double mutant *ab* due to mutation is $f_{ab} e^{-U_{ab}} \mu_{ab}^2$ and $2f_{Ab}e^{-U_{Ab}}\mu_{Ab}$, where the $e^{-U_k}$ terms are the Poisson probabilities that no deleterious mutations occured and express the load imposed by deleterious mutations. Together these transitions sum up to the probability that a new born individual is a double mutant $q_n$, assuming $U_{Ab}=U_{aB}, \mu_{Ab}=\mu_{aB}$:
$$
q_n = (1-\frac{U}{s})e^{-U_{ab}}\mu_{ab}^2 + 2 \frac{\mu}{s} e^{-U_{Ab}}\mu_{Ab}
$$

Denote the population size as *N* and assume that $\frac{s}{\mu} < N << (\frac{s}{\mu})^2$ so that drift due to small population size is unlikely and production of *AB* is not too common. For example, if $\mu=10^{-6}$ and $s=10^{-2}$ then $10^4 < N << 10^8$.

The condition on *N* guarantees that $Nq_n$ is very small. Therefore the probability that an *AB* individual would appear in a population without *AB* , $1-(1-q_n)^N$, can be approximated by $Nq_n$, using the Binomial series expansion:
$$
1-(1-q_n)^N = 1 - (1 + Nq_n + O(Nq_n^2)) = 
Nq_n + O(Nq_n^2) \Rightarrow
1-(1-q_n)^N \approx Nq_n
$$
The waiting time for the appearance of the double mutant *AB* follows a geometric distribution with expectation:
$$
E[T_1] \approx 1/Nq_n
$$

### Constant mutation rate

We define $\tau$ to be the multiplicate increase of the mutation rate by a mutator allele and denote $U_k=\tau U$ and $\mu_k=\tau \mu$ so that the mutation rate is constant:

$$
q_n = (1-\frac{\tau U}{s})e^{-\tau U}(\tau \mu)^2 + 2 \frac{\tau \mu}{s} e^{-\tau U}\tau \mu = \\\\
e^{-\tau U}\frac{(\tau \mu)^2}{s}(2+s-\tau U)
$$

Deriving this with respect to $\tau$, we get
$$
\frac{d q_n}{d \tau} = e^{-\tau U} \frac{\tau \mu^2}{s} (\tau U (\tau U -(s+5))+2s +4)
$$
Disregarding positive terms and approximating $s+5 \approx 5$ and $2s+4\approx4$, we get
$$
sign(\frac{d q_n}{d \tau} ) = sign(\tau U (\tau U -(s+5))+2s +4) \approx \\\\
sign(\tau U (\tau U -5)+4) = \\\\
sign( (\tau U)^2 -5 \tau U + 4)
$$
The last expression is positive for $\tau <1/U$ and $\tau > 4/U$. The maximum of $q_n$ is attained in the first root of the derivative, where $\tau = 1/U$:
$$
\max{q_n} = \frac{\mu^2}{U^2 s}e^{-1}(1+s)
$$
For example, for $U = 0.003, \mu=10^{-7}, s=0.01$, the probability of apperance of a double mutant in a mutator population is roughly $q_n = 4.13 \cdot 10^{-8}$.

### Stress-induced mutation

This is the case in which hypermutation is induced in individuals who are below **both** fitness peaks, that is, with fitness $\le 1-s$.

Denote $\mu_{ab} = \mu, U_{ab} = U$, $\mu_{Ab}= \tau \mu$ and $U_{Ab} = \tau U$, the probability of a appearance of a double mutant is:
$$
q_n = (1-\frac{U}{s})e^{-U}\mu^2 + 2 \frac{\mu}{s} e^{-\tau U}\tau \mu = \\\\
\frac{\mu^2}{s}e^{-U}(2 \tau e^{-U(\tau-1)} + s - U)
$$

Differentiating by $\tau$:
$$
\frac{d q_n}{d \tau} = 2\frac{\mu^2}{s}e^{-U\tau}(1-\tau U).
$$
Therefore the optimal mutation rate increase $\tau$ for decreasing the waiting time for the appearance of the double mutant $E[T_1]$ is $\tau=1/U$ which gives
$$
\max{q_n} = \frac{\mu^2}{s}e^{-U}(\frac{2}{U} e^{-1+U} +s -U)
$$

To take the same example as before, setting $U = 0.003, \mu=10^{-7}, s=0.01$ we get that $q_n = 2.45 \cdot 10^{-10}$.

Comparing the maximum appearance probability with constant mutation rate and stress-induced mutation rate,

$$
\frac{\frac{\mu^2}{U^2 s}e^{-1}(1+s)}{\frac{\mu^2}{s}e^{-U}(\frac{2}{U} e^{-1+U} +s -U)} = \\\\
\frac{1+s}{2U+s U^2 e^{1-U}+Ue^{1-U}} = \\\\
\frac{1+s}{2U+O(U^2)}
$$
Note that this advantage of CM depends on two components: higher mutation rate in *ab* individuals (hence the dependence on *s*), and greater genetic variation in the population at the mutation-selection balance which gives a higher frequency of single mutants *Ab* and *aB* (hence the negative dependence on *U*).

## Fixation time of the double mutant

Because the double mutant is at the global fitness peak (*H>0*), its fixation, if it doesn't go to extinction first, is a steady state and the population will converge toa new mutation-selection balance around the global fitness peak. Therefore, the fixation probability $\pi$ can be calculated as $\pi = 1-\epsilon$, where $\epsilon$ is the extinction probability.

Following [@Eshel1981], and assuming that the number of progeny of the double mutant follows a Poisson distribution with a mean of $\alpha$, we get that
$$
\pi = 2\frac{\alpha-1}{\alpha} + o(\alpha-1)
$$
Here, $\alpha$ is just the relative fitness of the double mutant divided by the population mean fitness at the time of its appearance:
$$
\alpha = \frac{ (1+sH) e^{-U_{AB}} }{ \bar{\omega} } 
$$
Since we are interested in the extinction probability, we are dealing with a very low frrequency of *AB* and therefore it doesn't contribute much to the mean fitness. The population mean fitness is that of the mutation selection balance before the adaptation appeared. 

### Stress-induced mutation

The fraction of individuals without deleterious mutations is $e^{-U}\approx 1-U$, but the fraction of individuals with *k* deleterious mtuations, including *Ab* and *aB*, is not simply $U^k e^{-U} /k!$, as shown in classical work [@Kimura1966; @Gordo2005], because the mutation rate is higher in individuals with deleterious mutations. 

We calculate the mean fitness by dividing the population first to mutation free individuals and individuals with one or more deleterious mutations. The fractions of these classes are $e^{-U/s}$ and $1-e^{-U/s}$. The number of additional deleterious mutations beyond the first is then Poisson distributed with a mean of $\tau U/s$:
$$
\bar{\omega} = e^{-U/s} + (1-e^{-U/s})\sum_{k \ge 0} {(1-s)^{k+1} e^{-\tau U/s}(\tau U/s)^k/k!} = \\\\
e^{-U/s} + (1-e^{-U/s})e^{-\tau U/s}(1-s)\sum_{k \ge 0} { ((1-s) \tau U/s)^k/k!} = \\\\
e^{-U/s} + (1-e^{-U/s})(1-s)e^{-s \tau U/s} \Rightarrow \\\\
\bar{\omega} = e^{-U/s} + (1-e^{-U/s})(1-s)e^{-\tau U}
$$

We plug that in the equation for $\alpha$, and assume that $U_{AB} = U$:
$$
\alpha = \frac{ (1+sH) e^{-U} }{ e^{-U} + (1-e^{-U})(1-s)e^{-s \tau U}} = \\\\
\frac{ 1+sH  }{ 1 + (1-e^{-U})(1-s)e^{(1-s \tau) U}} = \\\\
\frac{ 1+sH  }{ 1 + (1-e^{-U})(1-s)e^{(1-s \tau) U}}
$$

For the fixation probability to be >0 we need $\alpha>1$, so we need 
$$
H > \frac{1-s}{s}(1-e^{-U})e^{(1-s\tau)U}
$$
Now, subtituting $\alpha$ in the formula for the fixation probability:
$$
\pi \approx 2\frac{\frac{ 1+sH  }{ 1 + (1-e^{-U})(1-s)e^{(1-s \tau) U}}-1}{\frac{ 1+sH  }{ 1 + (1-e^{-U})(1-s)e^{(1-s \tau) U}}} = \\\\
2\frac{ 1+sH - 1 - (1-e^{-U})(1-s)e^{(1-s \tau) U} }{ 1+sH } = \\\\
2\frac{ sH - (1-e^{-U})(1-s)e^{(1-s \tau) U} }{ 1+sH }
$$
Differentiating this by $\tau$, we get
$$
\frac{d \pi}{d \tau} = \frac{2(1-s)s(e^U-1)Ue^{-s \tau U}}{1+s H}
$$
This is always positive (since $1>s>0, U>0$), and therefore the fixation probability always increases with $\tau$, although with diminishing returns, as $e^{-s\tau U}$ decreases when $\tau$ increases.

Setting $\tau=1/U$, the optimal mutation rates increase found in the previous section, we get 
$$
\pi \approx 2\frac{ sH - (1-e^{-U})(1-s)e^{U-s} }{ 1+sH }
$$
For our test case of $U=0.003, s=0.01$ and choosing $H=10$, this is equal to 0.1765.

## References
