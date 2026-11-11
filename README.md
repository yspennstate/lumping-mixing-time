# Lumping and the mixing time of a reversible chain: the dependence on the number of states

Yitzchak Shmalo, November 2026.

Lumping states of a reversible Markov chain (shorting vertices of the network together) does not increase the
maximal commute time, the average hitting time, the relaxation time or the Cheeger time constant. Aldous and Fill
(Open Problem 4.45 of *Reversible Markov Chains and Random Walks on Graphs*) asked whether the stopping-rule mixing
time tau_1^(2) of the lumped chain is at most an absolute constant times that of the original chain. Hermon and
Peres (Electron. J. Probab. 23, 2018) answered in the negative: for random walks on bounded-degree graphs the factor
can be of order log |V|, and no more.

This note treats general reversible chains with n states. Let K(n) be the largest possible ratio of tau_1^(2) of a
lumped chain to tau_1^(2) of the original chain. The paper proves

- floor((n-2)/2) <= K(n) <= n-1;
- tau_2 <= tau_1^(2) for every reversible chain (relaxation time against stopping-rule mixing time, with
  constant one), by a short martingale argument;
- lumping increases the mean hitting time of a set A by a factor of at most n - |A|, and a factor close to
  (n-2)/2 occurs.

The lower bounds come from a random walk on a weighted comb with 2N+2 states, for which the ratio tends to N.

## Files

- `lumping.pdf`, `lumping.tex`, `tail.tex`, `bib_main.tex`, `bib_own.tex`: the paper.
- `code/lumping_checks.py`: exact rational checks of every inequality used for the comb (N up to 8), of the chain of
  inequalities behind the upper bound and of the hitting-time bound on 400 random networks, and of the four-state
  example (`python code/lumping_checks.py data/lumping_checks.json`; standard library only).
- `data/lumping_checks.json`: its output.

An archived copy of the paper is at https://doi.org/10.5281/zenodo.24147769.

The paper is part of the Math Brain project: https://github.com/yspennstate/math-brain-public
