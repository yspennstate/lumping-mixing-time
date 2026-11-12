# Lumping and the mixing time of reversible chains with n states

Yitzchak Shmalo, November 2026.

Lumping states of a reversible Markov chain (gluing vertices of the network together) does not increase the
maximal commute time, the average hitting time, the relaxation time or the Cheeger time constant. Aldous and Fill
(Open Problem 4.45 of *Reversible Markov Chains and Random Walks on Graphs*) asked whether the stopping-rule mixing
time tau_1^(2) of the lumped chain is at most an absolute constant times that of the original chain. Hermon and
Peres (Electron. J. Probab. 23, 2018, paper 25) answered in the negative: for random walks on bounded-degree graphs
the factor can be of order log |V|, and no more.

This paper treats general reversible chains with n states. Let K(n) be the largest possible ratio of tau_1^(2) of a
lumped chain to tau_1^(2) of the original chain. The paper proves

- floor(n/2) <= K(n) <= max_q min(2^(n-q), q-1) <= max(1, n-2); hence K(n) = floor(n/2) for n <= 6 and n = 8,
  and 3 <= K(7) <= 4;
- lumping two states multiplies tau_1^(2) by less than 2, and 2 is the best constant;
- tau_2 <= tau_1^(2) for every reversible chain (relaxation time against stopping-rule mixing time, with
  constant one);
- lumping increases the mean hitting time of a set A by a factor of at most n - |A|, and a factor close to
  (n - |A|)/2 occurs.

The lower bounds come from a random walk on a weighted comb with 2N+4 states, for which the ratio tends to N+2.
If r successive lumpings of two states can multiply tau_1^(2) by at most r+1, then K(n) = floor(n/2) for all n;
this is open.

## Files

- `lumping.pdf`, `lumping.tex`, `tail.tex`, `notation.tex`, `bib_main.tex`, `bib_own.tex`: the paper.
- `code/merge_and_comb_checks.py`: exact rational checks of the formula for the lumped hitting-time matrix and of
  the two-state bound on 850 lumpings of 100 random networks, of the closed forms for the comb (N up to 3,
  t from 2 to 1000), and the table of bounds on K(n) (`python code/merge_and_comb_checks.py`; standard library only).
- `data/merge_and_comb_checks_out.txt`: its output.
- `code/lumping_checks.py`: exact rational checks of the hitting-time bound, of the chain of inequalities behind
  the upper bound and of the bound min(2^(n-q), q-1) on 400 random networks with random partitions
  (`python code/lumping_checks.py data/lumping_checks.json`).
- `data/lumping_checks.json`: its output.

An archived copy of the paper is at https://doi.org/10.5281/zenodo.24147769.

The paper is part of the Math Brain project: https://github.com/yspennstate/math-brain-public
