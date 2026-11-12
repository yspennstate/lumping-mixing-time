# Exact checks for "Lumping and the mixing time of a reversible chain".
#
# All arithmetic is rational (fractions.Fraction).  A chain is given by a symmetric
# matrix of weights W (loops on the diagonal); P = W / row sums.
#
#   python lumping_checks.py data/lumping_checks.json
#
# On random weighted chains with three to six states and random partitions it checks
#   - Proposition 4 (hitting times of sets) for every set of blocks and every start,
#   - tau12(hat X) <= tau0(hat X) <= tau0(X) <= (n-1) tau12(X),
#   - the upper bound of Theorem 1: tau12(hat X) <= min(2^(n-q), q-1) tau12(X), q blocks.
# Theorem 2 (one lumping of two states) and the comb are checked in merge_and_comb_checks.py.
import json
import random
import sys
from fractions import Fraction as F
from itertools import combinations


def solve(A, b):
    """Solve A x = b over the rationals (A square, invertible)."""
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = next(r for r in range(col, n) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        inv = 1 / M[col][col]
        M[col] = [v * inv for v in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [a - f * c for a, c in zip(M[r], M[col])]
    return [M[i][n] for i in range(n)]


def inverse(A):
    n = len(A)
    cols = [solve(A, [F(int(i == j)) for i in range(n)]) for j in range(n)]
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def weights(W):
    w = [sum(row) for row in W]
    tot = sum(w)
    return w, [x / tot for x in w]


def hit_set(W, A):
    """E_x T_A for every x, and the Green function of the chain killed on A."""
    n = len(W)
    w, _ = weights(W)
    out = [x for x in range(n) if x not in A]
    Q = [[F(int(x == y)) - W[x][y] / w[x] for y in out] for x in out]
    G = inverse(Q)
    E = [F(0)] * n
    green = {}
    for a, x in enumerate(out):
        E[x] = sum(G[a])
        green[x] = G[a][a]
    return E, green


def hitting_matrix(W):
    n = len(W)
    H = [[F(0)] * n for _ in range(n)]
    for j in range(n):
        E, _ = hit_set(W, {j})
        for i in range(n):
            H[i][j] = E[i]
    return H


def tau12(W):
    """max_{i,j} (E_i T_j - E_pi T_j), with a maximising pair."""
    n = len(W)
    _, pi = weights(W)
    H = hitting_matrix(W)
    best, arg = None, None
    for j in range(n):
        mean = sum(pi[k] * H[k][j] for k in range(n))
        for i in range(n):
            v = H[i][j] - mean
            if best is None or v > best:
                best, arg = v, (i, j)
    return best, arg, H


def tau0(W, H=None):
    _, pi = weights(W)
    H = H or hitting_matrix(W)
    n = len(W)
    return sum(pi[j] * H[0][j] for j in range(n))


def lump(W, blocks):
    m = len(blocks)
    return [[sum(W[x][y] for x in blocks[a] for y in blocks[b]) for b in range(m)] for a in range(m)]


# ---------------------------------------------------------------------------
# Random chains.

def random_chain(n, rng):
    W = [[F(0)] * n for _ in range(n)]
    order = list(range(n))
    rng.shuffle(order)
    for a, b in zip(order, order[1:]):            # a spanning path keeps the chain irreducible
        v = F(rng.randint(1, 9), rng.randint(1, 4))
        W[a][b] += v
        W[b][a] += v
    for a in range(n):
        for b in range(a, n):
            if rng.random() < 0.5:
                v = F(rng.randint(1, 30), rng.randint(1, 6))
                W[a][b] += v
                if a != b:
                    W[b][a] += v
    return W


def hitting_trials(trials, seed):
    rng = random.Random(seed)
    worst = F(0)
    worst_at = None
    count = 0
    ok = True
    ok_tau0 = True
    ok_upper = True
    for _ in range(trials):
        n = rng.randint(3, 6)
        W = random_chain(n, rng)
        m = rng.randint(2, n - 1)
        while True:
            lab = [rng.randrange(m) for _ in range(n)]
            if len(set(lab)) == m:
                break
        blocks = [[x for x in range(n) if lab[x] == a] for a in range(m)]
        V = lump(W, blocks)
        t_ind, _, Hh = tau12(V)
        t_or, _, H = tau12(W)
        if not (t_ind <= tau0(V, Hh) <= tau0(W, H) <= (n - 1) * t_or):
            ok_tau0 = False
        if not t_ind <= min(2 ** (n - m), m - 1) * t_or:
            ok_upper = False
        for size in range(1, m):
            for hatA in combinations(range(m), size):
                A = set(x for a in hatA for x in blocks[a])
                E, green = hit_set(W, A)
                Eh, _ = hit_set(V, set(hatA))
                s_green = sum(green.values())
                s_hit = sum(E)
                top = max(E)
                for a in range(m):
                    if a in hatA:
                        continue
                    count += 1
                    if not (Eh[a] <= s_green <= s_hit <= (n - len(A)) * top):
                        ok = False
                    q = Eh[a] / top
                    if q > worst:
                        worst, worst_at = q, (n, len(A))
    return {'trials': trials, 'seed': seed, 'inequalities_checked': count, 'all_hold': ok,
            'tau12_hat <= tau0_hat <= tau0 <= (n-1) tau12 in every trial': ok_tau0,
            'tau12_hat <= min(2^(n-q), q-1) tau12 in every trial': ok_upper,
            'largest hat E_a T_A / max_y E_y T_A': float(worst), 'attained_with_(n,|A|)': worst_at}


if __name__ == '__main__':
    out = {'hitting': hitting_trials(400, 20261112)}
    print('hitting', out['hitting'])
    if len(sys.argv) > 1:
        with open(sys.argv[1], 'w') as fh:
            json.dump(out, fh, indent=1)
