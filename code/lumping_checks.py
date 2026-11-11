# Exact checks for "Lumping and the mixing time of a reversible chain".
#
# All arithmetic is rational (fractions.Fraction).  A chain is given by a symmetric
# matrix of weights W (loops on the diagonal); P = W / row sums.
#
#   python lumping_checks.py data/lumping_checks.json
#
# Sections of the output:
#   comb      the comb chain: both sides of every inequality in the proof of the lower bound
#   hitting   the inequalities for hitting times of sets, and the chain of inequalities of the
#             upper bound, on random weighted chains and random partitions
#   small     a four-state chain lumped to three states
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
# The comb chain.  States z_0..z_N, y_0..y_{N-1}, c.

def comb(N, alpha, h, r, ell):
    names = ['z%d' % k for k in range(N + 1)] + ['y%d' % k for k in range(N)] + ['c']
    ix = {s: i for i, s in enumerate(names)}
    n = len(names)
    W = [[F(0)] * n for _ in range(n)]

    def edge(a, b, wt):
        W[ix[a]][ix[b]] += wt
        W[ix[b]][ix[a]] += wt

    for k in range(N):
        edge('z%d' % k, 'z%d' % (k + 1), r ** k)
        e = alpha * r ** k
        edge('y%d' % k, 'z%d' % (k + 1), e)
        W[ix['y%d' % k]][ix['y%d' % k]] += (h - 1) * e
    edge('z%d' % N, 'c', r ** N)
    W[ix['c']][ix['c']] += ell
    blocks = [[ix['z%d' % k], ix['y%d' % k]] for k in range(N)] + [[ix['z%d' % N]], [ix['c']]]
    return W, blocks, ix


def comb_row(N, t, heavy=True):
    alpha, h, r = F(t), F(t), F(1 + t ** 3)
    W0, _, ix = comb(N, alpha, h, r, F(0))
    rest = sum(sum(row) for row in W0) - r ** N      # total weight of the states other than c
    if heavy:
        ell = (N + 2) * t * t * rest - r ** N        # then w_c = (N+2) t^2 (w - w_c)
    else:
        ell = F(1)
    W, blocks, ix = comb(N, alpha, h, r, ell)
    n = len(W)
    w, pi = weights(W)
    delta = (2 + alpha + alpha * h) / (r - 1)
    c = ix['c']
    rho = w[c] / r ** N
    t_orig, arg_o, H = tau12(W)
    V = lump(W, blocks)
    t_ind, arg_i, Hh = tau12(V)
    max_to_c = max(H[i][c] for i in range(n))
    M = h + (N + 1) * (1 + delta)
    B = h + (N + 2) * (1 + delta)
    general = max(B, M + B - rho)
    low_hit = N * h * alpha / (1 + alpha)
    low_tau = pi[c] * low_hit
    hat_hit = Hh[0][len(V) - 1]
    checks = {
        'max_i E_i T_c <= h + (N+1)(1+delta)': max_to_c <= M,
        'tau12 <= max(B, M + B - rho)': t_orig <= general,
        'hat E_{b0} T_c >= N h alpha/(1+alpha)': hat_hit >= low_hit,
        'hat tau12 >= pi_c N h alpha/(1+alpha)': t_ind >= low_tau,
        'hat E_{b0} T_c <= (n-1) max_i E_i T_c': hat_hit <= (n - 1) * max_to_c,
        'hat tau12 <= (n-1) tau12': t_ind <= (n - 1) * t_orig,
    }
    if heavy:
        checks['rho >= h + (N+1)(1+delta)'] = rho >= M
        checks['tau12 <= h + (N+2)(1+delta)'] = t_orig <= B
    return {
        'N': N, 't': t, 'heavy_end': heavy, 'states': n, 'blocks': len(V), 'pi_c': float(pi[c]),
        'tau12_original': float(t_orig), 'tau12_induced': float(t_ind),
        'ratio': float(t_ind / t_orig), 'ratio_guaranteed': float(low_tau / general),
        'max_hit_c': float(max_to_c), 'hat_hit_c': float(hat_hit),
        'hit_ratio': float(hat_hit / max_to_c),
        'checks': checks,
    }


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
            'largest hat E_a T_A / max_y E_y T_A': float(worst), 'attained_with_(n,|A|)': worst_at}


# ---------------------------------------------------------------------------
# A four-state chain lumped to three states.

def small_example():
    W = [[F(0)] * 4 for _ in range(4)]
    for a, b, v in [(0, 2, 1), (1, 2, 1), (1, 3, 109), (2, 3, 54)]:
        W[a][b] = W[b][a] = F(v)
    W[1][1] = F(117)
    W[3][3] = F(10000)
    t_or, arg_o, _ = tau12(W)
    V = lump(W, [[0], [1], [2, 3]])
    t_in, arg_i, _ = tau12(V)
    return {'tau12_original': str(t_or), 'tau12_induced': str(t_in), 'ratio': str(t_in / t_or),
            'ratio_float': float(t_in / t_or), 'pair_original': arg_o, 'pair_induced': arg_i}


if __name__ == '__main__':
    out = {'comb': [], 'hitting': None, 'small': None}
    all_ok = True
    for N in (1, 2, 3, 4, 6, 8):
        for t in (4, 16, 64, 256):
            for heavy in (True, False):
                row = comb_row(N, t, heavy)
                all_ok = all_ok and all(row['checks'].values())
                out['comb'].append(row)
                print('N=%d t=%-4d heavy=%d  tau %.3f -> %.3f  ratio %.4f  guaranteed %.4f  hitting ratio %.4f  checks %s' % (
                    N, t, heavy, row['tau12_original'], row['tau12_induced'], row['ratio'],
                    row['ratio_guaranteed'], row['hit_ratio'], all(row['checks'].values())))
    out['hitting'] = hitting_trials(400, 20261112)
    print('hitting', out['hitting'])
    out['small'] = small_example()
    print('small', out['small'])
    out['all_comb_checks_hold'] = all_ok
    print('all comb checks hold:', all_ok)
    if len(sys.argv) > 1:
        with open(sys.argv[1], 'w') as fh:
            json.dump(out, fh, indent=1)
