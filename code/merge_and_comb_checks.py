"""Exact checks for Theorem 2 (lumping two states), the comb of Section 5 and the bounds on K(n).

Everything is computed in rational arithmetic. Mean hitting times are obtained by solving
the killed linear system for each target, and the mixing time is
    tmix = max_{i,j} (E_i T_j - E_pi T_j)
(Lovasz-Winkler). A loop w_xx counts once in the weight w_x of x; gluing two states turns
the edge between them into a loop carrying twice its weight.

Run:  python merge_and_comb_checks.py > ../data/merge_and_comb_checks_out.txt
"""
from fractions import Fraction as Q
from itertools import combinations
import random


def solve(A, b):
    n = len(b)
    a = [[Q(x) for x in row] + [Q(y)] for row, y in zip(A, b)]
    for j in range(n):
        k = next(k for k in range(j, n) if a[k][j])
        a[j], a[k] = a[k], a[j]
        z = a[j][j]
        a[j] = [x / z for x in a[j]]
        for k in range(n):
            if k != j and a[k][j]:
                z = a[k][j]
                a[k] = [x - z * y for x, y in zip(a[k], a[j])]
    return [row[-1] for row in a]


def analyse(w):
    """Return W, pi, H (mean hitting times), F (centred), kappa_i and tmix."""
    n = len(w)
    d = [sum(row) for row in w]
    W = sum(d)
    pi = [x / W for x in d]
    L = [[(d[i] - w[i][i] if i == j else Q(0)) - (w[i][j] if i != j else 0)
          for j in range(n)] for i in range(n)]
    H = [[Q(0)] * n for _ in range(n)]
    for j in range(n):
        ix = [i for i in range(n) if i != j]
        h = solve([[L[i][k] for k in ix] for i in ix], [d[i] for i in ix])
        for i, v in zip(ix, h):
            H[i][j] = v
    F = [[H[i][j] - sum(pi[k] * H[k][j] for k in range(n)) for j in range(n)]
         for i in range(n)]
    assert all(F[i][j] == F[j][i] for i in range(n) for j in range(n))
    kappa = [max(F[i]) for i in range(n)]
    return W, pi, H, F, kappa, max(kappa)


def lump(w, blocks):
    return [[sum(w[i][j] for i in a for j in b) for b in blocks] for a in blocks]


def harmonic(w, u, v):
    """h_x = P_x(T_u < T_v) and the effective resistance R between u and v."""
    n = len(w)
    S = [x for x in range(n) if x not in (u, v)]
    off = lambda x, y: w[x][y] if x != y else 0
    A = [[(sum(off(x, z) for z in range(n)) if x == y else -off(x, y)) for y in S] for x in S]
    b = [off(x, u) for x in S]
    hs = solve(A, b) if S else []
    h = [Q(0)] * n
    h[u] = Q(1)
    for x, val in zip(S, hs):
        h[x] = val
    energy = sum(off(x, y) * (h[x] - h[y]) ** 2 for x in range(n) for y in range(n)) / 2
    return h, 1 / energy


def random_network(n, rng):
    while True:
        w = [[Q(0)] * n for _ in range(n)]
        for x in range(n):
            for y in range(x, n):
                if rng.random() < 0.6:
                    val = Q(rng.choice([1, 2, 3, 5, 8, 13, 100, 1000]), rng.choice([1, 2, 7, 10]))
                    w[x][y] = w[y][x] = val
        seen, stack = {0}, [0]
        while stack:
            x = stack.pop()
            for y in range(n):
                if y != x and w[x][y] and y not in seen:
                    seen.add(y)
                    stack.append(y)
        if len(seen) == n:
            return w


def check_one_merge(count_per_size=25, seed=20261112):
    rng = random.Random(seed)
    nets = merges = entries = 0
    worst = Q(0)
    for n in range(3, 7):
        for _ in range(count_per_size):
            w = random_network(n, rng)
            W, pi, H, F, kappa, T = analyse(w)
            nets += 1
            for u, v in combinations(range(n), 2):
                blocks = [[x] for x in range(n) if x not in (u, v)] + [[u, v]]
                f = {x: k for k, b in enumerate(blocks) for x in b}
                _, _, _, Fh, _, Th = analyse(lump(w, blocks))
                h, R = harmonic(w, u, v)
                hbar = sum(pi[x] * h[x] for x in range(n))
                beta = [W * R * (h[x] - hbar) for x in range(n)]
                for i in range(n):
                    assert beta[i] == F[i][v] - F[i][u]
                    for j in range(n):
                        assert Fh[f[i]][f[j]] == F[i][j] + beta[i] * beta[j] / (W * R)
                        assert Fh[f[i]][f[j]] <= kappa[i] + kappa[j]
                        entries += 1
                assert Th < 2 * T
                worst = max(worst, Th / T)
                merges += 1
    print(f"one merge: {nets} random networks, {merges} merges, {entries} entries; "
          f"rank-one identity and F^ <= kappa_i + kappa_j hold everywhere; largest ratio {float(worst):.6f} < 2")


def comb(N, t):
    t = Q(t)
    r, al = t ** 3, t
    names = ["a", "b", "c"] + [f"z{k}" for k in range(N + 1)] + [f"y{k}" for k in range(N)]
    ix = {s: i for i, s in enumerate(names)}
    n = len(names)
    w = [[Q(0)] * n for _ in range(n)]

    def s(x, y, val):
        w[ix[x]][ix[y]] += val
        if x != y:
            w[ix[y]][ix[x]] += val

    s("a", "z0", Q(1))
    for k in range(N):
        s(f"z{k}", f"z{k+1}", r ** (k + 1))
        s(f"y{k}", f"z{k+1}", al * r ** (k + 1))
        s(f"y{k}", f"y{k}", (t - 1) * al * r ** (k + 1))
    s("b", f"z{N}", al * r ** (N + 1))
    s("c", f"z{N}", r ** (N + 2))
    s("a", "a", t - 1)
    s("b", "b", (t - 1) * al * r ** (N + 1))
    s("c", "c", (N + 2) * r ** (N + 4))
    blocks = [[ix["a"]], [ix["b"]], [ix[f"z{N}"], ix["c"]]] + \
             [[ix[f"z{k}"], ix[f"y{k}"]] for k in range(N)]
    return w, blocks, ix


def closed_forms(N, t):
    t = Q(t)
    r, al = t ** 3, t
    sigma = lambda k: sum(r ** l for l in range(1, k + 1))
    B = 2 + al * (t + 1)
    S = [t + 1 + r ** (k + 1) + B * sigma(k) for k in range(N)]
    SN = t + 1 + r ** (N + 2) + B * sigma(N) + (t + 1) * al * r ** (N + 1)
    W = (N + 2) * r ** (N + 4) + r ** (N + 2) + SN
    M = t + sum(S[k] / r ** (k + 1) for k in range(N)) + SN / r ** (N + 2)
    Shat = [S[k] + t * al * r ** (k + 1) for k in range(N)]
    Mhat = t + sum(Shat[k] / ((1 + al) * r ** (k + 1)) for k in range(N))
    EpiTc = (t ** 2 * (1 + al * sigma(N + 1)) + sum(S[k] ** 2 / r ** (k + 1) for k in range(N))
          + SN ** 2 / r ** (N + 2)) / W
    EpiTc_hat = (t ** 2 * (1 + al * r ** (N + 1)) + sum(Shat[k] ** 2 / ((1 + al) * r ** (k + 1)) for k in range(N))) / W
    cond = W / r ** (N + 2) >= M and (N == 0 or W / ((1 + al) * r ** N) >= Mhat - t)
    return M, Mhat, M - EpiTc, Mhat + t - EpiTc_hat, cond


def check_comb():
    print("comb: N, states, t, tmix(X), tmix(X^), ratio, Mhat/M (hitting ratio)")
    for N in (0, 1, 2, 3):
        for t in (2, 3, 10, 100, 1000):
            w, blocks, ix = comb(N, t)
            _, _, H, _, _, T = analyse(w)
            wh = lump(w, blocks)
            _, _, Hh, _, _, Th = analyse(wh)
            M, Mhat, Tc, Thc, cond = closed_forms(N, t)
            assert cond and T == Tc and Th == Thc
            c = ix["c"]
            assert M == H[ix["a"]][c] == max(H[x][c] for x in range(len(w)))
            assert Mhat == Hh[0][2]
            print(N, len(w), t, f"{float(T):.6f}", f"{float(Th):.6f}", f"{float(Th / T):.9f}",
                  f"{float(Mhat / M):.6f}")
    print("comb: closed forms equal the solver's values at every point")


def bounds_table():
    print("n  floor(n/2)  max_q min(2^(n-q), q-1)")
    for n in range(2, 13):
        up = max([1] + [min(2 ** (n - q), q - 1) for q in range(2, n)])
        print(n, n // 2, up)


if __name__ == "__main__":
    check_one_merge()
    check_comb()
    bounds_table()
