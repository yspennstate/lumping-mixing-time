# Exact check of the four-state family of Remark 7: w02=w12=1, w03=w13=t, w00=w11=t^2-1, w33=t^6, lump {2,3}.
# tmix is computed from max_{i,j}(E_i T_j - E_pi T_j) in rational arithmetic and compared with the closed forms.
# Usage: python four_state_family.py
from fractions import Fraction as Q


def solve(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]


def tmix(w):
    n = len(w)
    deg = [sum(row) for row in w]
    tot = sum(deg)
    pi = [d / tot for d in deg]
    P = [[w[i][j] / deg[i] for j in range(n)] for i in range(n)]
    H = [[Q(0)] * n for _ in range(n)]
    for j in range(n):
        idx = [i for i in range(n) if i != j]
        A = [[(1 if a == b else 0) - P[a][b] for b in idx] for a in idx]
        h = solve(A, [Q(1)] * len(idx))
        for k, i in enumerate(idx):
            H[i][j] = h[k]
    best, arg = None, None
    for j in range(n):
        Epi = sum(pi[i] * H[i][j] for i in range(n))
        for i in range(n):
            v = H[i][j] - Epi
            if best is None or v > best:
                best, arg = v, (i, j)
    return best, arg


def family(t):
    t = Q(t)
    w = [[Q(0)] * 4 for _ in range(4)]
    for a, b, x in [(0, 2, 1), (1, 2, 1), (0, 3, t), (1, 3, t)]:
        w[a][b] = w[b][a] = Q(x)
    w[0][0] = w[1][1] = t * t - 1
    w[3][3] = t ** 6
    # lump {2,3}: an internal edge counts in both orientations, so the loop weight is w22+w33+2*w23
    l = [[Q(0)] * 3 for _ in range(3)]
    for a in range(2):
        for b in range(2):
            l[a][b] = w[a][b]
        l[a][2] = l[2][a] = w[a][2] + w[a][3]
    l[2][2] = w[2][2] + w[3][3] + 2 * w[2][3]
    return w, l


ok = True
for t in [2, 3, 5, 10, 100, 1000, 10 ** 6]:
    w, l = family(t)
    a, arg_a = tmix(w)
    b, arg_b = tmix(l)
    T = Q(t)
    W = T ** 6 + 2 * T ** 2 + 4 * T + 2
    fa = (T + 1) * (T ** 6 + T ** 5 + 4 * T + 2) / W
    fb = 2 * T * (T ** 6 + T ** 2 + 3 * T + 2) / W
    good = a == fa and b == fb
    ok &= good
    print(f"t={t}: tmix(X)={float(a):.9f} at {arg_a}, tmix(Xhat)={float(b):.9f} at {arg_b}, "
          f"ratio={float(b / a):.12f}, closed forms {'agree' if good else 'DISAGREE'}")
assert ok
print("all closed forms agree")

# Symbolic check for every t >= 2: with s = t - 2, solve the hitting-time equations over the field Q(s), form
# D_ij = E_i T_j - E_pi T_j, and check that the claimed maximum minus every other D_ij is a quotient of two
# polynomials in s with nonnegative coefficients.
import sympy as sp
from sympy.polys.matrices import DomainMatrix

s = sp.symbols("s")
K = sp.QQ.frac_field(s)
t = K.convert(s) + 2


def sym_D(w):
    n = len(w)
    deg = [sum(row, K.zero) for row in w]
    tot = sum(deg, K.zero)
    H = [[K.zero] * n for _ in range(n)]
    for j in range(n):
        idx = [i for i in range(n) if i != j]
        A = DomainMatrix([[(K.one if a == b else K.zero) - w[a][b] / deg[a] for b in idx] for a in idx],
                         (n - 1, n - 1), K)
        h = A.lu_solve(DomainMatrix([[K.one] for _ in idx], (n - 1, 1), K))
        for k, i in enumerate(idx):
            H[i][j] = h[k, 0].element
    D = {}
    for j in range(n):
        Epi = sum((deg[i] / tot * H[i][j] for i in range(n)), K.zero)
        for i in range(n):
            D[(i, j)] = H[i][j] - Epi
    return D


def nonneg(x):
    pn, pd = sp.Poly(K.numer(x).as_expr(), s), sp.Poly(K.denom(x).as_expr(), s)
    if pd.LC() < 0:
        pn, pd = -pn, -pd
    return all(c >= 0 for c in pn.all_coeffs()) and all(c >= 0 for c in pd.all_coeffs())


Z = K.zero
w = [[Z] * 4 for _ in range(4)]
for a, b, x in [(0, 2, K.one), (1, 2, K.one), (0, 3, t), (1, 3, t)]:
    w[a][b] = w[b][a] = x
w[0][0] = w[1][1] = t ** 2 - 1
w[3][3] = t ** 6
l = [[w[a][b] for b in range(2)] + [w[a][2] + w[a][3]] for a in range(2)]
l.append([l[0][2], l[1][2], w[2][2] + w[3][3] + 2 * w[2][3]])
W = t ** 6 + 2 * t ** 2 + 4 * t + 2
for name, ww, top, closed in [("X", w, (3, 2), (t + 1) * (t ** 6 + t ** 5 + 4 * t + 2) / W),
                             ("Xhat", l, (1, 0), 2 * t * (t ** 6 + t ** 2 + 3 * t + 2) / W)]:
    D = sym_D(ww)
    assert D[top] == closed, name
    for key, val in D.items():
        assert nonneg(D[top] - val), (name, key)
    print(f"{name}: for every t >= 2 the maximum of E_iT_j - E_piT_j is at {top} and equals the closed form")
