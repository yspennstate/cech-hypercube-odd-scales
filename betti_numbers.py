"""Reduced Betti numbers over F_p (p = 1000003) of X_m = N(Q_m, 2m-3) for m = 3, 4, 5.

By Lemma 4.1 the faces of X_m are the sets whose complement contains an edge of Q_m, so X_m is the
Alexander dual of the independence complex Ind(Q_m) and H~_i(X_m) = H~^(N-i-3)(Ind(Q_m)), N = 2^m.
The script computes the Betti numbers of Ind(Q_m) and transfers them; for m = 3 it also computes X_3
directly. It prints, next to the result, the lower bound 2^h - 1 in degree N - N/2^h - 2 from Theorem 1.2
(for m < 6 only h = 2 occurs). Usage: python betti_numbers.py [m ...]   (default m = 3 4; m = 5 is slower)
"""
import itertools, sys
P = 1000003

def rank_mod(rows, ncols):
    piv = {}
    r = 0
    for row in rows:
        row = dict(row)
        while row:
            c = min(row)
            if c in piv:
                prow = piv[c]; f = row[c]
                for k, v in prow.items():
                    nv = (row.get(k, 0) - f * v) % P
                    if nv: row[k] = nv
                    else: row.pop(k, None)
            else:
                inv = pow(row[c], P - 2, P)
                piv[c] = {k: v * inv % P for k, v in row.items()}
                r += 1
                break
    return r

def betti(faces):
    """faces: set of frozensets incl. empty set. Reduced Betti numbers."""
    by = {}
    for f in faces: by.setdefault(len(f) - 1, []).append(tuple(sorted(f)))
    idx = {d: {f: i for i, f in enumerate(sorted(L))} for d, L in by.items()}
    top = max(by)
    rk = {}
    for d in range(0, top + 1):
        rows = []
        for f in idx[d]:
            row = {}
            for i in range(len(f)):
                g = f[:i] + f[i + 1:]
                row[idx[d - 1][g]] = (-1) ** i % P
            rows.append(row)
        rk[d] = rank_mod(rows, len(idx[d - 1]))
    out = {}
    for d in range(-1, top + 1):
        b = len(idx[d]) - rk.get(d, 0) - rk.get(d + 1, 0)
        if b: out[d] = b
    return out

def cube(m):
    V = list(range(1 << m)); E = [(u, u ^ (1 << j)) for u in V for j in range(m) if u < u ^ (1 << j)]
    return V, E

def ind_faces(m):
    V, E = cube(m); adj = {v: {v ^ (1 << j) for j in range(m)} for v in V}
    faces = set(); stack = [(frozenset(), 0)]
    while stack:
        s, start = stack.pop(); faces.add(s)
        for v in range(start, len(V)):
            if not (adj[v] & s): stack.append((s | {v}, v + 1))
    return faces

def x_faces(m):
    V, E = cube(m); N = len(V); faces = set()
    for mask in range(1 << N):
        comp = {v for v in V if not (mask >> v) & 1}
        if any(u in comp and w in comp for u, w in E):
            faces.add(frozenset(v for v in V if (mask >> v) & 1))
    return faces

for m in [int(a) for a in sys.argv[1:]] or [3, 4]:
    N = 1 << m
    bi = betti(ind_faces(m))
    dual = {N - j - 3: b for j, b in bi.items()}
    pred = {N - N // 2 ** h - 2: 2 ** h - 1 for h in range(2, 20) if 2 ** h - 1 <= m}
    print(f"m={m} N={N} Ind betti={bi} -> X betti via duality={dual} predicted_lower={pred}", flush=True)
    if m == 3:
        print(f"m=3 X direct betti={betti(x_faces(3))}", flush=True)
