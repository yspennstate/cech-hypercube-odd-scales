"""Checks the syndrome-fibre construction of Section 4 for m = 6 (h = 2) and m = 7, 8, 9 (h = 3).

For the parity-check matrix H of Section 4: every fibre S_s is independent and dominating in Q_m; the simple
columns span F_2^h; for every simple column t = H_j the subgraph induced on S_s and S_(s+t) is the perfect
matching x <-> x + e_j; and the vertex map f_(s,j) carries W_(s+t) bijectively onto W_s and every other W_u
onto a proper subset of W_s (Lemma 4.5). The sphere dimension printed is N - N/2^h - 2.
"""
import itertools, random
random.seed(1)

def run(m, h):
    cols = list(range(1, 2 ** h))
    c0 = cols[0]
    cols = cols + [c0] * (m - len(cols))
    assert len(cols) == m
    N = 1 << m
    syn = [0] * N
    for x in range(N):
        s = 0
        for j in range(m):
            if (x >> j) & 1: s ^= cols[j]
        syn[x] = s
    fib = {s: [x for x in range(N) if syn[x] == s] for s in range(2 ** h)}
    ok = True
    for s, F in fib.items():
        Fs = set(F)
        if any((x ^ (1 << j)) in Fs for x in F for j in range(m)): ok = False; print('not independent', s)
        if any(not any((y ^ (1 << j)) in Fs for j in range(m)) for y in range(N) if y not in Fs): ok = False; print('not dominating', s)
    uniq = [j for j in range(m) if cols.count(cols[j]) == 1]
    span = {0}
    for j in uniq: span |= {a ^ cols[j] for a in span}
    if len(span) != 2 ** h: ok = False; print('unique columns do not span')
    edges = [(u, u ^ (1 << j)) for u in range(N) for j in range(m) if u < u ^ (1 << j)]
    nmatch = 0
    for s in fib:
        for j in uniq:
            t = cols[j]
            P = set(fib[s]) | set(fib[s ^ t])
            ind = {(u, v) for u, v in edges if u in P and v in P}
            want = {(min(x, x ^ (1 << j)), max(x, x ^ (1 << j))) for x in fib[s]}
            if ind != want: ok = False; print('not a perfect matching', s, j)
            nmatch += 1
    # matrix step as set maps: f_s (child s, pairing across once-used column j to s+t) maps W_(s+t) onto W_s
    # bijectively and maps W_u to a proper subset of W_s for every other u
    setbad = 0
    V = set(range(N))
    for s in fib:
        for j in uniq:
            t = cols[j]
            f = {x: (x ^ (1 << j)) if syn[x] == s else x for x in V}
            Ws = V - set(fib[s])
            for u in fib:
                if u == s: continue
                img = {f[x] for x in V - set(fib[u])}
                if u == s ^ t:
                    if img != Ws or len(img) != len(V - set(fib[u])): setbad += 1
                elif not (img < Ws): setbad += 1
    bad = setbad
    print(f"m={m} h={h} N={N} fibres ok={ok} matchings checked={nmatch} sphere dim={N - N // 2 ** h - 2} set-map failures (bijection onto W_s from W_(s+t), proper image otherwise)={bad}")

run(7, 3)
run(8, 3)
run(9, 3)
run(6, 2)
