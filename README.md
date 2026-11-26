# Collapsibility numbers and high-degree homology of Čech complexes of hypercube graphs

Yitzchak Shmalo, Einstein Institute of Mathematics, The Hebrew University of Jerusalem.

Let N(Q_n, r) be the Čech complex of the hypercube graph Q_n at integer scale r. We show that its collapsibility number is 2^(k+1) - 1 at every even scale r = 2k and every n >= k + 1, the value proposed by Adams, Shukla and Singh (Čech complexes of hypercube graphs, Homology Homotopy Appl. 27 (2025) 83-105, doi:10.4310/HHA.2025.v27.n1.a6, arXiv:2212.05871). At odd scale r = 2m - 3 with m >= 3 we construct, for each h with 2 <= h <= floor(log2(m + 1)), a free direct summand of rank at least 2^h - 1 in the reduced integral homology of N(Q_n, r) in degree 2^m - 2^(m-h) - 2, for every n >= m. The summand for h = 3 lies in a degree outside the set proposed by Adams, Shukla and Singh for every odd r >= 11, and it forces the collapsibility number above the proposed odd-scale value 3 * 2^k - 1. Together with an upper bound of 2^(k+2) - 1, the collapsibility number at scale 2k + 1 divided by 2^(k+2) tends to 1, uniformly in n >= k + 2. For m >= 7 the complex Theta(Cube(m, 1)) = N(Q_m, 2m - 3) of Conant and Thistlethwaite is not homotopy equivalent to a wedge of spheres of a single dimension, which proves that case of the first part of their conjecture.

## Contents

- `paper/main.tex`, `paper/main.pdf`: the paper.
- `check_syndrome_fibres.py` checks the syndrome-fibre construction of Section 4 for m = 6, 7, 8, 9: the fibres are independent and dominating, the simple columns span, consecutive fibres induce perfect matchings, and the vertex maps of Lemma 4.5 behave as stated. Output: `fibres_output.txt`.
- `betti_numbers.py` computes the reduced Betti numbers over F_p of N(Q_m, 2m-3) for m = 3, 4, 5 through Alexander duality with the independence complex of Q_m. Output: `betti_output.txt`.
- `SHA256SUMS`: sha256 of each file above.

Both scripts use only the Python standard library (Python 3.8 or later):

    python check_syndrome_fibres.py
    python betti_numbers.py 3 4 5

The code is under the MIT License (`LICENSE`).
