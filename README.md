# Collapsibility numbers and high-degree homology of Čech complexes of hypercube graphs

Code accompanying the paper of the same title by Yitzchak Shmalo.

- `check_syndrome_fibres.py` checks the syndrome-fibre construction of Section 4 for m = 6, 7, 8, 9:
  the fibres are independent and dominating, the simple columns span, consecutive fibres induce perfect
  matchings, and the vertex maps of Lemma 4.5 behave as stated. Output: `fibres_output.txt`.
- `betti_numbers.py` computes the reduced Betti numbers over F_p of N(Q_m, 2m-3) for m = 3, 4, 5 through
  Alexander duality with the independence complex of Q_m. Output: `betti_output.txt`.

Both scripts use only the Python standard library (Python 3.8 or later):

    python check_syndrome_fibres.py
    python betti_numbers.py 3 4 5
