# Special Smith forms of the meander and chromatic-join matrices

Iaroslav Argunov · manuscript and computational supplement · 18 November 2026

Start with **meander_smith.pdf**. The main theorem gives special Smith equivalence over Z[X,Y] for the chromatic-join matrix. Joins are taken in the full lattice of partitions. The specialization described by [Miller and Stanton](https://doi.org/10.1007/s00605-017-1082-6) gives the meander matrix over Z[x]. The paper includes a direct integral construction and a Dyck-path proof that its quotient has determinant one.

## Reproduce the certificates

From the smith directory, run

~~~sh
python reproduce.py
~~~

Python 3.10 or later with the standard library suffices. This command verifies the manifest, checks 16 scalar identities over Z[z] and reconstructs and checks both polynomial equivalences at orders 6, 7 and 8. All coefficients are exact integers. A full replay with the generator, old explicit matrices and the auxiliary finite-ring experiments uses Python 3.11 or later and the versions in **requirements.txt**:

~~~sh
python -m pip install -r requirements.txt
python reproduce.py --all
~~~

Use plain Python with assertions enabled. The scripts write regenerated files in a temporary directory. The largest auxiliary finite-ring test forms a matrix of dimension 5720. It is separate from the standard-library polynomial verification.

## How the compact two-ring certificate works

Each **rowsN.json** contains the basis order and the sparse left multiplier L. A row entry specifies its column and pairs [degree, integer coefficient] in ascending degree. The one-variable checker forms G from connected components and divides LG entrywise by the predicted monic row factor. This reconstructs the integral quotient V uniquely. It checks that L is lower unitriangular by Dyck area.

For the two-colour lift, assign X degree difference +1 and Y degree difference −1. At fixed difference δ, a monomial xʲ after X=Y=x has exactly one possible lift with exponents (j+δ)/2 and (j−δ)/2. Every required exponent of both multipliers is checked. The program constructs the black and complementary white partitions and counts their two joins independently. It verifies both the sum and the difference of their exponents against the loop count and block grading at every position. Homogeneous specialization is injective for fixed δ, so these checks imply LJ=DV over the full polynomial ring. In both rings the degree of the product of the diagonal factors exhausts the maximum determinant degree and its leading coefficient at equal variables is one. Thus det V=1.

| Order | Dimension | Entries verified in each ring |
| --- | ---: | ---: |
| 6 | 132 | 17,424 |
| 7 | 429 | 184,041 |
| 8 | 1,430 | 2,044,900 |

**verify_rows.py** implements these tests without importing the projector construction. **construct_rows.py** and **check_bezout.py** generate L. **row_checks.json** records the independent verification results and hashes. At orders four and five, the separate **certificate** files store both entire transformations explicitly and have an independent SymPy verifier. **ADDITIONAL_CHECKS.md** explains the scalar weights and the 16 auxiliary tests of additive invariants in ramified rings.

To typeset the paper, run pdfLaTeX twice on **meander_smith.tex** or run Tectonic on that file. All included sections are supplied. SHA256SUMS records the checksum of every other file in this directory.

## Division-free right multiplier

The entry formula uses both block counts and the first prefix of selected arcs before a repeated loop. The proof keeps both possible outer colours and gives the spherical closure for each parity. direct_quotient.py computes V by tracing loops and accumulating univariate weight polynomials. It compares its result with LG coefficient by coefficient and checks the two-colour grading and the two joins. direct_checks.json records the resulting matrix fingerprints, numbers of nonzero polynomials and entries. The standard-library verification includes all orders 6–8. The full replay additionally reconstructs each small projector from its annihilation equations in both shadings and repeats the symbolic trace and denominator checks. See NORMALIZATION.md for the witness and trace_by_rescaling.md for a second proof at every size.

## Reuse

See the [license scope and source attribution](../NOTICE.md) for this paper, its code and its data.
