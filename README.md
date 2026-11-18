# Closed meanders and special Smith forms

Iaroslav Argunov · manuscripts dated 18 November 2026

| Manuscript | Computational material |
| --- | --- |
| [Local arch statistics and ratio limits for closed meanders](ratio/adjacent_arches.pdf) | [Methods and reproduction](ratio/README.md) · [Two-orientation counts](ratio/joint_counts.json) · [Weighted transfer check](ratio/weighted_check.py) |
| [Special Smith forms of the meander and chromatic-join matrices](smith/meander_smith.pdf) | [Methods and reproduction](smith/README.md) · [Explicit right multiplier](smith/direct_quotient.py) · [Polynomial certificates](smith/direct_checks.json) |

The first paper counts connected pairs of noncrossing matchings with weights on two orientations of adjacent arches. The consecutive ratio of the weighted sums converges for every pair of nonnegative weights, uniformly on compact subsets. The second paper constructs integral changes of basis for the chromatic-join and meander matrices. Each directory contains the manuscript source, input data, verification programs and expected outputs.

## Reproduce the calculations

From the repository root run

```sh
python ratio/reproduce.py
python smith/reproduce.py
```

These checks use Python 3.10 or later with its standard library. The ratio check regenerates the occupation tables and runs small independent enumerations. The Smith check verifies both polynomial rings at orders six, seven and eight.

For the full Smith calculation, use a Python 3.11 or newer environment compatible with the pinned NumPy release. The commands below regenerate the rows and repeat the projector and finite-ring checks.

```sh
python -m pip install -r smith/requirements.txt
python smith/reproduce.py --all
```

Run Python with assertions enabled. Each verifier checks the SHA-256 manifest for its directory before running. The [repository manifest](SHA256SUMS) covers all files except itself. To cite a result, cite the corresponding manuscript by its title and author.

## Reuse

Iaroslav Argunov licenses the original Python code under MIT and both manuscripts, authored documentation and independent certificates under CC BY 4.0. The OEIS input and the three tables derived from it are available under CC BY-SA 4.0 for the material covered by copyright or database rights. The [file-by-file scope and source credit](NOTICE.md) identify the applicable terms.
