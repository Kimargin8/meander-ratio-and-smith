# Local arch statistics and ratio limits for closed meanders

Iaroslav Argunov · manuscript and computational supplement · 18 November 2026

Read **adjacent_arches.pdf** for the definitions and proofs. The order is half the number of labelled crossings. The two planar matchings are ordered. The test sites start at even labels, and the final site does not wrap around the line.

The two marked orientations have a joint partition sum whose consecutive ratios converge for every nonnegative pair of activities, including the case where both orientations are forbidden. The result is locally uniform in the two activities. Classical transfer moves also give a sharp binomial shift and its ratio limit. The manuscript identifies the short deduction of the ordinary ratio limit from Jensen and Kesten and gives self-contained arguments for the weighted and shifted families. The constants are expressed in terms of the ordinary root growth constant.

## Reproduce the finite checks

From the ratio directory, run the following command there, with Python 3.10 or later and assertions enabled. The computation uses the standard library and arbitrary-precision integers.

~~~sh
python reproduce.py
~~~

The program verifies the file manifest, regenerates both tables in a temporary directory and replays an independent check of the surgery on bracket words. At orders 1 through 8 the latter checks all 2,248,355 ordered pairs before connectivity filtering, including every compatible or incompatible pattern prescription. For the explicit contraction/expansion round trips the upper limit is order 7. The data file records each range.

| File | Meaning |
| --- | --- |
| b005315.txt | Known connected counts through order 38, from [OEIS A005315](https://oeis.org/A005315/b005315.txt). The project retained the input on 15 November 2026. Index 0 is unused in the proofs. |
| joint_counts.json | One triangular array for each order 1–38. Entry at row i, column j counts diagrams with exactly i positive and j negative sites. Through order 8 the entire array is also computed directly from the ordered pairs. |
| longest_runs.json | Entry k of each row counts diagrams with longest positive run exactly k. The integer calculation also checks all binary words on at most 12 sites. |
| surgery_checks.json | Independent direct checks of the insertion, contraction, signed prescriptions, arbitrary disjoint exposed intervals, conditional law and both variants of the seam. This program does not read OEIS. |
| joint_counts.py, longest_runs.py, verify_surgery.py | Standalone implementations used to produce the respective tables. |
| SHA256SUMS | SHA-256 of every other file in the archive. |

To typeset the paper, run pdfLaTeX twice on **adjacent_arches.tex**, or run Tectonic on that file. The other TeX files are included sections. No separate bibliography processor is needed.

## Shifted counts and marked transfer

The coefficients of the shift by two are recorded in positivity_checks.json. Two computations use signed differences of the ordinary counts and a positive decoration of diagrams avoiding both arch signs. The short algebraic proof of this second computation is included in positive_composition.md. Through order eight, an independent ordered-pair oracle also computes the shifted class directly. Beyond that range the algebraic comparisons use the same ordinary counts. check_positivity.py tests all 741 differences through order 38 and rejects a modified avoidance input.

check_idle_moves.py starts from the two matchings and counts diagrams whose interior word avoids two specified two-letter blocks. check_return_words.py enumerates all odd cuts with at most sixteen ports and tests every return of two, three and four blocks found at the two-port state. See return_words.md for its geometric interpretation. weighted_check.py independently obtains the polynomials at orders 1–11 and checks the unmarked returns of the transfer with an extra letter of memory. In that code R moves a lower port to the upper half-plane. Reflect the names R and L when reading the manuscript's convention in which R moves an upper (left) port to the lower (right) side. Both record positive marks by the short upper arch on the first two vertices of a triple.

## Source data

The file [b005315.txt](b005315.txt) reproduces the table of orders 0–38 from *The On-Line Encyclopedia of Integer Sequences*, OEIS Foundation Inc. and its contributors. The snapshot was retrieved on 15 November 2026. Its bytes are unchanged. The entry credits Efthimios Alepis for the table, with earlier values from Iwan Jensen, Andrew Howroyd and Vladimir Osipov. The occupation tables and related scalar checks in the ratio directory use these published counts as inputs. [A005315](https://oeis.org/A005315) · [Source table](https://oeis.org/A005315/b005315.txt)

The OEIS material is distributed under the Creative Commons Attribution-ShareAlike 4.0 International license. That attribution and license apply to the source material wherever it is reproduced or adapted here. [OEIS End-User License Agreement](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement) · [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)

## Reuse

See the [license scope and source attribution](../NOTICE.md) for this paper, its code and its data.
