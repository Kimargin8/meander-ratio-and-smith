# Shaded projector denominators and spherical traces

This finite witness works in \(\mathbb Q(X,Y)\) for one through six strands and both choices of the leftmost region. It builds the Jones–Wenzl projector from its identity coefficient and its turnback equations. All subsequent trace and denominator comparisons are exact polynomial identities. The program reads no meander tables, Smith multipliers or projector coefficients.

## Reproduction

With Python 3.10 or later and SymPy 1.14.0 installed, run the following command in the extracted directory.

~~~sh
python check_colour_trace.py --max-strands 6 --write reproduced.json
~~~

The expected SHA-256 of the generated JSON is

~~~text
21bf11fb279b13bc58e4f18c2727831c78b4b3a0803949333732073026a01dd0
~~~

The file is encoded as UTF-8 with LF line endings. It contains no platform-dependent fields. A quick check through five strands can be run with **--max-strands 5**.

## Construction and normalization

Use clockwise boundary order, with the top row read left to right. The leftmost region has shade \(c\), its weight is \(w_c\), and the two face weights are \(X,Y\). For an adjacent generator whose cup starts at the first top port, a closed disk has the opposite weight \(w_{\bar c}\). Successive top columns alternate.

In the diagram basis set the straight-through coefficient to one. For each other diagram, impose the coefficient equation belonging to its first adjacent top cup in \(e_i f=0\). A fraction-free solve supplies a unique candidate. At every tested size the selected system is nonsingular. The program then checks *every* coefficient equation of \(e_i f=0\) and \(f e_i=0\), including all the unused equations. These conditions characterize the generic projector. Cancelling the polynomial gcd of all its numerators and its common denominator gives the actual monic least common denominator.

Set \(Q_{0,c}=0,\ Q_{1,c}=1\). For \(j\ge1\) put \(Q_{j+1,c}=\delta_j Q_{j,c}-Q_{j-1,c}\), where \(\delta_j=w_{\bar c}\) for odd \(j\) and \(\delta_j=w_c\) for even \(j\). In the partial trace on the right the scalar is \(Q_{m+1,c}/Q_{m,c}\) and the remaining projector has left shade \(c\). For the partial trace on the left, the scalar uses the original rightmost shade instead, and the remaining projector has left shade \(\bar c\). The two reductions are computed by deleting the joined port pair on every basis diagram and collecting coefficients. The program also compares the reduced common denominator with the product of primitive two-colour factors in [Hazi, Lemma 4.2 and the proof of Theorem A](https://arxiv.org/abs/2302.14476).

To check complete spherical closure, glue two disks with boundary intervals indexed clockwise. A chord with endpoints \(i,j\) identifies the intervals immediately before \(i\) and after \(j\), and those immediately after \(i\) and before \(j\). Make these identifications for both disks. Each resulting interval class is one face of the glued sphere. Count its shade and form the corresponding monomial, including every face. Apply this calculation to each term of the reconstructed projector and each external noncrossing matching. The only nonzero outcome is the through matching, whose weight is \(w_c Q_{m+1,c}\).

## Scope and controls

| Exact comparison through six strands | Number |
| --- | ---: |
| Common denominators, counting both shades | 12 |
| Coefficient equations for turnbacks on the two rows | 1,144 |
| Coefficients of left and right partial traces | 260 |
| Complete spherical evaluations against external matchings | 392 |

With \(X=2,Y=3\) and left shade \(X\), a one-strand projector gives spherical weight \(6\). Dropping the outer region gives \(3\). For two strands the projector is \(1-e_1/Y\). Its least denominator evaluates to \(3\), whereas exchanging its shade gives \(2\). The correct spherical weight is \(X(XY-1)=10\), whereas using the opposite outermost face gives \(Y(XY-1)=15\). The verifier rejects these three changes symbolically, before specialization. This certificate covers projector traces and their denominators in the stated range. The all-order assertions in the manuscript use Hazi's results and the proof of the normalization.
