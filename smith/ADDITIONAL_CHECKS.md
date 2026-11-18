# Integral row elimination and special fibres

The note integral_gluing.tex gives integer polynomial row operations for both special Smith forms. The scalar coefficients in this construction depend only on the maximum height of the row. The two checks below address different parts of the argument.

## Univariate Bézout weights

Run these commands from the directory containing the files.

~~~
python check_bezout.py --verify bezout_weights.json
python check_bezout.py 17 --write fresh_weights.json
~~~

Python 3.10 or later suffices. The distributed file covers \(N=k+1=2,\ldots,17\). For each record the key is the number \(r\), and its array encodes the coefficients of \(a_{N,r}(z)\) starting at degree zero. The verifier builds the primitive factors from the Chebyshev recurrence and checks the identity in \(\mathbb Z[z]\) by exact multiplication. The generator uses extended Euclidean division over the rationals and checks that all resulting coefficients belong to \(\mathbb Z\).

## Ramified finite rings

~~~
python check_ramified_fibres.py --orders 1,2,3,4,5,6,7 --write fresh_fibres.json
python check_ramified_fibres.py --orders 8 --write fresh_fibres_order8.json
~~~

These commands require NumPy. The archived records include the Python and NumPy versions. They record the additive Smith invariants of the regular representation of \(J/(XY)\) in the two algebras
\[
 (\mathbb Z/p^3)[u,v]/(u^2-p,v^2-p),\qquad p=2,3.
\]
The polynomial matrix is divided by \(XY\) before \(X=u,\ Y=v\) are substituted. The program enumerates all noncrossing set partitions and forms both lattice joins directly. Independent Dyck height counts and the two \(2\)-by-\(2\) blocks of multiplication by \(A+Buv\) give the predicted additive invariants. Direct enumeration of small linear maps checks the modular elimination routine. The largest check uses 1430 partitions and a \(5720\)-by-\(5720\) matrix over \(\mathbb Z/p^3\).

For each computation the four counts refer, in order, to valuations \(0,1,2,\ge3\). The SHA-256 input consists of the complete additive matrix in row-major order using four-byte unsigned little-endian entries, reduced to \(0,\ldots,p^3-1\). The basis blocks are \(1,u,v,uv\) in that order. The records verify additive invariants of these specializations. Equivalence over \(\mathbb Z[X,Y]\) follows from the mathematical proof.
