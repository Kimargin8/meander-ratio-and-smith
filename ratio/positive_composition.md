# An algebraic check of the sharp shift

Let \(M(z)\) count ordinary connected meanders, \(B(z)\) count those avoiding the positive arch orientation, and \(A(z)\) count those avoiding both orientations. Each series starts in degree one. The marked-run identity in the paper gives
\[
 A(z)=M\!\left(\frac{z(1-z)}{1+z}\right),\qquad
 B(z)=M\!\left(\frac{z}{1+z}\right).
\]
Set \(C(z)=M(z/(1+2z))\). For a direct verification of its positivity define \(h=z+q\) by
\[
 q=2z^3+(z+4z^2)q+(1+2z)q^2,\qquad q(0)=0 .
\]
This equation determines the coefficients successively as nonnegative integers. Rearrangement gives \((1+2z)h^2-(1+z)h+z=0\), hence \(h(1-h)/(1+h)=z/(1+2z)\). Substituting in the first marked-run identity proves \(C=A(h)\). This supplies the positive-composition check in the supplement.

Replacing \(z\) by \(z/(1-z)\) in \(C\) yields \(B\). If \(C_j=[z^j]C(z)\), then \(B_{N+1}=\sum_{j=0}^N\binom Nj C_{j+1}\). The root growth of \(B_n\) and the paper's binomial smoothing lemma now give its consecutive-ratio limit. The transfer argument in the paper also treats \(C\) itself, whose second coefficient vanishes.
