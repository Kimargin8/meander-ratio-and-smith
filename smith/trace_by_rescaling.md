# A second calculation of the spherical projector trace

Let a box have \(m\ge1\) ports on each of two rows. Label its boundary clockwise from the upper left. Call the shade of its leftmost region black, write \(W\) for its weight, and write \(V\) for the other shade. Thus \((W,V)\) is either \((X,Y)\) or \((Y,X)\). Closing the box on the right connects each top port to its own bottom port by disjoint exterior arcs. The exterior disk with this matching is denoted by \(I\). In the clockwise indexing its pairs are \(i\leftrightarrow 2m-1-i\). Drawing all the closure arcs to the left gives the same matching on the spherical complement disk; choosing an exterior point affects the planar evaluation.

**Proposition.** Put \(t^2=WV\) and let \(q_0=0,\ q_1=1,\ q_{j+1}=tq_j-q_{j-1}\). For the generic two-colour projector \(F_m\) with the indicated left shade, the complete spherical evaluation against \(I\), including its exterior face, is
\[
 S_m(W,V)=
 \begin{cases}
 W q_{m+1}(t),&m\ {\rm even},\\[2pt]
 t q_{m+1}(t),&m\ {\rm odd}.
 \end{cases}
\]
Against every other external matching the evaluation vanishes. In both cases \(S_m\in\mathbb Z[W,V]\), by the parity of \(q_{m+1}\).

**Proof.** Work generically over the extension \(\mathbb Q(W,V,t)\). For two disk diagrams \(P,Q\) let \(r(P)\) and \(r(Q)\) be their numbers of black regions before gluing. Their union has \(\ell\) circles. Let \(a,b\) count its black and white spherical faces. Then
\[
             a+b=\ell+1,\qquad a-b=r(P)+r(Q)-m-1.
\]
For the second identity, take one vertex in each black region of the two disks. In each such region draw disjoint spokes from that vertex to its incident black boundary intervals. After gluing, the graph has \(r(P)+r(Q)\) vertices, \(m\) edges, \(a\) connected components and \(b\) complementary regions. The graph thickens to the black portions of the sphere. The Euler formula for a possibly disconnected plane graph gives the assertion. The same region-count argument underlies the diagonal rescaling for chromatic joins and Temperley–Lieb matrices. [Chen and Przytycki](https://arxiv.org/abs/0806.0878)

Set \(\eta=W/t\), \(h(P)=\eta^{r(P)-r(I)}\), and \(\alpha=t\eta^{2r(I)-m-1}\). The previous identities give a termwise comparison with the usual unshaded pairing of loop weight \(t\):
\[
                       W^a V^b=\alpha h(P)h(Q)t^{\ell(P,Q)}.
\]
In the disk \(I\), the \(m\) parallel strands separate \(m+1\) successive strips. Starting at its left side these are alternately black and white, so \(r(I)=1+\lfloor m/2\rfloor\). Thus \(\alpha=W\) for even \(m\) and \(\alpha=t\) for odd \(m\).

Write the ordinary Jones–Wenzl projector as \(f_m=\sum_P c_P P\), with \(c_I=1\). Form the candidate \(\sum_P c_P h(P)^{-1}P\) in the shaded diagram basis. For every \(Q\ne I\), the outside matching has a same-row pairing: an all-through noncrossing pairing is forced to be \(I\). Among the pairs with endpoints on that row there is an adjacent pair, since the points it encloses can only be matched within the same interval. The resulting turnback annihilates \(f_m\). By the displayed pairing identity, the candidate vanishes against the same outside matchings. The genuine shaded projector \(F_m\) also vanishes against them by its turnback relations. If an adjacent boundary cap is presented in isolation, composing with the corresponding cup is injective on the diagram basis, so annihilation by the generator implies annihilation by the cap.

To justify uniqueness at this step, the ordinary Gram matrix \((t^{\ell(P,Q)})_{P,Q}\) is invertible over the generic field: each diagonal entry has degree \(m\) and every off-diagonal entry has smaller degree, so the determinant has monic leading term. The shaded matrix differs by invertible diagonal factors and the nonzero scalar \(\alpha\). Its columns with \(Q\ne I\) thus have a one-dimensional left annihilator. Both rows belong to that line and have coefficient \(1\) at \(I\), hence are equal. At \(Q=I\) their common spherical evaluation is \(\alpha\) times the ordinary, unnormalized full trace of \(f_m\).

For the ordinary positive loop weight \(t\), Wenzl recursion gives
\[
  f_m=P-\frac{q_{m-1}}{q_m}P e_{m-1}P,\qquad
  P=f_{m-1}\otimes 1,\qquad m\ge2 .
\]
In these conventions \(e_i^2=t e_i\). Closing the extra through strand shows \(\operatorname{tr}(P)=t\operatorname{tr}(f_{m-1})\). Cyclicity, \(P^2=P\) and the cap/cup at the right give \(\operatorname{tr}(P e_{m-1}P)=\operatorname{tr}(f_{m-1})\). Starting from \(\operatorname{tr}(f_1)=t\), induction therefore yields \(\operatorname{tr}(f_m)=q_{m+1}\). This proves the proposition. In the negative loop convention used by Morrison the sign of the recursion coefficient changes accordingly. [Morrison, *A formula for the Jones–Wenzl projections*, §2](https://arxiv.org/abs/1503.00384)

**Partial closure.** Closing the rightmost strand of \(F_m\) leaves an element annihilated by every surviving turnback. The same annihilator-line argument at size \(m-1\) shows this element is proportional to the projector with the same left shade. Closing the remaining strands determines the coefficient as \(S_m(W,V)/S_{m-1}(W,V)\), with \(S_0(W,V)=W\). A leftmost closure leaves opposite left shade, with coefficient \(S_m(W,V)/S_{m-1}(V,W)\). These equalities recover the coloured partial-trace scalars. The derivation uses unshaded Wenzl recursion and face counts; it permits weights to be inverted only during this normalization argument in the generic field. The integral divisibility in the manuscript uses the separate denominator lemma over \(\mathbb Z[X,Y]\). [Hazi, coloured projector and trace](https://arxiv.org/abs/2302.14476)
