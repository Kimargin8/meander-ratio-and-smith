# Return words at an odd meander cut

## Convention

At a vertical cut, list the upper ports from farthest to nearest to the river, followed by the lower ports from nearest to farthest. The arcs drawn to the left give a noncrossing perfect matching of this ordered list. A state records the entire matching and the position of the river. The move \(U\) inserts an initially paired port on each side of the river. The move \(R\) transfers the nearest lower port to the nearest upper position, and \(L\) does the reverse. The move \(D\) removes the two nearest ports and joins their partners. A transfer from an empty side is forbidden. A deletion requires both sides to be nonempty and the two removed ports to have different partners from each other. Deleting a paired pair would close a component prematurely.

Both side counts are odd after an odd number of elementary steps. Let \(e\) have one port on each side, paired together. A connected closed meander starts with a forced \(U\) into \(e\) and ends with the final closure from \(e\). All internal moves follow the preceding rules. On the odd cuts group the remaining letters into disjoint ordered pairs. In this convention the words \(LR\) and \(RL\) are each admissible identities on every state. Let \(Q\) count paths of one block with these two labels excluded and all other admissible labels retained, including multiplicities. The full two-step transfer is \(2I+Q\).

## Two explicit returns

Write \(b\) for the nearest old lower port, \(z\) for its partner and \(u_j,v_j\) for the nearest upper and lower ports introduced by the \(j\)-th \(U\). No restriction on \(z\) is needed.

* In \(URDL\), after \(UR\), the upper list starts at the river with \(v_1,u_1\) and the lower list starts with \(b\). The next \(D\) removes \(v_1,b\), which have different partners \(u_1,z\). This creates \(u_1-z\). The last \(L\) puts \(u_1\) in the original position of \(b\). Every other position and pair retains its initial value. The two blocks are \(UR,DL\).
* In \(UURDDL\), after \(UUR\), the upper list starts at the river with \(v_2,u_2,u_1\) and the lower list starts with \(v_1,b\). The first \(D\) joins \(v_2,v_1\), belonging to two different freshly inserted pairs. This creates \(u_2-u_1\). The second \(D\) joins \(u_2,b\). Their respective partners are \(u_1,z\), so this deletion also preserves admissibility, even if \(z\) is the nearest old upper port. The final \(L\) restores the original matching with \(u_1\) occupying the position of \(b\). The three blocks are \(UU,RD,DL\).

## A uniform bound from every return at the two-port state

**Proposition.** Put \(d_k=(Q^k)_{ee}\), including \(d_0=1\). Every word counted by \(d_k\) is an admissible return from every odd-cut state. Consequently
\[
Q^k\ \geq\ d_k I,\qquad
\min_s (Q^k)_{ss}=d_k .
\]
The inequality is entrywise on the infinite, unquotiented transfer graph.

**Proof.** Take any such word. Start at a general odd cut and mark its nearest upper and lower ports. Regard these as the two inputs to the word, kept separate at its left boundary. The remaining ports stay in place. Indeed, in the execution from \(e\), the number of local upper and lower ports is determined by the word, independently of their partners. Whenever a move consumes a port, there is an available local port on the required side. By induction the marked ports and their successors are nearer to the river than all untouched ports on their respective sides. Thus the word can be drawn in a rectangle with exactly two input and two output endpoints. These local height counts also apply when one side temporarily has zero local ports.

The rectangle has no closed component: such a component would still be present with the original input cup of \(e\) attached. Its two inputs cannot be paired with one another inside the rectangle either, since the input cup would then close a component during the run from \(e\). Hence each input goes through to an output. The rectangle is planar, and its boundary order makes the exchanged pairing crossing. The upper input therefore reaches the upper output and the lower input reaches the lower output. Restoring any original matching outside this rectangle reproduces the full starting matching and creates no closed component. Every prefix is admissible as well, since a closed component, once formed, persists in any continuation of the drawn diagram. Distinct labels remain distinct on every state. Taking the diagonal entry at \(e\) proves sharpness. \(\square\)

For \(k=2\) there are four words, \(LURD,RULD,ULDR,URDL\). For \(k=3\) there are ten. Thus \(Q^2\geq4I\) and \(Q^3\geq10I\). These constants concern multiplicities in the stated labelled transfer.

## Finite verification

The companion checker enumerates all noncrossing matchings on \(2m\) ordered ports and every odd river split, including states outside any specified finite meander calculation. Counts for \(m=1,\dots,8\) are \(1,4,15,56,210,792,3003,11440\), respectively. It first finds every return at \(e\) of two, three and four allowed blocks by exhaustive word enumeration. It then tests all \(4+10+70=84\) words separately on each of the 15,521 states, for 1,303,764 word–state tests in addition to the two single-block identity checks at each state. The complete words and output are in the JSON certificate.

| Number of internal blocks \(k\) | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| \(((2I+Q)^k)_{ee}=M_{k+1}\) | 1 | 2 | 8 | 42 | 262 | 1828 | 13820 | 110954 |
| \((Q^k)_{ee}\) | 1 | 0 | 4 | 10 | 70 | 376 | 2484 | 16578 |
| \(((I+Q)^k)_{ee}\) | 1 | 1 | 5 | 23 | 135 | 867 | 6051 | 44747 |

The last two rows also pass exact binomial inversion against the first row, with shifts \(-2\) and \(-1\), respectively. Reproduce with Python 3.9 or later, from this directory:

~~~
python check_return_words.py --write return_words.json
~~~
