"""Minimal cylinder-diagram (separatrix-diagram) calculator for TRANSLATION surfaces.

A9 Phase-1 step (ii-c) tool: given the horizontal-cylinder decomposition of a
translation surface as a pair of cyclic words per cylinder (bottom / top boundary
sequences of saddle connections), compute genus + stratum (singularity orders).

WHY custom (session 184, resolving the s183 tooling fork): the canonical
`surface_dynamics.CylinderDiagram` is Sage-only (import sage.all throughout, no pip
route; env has no Sage/conda, Python 3.13). But surface_dynamics would only VALIDATE an
already-assembled diagram -- our billiard tracer must produce the ribbon graph (cyclic
prong order + top/bottom SC assignment) regardless -- so the only thing it saves is the
genus/stratum checker below, which we already own a single-cylinder version of
(`o1_stratum.stratum_of_perm`, corner-UF, self-tested). This generalises that to the
MULTI-cylinder case, which is exactly the (ii-c) gate object: the orientation cover S_alpha
is an abelian (translation) surface, so this model applies directly to the gate.

MODEL (translation / abelian differential).  Saddle connections labelled 0..n-1, each
oriented rightward.  Each label appears exactly ONCE on some cylinder bottom and exactly
ONCE on some cylinder top (the segment is the top of the cylinder below it and the bottom
of the cylinder above it).  A cylinder = (bottom cyclic word, top cyclic word).
  b = permutation whose cycles are the bottom words;  t = ditto for tops.
  singularity permutation  sigma = t^{-1} o b ;  s = #cycles(sigma).
  Euler (cut each annulus along a vertical -> k disk faces, n+k edges, s vertices):
      chi = s - (n+k) + k = s - n  =>  g = (n - s + 2)/2 .
  A sigma-cycle of length L is a singularity of cone angle 2*pi*L, order = L - 1
      (L=1 => order 0 = regular marked point).  sum(order) = n - s = 2g-2.

Validated (see _selftest): torus H(0), H(2) 1-cyl, H(1,1) 1-cyl, 2-cyl torus H(0,0).

NOTE (scope): TRANSLATION surfaces only.  The half-translation BASE Q(P-2,Q-P-2,-1^Q)
needs signed/pole handling -- but its stratum is already independently confirmed
(s182-G1 Euler + prong-order = prongs-2), and the (ii-c) GATE is on the abelian cover, so
this suffices for the gate.  Do NOT feed a quadratic-differential diagram here unlabelled.
                                             f_leg_covariance.md A9 Phase-1 (ii-c); session 184
"""


def _perm_from_cycles(cycles, n):
    """Build a permutation (list p, p[x]=image) from a list of cyclic words on 0..n-1."""
    p = [None] * n
    seen = set()
    for cyc in cycles:
        k = len(cyc)
        for i in range(k):
            x = cyc[i]
            if x in seen:
                raise ValueError(f"label {x} appears twice among cycles")
            seen.add(x)
            p[x] = cyc[(i + 1) % k]
    if len(seen) != n or any(v is None for v in p):
        raise ValueError(f"cycles must be a permutation of 0..{n-1}; got labels {sorted(seen)}")
    return p


def _cycles_of(p):
    n = len(p)
    seen = [False] * n
    out = []
    for x in range(n):
        if seen[x]:
            continue
        cyc = []
        y = x
        while not seen[y]:
            seen[y] = True
            cyc.append(y)
            y = p[y]
        out.append(cyc)
    return out


class CylinderDiagram:
    """A translation-surface cylinder diagram: cylinders = (bottom word, top word)."""

    def __init__(self, cylinders):
        """cylinders: list of (bottom, top), each a sequence (cyclic word) of SC labels."""
        self.cylinders = [(tuple(b), tuple(t)) for (b, t) in cylinders]
        bot_labels = [x for (b, _) in self.cylinders for x in b]
        top_labels = [x for (_, t) in self.cylinders for x in t]
        if sorted(bot_labels) != sorted(top_labels):
            raise ValueError("each SC must appear once on a bottom and once on a top; "
                             f"bottoms={sorted(bot_labels)} tops={sorted(top_labels)}")
        self.n = len(bot_labels)
        # normalise labels to 0..n-1
        labels = sorted(set(bot_labels))
        if labels != list(range(self.n)):
            raise ValueError(f"labels must be 0..{self.n-1}, got {labels}")
        self.b = _perm_from_cycles([b for (b, _) in self.cylinders], self.n)
        self.t = _perm_from_cycles([t for (_, t) in self.cylinders], self.n)

    def _sigma(self):
        """singularity permutation sigma = t^{-1} o b."""
        tinv = [0] * self.n
        for x in range(self.n):
            tinv[self.t[x]] = x
        return [tinv[self.b[x]] for x in range(self.n)]

    def singularities(self):
        """list of cone-orders (one per singularity), including 0 for marked points."""
        return sorted((len(c) - 1 for c in _cycles_of(self._sigma())), reverse=True)

    def num_singularities(self):
        return len(_cycles_of(self._sigma()))

    def genus(self):
        s = self.num_singularities()
        g2 = self.n - s + 2
        return g2 // 2 if g2 % 2 == 0 else None

    def stratum(self):
        """(genus, [orders>0], n_marked_points).  H(orders) notation."""
        orders = self.singularities()
        pos = [o for o in orders if o > 0]
        return self.genus(), pos, sum(1 for o in orders if o == 0)

    def stratum_str(self):
        g, pos, m = self.stratum()
        inside = ",".join(str(o) for o in pos) if pos else "0"
        tail = f" +{m} marked" if m else ""
        return f"g={g} H({inside}){tail}  [sum(orders)={sum(pos)}, 2g-2={2*g-2 if g is not None else '?'}]"


def _selftest():
    cases = [
        ("torus H(0)          ", [((0,), (0,))],                       1, []),
        ("H(2) 1-cyl          ", [((0, 1, 2), (0, 2, 1))],             2, [2]),
        ("H(1,1) 1-cyl        ", [((0, 1, 2, 3), (0, 3, 2, 1))],       2, [1, 1]),
        ("2-cyl torus H(0,0)  ", [((0,), (1,)), ((1,), (0,))],         1, []),
    ]
    print("SELF-TESTS (cyl_diagram):")
    ok = True
    for name, cyls, exp_g, exp_orders in cases:
        cd = CylinderDiagram(cyls)
        g, pos, m = cd.stratum()
        good = (g == exp_g and pos == exp_orders)
        ok = ok and good
        # cross-check Euler: sum(orders) must equal 2g-2
        euler_ok = (sum(cd.singularities()) == 2 * g - 2)
        flag = "OK " if (good and euler_ok) else "FAIL"
        print(f"  [{flag}] {name} -> {cd.stratum_str()}"
              f"   (expected g={exp_g} orders={exp_orders})")
    print("ALL PASS" if ok else "*** FAILURES ***")
    return ok


if __name__ == "__main__":
    _selftest()
