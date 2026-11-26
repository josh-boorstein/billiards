"""[Deposit] Proposition A.1 and Appendix A: the exact certificates of the two finite partitions (8 and 10 beams), written as text. Paper: "Infinitely many perpendicular beams in right triangles at all but countably many irrational angles". Research name: probes/s975_cert_export.py, copied verbatim except imports (deposit names), paths (data/ under beam_note/; input: the two angles' rows of the research census record, shipped as data/lemma12_angles.json), the words naming a correspondence (removed) and this paragraph. Run through run_all.py.

s975 (queue s962-d): export the exact certificates of the two finite L1 partitions inside (π/4, π/3) — the counterexamples to
the STATEMENT of [Tro05] Lemma 12 ([TRO-L12-FALSE], general_triangles.md §s973 (3), audited s974) — as a self-contained text file
for the beam note's Remark 5.3.  Re-runs s973's exact certifier (probes/s973_exact_partition.py) and
writes, per angle: the minimal polynomial of σ = 4cos²α with its root, α in both conventions, each beam's full word and its exact
end points in K = ℚ(σ)(sin 2α), and the tiling check.  Not PRE-REGISTERED (an export).  PRIOR ART: audits/s963_certificates_cos2a_sqrt7.txt
(the format).  DATA: data/s973_coincidence_census_K20_certified.json (the angles).  → audits/s975_certificates_lemma12.txt
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = beam_note/)
import json
import math
import os
import sys

from mpmath import mp

import exact_partition as X                                             # noqa: E402

OUT = 'data/certificates_lemma12.txt'

HEAD = """Two irrational right triangles in which the perpendicular beam partition of a leg is FINITE (certificate data)
=====================================================================================================

Convention of Troubetzkoy, "Periodic billiard orbits in right triangles" (2005): the smaller acute angle is alpha, at the vertex A;
the beam is perpendicular to the leg L adjacent to A.  Both angles below have alpha in (pi/6, pi/4), the range of his Lemma 12.

Coordinates: right angle R = (0, 0), A = (0, 1), the third vertex O = (-tan(alpha), 0).  Sides: L = RA (the leg, length 1),
M = RO (the other leg), H = OA (the hypotenuse).  Heights on L are measured from R.  Write beta = pi/2 - alpha (the angle at O),
sigma = 4 cos^2(beta) = 4 sin^2(alpha), and t = sin(2 beta) = sin(2 alpha); then cos(2 beta) = sigma/2 - 1 and t^2 = 1 - (sigma/2 - 1)^2.

A trajectory leaves the point (0, s) of L perpendicularly, i.e. in the direction (-1, 0).  A BEAM is a maximal open interval of
starting heights s on which the sequence of sides hit, up to and including the first return to L at a right angle (before any
vertex), is defined and constant.  For each beam below: its word (the sides hit; the last letter is the perpendicular return to L,
after which the path retraces itself, so the orbit is periodic of period 2n) and its two end points, exactly, as elements
a + b*t with a, b in Q(sigma).  CHECK: along a fixed word the position of the k-th hit is affine in s, so the set of s realising the
word is an interval cut out by the conditions 0 <= u_k(s) <= 1 (u_k = position of the k-th hit along its side); its end points are
the ones listed.  The beams, left to right, TILE [0, 1]: each beam's left end EQUALS the previous beam's right end exactly, the
first starts at 0 and the last ends at 1.  So the first return to L is an interval exchange of finitely many intervals.  The
equalities and the side choices were checked in exact arithmetic in Q(sigma)(t) (signs from 120-digit evaluation, smallest
|value| tested 0.011).  Code: probes/s973_exact_partition.py (available on request).
"""


def export(name, coeffs, sig_guess):
    m = sum(X.sp.Integer(c) * X.S ** (len(coeffs) - 1 - i) for i, c in enumerate(coeffs))
    rts = [r for r in X.sp.Poly(m, X.S).nroots(n=130) if r.is_real and abs(r - sig_guess) < 1e-6]
    sig = mp.mpf(str(rts[0]))
    beta = mp.acos(mp.sqrt(sig) / 2)
    seeds, _ = X.float_seeds(float(2 * beta / mp.pi))
    SIDES = X.setup(m, sig)
    rows, cur, ok = [], X.E(0), True
    for s in sorted(seeds):
        w, lo, hi, _ = X.trace_branch(SIDES, X.E(X.Fraction(s).limit_denominator(10 ** 12)))
        ok &= (lo - cur).iszero()
        rows.append((w, lo, hi))
        cur = hi
    ok &= (cur - X.E(1)).iszero()
    assert ok, name
    al = mp.pi / 2 - beta
    poly = X.sp.sstr(X.sp.Poly(coeffs, X.sp.Symbol('x')).as_expr()).replace('**', '^')
    out = [f"\n{name}\n" + '-' * len(name),
           f"sigma = 4 sin^2(alpha) is the real root near {mp.nstr(sig, 15)} of   {poly}  (a unit: constant term +-1)",
           f"alpha = {mp.nstr(al, 12)} rad = {mp.nstr(al * 180 / mp.pi, 8)} deg   (alpha/pi irrational: the polynomial has a root in (-1, 0), while every conjugate of 2 + 2cos(2 pi r), r rational, lies in [0, 4])",
           f"number of beams: {len(rows)};  tiling of [0, 1] checked exactly: {ok}"]
    rename = {'L1': 'L', 'L2': 'M', 'H': 'H'}
    for i, (w, lo, hi) in enumerate(rows):
        word = ''.join(rename[x] for x in w)
        out.append(f"\nbeam {i + 1}:  ({mp.nstr(lo.num(), 12)}, {mp.nstr(hi.num(), 12)}),  n = {len(w)} reflections to the perpendicular return")
        out.append(f"  word: {word}")
        for tag, v in (('left end ', lo), ('right end', hi)):
            a, b = v.a.as_expr(), v.b.as_expr()
            sa = X.sp.sstr(X.sp.expand(a)).replace('S', 'sigma').replace('**', '^')
            sb = X.sp.sstr(X.sp.expand(b)).replace('S', 'sigma').replace('**', '^')
            out.append(f"  {tag} = {sa}" + ('' if b == 0 else f"  +  ({sb})*t"))
    return '\n'.join(out)


if __name__ == '__main__':
    mp.dps = 120
    cert = json.load(open('data/lemma12_angles.json'))
    pick = lambda xx: [c for c in cert if abs(c['x'] - xx) < 1e-6][0]
    parts = [HEAD]
    for name, xx in (('Triangle 1 (alpha ~ 40.53 deg)', 0.5497162364), ('Triangle 2 (alpha ~ 34.18 deg)', 0.6202381746)):
        c = pick(xx)
        parts.append(export(name, c['minpoly'], c['sigma1']))
    open(OUT, 'w').write('\n'.join(parts) + '\n')
    print(open(OUT).read()[:6000])
