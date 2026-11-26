"""[Deposit] Appendix A: the appendix's list of beams (end points and words), generated from the certificates and re-checked. Paper: "Infinitely many perpendicular beams in right triangles at all but countably many irrational angles". Research name: probes/s1011_beam_codes_tex.py, copied verbatim except imports (deposit names), paths (data/ under beam_note/), one sentence about correspondence (removed) and this paragraph. Run through run_all.py.

s1011 — LaTeX appendix of the beam codes at the two unit angles of the beam note (tex/beaminf, Appendix A).

Source: the certificates of the two Lemma 12 counterexamples (produced by probes/s973_exact_partition.py,
audited audits/s974_audit_tro_l12.md). This script only re-formats it, and re-checks what it prints:

  (1) exact: the end points chain (left end of beam k = right end of beam k-1), start at 0, end at 1, increase — in Q(sigma),
      reduced modulo the minimal polynomial;
  (2) exact: the beams pair into fans (Prop. "Fans" of the note): b_{2k-1} = (b_{2k-2} + b_{2k})/2;
  (3) the printed floats of each interval agree with the exact end points at the real root;
  (4) independent 50-digit trace: from the midpoint of each beam, the trajectory leaving L_1 horizontally hits exactly the
      printed word of sides and returns to L_1 perpendicularly at the last letter (no vertex within 1e-20 of a hit).

Letters as in the certificate (Troubetzkoy's): L = the leg L_1 = RA, M = the other leg L_2 = RO, H = the hypotenuse.
Output: data/scratch/s1011/appendix_codes.tex (pasted into tex/beaminf/beaminf.tex).
Run: .venv/bin/python3.13 probes/s1011_beam_codes_tex.py
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = beam_note/)
import re
from pathlib import Path

import mpmath as mp
import sympy as sp

SRC = Path("data/certificates_lemma12.txt")
OUT = Path("data/appendix_codes.tex")
s = sp.Symbol("sigma")
mp.mp.dps = 50


def parse(text):
    tris = []
    for block in re.split(r"\nTriangle \d", text)[1:]:
        poly = re.search(r"of\s+(.+?)\s+\(a unit", block).group(1).replace("^", "**")
        root = float(re.search(r"real root near ([\d.]+)", block).group(1))
        beams = []
        for m in re.finditer(r"beam (\d+):\s+\(([\d.]+), ([\d.]+)\),\s+n = (\d+) reflections.*?\n\s+word: (\w+)\n"
                             r"\s+left end\s+= (.+?)\n\s+right end = (.+?)\n", block):
            k, flo, fhi, n, word, le, re_ = m.groups()
            beams.append(dict(k=int(k), flo=float(flo), fhi=float(fhi), n=int(n), word=word,
                              le=sp.sympify(le.replace("^", "**"), locals={"sigma": s}),
                              re=sp.sympify(re_.replace("^", "**"), locals={"sigma": s})))
        tris.append(dict(poly=sp.Poly(sp.sympify(poly, locals={"x": s}).subs(sp.Symbol("x"), s), s), root=root, beams=beams))
    return tris


def red(e, P):
    """Reduce a polynomial in sigma (rational coefficients) modulo P."""
    return sp.Poly(sp.expand(e), s).rem(P).as_expr()


def tex_ql(e):
    """(numerator polynomial)/d, the remark's style."""
    e = sp.together(sp.expand(e))
    num, den = sp.fraction(e)
    num = sp.expand(num)
    body = sp.latex(num, order="lex").replace(r"\sigma", r"\sigma")
    if den == 1:
        return body
    if len(sp.Add.make_args(num)) == 1:
        return f"{body}/{den}"
    return f"({body})/{den}"


def word_tex(w):
    groups = [w[i:i + 10] for i in range(0, len(w), 10)]
    return r"\texttt{" + r"\,\allowbreak ".join(groups) + "}"


def trace(alpha_his, s0, nmax):
    """Leave (0, s0) horizontally to the left in the triangle R=(0,0), A=(0,1), O=(-tan a, 0); record sides hit
    until the trajectory meets L (x = 0) moving perpendicularly (direction (-1,0) reversed -> it arrives with dy = 0)."""
    tO = -mp.tan(alpha_his)
    x, y, dx, dy = mp.mpf(0), mp.mpf(s0), mp.mpf(-1), mp.mpf(0)
    word, mind = [], mp.mpf(1)
    for _ in range(nmax):
        cands = []
        if dx > 0:
            cands.append(((0 - x) / dx, "L"))
        if dy < 0:
            cands.append(((0 - y) / dy, "M"))
        # hypotenuse: x/tO + y = 1
        g = dx / tO + dy
        if g > 0:
            cands.append(((1 - x / tO - y) / g, "H"))
        t, side = min((c for c in cands if c[0] > mp.mpf(10) ** -40), key=lambda c: c[0])
        x, y = x + t * dx, y + t * dy
        word.append(side)
        # distance to the three vertices
        for vx, vy in ((0, 0), (0, 1), (tO, 0)):
            mind = min(mind, mp.sqrt((x - vx) ** 2 + (y - vy) ** 2))
        if side == "L":
            if abs(dy) < mp.mpf(10) ** -30:
                return "".join(word), mind
            dx = -dx
        elif side == "M":
            dy = -dy
        else:  # reflect in the hypotenuse, direction vector (-tO, 1)/|.|
            ux, uy = -tO, mp.mpf(1)
            nrm = mp.sqrt(ux * ux + uy * uy)
            ux, uy = ux / nrm, uy / nrm
            d = dx * ux + dy * uy
            dx, dy = 2 * d * ux - dx, 2 * d * uy - dy
    return "".join(word) + "...", mind


def main():
    tris = parse(SRC.read_text())
    assert len(tris) == 2 and [len(t["beams"]) for t in tris] == [8, 10], [len(t["beams"]) for t in tris]
    out = []
    for ti, T in enumerate(tris):
        P, bs = T["poly"], T["beams"]
        root = mp.findroot(lambda z: sp.lambdify(s, P.as_expr(), "mpmath")(z), T["root"])
        alpha_his = mp.asin(mp.sqrt(root) / 2)  # sigma = 4 sin^2(alpha_his)
        ends = [bs[0]["le"]] + [b["re"] for b in bs]
        # (1) chaining and range, exact
        assert red(ends[0], P) == 0 and red(ends[-1] - 1, P) == 0
        for k in range(1, len(bs)):
            assert red(bs[k]["le"] - bs[k - 1]["re"], P) == 0, (ti, k)
        vals = [mp.mpf(sp.lambdify(s, e, "mpmath")(root)) for e in ends]
        assert all(vals[i] < vals[i + 1] for i in range(len(vals) - 1))
        # (2) fans
        for k in range(1, len(bs) // 2 + 1):
            assert red(2 * ends[2 * k - 1] - ends[2 * k - 2] - ends[2 * k], P) == 0, (ti, k)
        # (3) printed floats
        for b, lo, hi in zip(bs, vals[:-1], vals[1:]):
            assert abs(lo - b["flo"]) < 1e-9 and abs(hi - b["fhi"]) < 1e-9
        # (4) independent trace
        for b, lo, hi in zip(bs, vals[:-1], vals[1:]):
            w, mind = trace(alpha_his, (lo + hi) / 2, 400)
            assert w == b["word"] and len(w) == b["n"], (ti, b["k"], w[:40], b["word"][:40])
            assert mind > mp.mpf(10) ** -20
        print(f"triangle {ti + 1}: {len(bs)} beams, chain/fans/floats/trace all OK; alpha_ours = "
              f"{mp.nstr(mp.pi / 2 - alpha_his, 12)}")
        # LaTeX
        polytex = sp.latex(P.as_expr(), order="lex")
        head = ("First" if ti == 0 else "Second") + " angle:"
        out.append(r"\medskip\noindent\emph{" + head + "} "
                   + rf"$\sigma = {mp.nstr(root, 10)}\ldots$ is the root of ${polytex} = 0$, "
                   + rf"$\alpha = {mp.nstr(mp.pi / 2 - alpha_his, 9)}\ldots$, $n(\alpha) = {len(bs)}$. End points:"
                   + "\n\\[\\begin{array}{ll}")
        rows = [rf"b_{{{i}}} = {tex_ql(red(e, P))} & ({mp.nstr(v, 7)})" for i, (e, v) in enumerate(zip(ends, vals))]
        out.append(" \\\\\n".join(rows) + "\n\\end{array}\\]")
        out.append(r"\begin{itemize}\setlength\itemsep{1pt}\small\raggedright")
        for b in bs:
            out.append(rf"\item[{b['k']}.] $(b_{{{b['k'] - 1}}}, b_{{{b['k']}}})$, ${b['n']}$ reflections: {word_tex(b['word'])}")
        out.append(r"\end{itemize}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(out) + "\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
