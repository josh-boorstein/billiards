"""s302_model.py -- EXPLORATORY: a pure-arithmetic simulator of the a2=1 fan chain,
validated cell-by-cell against the exact ring.

THE MODEL (derived s302, from the sec.s255 exit rule + the point-reflection escape).
A fan is a pinwheel of triangle copies about its O-corner at height c, boundary vertices
at signed angles Phi + eps*i*alpha (i = 1,2,3,...), radius cot(alpha) at EVEN i (an R
vertex) and csc(alpha) at ODD i (an A vertex), so the vertex heights are
        y_i = c + rho_i * sin(Phi + eps*i*alpha).
The fan runs from i = j0 and EXITS at the first i whose vertex is on the ray's
centre-side (sec.s255's convexity dichotomy).  The crossed L1 side is [v_{j*-1}, v_{j*}]
and its R endpoint is the vertex at the largest EVEN index <= j*.  Because the triangle's
right angle sits at R, OR is PERPENDICULAR to L1, so reflecting in L1 is the POINT
reflection through R:
        c' = 2 y_R - c,   Phi' = (Phi + eps*jR*alpha) + pi,   eps' = -eps,
and the next fan starts at j0' = 1 + (j* mod 2) (letter-parity continuity).
Letters: index i even -> 'H', odd -> '2'; an L1 letter '.' between fans.

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
     .venv/bin/python3.13 -u probes/s302_model.py [validate|scan]
"""
import json
import os
import sys
from math import gcd

import mpmath as mp

mp.mp.dps = 40


# --------------------------------------------------------------- basic constants

def base(P, r, a1):
    Q = P * a1 + r
    th = mp.pi / (2 * Q)
    al = P * th
    cot = mp.cos(al) / mp.sin(al)
    csc = 1 / mp.sin(al)
    psi = 2 * r * th
    H = cot * mp.sin(psi)
    sP = 2 * mp.cos(al) ** 2
    Hp = mp.sin((2 * r - P) * th) / mp.sin(al)
    return dict(P=P, r=r, a1=a1, Q=Q, th=th, al=al, cot=cot, csc=csc, psi=psi,
                H=H, sP=sP, Hp=Hp, S=2 * H - 1, D=sP - H)


def rho(B, i):
    return B['cot'] if i % 2 == 0 else B['csc']


# ------------------------------------------------------------------ the simulator

def simulate(P, r, a1, s, max_fans=400):
    """Fan chain of the ray at height s.  Returns dict with fans, word, headon."""
    B = base(P, r, a1)
    al, cot = B['al'], B['cot']
    s = mp.mpf(s)
    c, Phi, eps, j0 = mp.mpf(0), mp.mpf(0), 1, 2
    depth = 0                      # depth of the last letter written
    fans, letters = [], []
    headon = None
    for f in range(max_fans):
        # find the exit index
        jstar = None
        for i in range(j0, 4 * a1 + 6):
            y = c + rho(B, i) * mp.sin(Phi + eps * i * al)
            if (y < s) if (c < s) else (y > s):
                jstar = i
                break
        if jstar is None:
            return dict(ok=False, why='no-exit', fans=fans)
        jR = jstar - (jstar % 2)                       # largest even <= jstar
        Om = Phi + eps * jR * al
        yR = c + cot * mp.sin(Om)
        blk = ''.join('H' if i % 2 == 0 else '2' for i in range(j0, jstar + 1))
        entry = depth if f else 0                      # O-split depth (fan 1 has none)
        fans.append(dict(f=f + 1, c=c, Phi=Phi, eps=eps, j0=j0, jstar=jstar,
                         jR=jR, blk=blk, nH=blk.count('H'), entry=entry))
        letters.append(blk)
        depth = entry + len(blk)
        # HEAD-ON: the R VERTEX sits at the centre height, y_R = c -- equivalently
        # sin(Omega) = 0, equivalently c' = c, which is how [OPS-025] states it.
        # ⚠ AMENDED s341 ([WFLOOR-115]): this used to test the EXIT vertex,
        # sin(Phi + eps*jstar*al).  Retrace happens when the ray meets L1
        # PERPENDICULARLY and OR is perpendicular to L1, so the condition belongs to the
        # R vertex jR; the two coincide exactly when jstar is EVEN, and at odd P the
        # phi-subgroup forces that, so the old form was right on everything odd-P ever
        # measured.  At even P (hence Q odd) an ODD solution exists and puts an A vertex
        # at the centre height WITHOUT making L1 vertical -- a false head-on.
        if abs(mp.sin(Om)) < mp.mpf(10) ** -25:
            headon = depth
            break
        letters.append('.')
        depth += 1
        c, Phi, eps = 2 * yR - c, Om + mp.pi, -eps
        j0 = 1 + (jstar % 2)
    fwd = ''.join(letters)
    if headon is None:
        return dict(ok=False, why='no-headon', fans=fans, fwd=fwd)
    # retrace: the head-on L1 is the palindrome CENTRE, so |word| = 2*headon + 1
    word = fwd + '.' + fwd[::-1]
    return dict(ok=True, fans=fans, word=word, headon=headon,
                closure=2 * headon + 1, nfan=2 * len(fans),
                nH=[f['nH'] for f in fans], B=B)


def nu_max(w):
    nu, mx = 0, 0
    for j, ch in enumerate(w, start=1):
        if ch == 'H':
            nu += (-1) ** j
        mx = max(mx, abs(nu))
    return mx


def skeleton_counts(a1, k):
    return [a1] + [a1 - 1] * (2 * k) + [a1]


# ------------------------------------------------------------------- validation

def validate(targets, max_depth=40000):
    from orphan_region import find_orphan_ring
    from s242_m0fan import trace_branch
    out = []
    for (P, r, a1) in targets:
        Q = P * a1 + r
        if gcd(P, Q) != 1:
            continue
        B = base(P, r, a1)
        if B['H'] < 1:
            continue
        leaves = sorted(find_orphan_ring(P, Q, 0.0, 1.0, max_depth=max_depth,
                                         return_leaves=True, leaf_detail=True),
                        key=lambda lf: lf['lo'])
        cert = all(lf.get('reason') == 'closed' for lf in leaves)
        above = [lf for lf in leaves if lf['lo'] >= float(B['Hp']) - 1e-9]
        okw = okl = n = 0
        bad = []
        for lf in above:
            s = 0.5 * (lf['lo'] + lf['hi'])
            br = trace_branch(P, r, a1, s, 40000)
            sim = simulate(P, r, a1, s)
            n += 1
            if sim['ok'] and sim['word'] == br['word']:
                okw += 1
            else:
                bad.append((round(s, 8), len(br['word']),
                            len(sim.get('word', '')) if sim['ok'] else sim['why']))
            if sim['ok'] and sim['closure'] == br['closure']:
                okl += 1
        row = dict(P=P, r=r, a1=a1, Q=Q, cert=bool(cert), n_above=n,
                   word_match=okw, len_match=okl, bad=bad[:4],
                   H=float(B['H']), sP=float(B['sP']), D=float(B['D']),
                   S=float(B['S']), Hp=float(B['Hp']))
        out.append(row)
        print(f"  {P:3d}/{Q:3d} a1={a1} cert={int(cert)} above={n:2d} "
              f"word {okw}/{n}  closure {okl}/{n}  D={float(B['D']):+.5f}"
              + (f"  BAD {bad[:2]}" if bad else ""))
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'validate'
    if mode == 'validate':
        tg = []
        for P in range(3, 12):
            for r in range(P // 2 + 1, P):
                for a1 in range(2, 7):
                    Q = P * a1 + r
                    if Q > 50 or gcd(P, r) != 1 or gcd(P, Q) != 1:
                        continue
                    tg.append((P, r, a1))
        print(f"=== MODEL vs RING, {len(tg)} candidate centres ===")
        rows = validate(tg)
        tot = sum(x['n_above'] for x in rows if x['cert'])
        okw = sum(x['word_match'] for x in rows if x['cert'])
        okl = sum(x['len_match'] for x in rows if x['cert'])
        print(f"\nCERTIFIED: word {okw}/{tot}   closure {okl}/{tot}   "
              f"({sum(1 for x in rows if x['cert'])} centres)")
        os.makedirs('data', exist_ok=True)
        json.dump(rows, open('data/s302_model_validate.json', 'w'))


if __name__ == '__main__':
    main()
