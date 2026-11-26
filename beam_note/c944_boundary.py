"""[Deposit] Remark 5.3: the period-944 non-perpendicular cylinder and the adjacent period-7338 perpendicular cylinder, certified exactly, and their common boundary. Paper: "Infinitely many perpendicular beams in right triangles at all but countably many irrational angles". Research name: probes/s963_c944_boundary.py, copied verbatim except imports (deposit names), paths (data/ under beam_note/); the cylinder record is read from ../mirror_A/tro05_record.json and this paragraph. Run through run_all.py.

s963: what lies across the boundary of Mirror A's period-944 non-perpendicular cylinder (Remark 6.10)?

QUESTION.  [Tro05] Lemma 12's proof (and Thm 3's, which reuses it) claims that every periodic loop with perpendicular
periodic orbits on one side has them on BOTH sides (s962 audits: the step is unproved for links not meeting the leg,
gap G2, [LIT-060]).  Mirror A Thm 6.1 / Rem 6.10 put a non-perpendicular periodic cylinder C (period 944) in the
perpendicular invariant surface at cos 2a = -(sqrt 7 - 1)/4 (the record's own frame: class perp OB, a at A, so the
leg OB is adjacent to the SMALLER angle pi/2 - a ~ 32.85 deg, inside Troubetzkoy's range).  If a perpendicular
cylinder P sits directly across a boundary loop of C, that loop has perpendicular periodic orbits on one side and the
non-perpendicular C on the other: the claimed step is then FALSE at an angle in his range, not merely unproved.

METHOD.  (1) certify C from its word (deposit `lib/mirror.py`, exact in Q(2 cos a)); control: the midpoint of its
launch interval on side B, traced at 50 digits, closes with the 944-word and no perpendicular hit.  (2) for each end of
the launch interval, start leaves at distance eps = 1e-6, 1e-12, 1e-20 outside it, same direction, trace forward and
backward (cap 20 000 bounces each) until a perpendicular hit on OB, a vertex, a closing, or the cap.  (3) if both
directions hit OB perpendicularly with the same word for every eps, build the perpendicular orbit's periodic word,
rotate it to start at the launch point, certify it exactly, and test EXACT adjacency in K: hi' = lambda lo or
lo' = lambda hi, lambda = beta'/beta (real, in K).
Not pre-registered (a structural check on one certified instance). PRIOR ART: rulings --grep 944 / Troubetzkoy
(-> [LIT-060], [LIT-038]); mirror_classes.md s685 (7) "the rest of the leg" (perp-vs-nonperp shares, not adjacency);
deposit mirror_A/s6_4_troubetzkoy.py (certifies C only). DATA: data/s963_c944_boundary.json (written here).
Run: PYTHONPATH=.:engine:probes .venv/bin/python3.13 probes/s963_c944_boundary.py
"""
import _setup  # noqa: F401  (deposit: this directory and ../lib on the path, cwd = beam_note/)
import json
import sys

import mpmath as mp

import mirror as m                                          # noqa: E402

DPS = 50
CAP = 20000
REC = '../mirror_A/tro05_record.json'
OUT = 'data/c944_boundary.json'
SIDE = m.Triangle.SIDE                                     # {'A': ('O','A'), 'B': ('O','B'), 'H': ('A','B')}


def frame(y):
    a = mp.acos(mp.mpf(y)) / 2
    return a, {'O': mp.mpc(0), 'A': mp.mpc(1), 'B': mp.mpc(0, mp.tan(a))}


def reflect(d, e):
    e = e / abs(e)
    return e * e * mp.conj(d)


def trace(V, p, d, last, cap, p0=None, d0=None):
    """Billiard trace from p (on side `last`) in direction d. Stops at: a perpendicular hit on B ('PERP'),
    a vertex ('VERTEX'), a return to (p0, d0) ('CLOSED'), or the cap. Returns (status, word, info)."""
    tol = mp.mpf(10) ** (-DPS // 2)
    word, minclr = [], mp.inf
    for n in range(cap):
        best = None
        for ch, (i, j) in SIDE.items():
            if ch == last:
                continue
            P, Q = V[i], V[j]
            e = Q - P
            den = mp.im(mp.conj(d) * e)
            if den == 0:
                continue
            t = mp.im(mp.conj(e) * (P - p)) / mp.im(mp.conj(e) * d)
            u = mp.im(mp.conj(d) * (p - P)) / den
            if t > tol and -tol <= u <= 1 + tol and (best is None or t < best[0]):
                best = (t, ch, e, u)
        if best is None:
            return 'LOST', word, {}
        t, ch, e, u = best
        p = p + t * d
        word.append(ch)
        clr = min(u, 1 - u) * abs(e)
        minclr = min(minclr, clr)
        if clr < mp.mpf(10) ** -30:
            return 'VERTEX', word, {'n': n + 1}
        if ch == 'B' and abs(mp.re(mp.conj(d) * e)) < tol * abs(e):     # d normal to OB
            return 'PERP', word, {'n': n + 1, 'height': float(mp.im(p)), 'minclr': float(minclr)}
        d = reflect(d, e)
        last = ch
        if p0 is not None and ch == 'B' and abs(p - p0) < tol and abs(d - d0) < tol:
            return 'CLOSED', word, {'n': n + 1, 'minclr': float(minclr)}
    return 'CAP', word, {'n': cap}


def main():
    mp.mp.dps = DPS
    r = json.load(open(REC))['rows'][0]
    y = float(r['y'])
    K = m.t_field(r['minpoly'], y)
    c = m.certify(K, r['word'])
    assert c['ok'] and not c['perpendicular'] and c['period'] == 944, c.get('why')
    l1, l2 = m.launch_interval(c, 'B')
    beta = c['beta']
    th = mp.mpf(mp.arg(beta.complex(40)))
    yx = (1 - mp.sqrt(7)) / 4                                  # exact root of 8y^2 - 4y - 3 (the record's y)
    assert abs(yx - mp.mpf(r['y'])) < 1e-15 and abs(8 * yx**2 - 4 * yx - 3) < mp.mpf(10) ** (-DPS + 2)
    a, V = frame(yx)
    L1, L2 = mp.mpf(l1.mp(40)), mp.mpf(l2.mp(40))
    print(f"C certified: period 944, class {c['cls']}, launch interval on OB = ({mp.nstr(L1, 12)}, {mp.nstr(L2, 12)}), "
          f"direction {mp.nstr(th, 10)}, width {c['width']:.3e}")
    d0 = mp.expj(th)
    out = {'y': r['y'], 'launch': [str(L1), str(L2)], 'theta': str(th), 'rows': []}

    # control: midpoint closes with the 944-word, no perpendicular hit
    pm = V['O'] + (L1 + L2) / 2 * (V['B'] - V['O'])
    st, w, info = trace(V, pm, d0, 'B', 2000, p0=pm, d0=d0)
    ok = st == 'CLOSED' and ''.join(w) == r['word']
    print(f"control (midpoint): {st} after {info.get('n')} bounces, word == record: {ok}")
    assert ok
    out['control'] = {'status': st, 'n': info.get('n'), 'word_matches': ok}

    for end, base, sgn in (('l1', L1, -1), ('l2', L2, +1)):
        for eps in ('1e-6', '1e-12', '1e-20'):
            lam = base + sgn * mp.mpf(eps)
            row = {'end': end, 'eps': eps, 'lam': str(lam)}
            if not (0 < lam < 1):
                row['status'] = 'outside side'
                out['rows'].append(row); print(row); continue
            p = V['O'] + lam * (V['B'] - V['O'])
            sf, wf, inf_f = trace(V, p, d0, 'B', CAP, p0=p, d0=d0)
            db = reflect(-d0, V['B'] - V['O'])                     # backward: bounce at the start point first
            sb, wb, inf_b = trace(V, p, db, 'B', CAP)
            row.update(fwd=sf, fwd_n=inf_f.get('n'), fwd_h=inf_f.get('height'),
                       bwd=sb, bwd_n=inf_b.get('n'), bwd_h=inf_b.get('height'),
                       fwd_word=''.join(wf) if sf in ('PERP', 'CLOSED') else None, bwd_word=''.join(wb) if sb == 'PERP' else None)
            out['rows'].append(row)
            print({k: v for k, v in row.items() if not k.endswith('_word')})

    # exact adjacency for each end whose three eps agree on a perpendicular word
    out['adjacency'] = {}
    for end in ('l1', 'l2'):
        rows = [x for x in out['rows'] if x['end'] == end and x.get('fwd') == 'PERP' and x.get('bwd') == 'PERP']
        if len(rows) < 3 or len({(x['fwd_word'], x['bwd_word']) for x in rows}) != 1:
            out['adjacency'][end] = 'no stable perpendicular word'
            continue
        F, Bk = list(rows[0]['fwd_word']), list(rows[0]['bwd_word'])
        Lw = Bk[::-1] + ['B'] + F                                   # P1 ... start bounce ... P2
        Wp = Lw + Lw[::-1][1:-1]                                    # the perpendicular orbit, one period
        k = len(Bk) + 1
        Wr = ''.join(Wp[k:] + Wp[:k])                               # rotated: T_0 segment = the launch segment
        cp = m.certify(K, Wr)
        res = {'period': len(Wr), 'ok': cp['ok'], 'perpendicular': cp.get('perpendicular'), 'why': cp.get('why')}
        if cp['ok']:
            lam_ = (cp['beta'] * beta.conj()).re() / beta.norm()
            adj_lo = (cp['hi'] - lam_ * c['lo']).iszero()
            adj_hi = (cp['lo'] - lam_ * c['hi']).iszero()
            pl1, pl2 = m.launch_interval(cp, 'B')
            res.update(lam_positive=lam_.sign() > 0, shares_C_lo=adj_lo, shares_C_hi=adj_hi,
                       launch=[float(pl1.mp(20)), float(pl2.mp(20))], width=cp['width'])
        out['adjacency'][end] = res
        print(end, res)
    # the other side: a closed non-perpendicular word -> certify and test adjacency to C's hi (or lo)
    for end in ('l1', 'l2'):
        rows = [x for x in out['rows'] if x['end'] == end and x.get('fwd') == 'CLOSED']
        if len(rows) < 3 or len({x['fwd_word'] for x in rows}) != 1:
            continue
        Wc = rows[0]['fwd_word']                                     # starts after the launch, ends with 'B'
        cc = m.certify(K, Wc)
        res = {'period': len(Wc), 'same_as_C': Wc == r['word'], 'ok': cc['ok'],
               'perpendicular': cc.get('perpendicular'), 'cls': cc.get('cls'), 'why': cc.get('why')}
        if cc['ok']:
            lam_ = (cc['beta'] * beta.conj()).re() / beta.norm()
            res.update(lam=float(lam_.mp(20)), shares_C_hi=(cc['lo'] - lam_ * c['hi']).iszero(),
                       shares_C_lo=(cc['hi'] - lam_ * c['lo']).iszero(), width=cc['width'],
                       launch=[float(v.mp(20)) for v in m.launch_interval(cc, 'B')])
        res['cyclic_rotation_of_C'] = Wc in r['word'] * 2
        res['cyclic_rotation_of_reversed_C'] = Wc in r['word'][::-1] * 2      # C's time-reversal
        out['adjacency'][end + '_nonperp'] = res
        print(end, 'non-perp neighbour:', res)
    json.dump(out, open(OUT, 'w'), indent=1, default=str)
    print('wrote', OUT)


if __name__ == '__main__':
    main()
