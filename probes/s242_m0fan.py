"""s242_m0fan.py -- DE-RISK for [M0-FAN]: measure the m=0 rung (the CAP/MIR J-pair) the
way s239 measured the m=1 rung (the beat), and test whether the s239 SIX-fan machine
degenerates to a TWO-fan machine there.

WHY THESE PREDICTIONS (stated before running; they are forced by two already-recorded laws,
not guessed).  s239 PROVED the beat branch is 6 O-fans, head-on at depth 12a1-8, closure
(= word length) at 24a1-15 = 2*(12a1-8)+1.  s241 measured the generation ladder
    wlen(m) = 8(2m+1)a1 + 1 - 16 m^2,   #L1(m) = 8m + 3.
Inverting the retrace rule wlen = 2*D_headon + 1 on s241's law gives
    D_headon(m) = 4(2m+1) a1 - 8 m^2,
which at m=1 is EXACTLY s239's proved 12a1-8 (a consistency check of the ladder law against
a proof, not an assumption).  Likewise #L1(m) = 8m+3 gives fans-before-head-on
    = (8m+3-1)/2 + 1 = 4m+2 = 2(2m+1),
which at m=1 is exactly s239's 6.  Reading both at m=0 is the pre-registration below.

PRE-REGISTERED HYPOTHESES (scope: a2=2 cone, 2r<P<3r, gcd(P,r)=1, a3>=3, a1>=5; the
a3=2 family (5,2) and a1=4 are carried as KNOWN-EXCEPTION controls, s241 (3)).

  M0-A   the CAP branch (ray in (H,1)) is exactly 2 O-fans before the head-on.
  M0-B   its head-on is at depth 4*a1.
  M0-C   it closes at depth 8*a1 + 1  (= the word length, s241's m=0 rung).
  M0-D   its word has exactly 3 L1 letters and is a palindrome.
  M0-E   MIR (ray in (S,H)) has the same length and the reversed word  [s230 (Mirror);
         a CONSISTENCY control on the setup, since (Mirror) is already PROVED].
  M0-F   no split of the CAP branch lies in (H,1)  [ALREADY PROVED, s230 Thm (Top) via
         Thm (D)/(L): the capped cell is a single cell.  So this is a SETUP CONTROL --
         a failure here indicts the probe, not the mathematics.]
  M0-G   CONTROL, s239-proved: the beat ray in (B1,S) gives 6 fans, head-on 12a1-8,
         closure 24a1-15, #L1 = 11.

  ** NOT pre-registered as a prediction, MEASURED and reported raw: whether the CAP word
     literally equals the sec.15.7/sec.16.4 skeleton Sigma(a1).  s241 (2) identifies the
     m=0 rung with Sigma(a1) BY LENGTH ALONE, but sec.16.4 explicitly scopes the clean
     single-Sigma(a1) fan to a2>=3 and warns that a2=2 adds sec.13 cot-undercut peak
     words -- and the whole ladder was measured AT a2=2.  So the word identity is an open
     sub-question, not an input.  The probe prints the words and the L1 positions. **

Run: PYTHONPATH=.:engine:archive/scripts_2026-07:archive/scripts_2026-08:probes \
     .venv/bin/python3.13 -u probes/s242_m0fan.py > logs/s242_m0fan.log 2>&1
"""
import json
import math
import sys

from orphan_region import (SIDES_AFTER, SPLIT_VERTEX, SIDE_ENDPOINTS,
                           _hit_side_ring, _exact_vertical)
from exact_state_ring import RingTriangleState
from ring_cyc import RingContext
from s238_fan_anatomy import landmarks

LET = {'H': 'H', 'L2': '2', 'L1': '.'}


def phase_index(st, vtx, Q, cot):
    """Exact direction of `vtx` from the fan centre, in units pi/Q ('R': height c+s_n)
    or pi/2Q ('A': height c+a_n).  atan2, so no n -> 2Q-n aliasing."""
    if vtx == 'R':
        dx, dy = (st.Rx_c - st.Ox_c).to_float() * cot, (st.Ry - st.Oy).to_float()
        return math.atan2(dy, dx) * Q / math.pi
    dx, dy = (st.Ax_c - st.Ox_c).to_float() * cot, (st.Ay - st.Oy).to_float()
    return math.atan2(dy, dx) * 2 * Q / math.pi


def _closed(st, tol=1e-6):
    dx, dy = (st.Rx_c - st.Ox_c).to_float(), (st.Ry - st.Oy).to_float()
    dax, day = (st.Ax_c - st.Ox_c).to_float(), (st.Ay - st.Oy).to_float()
    return ((abs(dx - 1) < tol and abs(dy) < tol
             and abs(dax - 1) < tol and abs(day - 1) < tol)
            or (abs(dx + 1) < tol and abs(dy) < tol
                and abs(dax + 1) < tol and abs(day - 1) < tol))


def trace_branch(P, r, a1, s, max_steps):
    """Full itinerary of the ray at height s: word, fans, head-on depth, closure."""
    Q = P * a1 + r
    ctx = RingContext(P, Q)
    cot = 1.0 / math.tan(ctx.alpha)
    st = RingTriangleState(ctx)
    st.reflect('H')
    prev, depth = 'H', 1
    word = ['H']
    headon, closure, ev_disagree = None, None, []
    fans = [dict(entry=1, c=st.Oy.to_float(), splits=[])]
    while depth < max_steps:
        if depth > 3 and _closed(st):
            closure = depth
            break
        cand = SIDES_AFTER[prev]
        vtx = SPLIT_VERTEX[frozenset(cand)]
        y = st.vertex_y(vtx).to_float()
        n = None if vtx == 'O' else phase_index(st, vtx, Q, cot)
        fans[-1]['splits'].append((depth, vtx, n, y))
        side = _hit_side_ring(st, cand, s, cot)
        if side is None:
            break
        # HEAD-ON detection.  NOTE (s242): orphan_region._exact_vertical returns False on
        # walls whose endpoint x-gap is exactly 0.0 in float -- the ring difference is a
        # NON-CANONICAL representative of zero, so .is_zero() misses it.  We therefore key
        # off the float gap and CROSS-CHECK the result against closure == 2*headon+1
        # (s239's retrace rule), which is an independent integer identity.
        p_, q_ = SIDE_ENDPOINTS[side]
        vgap = abs((st.vertex_xc(p_) - st.vertex_xc(q_)).to_float())
        if vgap == 0.0 and headon is None:
            headon = depth
            ev_disagree.append((depth, side, bool(_exact_vertical(st, side))))
        st.reflect(side)
        word.append(LET[side])
        prev, depth = side, depth + 1
        if side == 'L1':
            fans[-1]['exit'] = depth - 1
            fans[-1]['rho'] = st.Ry.to_float()
            fans.append(dict(entry=depth, c=st.Oy.to_float(), splits=[]))
    return dict(word=''.join(word), fans=fans, headon=headon, closure=closure,
                Q=Q, ev_disagree=ev_disagree)


def sigma(a1):
    """The sec.15.7 / sec.16.4 leading-digit skeleton Sigma(a1), CAP orientation.
    Blocks 2a1-1, 2a1, 2a1, 2a1-1 separated by L1 ('.').  Instantiated at a1=3 this is
    'H2H2H.2H2H2H.H2H2H2.H2H2H' = sec.16.4's Sigma(3) character-for-character; the
    block-swapped orientation below is sec.16.4's a1=2 example, i.e. the J-partner."""
    return '.'.join(['H2' * (a1 - 1) + 'H', '2H' * a1, 'H2' * a1, 'H2' * (a1 - 1) + 'H'])


def sigma_J(a1):
    """The J-partner orientation: blocks 2a1, 2a1-1, 2a1-1, 2a1."""
    return '.'.join(['H2' * a1, 'H2' * (a1 - 1) + 'H',
                     'H2' * (a1 - 1) + 'H', '2H' * a1])


def analyse(P, r, a1):
    H, S, B1 = landmarks(P, r, a1)
    out = dict(P=P, r=r, a1=a1, r2=P - 2 * r, a3=r // (P - 2 * r),
               H=H, S=S, B1=B1)
    rays = dict(CAP=0.5 * (H + 1.0), MIR=0.5 * (S + H), BEAT=0.5 * (B1 + S))
    for nm, s in rays.items():
        cap = (40 * a1 + 80) if nm == 'BEAT' else (20 * a1 + 60)
        t = trace_branch(P, r, a1, s, cap)
        w = t['word']
        # M0-F control: does any split of this branch land strictly inside its own cell?
        lo, hi = (H, 1.0) if nm == 'CAP' else ((S, H) if nm == 'MIR' else (B1, S))
        intr = [(d, v, round(y, 12))
                for f in t['fans'] for (d, v, n, y) in f['splits']
                if lo + 1e-12 < y < hi - 1e-12]
        out[nm] = dict(word=w, wlen=len(w), nL1=w.count('.'),
                       L1pos=[i + 1 for i, ch in enumerate(w) if ch == '.'],
                       nfan_total=len(t['fans']),
                       headon=t['headon'], closure=t['closure'],
                       palindrome=(w == w[::-1]), intruders=intr,
                       centres=[round(f['c'], 12) for f in t['fans']],
                       exits=[f.get('exit') for f in t['fans']],
                       ev_disagree=t['ev_disagree'],
                       is_sigma=(w == sigma(a1)), is_sigmaJ=(w == sigma_J(a1)),
                       retrace_ok=(t['headon'] is not None
                                   and t['closure'] == 2 * t['headon'] + 1))
        # fans strictly before the head-on
        if t['headon'] is not None:
            out[nm]['nfan_pre'] = sum(1 for f in t['fans'] if f['entry'] <= t['headon'])
    return out


FAMILIES = [(7, 3), (9, 4), (11, 5), (13, 6), (16, 7), (23, 10),
            (5, 2), (12, 5), (19, 8)]          # last three: a3=2 controls (out of scope)
A1S = [4, 5, 6, 7, 8, 10]


def main():
    rows, tally = [], {}
    for (P, r) in FAMILIES:
        for a1 in A1S:
            if math.gcd(P, r) != 1 or not (2 * r < P < 3 * r):
                continue
            try:
                rows.append(analyse(P, r, a1))
            except Exception as e:                       # noqa: BLE001
                rows.append(dict(P=P, r=r, a1=a1, error=repr(e)))
    inscope = [x for x in rows if 'error' not in x and x['a3'] >= 3 and x['a1'] >= 5]

    def score(name, pred):
        ok = sum(1 for x in inscope if pred(x))
        tally[name] = f"{ok}/{len(inscope)}"

    score('M0-A  2 fans pre-head-on', lambda x: x['CAP'].get('nfan_pre') == 2)
    score('M0-B  head-on = 4a1', lambda x: x['CAP']['headon'] == 4 * x['a1'])
    score('M0-C  closure = 8a1+1', lambda x: x['CAP']['closure'] == 8 * x['a1'] + 1)
    score('M0-D  #L1 = 3 & palindrome',
          lambda x: x['CAP']['nL1'] == 3 and x['CAP']['palindrome'])
    score('M0-E  MIR = reverse(CAP)  [MISS]',
          lambda x: x['MIR']['word'] == x['CAP']['word'][::-1])
    score('M0-E2 J-pair: equal width+wlen',
          lambda x: x['MIR']['wlen'] == x['CAP']['wlen']
          and abs((1.0 - x['H']) - (x['H'] - x['S'])) < 1e-12)
    score('M0-F  no intruder in (H,1)', lambda x: not x['CAP']['intruders'])
    score('M0-G  beat 6 fans/12a1-8/24a1-15/11',
          lambda x: (x['BEAT'].get('nfan_pre') == 6
                     and x['BEAT']['headon'] == 12 * x['a1'] - 8
                     and x['BEAT']['closure'] == 24 * x['a1'] - 15
                     and x['BEAT']['nL1'] == 11))
    # POST-HOC on 6 rows, then tested on all 30 (labelled as such, not a blind test):
    score('M0-H  CAP word == Sigma(a1)', lambda x: x['CAP']['is_sigma'])
    score('M0-I  MIR word == Sigma_J(a1)', lambda x: x['MIR']['is_sigmaJ'])
    score('M0-J  closure == 2*headon+1', lambda x: x['CAP']['retrace_ok'])
    score('M0-K  beat closure == 2*headon+1', lambda x: x['BEAT']['retrace_ok'])

    print(f"IN-SCOPE ROWS (a3>=3, a1>=5): {len(inscope)}   ALL ROWS: {len(rows)}")
    print("\n=== PRE-REGISTERED SCORES ===")
    for k, v in tally.items():
        print(f"  {k:34s} {v}")

    print("\n=== PER-ROW (CAP) ===")
    print(f"  {'P':>3} {'r':>3} {'a1':>3} {'a3':>3} | {'wlen':>5} {'8a1+1':>6} "
          f"{'ho':>5} {'4a1':>5} {'#L1':>4} {'fanpre':>6} {'pal':>4} {'L1 positions'}")
    for x in rows:
        if 'error' in x:
            print(f"  {x['P']:3d} {x['r']:3d} {x['a1']:3d}  ERROR {x['error']}")
            continue
        c = x['CAP']
        mk = '' if (x['a3'] >= 3 and x['a1'] >= 5) else '  (out of scope)'
        print(f"  {x['P']:3d} {x['r']:3d} {x['a1']:3d} {x['a3']:3d} | {c['wlen']:5d} "
              f"{8*x['a1']+1:6d} {str(c['headon']):>5} {4*x['a1']:5d} {c['nL1']:4d} "
              f"{str(c.get('nfan_pre')):>6} {str(c['palindrome']):>5} {c['L1pos']}{mk}")

    print("\n=== RAW CAP WORDS (for the Sigma(a1) question -- measured, not predicted) ===")
    for x in rows:
        if 'error' in x:
            continue
        print(f"  (P,r,a1)=({x['P']},{x['r']},{x['a1']}) a3={x['a3']} "
              f"len={x['CAP']['wlen']}\n      {x['CAP']['word']}")

    ev = [(x['P'], x['r'], x['a1'], x['CAP']['ev_disagree']) for x in rows
          if 'error' not in x and x['CAP']['ev_disagree']
          and not all(f for (_, _, f) in x['CAP']['ev_disagree'])]
    print(f"\n=== _exact_vertical FALSE NEGATIVES (float gap 0.0, is_zero() False): "
          f"{len(ev)}/{len([x for x in rows if 'error' not in x])} rows ===")
    for e in ev[:4]:
        print(f"    {e}")

    print("\n=== FAN CENTRES (CAP) vs landmark 2H ===")
    for x in rows:
        if 'error' in x:
            continue
        print(f"  ({x['P']},{x['r']},{x['a1']}) H={x['H']:.9f} 2H={2*x['H']:.9f} "
              f"centres={[round(v, 9) for v in x['CAP']['centres']]} "
              f"exits={x['CAP']['exits']}")

    with open('data/s242_m0fan.json', 'w') as fh:
        json.dump(dict(rows=rows, tally=tally), fh, indent=1)
    print("\nwrote data/s242_m0fan.json")


if __name__ == '__main__':
    sys.exit(main())
