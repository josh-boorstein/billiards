"""Section 5.5, Table 2: the census of irrational right triangles carrying certified
non-perpendicular periodic cylinders, by class, and the overlaps stated under the table.

Each row of `census_record.json` is a (scanned y, word) point of one of the research searches
(the two leg boxes and the hypotenuse census; see the record's note).  For every row this
script certifies the word EXACTLY at its algebraic y (`lib/mirror.py`): a non-perpendicular,
primitive periodic cylinder whose direction lies in a mirror class -- the class is read off
the certificate, not off the record -- and checks a/pi irrational.  It then normalises to the
paper's convention (a triangle is {y, -y}, listed at y > 0 where OB is the shorter leg) and
counts triangles per class.

Statements checked:
  * every row certifies: a non-perpendicular cylinder in a mirror direction, at an irrational
    angle;
  * Table 2: perp OB 1 024, perp OA 428, perp H 182, par H 215 triangles (the range of degrees
    of cos 2a per class is printed, not checked: the leg classes include four rows of the second
    leg box whose algebraic point has degree 8, 16, 19 or 20, outside that box's degree bound,
    found because they lie within 1e-12 of a box root);
  * the overlaps: 8 triangles in both leg classes, 9 in both hypotenuse classes, 138
    hypotenuse-class triangles also in a leg class, and cos 2a = (1 + sqrt 41)/8 in all four.

Verdicts are written to `census_verdicts.json` as they are reached (resumable: a rerun skips
rows already decided).  By default the script re-certifies a stratified sample of rows and
reads the rest from the verdict file; `--full` re-certifies every row (about half an hour on
one core; `--shard i n` splits it, `--merge` only counts).
"""
import json
import math
import os
import random
import sys
import time

import sympy as sp

from common import HERE, FULL, check, done, load, normalise
import mirror as m

Y = sp.Symbol('y')
VERD = os.path.join(HERE, 'census_verdicts.json')


def key(r):
    return f"{r['src']}|{r['y']}|{r['period']}"


def certify_row(r):
    K = m.t_field(r['minpoly'], float(r['y']))
    c = m.certify(K, r['word'])
    return dict(ok=bool(c['ok']), cls=c.get('cls'), e=c.get('e'), perp=c.get('perpendicular'),
                primitive=c.get('primitive'), irrational=m.angle_is_irrational(r['minpoly']),
                deg=sp.Poly(sp.sympify(r['minpoly']), Y).degree(), width=c.get('width'),
                why=c.get('why', ''))


def load_verdicts(path=VERD):
    return json.load(open(path)) if os.path.exists(path) else {}


def run(rows, path, label):
    V = load_verdicts(path)
    t0 = time.time()
    for i, r in enumerate(rows):
        k = key(r)
        if k in V:
            continue
        V[k] = certify_row(r)
        if i % 25 == 0:
            json.dump(V, open(path, 'w'))
            print(f'  {label}: {i + 1}/{len(rows)}  {time.time() - t0:.0f}s', flush=True)
    json.dump(V, open(path, 'w'))
    return V


def main(argv):
    R = load('census_record.json')['rows']
    if '--shard' in argv:
        i, n = int(argv[argv.index('--shard') + 1]), int(argv[argv.index('--shard') + 2])
        run(R[i::n], VERD.replace('.json', f'_shard{i}.json'), f'shard {i}/{n}')
        return
    if '--merge' in argv or FULL:
        V = {}
        for f in sorted(os.listdir(HERE)):
            if f.startswith('census_verdicts_shard'):
                V.update(json.load(open(os.path.join(HERE, f))))
        if os.path.exists(VERD):
            V.update(load_verdicts())
        json.dump(V, open(VERD, 'w'))
        if FULL:
            V = run(R, VERD, 'full')
        fresh = None
    else:
        V = load_verdicts()
        rng = random.Random(734)
        by_src = {}
        for r in R:
            by_src.setdefault(r['src'], []).append(r)
        sample = [r for s in sorted(by_src) for r in rng.sample(by_src[s], min(15, len(by_src[s])))]
        sample += [r for r in R if 'collision' in r]
        fresh = {}
        for r in sample:
            fresh[key(r)] = certify_row(r)
        agree = all(V.get(k) is None or (V[k]['ok'], V[k]['cls']) == (v['ok'], v['cls'])
                    for k, v in fresh.items())
        check(f'{len(fresh)} rows re-certified from scratch (15 per source, and the 4 root '
              f'collisions) agree with the verdict record', agree and all(v['ok'] for v in fresh.values()))
        V.update(fresh)

    missing = [key(r) for r in R if key(r) not in V]
    check(f'every one of the {len(R)} census rows has a verdict', not missing)
    bad = [k for k, v in V.items() if not (v['ok'] and not v['perp'] and v['primitive'] and v['irrational'])]
    check(f'every row is a primitive non-perpendicular periodic cylinder in a mirror direction '
          f'at an irrational angle ({len(bad)} exceptions)', not bad)

    tri = {}                                        # (|y| rounded, class) -> degree
    for r in R:
        v = V.get(key(r))
        if not v or not v['ok']:
            continue
        ay, cls = normalise(float(r['y']), v['cls'])
        tri[(round(ay, 9), cls)] = v['deg']
    count = lambda c: sorted(d for (ay, cc), d in tri.items() if cc == c)
    want = {'perp OB': 1024, 'perp OA': 428, 'perp H': 182, 'par H': 215}
    for c, n in want.items():
        ds = count(c)
        print(f'    {c:8s} {len(ds):5d} triangles, degree of cos 2a {ds[0] if ds else "-"}..{ds[-1] if ds else "-"}')
        check(f'Table 2: {c} carries {n} triangles', len(ds) == n)
    by = {}
    for (ay, c) in tri:
        by.setdefault(ay, set()).add(c)
    legs = sum(1 for s in by.values() if {'perp OB', 'perp OA'} <= s)
    hyps = sum(1 for s in by.values() if {'perp H', 'par H'} <= s)
    mixed = sum(1 for s in by.values() if s & {'perp H', 'par H'} and s & {'perp OB', 'perp OA'})
    four = [ay for ay, s in by.items() if len(s) == 4]
    print(f'    both leg classes {legs}; both hypotenuse classes {hyps}; hypotenuse and leg {mixed}; '
          f'all four at |y| = {[round(x, 6) for x in four]}')
    check('8 triangles carry cylinders in both leg classes', legs == 8)
    check('9 triangles carry cylinders in both hypotenuse classes', hyps == 9)
    check('138 hypotenuse-class triangles also carry one in a leg class', mixed == 138)
    check('the triangle cos 2a = (1 + sqrt 41)/8 carries all four, and it is the only one',
          len(four) == 1 and abs(four[0] - (1 + math.sqrt(41)) / 8) < 1e-9)
    done()


if __name__ == '__main__':
    main(sys.argv[1:])
