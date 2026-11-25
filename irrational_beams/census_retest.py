"""[Deposit] Section 5: the floating-point finiteness retest of a fresh sweep slice (run_all.py --sweep). Paper: "Finite perpendicular-beam partitions at irrational angles". Research name: data/scratch/s977/retest.py (the research run's driver, run after each slice), copied verbatim except imports (deposit names), this paragraph and the _setup import."""
import json, sys
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = irrational_beams/)
P_ARG = sys.argv[1]; sys.argv = sys.argv[:1]
from multiprocessing import Pool
import coincidence_census as C0
p = P_ARG
d = json.load(open(p))
cand = [r for r in d['rows'] if r.get('screen_finite')]
with Pool(2) as P:
    res = P.map(C0.test, [(r['alpha'], ['conn']) for r in cand], chunksize=2)
for r, t in zip(cand, res):
    r['finite'] = bool(t['finite']); r['n6'] = t['n6']; r['n9'] = t.get('n9'); r['c9'] = t.get('c9')
    if t['finite']:
        r['branches'] = t['branches']
json.dump(d, open(p, 'w'), default=str)
fin = [r for r in cand if r['finite']]
print(len(cand), 'screened ->', len(fin), 'finite (C0.test)')
for r in fin:
    print('  x=%.10f n=%s A@%d O@%d' % (r['x'], r['n9'], r['A'][1], r['O'][1]))
