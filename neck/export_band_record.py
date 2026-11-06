"""Write band_record.json from the research repository's flow-decomposition sweeps.

  python export_band_record.py <sweep.json> <census_H.jsonl>
The sweeps ran `lib/flow.py`'s decomposition (`n_flow`) over the band: every coprime centre with
2 <= Q <= 60 in the direction perpendicular to L1 (the sweep's part A, 1101 centres), and in the
direction perpendicular to H at every even Q <= 60 (373 centres) and at odd Q <= 31 (212, the
control).  Kept per centre: the branch count, the numbers of cylinders and of components, and
whether the cylinders tile every edge of the beam's side exactly.  Every script that reads the
record also recomputes a slice of it from scratch and compares; `--full` recomputes all of it."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RECORD = os.path.join(HERE, 'band_record.json')


def main():
    sweep = json.load(open(sys.argv[1]))
    L1 = {f'{v["P"]}/{v["Q"]}': {'n': v['n'], 'ncyl': v['ncyl'], 'ncomp': v['ncomp'],
                                 'tiled': bool(v['tiled'])}
          for v in sweep.values() if v['part'] == 'A'}
    H = {}
    for line in open(sys.argv[2]):
        v = json.loads(line)
        H[f'{v["P"]}/{v["Q"]}'] = {'n': v['n_H'], 'ncyl': v['ncyl'], 'ncomp': v['ncomp'],
                                   'tiled': bool(v['tiled'])}
    json.dump({'L1': L1, 'H': H}, open(RECORD, 'w'), indent=0, sort_keys=True)
    print(f'{len(L1)} centres perpendicular to L1, {len(H)} perpendicular to H -> {RECORD}')


if __name__ == '__main__':
    main()
