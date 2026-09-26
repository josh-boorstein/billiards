"""Write partition_record.json from a census run: complete centres only, no walk state.

  python export_record.py <census.json>
The census itself is `probes/s621_partition_census.py` in the research repository, which runs
`lib/partition.py` over every centre the paper's ranges need, on a resumable step ladder.  Every
script here that reads the record also recomputes a slice of it from scratch and compares."""
import json
import sys

from common import RECORD

KEEP = ('P', 'Q', 'n', 'orphans', 'orphan_head_on', 'J_ok', 'tiles', 'bounds', 'widths',
        'prefix_len', 'steps')


def main():
    src = json.load(open(sys.argv[1]))
    out = {k: {f: v[f] for f in KEEP if f in v} for k, v in src.items() if v['open'] == 0}
    json.dump(out, open(RECORD, 'w'), indent=0, sort_keys=True)
    print(f'{len(out)} complete centres of {len(src)} -> {RECORD}')


if __name__ == '__main__':
    main()
