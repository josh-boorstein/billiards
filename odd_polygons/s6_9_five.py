"""[Deposit] Remark 6.9: (B) and (C) for the five acyclic trapezoids without a π/N corner (3102 families). Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s932_five_acyclic.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s932 audit check, promoted to a keeper in s958 (queue `s955-f`): Theorem 6.4's finite hypotheses (B) and (C) on the five acyclic
trapezoids WITHOUT a π/N corner (paper 2, Remark 6.9 `cor:trap` and the paragraph after it: "(A) and (C) are the output of a finite
exact computation … (3102 of them for (3,4,9,10)/13)").

`probes/s932_transfer_census.py` gives (A) (acyclic transfer graph) for all 78 vectors and enumerates wall families only for chains
≤ 6; (3,4,9,10)/13 has chain 7 and was "checked separately in-session" — this is that check, and it covers all five.
(B): every stage-0 cap path ends at an exit or perpendicularly (`stage0` kinds ⊆ {exit, perp}). (C): every wall family enumerated by
`s932_multistage_walls.families` has W ≠ 0 (`W0 = 0`); `families` RAISES if the transit depth reaches maxT, so no family is capped.
s932 audit record (`data/scratch/s932/audit/check_BC.json`): families 846, 3102, 3078, 230, 1760 for (2,5,8,11)/13, (3,4,9,10)/13,
(3,5,10,12)/15, (2,3,6,7)/9, (3,4,5,6)/9; W0 = 0 on all; max transits 6, 8, 6, 4, 6.
Copied verbatim from `data/scratch/s932/audit/check_BC.py` except this docstring and the output (path, verdict line).
PRIOR ART: §s932 (1); `audits/s932_*` (the audit that wrote it). DATA: the audit record above.
Usage: python probes/s932_five_acyclic.py → data/s932_five_acyclic.json
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import sys, json, collections
import s6_10_walls as TW, multistage_walls as M
from poly_exact import Fam
rows=[(13,[2,5,8,11]),(13,[3,4,9,10]),(15,[3,5,10,12]),(9,[2,3,6,7]),(9,[3,4,5,6])]
out=[]
for N,a in rows:
    TW.N=N; G=TW.Geo(a); E=Fam.closing(a,N,[(1,0),(0,1)])
    st=M.stage0(G,a,N,E.slots())
    F=M.families(G,st,maxT=30)
    bpk=collections.Counter()
    for (X,d),pcs in G._maps.items():
        for p in pcs:
            if 'bp' in p: bpk[(p['bp']['kind'],p['bp'].get('vertex'))]+=1
    r=dict(N=N,a=a,slots=len(st),stage0=collections.Counter(o['kind'] for o in st.values()),fams=len(F),
           W0=sum(f['W'].iszero() for f in F),maxtransits=max(f['transits'] for f in F),bp_kinds={str(k):v for k,v in bpk.items()})
    print(r,flush=True); out.append({k:(dict(v) if isinstance(v,collections.Counter) else v) for k,v in r.items()})
ok = all(r['W0'] == 0 and set(r['stage0']) <= {'exit', 'perp'} for r in out)
json.dump({'rows': out, 'all_ok': ok}, open('data/s932_five_acyclic.json', 'w'), indent=1, default=str)
print('OK ' if ok else 'FAIL', '(B) no stage-0 cap path ends at a vertex, (C) no wall family with W = 0, on all five', [r['fams'] for r in out])
