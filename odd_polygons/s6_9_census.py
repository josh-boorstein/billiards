"""[Deposit] Remark 6.9: the 78 trapezoid vectors, 32 acyclic, (B) and (C) for chains ≤ 6. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s932_transfer_census.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s932 keeper (queue `s931-c`, `s931-b`): the REACH of [POLY-TRAP-ACYC] — exact cap transfer graph of every primitive trapezoid
(a0, a1, N − a1, N − a0)/N, odd N ≤ NMAX, both base angles acute; for the acyclic ones the theorem's two finite hypotheses
(no stage-0 cap path ends at a vertex; every wall family has W ≠ 0 — families enumerated by `s932_multistage_walls.families`
when chain ≤ 6; larger chains checked separately in-session: (3,4,9,10)/13, 3102 families, W ≠ 0).
Result (s932, N ≤ 15): 78 angle vectors up to mirror, 32 acyclic (27 with a π/N corner, 5 without: (2,3,6,7)/9, (3,4,5,6)/9,
(2,5,8,11)/13, (3,4,9,10)/13, (3,5,10,12)/15); every acyclic one passes both hypotheses.
PRIOR ART: `probes/s931_transfer_graph.py` (nine families; the graph definition is copied); §s930 (3) ⊕ audit S28 (the
discrete-change-set principle).  DATA: none before; writes data/s932_transfer_census.json.
Usage: python probes/s932_transfer_census.py NMAX
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import sys, json, time, math
import s6_10_walls as TW
import multistage_walls as M
from poly_exact import Fam
out=[]
for N in range(5, int(sys.argv[1])+1, 2):
    for a0 in range(1,(N+1)//2):
        for a1 in range(1,(N+1)//2):
            a=[a0,a1,N-a1,N-a0]
            if math.gcd(math.gcd(a0,a1),N)!=1: continue
            t=time.time(); TW.N=N; G=TW.Geo(a); nodes={}
            for d in range(1,4*N,2):
                cs=G.sg(G.cos(d))
                if cs==0: continue
                X='R' if cs>0 else 'L'; outs=set()
                for p in G.cap_map(X,d):
                    if p['kind']=='exit':
                        d1=p['d']; Y='L' if X=='R' else 'R'; outs|={(Y,d1),(Y,(-d1)%(4*N))}
                nodes[(X,d)]=outs
            state,longest,cyc={},{},[False]
            def dfs(v):
                state[v]=1; best=0
                for w in nodes.get(v,()):
                    if state.get(w)==1: cyc[0]=True; continue
                    if state.get(w)!=2: dfs(w)
                    best=max(best,1+longest.get(w,0))
                state[v]=2; longest[v]=best
            for v in nodes:
                if v not in state: dfs(v)
            row=dict(N=N,a=a,cyclic=cyc[0],chain=None if cyc[0] else max(longest.values()))
            if not cyc[0]:
                E=Fam.closing(a,N,[(1,0),(0,1)]); st=M.stage0(G,a,N,E.slots())
                row['stage0_vertex']=sum(o['kind']=='vertex' for o in st.values())
                if row['chain']<=6:
                    F=M.families(G,st); row['families']=len(F); row['W0']=sum(f['W'].iszero() for f in F)
            row['sec']=round(time.time()-t,1); out.append(row)
            print(row, flush=True)
json.dump(out,open('data/s932_transfer_census.json','w'),indent=1)
