"""[Deposit] Proposition 4.4: the rhombus (3,2,3,2)π/5 has free fraction exactly φ/5. Paper: "Perpendicular billiards in rational polygons of odd denominator". Research name: probes/s930_rhombus_exact.py, copied verbatim except imports (deposit names), paths (data/ and figs/ under odd_polygons/) and this paragraph. Run through run_all.py; the usage lines below are the research repo's. Bracketed IDs and section marks (§sNNN) refer to the project's internal ledgers.

s930 keeper (queue `s928-e`): the rhombus (3,2,3,2)/5 — its free area is EXACTLY φ/5 of S, and R = 1 − φ/5 exactly.

Not pre-registered (an exact evaluation of a measured constant, §s928 (4)(b): 1 − R = φ/5 to 25 digits).
Side c = sin(3π/10) = φ/2 (keeper convention of s928_kites.py). Two computations, both exact (sympy):
  (A) the free cylinder: the period-4 billiard orbit (sides 1, 3, 2, 0) in direction index 3 (54° to side 0); its band
      on side 0 is s ∈ (1/φ², 1) (lower edge: the orbit through V3, s = 1/φ² read off the scan and checked exactly by the audit);
      λ = φ² sin72° (= the mirror-couple loop length), h = (1 − 1/φ²)·c·sin54° = φ/4; free fraction 2λh/area(S) = φ/5.
      By hand: opposite sides are parallel, so σ1σ3 and σ2σ0 are translations by 2H·n₁, 2H·n₀ (H = c sin72° the
      rhombus height, n₀, n₁ the normals at 72°), λ = 2H|n₀ + n₁| = 2Hφ = φ² sin72°.
  (B) the invariant cylinders: every branch on the four perpendicular arcs, words read off a 60-digit float walk,
      values exact along the words: breakpoints at s = 1 − 1/(2φ) (sides 0, 2) and 1/(2φ) (sides 1, 3), λ ∈ {H, 2φH};
      R = (2/5)[(3 − φ)/2 + 2φ·1/(2φ)] = 1 − φ/5. (A) + (B) = 1 exactly: no other free set.
PRIOR ART: §s928 (4)(b) (the measured constant, 20 digits; audit 1-5 to 25)
AUDIT: audits/s930_audit_kite.md (W: band-edge / return-map identities printed, not asserted — checked exactly by the audit;
  a capped prong would be silently dropped — none here).; `rulings.py --grep` *phi/5* → none.
DATA: data/s928_kites.json (the float R). Output: printed (exact identities). ~1 min.
"""
import _setup  # noqa: F401  (deposit: this directory on the path, cwd = odd_polygons/)
import sympy as sp
phi=(1+sp.sqrt(5))/2
c=phi/2
def u(deg): return sp.Matrix([sp.cos(sp.pi*deg/180), sp.sin(sp.pi*deg/180)])
V=[sp.Matrix([0,0])]; V.append(V[0]+c*u(0)); V.append(V[1]+c*u(108)); V.append(V[0]+c*u(108))
sides=[(V[0],V[1]),(V[1],V[2]),(V[2],V[3]),(V[3],V[0])]
s=sp.symbols('s')
P=V[0]+s*(V[1]-V[0]); d=u(54)   # dir index 3 = 3*pi/10
lam=0; word=[1,3,2,0]
for j in word:
    A,B=sides[j]; e=B-A
    t,r=sp.symbols('t r')
    sol=sp.solve(list(P+t*d-(A+r*e)),[t,r],dict=True)[0]
    tt=sp.simplify(sol[t]); lam+=tt; P=sp.simplify(P+tt*d)
    n=sp.Matrix([-e[1],e[0]])/sp.sqrt(e.dot(e)); d=sp.simplify(d-2*d.dot(n)*n)
    print('hit side',j,'r =',sp.nsimplify(sp.simplify(sol[r])), '  r(s=1/phi^2)=',sp.simplify(sol[r].subs(s,1/phi**2)))
lam=sp.simplify(lam); print('lambda =',lam,'=',sp.N(lam,25))
print("A:",'return point s-coord:',sp.simplify(P[0]/c), ' dir back to 54deg:',sp.simplify(d-u(54)))
h=(1-1/phi**2)*c*sp.sin(sp.pi*54/180)
area_S=10*c**2*sp.sin(sp.pi*72/180)
frac=sp.simplify(2*lam*h/area_S)
print("A:",'h =',sp.simplify(h),' free fraction - phi/5: minimal_polynomial ->',sp.minimal_polynomial(frac-phi/5,sp.Symbol('x')),' value',sp.N(frac,25))
print("A:",'lambda - phi^2 sin72: minimal_polynomial ->',sp.minimal_polynomial(lam-phi**2*sp.sin(2*sp.pi/5),sp.Symbol('x')))

import mpmath as mp
from poly_slots import Poly
mp.mp.dps=60; N=5
P=Poly([3,2,3,2],N,lengths=[mp.sin(3*mp.pi/10)]*4)
c=phi/2
def U(d): return sp.Matrix([sp.cos(sp.pi*d/(2*N)), sp.sin(sp.pi*d/(2*N))])
V=[sp.Matrix([0,0])]
for j in range(3): V.append(V[-1]+c*U(P.idx[j]))
k=4; m=P.m
def fwalk(X,d,cap=500):
    w=[]
    for _ in range(cap):
        h=P.first_hit(X,d,None)
        if h[0]=='vertex': return ('vertex',w+[('v',h[1])])
        _,j,X,s=h; w.append(j)
        if (d-m[j]-N)%(2*N)==0: return ('perp',w)
        d=(2*m[j]-d)%(4*N)
    return ('cap',w)
def ewalk(X,d,word):
    """exact: follow the side word from exact point X; return (end point, total length, side param of last hit)"""
    lam=0
    for j in word:
        if isinstance(j,tuple): # ends at vertex
            Vt=V[j[1]]; lam+=sp.sqrt(((Vt-X).T*(Vt-X))[0]); return Vt,sp.simplify(lam),None
        A,B=V[j],V[(j+1)%k]; e=B-A; ud=U(d)
        t,r=sp.symbols('t r'); sol=sp.solve(list(X+t*ud-(A+r*e)),[t,r],dict=True)[0]
        lam+=sol[t]; X=A+sol[r]*e; last=(j,sp.simplify(sol[r]))
        d=(2*m[j]-d)%(4*N)
    return X,sp.simplify(lam),last
area=sp.simplify(c**2*sp.sin(2*sp.pi/5))
# breakpoints on side j from prongs that land perpendicularly
cls=P.classes()[0]
bps={j:[sp.Integer(0),sp.Integer(1)] for j in range(k)}
for (i,d) in P.slots(cls):
    st,w=fwalk(P.V[i],d)
    if st=='perp':
        _,_,last=ewalk(V[i],d,w); bps[last[0]].append(last[1])
S=0
for j in range(k):
    pts=sorted(set(bps[j]),key=lambda x:float(x))
    d0=(P.idx[j]+N)%(4*N); Lj=c
    for a,b in zip(pts,pts[1:]):
        if abs(float(b-a))<1e-30: continue
        mid=(a+b)/2; Xf=(P.V[j][0]+mp.mpf(float(mid))*(P.V[(j+1)%k][0]-P.V[j][0]), P.V[j][1]+mp.mpf(float(mid))*(P.V[(j+1)%k][1]-P.V[j][1]))
        st,w=fwalk(Xf,d0); assert st=='perp',(j,st)
        Xe=V[j]+mid*(V[(j+1)%k]-V[j]); _,lam,_=ewalk(Xe,d0,w)
        S+=(b-a)*Lj*lam
        print(f'side {j} branch [{float(a):.6f},{float(b):.6f}] word {w} lam={float(lam):.10f}')
R=S/(2*N*area)
diff=sp.nsimplify(sp.N(R-(1-phi/5),80),[sp.sqrt(5)])
print('R =',sp.N(R,40),'  1-phi/5 =',sp.N(1-phi/5,40))
print('R - (1 - phi/5): simplify ->',sp.simplify(R-(1-phi/5)), '| minimal_polynomial ->', sp.minimal_polynomial(sp.simplify(R-(1-phi/5)),sp.Symbol('x')))
