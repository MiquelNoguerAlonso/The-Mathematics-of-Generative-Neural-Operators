"""Compute the exact worst-case energy gain and construct an attaining vector."""
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def calculate(L):
    g=1.;eps=np.ones(1)
    for ell in L:
        matrix=np.array([[(1+ell**2)*g,ell*np.sqrt(g)],[ell*np.sqrt(g),1.]])
        values,vectors=np.linalg.eigh(matrix)
        v=abs(vectors[:,-1])
        eps=np.r_[eps*v[0],v[1]]
        g=float(values[-1])
    d=eps[0]
    for ell,e in zip(L,eps[1:]):d=np.hypot(d,e+ell*d)
    assert abs(np.dot(eps,eps)-1)<1e-12
    assert abs(d*d-g)<=2e-12*max(1,g)
    return g,eps,d

rows=[]
for p in [1.,.5,.25]:
    for m in [8,32,128]:
        L=.5/np.arange(1,m+1,dtype=float)**p
        g,eps,d=calculate(L)
        rows.append({"p":p,"last_band":m,"optimal_gain":float(np.sqrt(g)),
            "exponential_gain":float(np.sqrt(2)*np.exp(np.dot(L,L))),
            "attainment_relative_residual":float(abs(d*d-g)/g),
            "attaining_errors":eps.tolist()})
# Independent check on many nonnegative error vectors at a fixed sensitivity schedule.
rng=np.random.default_rng(302841)
L=.5/np.arange(1,33,dtype=float)
g,_,_=calculate(L)
eps=rng.exponential(size=(10000,33));eps/=np.linalg.norm(eps,axis=1)[:,None]
d=eps[:,0]
for j,ell in enumerate(L):d=np.hypot(d,eps[:,j+1]+ell*d)
assert np.max(d*d)<=g+1e-12
report={"rows":rows,"random_checks":10000,"random_max_squared_gain":float(np.max(d*d)),
        "reference_squared_gain":g,"seed":302841}
(ROOT/'results/optimal_gain.json').write_text(json.dumps(report,indent=2)+'\n')
lines=[r'\begin{tabular}{rrrr}',r'\toprule',r'$p$ & Last band $m$ & Optimal gain & Exponential gain \\',r'\midrule']
for r in rows:
    lines.append(f"{r['p']:.2f} & {r['last_band']} & {r['optimal_gain']:.4f} & {r['exponential_gain']:.4f} "+r'\\')
lines.extend([r'\bottomrule',r'\end{tabular}'])
(ROOT/'tables/optimal_gain.tex').write_text('\n'.join(lines)+'\n')
print('PASS: exact gain attained at all 9 schedules/levels; 10000 independent random-vector inequalities.')
