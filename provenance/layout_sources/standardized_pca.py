from pathlib import Path
from fractions import Fraction as F
import json, sys
import numpy as np

B = Path(sys.argv[1] if len(sys.argv)>1 else 'outputs/standardized_pca_python')
B.mkdir(parents=True, exist_ok=True)
G = np.array([[1,1,1,0,0],[0,1,2,1,2],[2,1,1,0,1],
              [0,0,1,2,2],[2,1,1,0,0],[0,0,1,1,1],[2,2,1,1,0]])
M, N = G.shape
p = G.mean(axis=1) / 2
sd = np.sqrt(2*p*(1-p))
X = (G-2*p[:,None]) / sd[:,None]
psi = X.T @ X / M
# The matrix entries also admit a completely independent rational calculation.
pf = [F(int(row.sum()),2*N) for row in G]
exact = [[sum((F(int(G[i,j]))-2*pf[i])*(F(int(G[i,k]))-2*pf[i]) /
              (2*pf[i]*(1-pf[i])) for i in range(M))/M for k in range(N)] for j in range(N)]
assert np.max(abs(psi-np.array(exact,dtype=float))) < 1e-14
d, V = np.linalg.eigh(psi)
d, V = d[::-1], V[:,::-1]
for k in range(N):
    if V[np.argmax(abs(V[:,k])),k] < 0: V[:,k] *= -1
assert np.max(abs(X.sum(axis=1))) < 1e-14
assert np.max(abs(psi.sum(axis=1))) < 1e-14
assert np.max(abs(psi @ V - V*d)) < 1e-13
assert np.max(abs(V.T@V-np.eye(N))) < 1e-13
assert d[-1] > -1e-13
score = V[:,:2] * np.sqrt(M*d[:2])
theta = .52
axes = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
small = dict(G=G.tolist(),X=X.tolist(),psi=psi.tolist(),V=V.tolist(),d=d.tolist(),
             score=score.tolist(),geometry=(score@axes.T).tolist(),axes=axes.tolist(),
             means=G.mean(axis=1).tolist(),frequencies=p.tolist(),scales=sd.tolist(),
             explained=(d/d.sum()).tolist(),exact=[[str(v) for v in row] for row in exact],
             mean_diagonal=float(np.trace(psi)/N))
(B/'small.json').write_text(json.dumps(small,indent=2))
