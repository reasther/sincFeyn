"""Count the lobes of an integrand from a SINCQED_DUMP file: the connected pieces (3^D lattice
neighbours) of {|f| >= t M}, {f >= t M} and {f <= -t M}, M = max |f|, for t = 0.5 ... 0.01.

    SINCQED_DUMP=d bin/sincqed run 1.5 65 1e-3 0 g4_133      # top = 0: the whole grow scan at step h
    python3 tools/lobes.py d_g4_133.pts

Prints per graph the ratio of negative to positive mass and the piece counts |f|/f+/f- per level;
"(cov x)" marks a level where only a fraction x of the points has all 2D axis neighbours in the
dump (lower eps to sample it fully).  A split that disappears at half the step is a narrow ridge
falling between lattice points, not a second peak.  Keys use base 67, so D <= 10.  Needs numpy, scipy.
"""
import numpy as np, sys, itertools
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
def load(fn):
    f=open(fn,'rb'); hdr=f.readline().split(); D=int(hdr[2])
    dt=np.dtype([('l','u1'),('k','<i2',(D,)),('v','<f8')])
    a=np.fromfile(f,dtype=dt); return hdr[1].decode(),D,a
LEV=(0.5,0.3,0.1,0.03,0.01)
def comps(k,v,D):
    B=67; key=np.zeros(len(k),np.int64)
    for i in range(D): key=key*B+(k[:,i]+33)
    o=np.argsort(key); key=key[o]; v=v[o]
    n=len(v); lab=np.arange(n); I=[];J=[]; ne=0
    def merge(lab,I,J):
        I=np.concatenate(I+[np.arange(n)]); J=np.concatenate(J+[lab])
        g=coo_matrix((np.ones(len(I),np.int8),(I,J)),shape=(n,n))
        return connected_components(g,directed=False)[1]
    for d in itertools.product((-1,0,1),repeat=D):
        if d<=tuple([0]*D): continue      # half the offsets
        dk=0
        for x in d: dk=dk*B+x
        nk=key+dk; j=np.searchsorted(key,nk); j[j>=len(key)]=len(key)-1
        f=np.where(key[j]==nk)[0]; I.append(f); J.append(j[f]); ne+=len(f)
        if ne>2e7: lab=merge(lab,I,J); I=[];J=[]; ne=0
    lab=merge(lab,I,J); nc=lab.max()+1
    mx=sorted((np.abs(v[lab==c]).max() for c in range(nc)),reverse=True)
    return nc,mx
def analyse(fn):
    name,D,a=load(fn)
    a=a[(a['l']&0x7f)==(a['l']&0x7f).min()]
    k=a['k'].astype(np.int64); v=a['v']
    key=np.zeros(len(k),np.int64)
    for i in range(D): key=key*67+(k[:,i]+33)
    _,u=np.unique(key,return_index=True); k=k[u]; v=v[u]; key=key[u]
    M=np.abs(v).max(); r={}
    o=np.argsort(key); key=key[o]; k=k[o]; v=v[o]
    cov={}
    for t in LEV:
        sel=np.where(np.abs(v)>=t*M)[0]; ok=np.ones(len(sel),bool)
        for i in range(D):
            for sg in (-1,1):
                nk=key[sel]+sg*67**(D-1-i); j=np.searchsorted(key,nk); j[j>=len(key)]=len(key)-1; ok&=key[j]==nk
        cov[t]=ok.mean()
    for t in LEV:
        row=[]
        for s in (np.abs(v)>=t*M, v>=t*M, v<=-t*M):
            if s.sum()==0: row.append((0,[])); continue
            row.append(comps(k[s],v[s],D))
        r[t]=row
    pos=v[v>0].sum(); neg=-v[v<0].sum()
    return name,D,r,neg/pos,cov
for fn in sys.argv[1:]:
    name,D,r,q,cov=analyse(fn)
    s=" ".join("%g:%d/%d/%d%s"%(t,r[t][0][0],r[t][1][0],r[t][2][0],"" if cov[t]>0.999 else "(cov %.2f)"%cov[t]) for t in LEV)
    sec=[]
    for t in LEV:
        if r[t][0][0]>1: sec.append("2nd |f| piece at t=%g peaks %.2f"%(t,r[t][0][1][1]/r[t][0][1][0])); break
    print("%-7s D=%d neg/pos %-7.3g  pieces |f|/f+/f- at t*max  %s   %s"%(name,D,q,s," ".join(sec)),flush=True)
