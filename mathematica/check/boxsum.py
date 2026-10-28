import numpy as np, itertools, sys, time, math
from gen_pw import load
proj=lambda v,n: v**-n
schw=lambda v,n: np.exp(-v)/(1.0 if n<3 else float(math.factorial(n-1)))
def sincsum(F, D, h, K, chunk_axes=2):
    k=np.arange(-K,K+1)*h
    tot=0.0
    rest=D-chunk_axes
    grid=np.meshgrid(*([k]*rest), indexing='ij'); grid=[g.ravel() for g in grid]
    for head in itertools.product(k, repeat=chunk_axes):
        z=[np.full(grid[0].shape, c) for c in head]+grid
        s=[np.exp(t) for t in z]
        w=np.exp(sum(z))*h**D
        tot+=np.sum(F(s)*w)
    return tot
def forms(g):
    k=load(g); pt=k['pt']
    one=lambda s: np.ones_like(s[0])
    if g=='g2_01':  # x4 photon fixed; delta form keeps x1, x3 separate
        d5=lambda s: pt([s[0],s[1],s[2],s[3],one(s),s[4]],proj)/256
    else:
        d5=lambda s: pt([s[0],s[1],s[2],s[3],one(s),s[4]],proj)/256
    a6=lambda s: pt(s,schw)/256
    gl=lambda s: k['f'](s,proj)
    return a6,d5,gl
exact={'g2_00':-0.0640193/2,'g2_01':-0.5899758/2,'g2_02':-0.4676454,'g2_03':0.7774780}
if __name__=='__main__':
    g=sys.argv[1]; form=sys.argv[2]; h=float(sys.argv[3]); K=int(sys.argv[4])
    a6,d5,gl=forms(g)
    F,D={'all':(a6,6),'delta':(d5,5),'gl':(gl,3 if g=='g2_01' else 4)}[form]
    t=time.time(); v=sincsum(F,D,h,K)
    print(g,form,h,K,'%.10f'%v,'diff %.2e'%(v-exact[g]),'%.0fs'%(time.time()-t))
