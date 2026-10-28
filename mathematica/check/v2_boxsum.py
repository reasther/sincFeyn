import numpy as np, math, sys
ex=119/36-math.pi**2/3
def kern(x0,x1,x2,x3,x4,x5,schw):
    a,b=x3,x5; S=a+b; c=a*b/S; om=-2*a*b/S**4; X2=x2+x4
    def base(y2):
        U=x0+x1+y2; V=(x0+x1)**2/U; cp0=(x1+y2)/U; cpp0=-x1/U; cp1=-x0/U; cpp1=(x0+y2)/U
        T=cp0+cp1*(1-cp0)+cpp0+cpp1*(1-2*cp0-cpp0); W=V+S
        return -64*T/U**2*(np.exp(-W) if schw else 1/W)
    return -4*om*(base(X2+c)-base(X2))/X2/256
def run(form,h,K):
    z=h*np.arange(-K,K+1); e=np.exp(z); tot=0.0
    if form=='DeltaGL':
        g=[(0.5-math.sqrt(15)/10,5/18),(0.5,4/9),(0.5+math.sqrt(15)/10,5/18)]
        S,A,B=np.meshgrid(e,e,e,indexing='ij'); acc=0
        for u,w in g: acc=acc+w*kern(S*u,S*(1-u),0.5,A,0.5,B,False)
        return h**3*np.sum(S*A*B*S*acc)
    if form=='Delta':
        X1,X2,X3,X5=np.meshgrid(e,e,e,e,indexing='ij')
        for x0 in e: tot+=h**5*np.sum(x0*X1*X2*X3*X5*kern(x0,X1,X2,X3,1.0,X5,False))
        return tot
    if form=='All':
        X2,X3,X4,X5=np.meshgrid(e,e,e,e,indexing='ij')
        for x0 in e:
            for x1 in e: tot+=h**6*np.sum(x0*x1*X2*X3*X4*X5*kern(x0,x1,X2,X3,X4,X5,True))
        return tot
form,h,K=sys.argv[1],float(sys.argv[2]),int(sys.argv[3])
v=run(form,h,K); print(form,h,K,'%.10f'%v,'%+.1e'%(v-ex))
