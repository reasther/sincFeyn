import sympy as sp
from sympy.polys.matrices import DomainMatrix
def symbolic(A,B,m2):
    n=len(A); L=len(A[0])
    s=sp.symbols('s0:%d'%n, positive=True)
    M=sp.Matrix(L,L,lambda l,k: sum(s[i]*A[i][l]*A[i][k] for i in range(n)))
    # det and adjugate over ZZ[s]: Matrix.det()/adjugate() (Bareiss with dotprodsimp) stall at five loops
    R=sp.ZZ[s]; dM=DomainMatrix.from_list_sympy(L,L,M.tolist()).convert_to(R)
    adjD,detD=dM.adj_det()
    U=R.to_sympy(detD); adj=adjD.to_Matrix()
    J=[[sum(s[i]*A[i][l]*B[i][c] for i in range(n)) for c in range(2)] for l in range(L)]
    Js=[J[l][0]+J[l][1] for l in range(L)]
    F=sp.expand(sum(Js[l]*adj[l,k]*Js[k] for l in range(L) for k in range(L)) - U*sum(s[i]*((B[i][0]+B[i][1])**2-m2[i]) for i in range(n)))
    X=[[sum(adj[l,k]*J[k][c] for k in range(L)) for c in range(2)] for l in range(L)]   # U * M^-1 J
    coefU=[[sp.expand(U*B[i][c]-sum(A[i][l]*X[l][c] for l in range(L))) for c in range(2)] for i in range(n)]
    CU=[[sp.expand(sum(A[i][l]*adj[l,k]*A[j][k] for l in range(L) for k in range(L))) for j in range(n)] for i in range(n)]
    return s,U,F,coefU,CU
if __name__=='__main__':
    A=[[-1,0],[-1,-1],[-1,-1],[0,-1],[1,0],[0,1]]
    B=[[1,0],[1,0],[0,1],[0,1],[0,0],[0,0]]
    m2=[1,1,1,1,0,0]
    s,U,F,coefU,CU=symbolic(A,B,m2)
    print('U =',U); print('F =',sp.factor(F))
