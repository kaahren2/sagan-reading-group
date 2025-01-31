from sage.all import *

def transition_matrix(n):
    # 1's on diagonal because we can transition
    # to ourselves
    M  = zero_matrix(QQ, factorial(n))
    P  = SymmetricGroup(list(range(n)))
    Pl = list(P)
    Pi = {p:i for i,p in enumerate(Pl)}
    for i in range(n):
        for j in range(n):
            if j == i:
                p = P.identity()
            else:
                p = P((i,j))
            for qi,q in enumerate(Pl):
                pqi = Pi[p*q]
                M[pqi,qi] += 1
    for ci,c in enumerate(M.columns()):
        M[:,ci] /= sum(c)
    return M, Pl, Pi

M, Pl, Pi = transition_matrix(4)
print(M.str())
print(Pi)