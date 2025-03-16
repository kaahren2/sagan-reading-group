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

n = 4
M, Pl, Pi = transition_matrix(n)
print("Here is the transition matrix for the random walk")
print(M.str())

def prange(n):
    yield from (i+1 for i in range(n))

def enumerate_tableaux(lam):
    n = sum(lam)
    for p in Permutations(n):
        s = 0
        yield tuple([tuple(p[s:(s:=s+ell)]) for ell in lam])

def tableau_to_tabloid(t):
    return tuple(frozenset(r) for r in t)

def list_tabloidx(lam):
    return sorted(set(
        tableau_to_tabloid(t) for t in enumerate_tableaux(lam)
    ))

def remove_a_box(lam):
    """
    iterate over all possible ways of removing
    one box from lambda
    """
    for i,r in enumerate(lam):
        if (i+1 == len(lam) or lam[i+1] < lam[i]):
            if lam[i] == 1:
                yield i, lam[:i]
            else:
                yield i, lam[:i] + (lam[i]-1,) + lam[i+1:]

def enumerate_SYT(lam, cache={}):
    if lam in cache:
        yield from cache[lam]
    elif len(lam) == 1: # if lambda is a single row...
        ans = (tuple(prange(lam[0])),)
        cache[lam] = [ans]
        yield ans
    elif all(li == 1 for li in lam): # lambda is a single column
        ans = tuple((i,) for i in prange(len(lam)))
        cache[lam] = [ans]
        yield ans
    else:
        n = sum(lam)
        cache[lam] = []
        for i,mu in remove_a_box(lam):
            for t in enumerate_SYT(mu):
                if len(t) <= i:
                    ans = t + ((n,),)
                else:
                    ans = t[:i] + (t[i] + (n,),) + t[i+1:]
                cache[lam].append(ans)
                yield ans

def act_on_thing(p, t):
    """
    t should be a tuple of iterables of inputs to p
    """
    tt  = type(t[0])
    ans = tuple(
        tt(p(i) for i in r) for r in t
    )
    return ans

def shape(t):
    return tuple(len(r) for r in t)

def enumerate_columns(t):
    for j in range(len(t[0])):
        yield tuple(t[i][j] for i in range(len(t)) if j < len(t[i]))

def involutions_to_generate_Sn(s):
    return [[(s[0], s[i])] for i in range(1, len(s))]

def column_stabilizer(t):
    col_gens = [g for C in enumerate_columns(t) for g in involutions_to_generate_Sn(C)]
    return PermutationGroup(col_gens)

def polytabloid(t):
    CS  = column_stabilizer(t)
    tt  = tableau_to_tabloid(t)
    ans = {act_on_thing(p, tt): p.sign() for p in CS}
    return ans

def act_on_polytabloid(p, pt):
    return {act_on_thing(p, t):v for t,v in pt.items()}

def tabloid_basis(n, bases={}):
    if n not in bases:
        bases[n] = {'basis': []}
        for sh in Partitions(n):
            bases[n]['basis'].extend(list_tabloidx(tuple(sh)))
        bases[n]['ttoi'] = {b:i for i,b in enumerate(bases[n]['basis'])}
    return bases[n]

def polytabloid_to_vector(pt):
    if not isinstance(pt, dict):
        pt = {pt:1}
    sh   = shape(next(iter(pt)))
    ttoi = tabloid_basis(sum(sh))['ttoi']
    ans  = vector(QQ, len(ttoi))
    for t,v in pt.items():
        ans[ttoi[t]] = v
    return ans

def tabloid_action_matrix(p):
    if not isinstance(p, Permutation):
        p = Permutation(p)
    B    = tabloid_basis(len(p))
    cols = [polytabloid_to_vector({act_on_thing(p, b):1}) for b in B['basis']]
    M    = matrix(QQ, cols).transpose()
    return M

def specht_basis_as_polytabloids(n):
    ans = []
    for sh in Partitions(n):
        syt  = enumerate_SYT(tuple(sh))
        sytB = [polytabloid(t) for t in syt]
        ans.append(sytB)
    return ans

def specht_basis_as_vectors(n):
    SB = specht_basis_as_polytabloids(n)
    ans = [
        [polytabloid_to_vector(pt) for pt in sb] for sb in SB
    ]
    return ans

def specht_basis_column_matrix(n, cache={}):
    if n in cache:
        return cache[n]
    cols = [v for subrep in specht_basis_as_vectors(n) for v in subrep]
    B    = matrix(QQ, cols).transpose()
    cache[n] = B
    return B

def specht_bases_as_column_matrices(n, cache={}):
    if n in cache:
        return cache[n]
    subreps = specht_basis_as_vectors(n)
    Bs = [matrix(QQ, [v for v in subrep]).transpose() for subrep in subreps]
    cache[n] = Bs
    return Bs

def action_on_specht_basis(p, n):
    M    = tabloid_action_matrix(p)
    SB   = specht_basis_column_matrix(n)
    ans  = SB.solve_right(M*SB)
    return ans

def actions_on_specht_bases(p, n):
    M = tabloid_action_matrix(p)
    SBs = specht_bases_as_column_matrices(n)
    acts = [sb.solve_right(M*sb) for sb in SBs]
    ans = {}
    for lam, sb, act in zip(Partitions(n), SBs, acts):
        ans[lam] = {
            'basis': sb,
            'matrix': act
        }
    return ans

Ps = Partitions(n)
SGA = SymmetricGroupAlgebra(CC, n)
for lam in Ps:
    print('---------------')
    print(lam)
    print(lam.ferrers_diagram())
    sm = SGA.specht_module(lam)
    print(sm)
    print(sm.basis())
    # for g in SGA.group():
    #     print(g)

print(actions_on_specht_bases(((1,2),(3,4)), 4))

def shuffle_density(p):
    """
    Assign mass 1/n to the identity
    and mass 2/n*n to the transpositions
    """
    from collections import Counter
    n = len(p)
    ct = p.cycle_tuples()
    ctl = Counter(len(t) for t in ct)
    if len(ctl) == 0 or set(ctl) == set([1]):
        return 1/n
    elif ctl[2] == 1 and set(ctl).issubset([1,2]):
        return 2/(n*n)
    return 0

def FT(p_dens, X, n):
    summands = [
        (p_dens(p), X(p)) for p in Permutations(n)
    ]
    return (sum(a*b for a,b in summands), summands)

n = 4
for lam in Partitions(4):
    print(f"{lam} ************** ")
    just_X = lambda p:actions_on_specht_bases(p,n)[lam]['matrix']
    print(FT(shuffle_density, just_X, n))