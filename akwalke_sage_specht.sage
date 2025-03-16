from sage.all import *
from sage.combinat.specht_module import polytabloid
SGA = SymmetricGroupAlgebra(QQ, 4)
dtov = lambda d,M:sum(v*M(x) for x,v in d.items())
TM  = SGA.tabloid_module([2,2])

def sym_act_tab(p, t):
    return tuple(frozenset(p(i) for i in s) for s in t)

def sym_act_tab_mod(p, tt):
    ttp = tt.parent()
    return sum(v*ttp(sym_act_tab(p, t)) for t,v in zip(tt.support(), tt.coefficients()))

t   = StandardTableau([[1,2],[3,4]])
et  = polytabloid(t)

print(et)
et = dtov(et,TM)

ST = StandardTableaux([2,2])
basis = [dtov(polytabloid(t), TM) for t in ST]

p = Permutation(((1,2),(3,4)))
print("This is t")
print(t)
print("This is e_t")
print(et)
print("This is p(et)")
pet = sym_act_tab_mod(p, et)
print(pet)

print("This is a basis for S^(2,2)")
print(basis)

p_basis = [sym_act_tab_mod(p, b) for b in basis]
m_basis = matrix(QQ, [b.to_vector() for b in basis])
m_p_basis = matrix(QQ, [b.to_vector() for b in p_basis])

print(f"This is a matrix for p = {p} acting on S^(2,2)")
p_matrix = m_basis.solve_left(m_p_basis)
print(p_matrix)

p2 = Permutation(((1,2,3,4)))
p2_basis = [sym_act_tab_mod(p2, b) for b in basis]
m_basis = matrix(QQ, [b.to_vector() for b in basis])
m_p2_basis = matrix(QQ, [b.to_vector() for b in p2_basis])

print(f"This is a matrix for p2 = {p2} acting on S^(2,2)")
p2_matrix = m_basis.solve_left(m_p2_basis)
print(p2_matrix)

