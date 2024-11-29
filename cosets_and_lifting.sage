

def partition_to_product_subgroup(p):
  n  = sum(p)
  #Sn = SymmetricGroup(domain=range(1,n+1))
  SGs = []
  offset = 1
  for i in p:
    SGs.append(SymmetricGroup(domain=[offset+j for j in range(i)]))
    offset += i
  return direct_product_permgroups(SGs)

def induce_trivial_rep(G, S):
  t    = transversal(G, S)
  tinv = [ti.inverse() for ti in t]
  ind  = len(t)
  def X(g):
    return matrix(QQ, [[(1 if tinv[i]*g*t[j] in S else 0) for j in range(ind)] for i in range(ind)])
  return X, t

def transversal(G, S):
  cosets = G.cosets(S, side='left')
  return [c[0] for c in cosets]

S5 = SymmetricGroup(domain=[1,2,3,4,5])
print(S5)
print("Partition [4,1]")
dp41 = partition_to_product_subgroup([4,1])
X41, t41 = induce_trivial_rep(S5, dp41)
print("product: ", dp41)
print("X((1,2,3))")
print(X41(S5((1,2,3))))

print("Partition [3,2]")
dp32 = partition_to_product_subgroup([3,2])
X32, t32 = induce_trivial_rep(S5, dp32)
print("product:", dp32)
print("X(e)")
print(X32(S5.identity()))
print("X((1,2,3))")
print(X32(S5((1,2,3))))

