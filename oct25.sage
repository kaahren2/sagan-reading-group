#Sagan sections 1.1 to 1.3, doodles

S4 = SymmetricGroup(4)
print(S4)
print()

S4_elements = [g for g in S4]
print(f"S4 has {len(S4_elements)} elements and they are {S4_elements}")
print()

#how many conjugacy classes of S4? Cycle type determines conjugacy class, so I expect the same as the number of integer partitions of 4. I will see ins momentarily, because Sage helpful prints the cycle type!

print(S4.conjugacy_classes())
print()

#That's a little overwhelming. Maybe just one representative from each class:

print(S4.conjugacy_classes_representatives())
print()

#I can compute which cycle class each element lives in

for g in S4_elements:
    print(g.conjugacy_class())

print()
#ok, Sagan says (order of centralizer of element)*(order of elements conjugacy class)=order of group. Let's verify. 

for g in S4_elements:
    assert(S4.centralizer(g).order() * g.conjugacy_class().cardinality() == S4.order())

#Let's write down the defining representation
S4_matrices = [f"{g.matrix()}     {g}" for g in S4_elements]
print("\n\n".join(S4_matrices))

