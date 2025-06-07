def from_lists(L):
    """
    Map entry j in list i to answer[(i,j)]
    """
    return {(i,j):k for i,r in enumerate(L) for j,k in enumerate(r)}

def is_standard(T):
    """
    T is standard if for all i,j we have
    T[(i,j)] is less than T[(i+1,j)] and T[(i,j+1)]
    """
    one = all(((i+1,j) not in T or T[(i,j)] < T[(i+1,j)]) for (i,j) in T)
    two = all(((i,j+1) not in T or T[(i,j)] < T[(i,j+1)]) for (i,j) in T)
    return one and two

def up_rows_left_cols(T):
    """
    Order (i,j) in T from right (higher) to left columns and
    bottom (higher) to top rows lexicographically
    """
    return sorted(T, key=lambda p:(-p[1],-p[0]))

def maximal_position_with_standard_subtableau(T):
    """
    Pick the largest (as sorted in up_rows_left_cols)
    position c such that T^{<c} is standard,
    or None if the entire T is standard
    """
    poses = up_rows_left_cols(T)
    for (i,j) in poses:
        if (((i+1,j) in T and T[(i,j)] > T[(i+1,j)]) or
            ((i,j+1) in T and T[(i,j)] > T[(i,j+1)])):
            return (i,j)
    return None

def slide_pos(T, pos):
    """
    While T^{<=pos} is not standard, swap pos and min of T[(i+1,j)], T[(i,j+1)]
    Assumes that T^{<=pos} is not standard initially
    return the path
    """
    TT = T.copy()
    i,j = pos
    path = [(i,j)]
    while True:
        print(TT)
        nbrs    = [(i+1,j), (i,j+1)]
        options = [(TT[(k,l)], (k,l)) for (k,l) in nbrs if (k,l) in TT and TT[(k,l)] < TT[(i,j)]]
        print((i,j), nbrs, options)
        if len(options) == 0:
            break
        _, (ni,nj) = min(options)
        path.append((ni,nj))
        TT[(i,j)], TT[(ni,nj)], (i,j) = TT[(ni,nj)], TT[(i,j)], (ni,nj)
        
    return path



            
