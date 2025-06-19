import itertools

def tableau_from_lists(L):
    """
    Map entry j in list i to answer[(i,j)]
    """
    return {(i,j):k for i,r in enumerate(L) for j,k in enumerate(r)}

def all_tableau_of_shape(s):
    n = sum(s)
    P = itertools.permutations(range(1, n+1))
    for p in P:
        T = {}
        ptr = 0
        for i in range(len(shape)):
            for j in range(shape[i]):
                T[(i,j)] = p[ptr]
                ptr += 1
        yield T

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

def modified_forward_slide(T, pos):
    """
    While T^{<=pos} is not standard, swap pos and min of T[(i+1,j)], T[(i,j+1)]
    Assumes that T^{<=pos} is not standard initially
    return the final tableau and the path
    """
    TT = T.copy()
    i,j = pos
    path = [(i,j)]
    while True:
        # print(TT)
        nbrs    = [(i+1,j), (i,j+1)]
        options = [(TT[(k,l)], (k,l)) for (k,l) in nbrs if (k,l) in TT and TT[(k,l)] < TT[(i,j)]]
        # print((i,j), nbrs, options)
        if len(options) == 0:
            break
        _, (ni,nj) = min(options)
        path.append((ni,nj))
        TT[(i,j)], TT[(ni,nj)], (i,j) = TT[(ni,nj)], TT[(i,j)], (ni,nj)
    return TT, path

def funky_J(J, path):
    """
    Given the current J_{k-1} and next slide path,
    return the next J_k
    """
    JJ = J.copy()
    (i0,j0), *_, (i1,j1) = path
    HJ = [(h,j) for (h,j) in JJ if j == j0]
    for (h,j) in HJ:
        if i0 <= h and h < i1:
            JJ[(h,j)] = J[(h+1,j)]-1
        elif h == i1:
            JJ[(h,j)] = j1 - j0
    return JJ

def NPS(T):
    """
    Do the bijection tableau T -> pair (P, J) with 
    P standard, J hook, and all of the same shape
    """
    P_history = [T.copy()]
    J_history = [{(i,j):0 for (i,j) in T}]
    pos_history = []
    while (pos := maximal_position_with_standard_subtableau(P_history[-1])) is not None:
        pos_history.append(pos)
        next_P, path = modified_forward_slide(P_history[-1], pos)
        next_J = funky_J(J_history[-1], path)
        P_history.append(next_P)
        J_history.append(next_J)
    return P_history, J_history, pos_history

def SPN_Ck(i0, j0, J, nrows, ncols):
    return set((ii,j0 + Jiij0) for ii in range(i0, nrows) if (Jiij0 := J[(ii,j0)]) >= 0)

def SPN_modified_backward_slide(T, i0, j0, i, j):
    T     = T.copy()
    path  = [(i,j)]
    steps = []
    # print(f"Sliding pos {(i0,j0)}, {(i,j)} in:\n{display_tableau(T)}")
    Tget = lambda k,l: (T.get((k,l),-1) if l >= j0 else -1)
    while i != i0 or j != j0:
        cp   = ((i-1, j) if Tget(i-1,j) > Tget(i,j-1) else (i,j-1))
        step = ('N' if cp == (i-1, j) else 'W')
        T[(i,j)], T[cp], (i, j) = T[cp], T[(i,j)], cp
        path.append((i,j))
        steps.append(step)
    return T, path, ''.join(steps[::-1])

def SPN(P, J):
    T_history = [T:=P.copy()]
    J_history = [J:=J.copy()]
    nrows = max(i for (i,j) in T) + 1
    ncols = max(j for (i,j) in T) + 1
    ordered_cells = sorted(P, key=lambda p:(-p[1], -p[0]), reverse=True)
    for (i0,j0) in ordered_cells:
        # print(f"Current T:\n{display_tableau(T)}")
        # print(f"Current J:\n{display_tableau(J)}")
        Ck = SPN_Ck(i0, j0, J, nrows, ncols)
        # print(f"Doing position {(i0,j0)} with Ck {Ck}")
        ck_data = {(i,j):SPN_modified_backward_slide(T, i0, j0, i, j) for (i,j) in Ck}
        max_path_len = max(len(path) for _,_,path in ck_data.values())
        for (i,j),(_,_,steps) in ck_data.items():
            ck_data[(i,j)] += (f"{steps:O<{max_path_len}}",)
        # print("Options:", ck_data)
        (im,jm) = max(ck_data, key=lambda p:ck_data[p][-1])
        # print(f"Found max ck {(im,jm)} val {ck_data[(im,jm)]}")
        new_J = J.copy()
        for h in range(i0+1, im+1):
            new_J[(h,j0)] = J[(h-1,j0)] + 1
        new_J[(i0,j0)] = 0
        T_history.append(T := ck_data[(im,jm)][0])
        J_history.append(J := new_J)
    return T_history, J_history
        



def n_ascii(s):
    """
    Return the number of ascii characters in the string s
    """
    return len([x for x in s if x.isascii()])

def display_tableau(T, special_pos=None):
    """
    Return a string which displays a tableau
    """
    import unicodedata
    n_rows    = max(i for (i,j) in T) + 1
    n_cols    = max(j for (i,j) in T) + 1
    flair     = lambda pos,v: (f'{v}\u0332' if pos == special_pos else str(v))
    Tstrs     = {x:flair(x,v) for (x,v) in T.items()}
    col_width = max(n_ascii(s) for s in Tstrs.values())
    ans = [
        ' '.join(f"{Tstrs.get((i,j), ''):>{col_width}}" for j in range(n_cols))
        for i in range(n_rows)
    ]
    return ans

def display_tableau_list(L, special_poses=None):
    """
    Return a string which displays a list of tableaus horizontally next to each other
    """
    if special_poses is not None:
        blocks = [display_tableau(T, sp) for T,sp in zip(L, special_poses)]
    else:
        blocks = [display_tableau(T) for T in L]
    height = max(len(b) for b in blocks)
    for b in blocks:
        if len(b) < height:
            pad = n_ascii(b[0])+' '
            b.extend((height-len(b))*[pad])
    ans = []
    for i in range(height):
        line = '   '.join(b[i] for b in blocks)
        ans.append(line)
    ans = '\n'.join(ans) + '\n'
    return ans




T = tableau_from_lists(
    [
        [6, 2],
        [4, 3],
        [5, 1]
    ]
)

PH, JH, poses = NPS(T)

TH, JJH = SPN(PH[-1], JH[-1])

print(display_tableau_list(PH, poses+[None]))
print()
print(display_tableau_list(JH))
print()

print(display_tableau_list(TH[::-1]))
print()
print(display_tableau_list(JJH[::-1]))

shapes = [(2,2,2)]
for shape in shapes:
    for T in all_tableau_of_shape(shape):
        TH, JH, poses = NPS(T)
        TTH, JJH = SPN(TH[-1], JH[-1])
        assert TTH[-1] == TH[0]
        assert all(v==0 for v in JJH[-1].values())
    print(f"Checked all tableau of shape {shape}")
        

            
