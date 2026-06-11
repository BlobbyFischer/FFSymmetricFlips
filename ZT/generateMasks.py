## Given three maps U, V, W which map the indices 1-8 to \pm 1, we get the group action given by
## a_ij -> U(i)/V(j) a_ij, b_jk -> V(j)/W(k) b_jk, c_ki -> W(k)/U(i) c_ki
## We hope to search for matrix multiplication tensor decompositions that are made from orbits of this group action

# U, V, W should be stored as dictionaries (e.g U = {0:1, 1:-1, ...})

def is_correct_map(U):
    for i in range(8):
        if U[i] not in [-1,1]:
            return False
    return True

def generate_masks(U,V,W):
    assert is_correct_map(U) and is_correct_map(V) and is_correct_map(W) # maybe this line is overkill, but oh well
    mask1 = 0
    mask2 = 0
    mask3 = 0
    for i in range(8):
        for j in range(8):
            if U[i] ^ V[j]:
                mask1 += (1 << (8*i+j))
            if V[i] ^ W[j]:
                mask2 += (1 << (8*i+j))
            if W[i] ^ U[j]:
                mask3 += (1 << (8*i+j))
    return hex(mask1), hex(mask2), hex(mask3)

if __name__ == "__main__":
    ## these are the exact two that make up the C2 x C2 symmetry in Smirnov's 336
    U1 = dict()
    V1 = dict()
    W1 = dict()
    U2 = dict()
    V2 = dict()
    W2 = dict()
    for i in range(8):
        U1[i] = 1
        V1[i] = 1
        W1[i] = 1
        U2[i] = 1
        V2[i] = 1
        W2[i] = 1

    ## Now they are initialised, we put in the sign flips

    U1[1] = -1
    U1[2] = -1

    V1[1] = -1

    W1[3] = -1
    W1[5] = -1

    ## Onto the second one

    U2[1] = -1

    V2[1] = -1
    V2[2] = -1

    W2[1] = -1
    W2[2] = -1

    print(generate_masks(U1,V1,W1))
    print(generate_masks(U2,V2,W2))

    ## How could we extend this to 8x8x8?