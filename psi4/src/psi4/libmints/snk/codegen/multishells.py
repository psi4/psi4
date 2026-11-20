#
#   For multishell-specific information
#

MULTISHELLS = {
    "S": (0,),
    "P": (1,),
    "D": (2,), 
    "F": (3,), 
    "G": (4,), 
    "H": (5,), 
    "I": (6,), 
    "J": (7,), 
    "K": (8,), 
    "L": (9,), 
    "M": (10,), 
    "N": (11,), 
    "O": (12,)
}
MULTISHELLS_LMAX = max(max(ms) for ms in MULTISHELLS.values())

def num_cartesian_integrals(vLA, vLB):
    nci = 0
    for LA in vLA:
        for LB in vLB:
            nci += (LA+1)*(LA+2)*(LB+1)*(LB+2)//4
    return nci

def get_hrr_blist(vLA, vLB, same_center):
    """Base integrals of the HRR: (lx,ly,lz,0,0,0,LA,LB) or (0,0,0,lx,ly,lz,LA,LB) -> index."""
    hrr_blist = {} 
    idx = 0
    if same_center:
        for LA in vLA:
            for LB in vLB:
                L = LA+LB
                for lx in range(L, -1, -1):
                    for ly in range(L-lx, -1, -1):
                        lz = L-lx-ly
                        lqns = (lx,ly,lz,0,0,0,LA,LB) if (LA>=LB) else (0,0,0,lx,ly,lz,LA,LB) # need LA,LB to keep keys unique
                        hrr_blist[lqns] = idx
                        idx += 1
    else:
        for LA in vLA:
            for LB in vLB:
                for L in range(max(LA,LB), LA+LB+1):
                    for lx in range(L, -1, -1):
                        for ly in range(L-lx, -1, -1):
                            lz = L-lx-ly
                            lqns = (lx,ly,lz,0,0,0,LA,LB) if (LA>=LB) else (0,0,0,lx,ly,lz,LA,LB) # need LA,LB to keep keys unique
                            hrr_blist[lqns] = idx
                            idx += 1
    return hrr_blist
