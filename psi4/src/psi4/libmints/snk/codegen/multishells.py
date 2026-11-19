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
