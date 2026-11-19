#
#   For unrolling the basic recurrence relations
#

from multishells import num_cartesian_integrals
from sympy import Function, Symbol, symbols, S

#######################################################################################################
#                                                                                                     #
#   Vertical Recurrence Relation (VRR):                                                               #
#                                                                                                     #
#   [(i+1)jklmn]^[N] = (PA)_x*[ijklmn]^[N]   + (1/2p)*(i*[(i-1)jklmn]^[N]   + j*[i(j-1)klmn]^[N]  )   #
#                    - (PG)_x*[ijklmn]^[N+1] + (1/2p)*(i*[(i-1)jklmn]^[N+1] + j*[i(j-1)klmn]^[N+1])   #
#                                                                                                     #
#   Begin : [000000]^[N]                                                                              #
#   Target: [(lax)(lbx)(lay)(lby)(laz)(lbz)]^[0]                                                      #
#   Note: for same-center integrals A = B, (PA) = 0 and the first term in the VRR vanishes.           #
#                                                                                                     #
#######################################################################################################

def check_VRR_input(args):
    if len(args) != 13:
        raise ValueError(f"recurrences.VRR() expects 13 arguments but {len(args)} were given.")
    i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13 = args
    if i1.is_integer is False or \
       i2.is_integer is False or \
       i3.is_integer is False or \
       i4.is_integer is False or \
       i5.is_integer is False or \
       i6.is_integer is False or \
       i7.is_integer is False:
        raise TypeError("recurrences.VRR() expects integer Lqn inputs.")
    elif (i7<0) is S.true:
        raise ValueError("recurrences.VRR() expects N >= 0.")

class VRR(Function):
    @classmethod
    def eval(cls, *args):
        tpi = Symbol('tpi') # 1 / (2p) where p = a + b is the combined Gaussian exponent
        PGx, PGy, PGz = symbols('PGx[g] PGy[g] PGz[g]') # (Pi - Gi) where P is the combined Gaussian center, G is a gridpoint
        
        check_VRR_input(args)
        i,j,k,l,m,n,N,PAx,PAy,PAz,PBx,PBy,PBz = args
        if ((i<0) is S.true) or ((j<0) is S.true) or \
           ((k<0) is S.true) or ((l<0) is S.true) or \
           ((m<0) is S.true) or ((n<0) is S.true):
            return 0
        elif (i==0) and (j==0) and \
             (k==0) and (l==0) and \
             (m==0) and (n==0):
            return Symbol(f"scaled_boys_{N}[g]")
        elif (i>0) is S.true:
            return PAx*VRR(i-1,j,k,l,m,n,N,PAx,PAy,PAz,PBx,PBy,PBz) - PGx*VRR(i-1,j,k,l,m,n,N+1,PAx,PAy,PAz,PBx,PBy,PBz) + tpi*(
                (i-1)*(VRR(i-2,j  ,k,l,m,n,N,PAx,PAy,PAz,PBx,PBy,PBz) - VRR(i-2,j  ,k,l,m,n,N+1,PAx,PAy,PAz,PBx,PBy,PBz)) + \
                j*(VRR(i-1,j-1,k,l,m,n,N,PAx,PAy,PAz,PBx,PBy,PBz) - VRR(i-1,j-1,k,l,m,n,N+1,PAx,PAy,PAz,PBx,PBy,PBz)) 
            )
        elif (j>0) is S.true:
            return PBx*VRR(i,j-1,k,l,m,n,N,PAx,PAy,PAz,PBx,PBy,PBz) - PGx*VRR(i,j-1,k,l,m,n,N+1,PAx,PAy,PAz,PBx,PBy,PBz) + tpi*(
                i*(VRR(i-1,j-1,k,l,m,n,N,PAx,PAy,PAz,PBx,PBy,PBz) - VRR(i-1,j-1,k,l,m,n,N+1,PAx,PAy,PAz,PBx,PBy,PBz)) + \
                (j-1)*(VRR(i  ,j-2,k,l,m,n,N,PAx,PAy,PAz,PBx,PBy,PBz) - VRR(i  ,j-2,k,l,m,n,N+1,PAx,PAy,PAz,PBx,PBy,PBz)) 
            )
        elif (k>0) is S.true:
            return PAy*VRR(i,j,k-1,l,m,n,N,PAx,PAy,PAz,PBx,PBy,PBz) - PGy*VRR(i,j,k-1,l,m,n,N+1,PAx,PAy,PAz,PBx,PBy,PBz) + tpi*(
                (k-1)*(VRR(i,j,k-2,l  ,m,n,N,PAx,PAy,PAz,PBx,PBy,PBz) - VRR(i,j,k-2,l  ,m,n,N+1,PAx,PAy,PAz,PBx,PBy,PBz)) + \
                l*(VRR(i,j,k-1,l-1,m,n,N,PAx,PAy,PAz,PBx,PBy,PBz) - VRR(i,j,k-1,l-1,m,n,N+1,PAx,PAy,PAz,PBx,PBy,PBz)) 
            )
        elif (l>0) is S.true:
            return PBy*VRR(i,j,k,l-1,m,n,N,PAx,PAy,PAz,PBx,PBy,PBz) - PGy*VRR(i,j,k,l-1,m,n,N+1,PAx,PAy,PAz,PBx,PBy,PBz) + tpi*(
                k*(VRR(i,j,k-1,l-1,m,n,N,PAx,PAy,PAz,PBx,PBy,PBz) - VRR(i,j,k-1,l-1,m,n,N+1,PAx,PAy,PAz,PBx,PBy,PBz)) + \
                (l-1)*(VRR(i,j,k  ,l-2,m,n,N,PAx,PAy,PAz,PBx,PBy,PBz) - VRR(i,j,k  ,l-2,m,n,N+1,PAx,PAy,PAz,PBx,PBy,PBz)) 
            )
        elif (m>0) is S.true:
            return PAz*VRR(i,j,k,l,m-1,n,N,PAx,PAy,PAz,PBx,PBy,PBz) - PGz*VRR(i,j,k,l,m-1,n,N+1,PAx,PAy,PAz,PBx,PBy,PBz) + tpi*(
                (m-1)*(VRR(i,j,k,l,m-2,n  ,N,PAx,PAy,PAz,PBx,PBy,PBz) - VRR(i,j,k,l,m-2,n  ,N+1,PAx,PAy,PAz,PBx,PBy,PBz)) + \
                n*(VRR(i,j,k,l,m-1,n-1,N,PAx,PAy,PAz,PBx,PBy,PBz) - VRR(i,j,k,l,m-1,n-1,N+1,PAx,PAy,PAz,PBx,PBy,PBz)) 
            )
        elif (n>0) is S.true:
            return PBz*VRR(i,j,k,l,m,n-1,N,PAx,PAy,PAz,PBx,PBy,PBz) - PGz*VRR(i,j,k,l,m,n-1,N+1,PAx,PAy,PAz,PBx,PBy,PBz) + tpi*(
                m*(VRR(i,j,k,l,m-1,n-1,N,PAx,PAy,PAz,PBx,PBy,PBz) - VRR(i,j,k,l,m-1,n-1,N+1,PAx,PAy,PAz,PBx,PBy,PBz)) + \
                (n-1)*(VRR(i,j,k,l,m  ,n-2,N,PAx,PAy,PAz,PBx,PBy,PBz) - VRR(i,j,k,l,m  ,n-2,N+1,PAx,PAy,PAz,PBx,PBy,PBz)) 
            )
        else:
            raise RuntimeError("recurrences.VRR(): bad recurrence, should not be here!")

def do_VRR(lax,lay,laz,lbx,lby,lbz,same_center):
    PAx, PAy, PAz, PBx, PBy, PBz = (0,0,0,0,0,0) if same_center else symbols('PAx PAy PAz PBx PBy PBz')
    return VRR(lax,lbx,lay,lby,laz,lbz,0,PAx,PAy,PAz,PBx,PBy,PBz)

def do_VRR_multishell(vLA,vLB,same_center,verbose=False):
    vrr_exprs = []
    nci = num_cartesian_integrals(vLA, vLB)
    ic = 0
    
    if verbose:
        print(f"Evaluating VRR ({nci} final integrals).")
    
    for LA in vLA:
        for lax in range(LA,-1,-1):
            for lay in range(LA-lax,-1,-1):
                laz = LA-lax-lay
                for LB in vLB:
                    for lbx in range(LB,-1,-1):
                        for lby in range(LB-lbx,-1,-1):
                            lbz = LB-lbx-lby
                            vrr_exprs.append(do_VRR(lax,lay,laz,lbx,lby,lbz,same_center))
                            ic += 1
                            if verbose:
                                print(f"{ic}/{nci} : Evaluated [{lax},{lay},{laz}|{lbx},{lby},{lbz}]")
    if (ic != nci):
        raise RuntimeError(f"recurrences.do_VRR_multishell(): Expected {nci} integrals, made {ic} integrals!")
    return vrr_exprs

################################################################################################
#                                                                                              #
#   Horizontal Recurrence Relation (HRR):                                                      #
#                                                                                              #
#   [i(j+1)klmn] = [(i+1)jklmn] + (AB)_x*[ijklmn] ("left-to-right," ltr)                       #
#   [(i+1)jklmn] = [i(j+1)klmn] - (AB)_x*[ijklmn] ("right-to-left," rtl)                       #
#                                                                                              #
#   Begin : [(lax+lbx)0(lay+lby)0(laz+lbz)0]^[0] (ltr)                                         #
#           [0(lax+lbx)0(lay+lby)0(laz+lbz)]^[0] (rtl)                                         #
#                                                                                              #
#   Target: [(lax)(lbx)(lay)(lby)(laz)(lbz)]^[0]                                               #
#   Note: for same-center integrals A = B, (AB) = 0 and the second term in the HRR vanishes.   #
#                                                                                              #
################################################################################################

def get_hrr_blist(vLA, vLB, same_center):
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

def check_HRR_input(args):
    if len(args) != 13:
        raise ValueError(f"recurrences.HRR_*() expects 10 arguments but {len(args)} were given.")
    i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13 = args
    if i1.is_integer is False or \
       i2.is_integer is False or \
       i3.is_integer is False or \
       i4.is_integer is False or \
       i5.is_integer is False or \
       i6.is_integer is False or \
       i11.is_integer is False or \
       i12.is_integer is False:
        raise TypeError("recurrences.HRR_*() expects integer Lqn inputs.")

class HRR_ltr(Function):
    @classmethod
    def eval(cls, *args):
        check_HRR_input(args)
        i,j,k,l,m,n,ABx,ABy,ABz,hrr_blist,LA,LB,sc = args
        if sc:
            if ((i<0) is S.true) or ((j<0) is S.true) or \
               ((k<0) is S.true) or ((l<0) is S.true) or \
               ((m<0) is S.true) or ((n<0) is S.true):
                return 0
            elif (j==0) and (l==0) and (n==0):
                return Symbol(f"base_integral_{hrr_blist[(i,k,m,0,0,0,LA,LB)]}[g]")
            elif (j>0) is S.true:
                return HRR_ltr(i+1,j-1,k,l,m,n,ABx,ABy,ABz,hrr_blist,LA,LB,sc)
            elif (l>0) is S.true:
                return HRR_ltr(i,j,k+1,l-1,m,n,ABx,ABy,ABz,hrr_blist,LA,LB,sc)
            elif (n>0) is S.true:
                return HRR_ltr(i,j,k,l,m+1,n-1,ABx,ABy,ABz,hrr_blist,LA,LB,sc)
            else:
                raise RuntimeError("recurrences.HRR_ltr(): bad recurrence, should not be here!")
        else:
            if ((i<0) is S.true) or ((j<0) is S.true) or \
               ((k<0) is S.true) or ((l<0) is S.true) or \
               ((m<0) is S.true) or ((n<0) is S.true):
                return 0
            elif (j==0) and (l==0) and (n==0):
                return Symbol(f"base_integral_{hrr_blist[(i,k,m,0,0,0,LA,LB)]}[g]")
            elif (j>0) is S.true:
                return HRR_ltr(i+1,j-1,k,l,m,n,ABx,ABy,ABz,hrr_blist,LA,LB,sc) + ABx*HRR_ltr(i,j-1,k,l,m,n,ABx,ABy,ABz,hrr_blist,LA,LB,sc)
            elif (l>0) is S.true:
                return HRR_ltr(i,j,k+1,l-1,m,n,ABx,ABy,ABz,hrr_blist,LA,LB,sc) + ABy*HRR_ltr(i,j,k,l-1,m,n,ABx,ABy,ABz,hrr_blist,LA,LB,sc)
            elif (n>0) is S.true:
                return HRR_ltr(i,j,k,l,m+1,n-1,ABx,ABy,ABz,hrr_blist,LA,LB,sc) + ABz*HRR_ltr(i,j,k,l,m,n-1,ABx,ABy,ABz,hrr_blist,LA,LB,sc)
            else:
                raise RuntimeError("recurrences.HRR_ltr(): bad recurrence, should not be here!")

class HRR_rtl(Function):
    @classmethod
    def eval(cls, *args):
        check_HRR_input(args)
        i,j,k,l,m,n,ABx,ABy,ABz,hrr_blist,LA,LB,sc = args
        if sc:
            if ((i<0) is S.true) or ((j<0) is S.true) or \
               ((k<0) is S.true) or ((l<0) is S.true) or \
               ((m<0) is S.true) or ((n<0) is S.true):
                return 0
            elif (i==0) and (k==0) and (m==0):
                return Symbol(f"base_integral_{hrr_blist[(0,0,0,j,l,n,LA,LB)]}[g]")
            elif (i>0) is S.true:
                return HRR_rtl(i-1,j+1,k,l,m,n,ABx,ABy,ABz,hrr_blist,LA,LB,sc)
            elif (k>0) is S.true:
                return HRR_rtl(i,j,k-1,l+1,m,n,ABx,ABy,ABz,hrr_blist,LA,LB,sc)
            elif (m>0) is S.true:
                return HRR_rtl(i,j,k,l,m-1,n+1,ABx,ABy,ABz,hrr_blist,LA,LB,sc)
            else:
                raise RuntimeError("recurrences.HRR_ltr(): bad recurrence, should not be here!")
        else:
            if ((i<0) is S.true) or ((j<0) is S.true) or \
               ((k<0) is S.true) or ((l<0) is S.true) or \
               ((m<0) is S.true) or ((n<0) is S.true):
                return 0
            elif (i==0) and (k==0) and (m==0):
                return Symbol(f"base_integral_{hrr_blist[(0,0,0,j,l,n,LA,LB)]}[g]")
            elif (i>0) is S.true:
                return HRR_rtl(i-1,j+1,k,l,m,n,ABx,ABy,ABz,hrr_blist,LA,LB,sc) - ABx*HRR_rtl(i-1,j,k,l,m,n,ABx,ABy,ABz,hrr_blist,LA,LB,sc)
            elif (k>0) is S.true:
                return HRR_rtl(i,j,k-1,l+1,m,n,ABx,ABy,ABz,hrr_blist,LA,LB,sc) - ABy*HRR_rtl(i,j,k-1,l,m,n,ABx,ABy,ABz,hrr_blist,LA,LB,sc)
            elif (m>0) is S.true:
                return HRR_rtl(i,j,k,l,m-1,n+1,ABx,ABy,ABz,hrr_blist,LA,LB,sc) - ABz*HRR_rtl(i,j,k,l,m-1,n,ABx,ABy,ABz,hrr_blist,LA,LB,sc)
            else:
                raise RuntimeError("recurrences.HRR_ltr(): bad recurrence, should not be here!")

def do_HRR(lax,lay,laz,lbx,lby,lbz,same_center,hrr_blist):
    LA = lax+lay+laz
    LB = lbx+lby+lbz
    ABx, ABy, ABz = symbols('ABx ABy ABz')
    return HRR_ltr(lax,lbx,lay,lby,laz,lbz,ABx,ABy,ABz,hrr_blist,LA,LB,same_center) if (LA>=LB) else HRR_rtl(lax,lbx,lay,lby,laz,lbz,ABx,ABy,ABz,hrr_blist,LA,LB,same_center)

def do_HRR_multishell(vLA,vLB,same_center,hrr_blist,verbose=False):
    hrr_bi_exprs = []
    ncbi = len(hrr_blist)
    bic = 0
    
    if verbose:
        print(f"Evaluating HRR base integrals via VRR ({ncbi} base integrals).")
    
    for lqns, idx in hrr_blist.items():
        lax, lay, laz, lbx, lby, lbz = lqns[0], lqns[1], lqns[2], lqns[3], lqns[4], lqns[5]
        hrr_bi_exprs.append(do_VRR(lax,lay,laz,lbx,lby,lbz,same_center))
        bic += 1
        if verbose:
            print(f"{bic}/{ncbi} : Evaluated [{lax},{lay},{laz}|{lbx},{lby},{lbz}]")
    
    hrr_exprs= []
    nci = num_cartesian_integrals(vLA, vLB)
    ic = 0
    
    if verbose:
        print(f"Evaluating HRR ({nci} final integrals).")
    
    for LA in vLA:
        for lax in range(LA,-1,-1):
            for lay in range(LA-lax,-1,-1):
                laz = LA-lax-lay
                for LB in vLB:
                    for lbx in range(LB,-1,-1):
                        for lby in range(LB-lbx,-1,-1):
                            lbz = LB-lbx-lby
                            hrr_exprs.append(do_HRR(lax,lay,laz,lbx,lby,lbz,same_center,hrr_blist))
                            ic += 1
                            if verbose:
                                print(f"{ic}/{nci} : Evaluated ({lax},{lay},{laz}|{lbx},{lby},{lbz})")
    if (ic != nci):
        raise RuntimeError(f"recurrences.do_HRR_multishell(): Expected {nci} integrals, made {ic} integrals!")
    return hrr_bi_exprs, hrr_exprs

