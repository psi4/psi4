C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
         SUBROUTINE  OED__OVL3C_SET_IJK
     +
     +                    ( NCGTOA,NCGTOB,NCGTOC,
     +                      NPGTOA,NPGTOB,NPGTOC,
     +                      SHELLA,SHELLB,SHELLC,
     +                      EQUALAB,EQUALAC,EQUALBC,
     +                      XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,
     +                      RNABSQ,RNACSQ,RNBCSQ,
     +                      INDEXA,INDEXB,INDEXC,
     +                      LEXPA,LEXPB,LEXPC,
     +                      LCCA,LCCB,LCCC,
     +                      LCCSEGA,LCCSEGB,LCCSEGC,
     +
     +                             NCGTOI,NCGTOJ,NCGTOK,NCGTOIJ,
     +                             NPGTOI,NPGTOJ,NPGTOK,NPGTOIJ,
     +                             SHELLI,SHELLJ,SHELLK,
     +                             EQUALIJ,EQUALIK,EQUALJK,
     +                             XI,YI,ZI,XJ,YJ,ZJ,XK,YK,ZK,
     +                             RNIJSQ,RNIKSQ,RNJKSQ,
     +                             SWAPAC,SWAPBC,SWAPRS,
     +                             INDEXI,INDEXJ,INDEXK,
     +                             LEXPI,LEXPJ,LEXPK,
     +                             LCCI,LCCJ,LCCK,
     +                             LCCSEGI,LCCSEGJ,LCCSEGK )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__OVL3C_SET_IJK
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation defines the abc -> ijk map for optimum
C                performance of the contraction steps.
C
C                The contraction step from primitive to contracted
C                [F00] 3-center overlap integrals will be performed
C                in two steps in the following order:
C
C                        1) contract over a single index k
C                        2) contract over index pairs ij.
C
C                Only the last contraction step over the index pairs
C                can accomodate any restriction i>=j due to shell
C                equalities. Since shell equalities can arise between
C                all possible shell pairs ab,ac and bc, a remapping
C                has to be performed: a,b,c -> i,j,k to adjust for
C                the contraction steps. Two criteria are used for
C                the mapping:
C
C                  If any two shells xy are equal then:
C
C                              xy
C                              ab,c -> ij,k    (no swap)
C                              ac,b -> ij,k    (swap b <-> c)
C                              cb,a -> ij,k    (swap a <-> c)
C                  else
C                     if # of 'a' prim >= b,c:  a -> k   (swap a <-> c)
C                     if # of 'b' prim >= a,c:  b -> k   (swap b <-> c)
C                     if # of 'c' prim >= a,b:  c -> k   (no swap)
C                  endif
C
C                Note, that this scheme does not consider the purely
C                atomic case where all shells are equal, which could
C                in principle use the index restrictions i>=j>=k.
C                The reason lies in the way the contraction step is
C                splitted into the two parts mentioned above. However
C                these cases present only a minute fraction during
C                an electronic structure calculation, hence efficieny
C                of the package will not be affected.
C
C                In case any two shells are equal, the corresponding
C                i>=j index restriction is transmitted via setting
C                EQUALIJ = .TRUE.
C
C
C                  Input (x = A,B and C):
C
C                   NCGTOx        =  # of contractions for csh x
C                   NPGTOx        =  # of primitives per contraction
C                                    for csh x
C                   SHELLx        =  the shell type for csh x
C                   EQUALxy       =  indicates, if csh x and csh y are
C                                    considered to be equal for the
C                                    pairs AB,AC and BC
C                   Xy,Yy,Zy      =  the x,y,z-coordinates for centers
C                                    y = A,B and C
C                   RNxySQ        =  square of the magnitude of the
C                                    distance between centers xy = AB,AC
C                                    and BC
C                   INDEXx        =  index A,B,C -> 1,2,3 map
C                   LEXPx         =  pointers to locate appropriate
C                                    section of the exponent array
C                                    corresponding to csh x
C                   LCCx          =  pointers to locate appropriate
C                                    section of the contraction coeff
C                                    array corresponding to csh x
C                   LCCSEGx       =  pointers to locate appropriate
C                                    section of the lowest and highest
C                                    primitive index array defining
C                                    segmented contraction boundaries
C                                    for csh x
C
C                  Output (x = I,J and K):
C
C                   NCGTOx        =  # of contractions for csh x
C                   NCGTOIJ       =  # of contractions to handle for
C                                    the csh pair I,J
C                   NPGTOx        =  # of primitives per contraction
C                                    for csh x
C                   NPGTOIJ       =  # of primitive pairs to handle for
C                                    contraction csh pair I,J
C                   SHELLx        =  the shell type for csh x
C                   EQUALxy       =  indicates, if csh x and csh y are
C                                    considered to be equal for the
C                                    pairs IJ,IK and JK
C                   Xy,Yy,Zy      =  the x,y,z-coordinates for centers
C                                    y = I,J and K
C                   RNxySQ        =  square of the magnitude of the
C                                    distance between centers xy = IJ,IK
C                                    and JK
C                   SWAPxy        =  is .true., if a swap x <-> y has
C                                    been performed for xy = AC and BC
C                   SWAPRS        =  is set .true. if the contraction
C                                    order of the primitives pair IJ
C                                    will be performed in reverse order
C                                    JI for efficiency reasons
C                   INDEXx        =  index I,J,K -> 1,2,3 map
C                   LEXPx         =  pointers to locate appropriate
C                                    section of the exponent array
C                                    corresponding to csh x
C                   LCCx          =  pointers to locate appropriate
C                                    section of the contraction coeff
C                                    array corresponding to csh x
C                   LCCSEGx       =  pointers to locate appropriate
C                                    section of the lowest and highest
C                                    primitive index array defining
C                                    segmented contraction boundaries
C                                    for csh x
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT  NONE

         LOGICAL   EQUALAB,EQUALAC,EQUALBC
         LOGICAL   EQUALIJ,EQUALIK,EQUALJK
         LOGICAL   NPAMAX,NPBMAX
         LOGICAL   SWAPAC,SWAPBC,SWAPRS

         INTEGER   INDEXA,INDEXB,INDEXC
         INTEGER   INDEXI,INDEXJ,INDEXK
         INTEGER   LCCA,LCCB,LCCC
         INTEGER   LCCI,LCCJ,LCCK
         INTEGER   LCCSEGA,LCCSEGB,LCCSEGC
         INTEGER   LCCSEGI,LCCSEGJ,LCCSEGK
         INTEGER   LEXPA,LEXPB,LEXPC
         INTEGER   LEXPI,LEXPJ,LEXPK
         INTEGER   NCGTOA,NCGTOB,NCGTOC
         INTEGER   NCGTOI,NCGTOJ,NCGTOK,NCGTOIJ
         INTEGER   NPGTOA,NPGTOB,NPGTOC
         INTEGER   NPGTOI,NPGTOJ,NPGTOK,NPGTOIJ
         INTEGER   SHELLA,SHELLB,SHELLC
         INTEGER   SHELLI,SHELLJ,SHELLK

         DOUBLE PRECISION  RNABSQ,RNACSQ,RNBCSQ
         DOUBLE PRECISION  RNIJSQ,RNIKSQ,RNJKSQ
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB,XC,YC,ZC
         DOUBLE PRECISION  XI,YI,ZI,XJ,YJ,ZJ,XK,YK,ZK
C
C
C------------------------------------------------------------------------
C
C
C             ...determine primitives index restrictions and their
C                index order. Remap associated data: a,b,c -> i,j,k.
C
C
         NPAMAX = (NPGTOA.GE.NPGTOB) .AND. (NPGTOA.GE.NPGTOC)
         NPBMAX = (NPGTOB.GE.NPGTOA) .AND. (NPGTOB.GE.NPGTOC)

         IF (EQUALAB) THEN
             SWAPAC = .FALSE.
             SWAPBC = .FALSE.
             NPGTOIJ = NPGTOA * (NPGTOA + 1) / 2
             NCGTOIJ = NCGTOA * (NCGTOA + 1) / 2
         ELSE IF (EQUALAC) THEN
             SWAPAC = .FALSE.
             SWAPBC = .TRUE.
             NPGTOIJ = NPGTOA * (NPGTOA + 1) / 2
             NCGTOIJ = NCGTOA * (NCGTOA + 1) / 2
         ELSE IF (EQUALBC) THEN
             SWAPAC = .TRUE.
             SWAPBC = .FALSE.
             NPGTOIJ = NPGTOB * (NPGTOB + 1) / 2
             NCGTOIJ = NCGTOB * (NCGTOB + 1) / 2
         ELSE
             IF (NPAMAX) THEN
                 SWAPAC = .TRUE.
                 SWAPBC = .FALSE.
                 NPGTOIJ = NPGTOB * NPGTOC
                 NCGTOIJ = NCGTOB * NCGTOC
             ELSE IF (NPBMAX) THEN
                 SWAPAC = .FALSE.
                 SWAPBC = .TRUE.
                 NPGTOIJ = NPGTOA * NPGTOC
                 NCGTOIJ = NCGTOA * NCGTOC
             ELSE
                 SWAPAC = .FALSE.
                 SWAPBC = .FALSE.
                 NPGTOIJ = NPGTOA * NPGTOB
                 NCGTOIJ = NCGTOA * NCGTOB
             END IF
         END IF

         IF (SWAPAC) THEN
             XI = XC
             XJ = XB
             XK = XA
             YI = YC
             YJ = YB
             YK = YA
             ZI = ZC
             ZJ = ZB
             ZK = ZA
             EQUALIJ = EQUALBC
             EQUALIK = EQUALAC
             EQUALJK = EQUALAB
             LCCI = LCCC
             LCCJ = LCCB
             LCCK = LCCA
             LCCSEGI = LCCSEGC
             LCCSEGJ = LCCSEGB
             LCCSEGK = LCCSEGA
             LEXPI = LEXPC
             LEXPJ = LEXPB
             LEXPK = LEXPA
             INDEXI = INDEXC
             INDEXJ = INDEXB
             INDEXK = INDEXA
             NCGTOI = NCGTOC
             NCGTOJ = NCGTOB
             NCGTOK = NCGTOA
             NPGTOI = NPGTOC
             NPGTOJ = NPGTOB
             NPGTOK = NPGTOA
             SHELLI = SHELLC
             SHELLJ = SHELLB
             SHELLK = SHELLA
             RNIJSQ = RNBCSQ
             RNIKSQ = RNACSQ
             RNJKSQ = RNABSQ
         ELSE IF (SWAPBC) THEN
             XI = XA
             XJ = XC
             XK = XB
             YI = YA
             YJ = YC
             YK = YB
             ZI = ZA
             ZJ = ZC
             ZK = ZB
             EQUALIJ = EQUALAC
             EQUALIK = EQUALAB
             EQUALJK = EQUALBC
             LCCI = LCCA
             LCCJ = LCCC
             LCCK = LCCB
             LCCSEGI = LCCSEGA
             LCCSEGJ = LCCSEGC
             LCCSEGK = LCCSEGB
             LEXPI = LEXPA
             LEXPJ = LEXPC
             LEXPK = LEXPB
             INDEXI = INDEXA
             INDEXJ = INDEXC
             INDEXK = INDEXB
             NCGTOI = NCGTOA
             NCGTOJ = NCGTOC
             NCGTOK = NCGTOB
             NPGTOI = NPGTOA
             NPGTOJ = NPGTOC
             NPGTOK = NPGTOB
             SHELLI = SHELLA
             SHELLJ = SHELLC
             SHELLK = SHELLB
             RNIJSQ = RNACSQ
             RNIKSQ = RNABSQ
             RNJKSQ = RNBCSQ
         ELSE
             XI = XA
             XJ = XB
             XK = XC
             YI = YA
             YJ = YB
             YK = YC
             ZI = ZA
             ZJ = ZB
             ZK = ZC
             EQUALIJ = EQUALAB
             EQUALIK = EQUALAC
             EQUALJK = EQUALBC
             LCCI = LCCA
             LCCJ = LCCB
             LCCK = LCCC
             LCCSEGI = LCCSEGA
             LCCSEGJ = LCCSEGB
             LCCSEGK = LCCSEGC
             LEXPI = LEXPA
             LEXPJ = LEXPB
             LEXPK = LEXPC
             INDEXI = INDEXA
             INDEXJ = INDEXB
             INDEXK = INDEXC
             NCGTOI = NCGTOA
             NCGTOJ = NCGTOB
             NCGTOK = NCGTOC
             NPGTOI = NPGTOA
             NPGTOJ = NPGTOB
             NPGTOK = NPGTOC
             SHELLI = SHELLA
             SHELLJ = SHELLB
             SHELLK = SHELLC
             RNIJSQ = RNABSQ
             RNIKSQ = RNACSQ
             RNJKSQ = RNBCSQ
         END IF
C
C
C             ...set the control variable to determine contraction
C                order between csh I and J
C
C
         SWAPRS = NPGTOI .GT. NPGTOJ
C
C
C             ...ready!
C
C
         RETURN
         END
