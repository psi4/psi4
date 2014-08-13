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
         SUBROUTINE  OED__OVL3C_F00_PCGTO_BLOCK
     +
     +                    ( NBATCH,NINT1D,
     +                      ATOMIC,
     +                      MIJ,MK,MIJK,
     +                      NIJ,NIJBEG,NIJEND,NK,NKBEG,NKEND,
     +                      NPGTOI,NPGTOJ,NPGTOK,
     +                      NXYZFT,NXYZQ,
     +                      SHELLA,SHELLQ,
     +                      XA,YA,ZA,XI,YI,ZI,XJ,YJ,ZJ,XK,YK,ZK,
     +                      ALPHAI,ALPHAJ,ALPHAK,
     +                      PRIMI,PRIMJ,PRIMK,
     +                      NORMI,NORMJ,NORMK,
     +                      RHOIJK,
     +                      QAX,QAY,QAZ,QINVHF,SCALE,
     +                      INT1DX,INT1DY,INT1DZ,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__OVL3C_F00_PCGTO_BLOCK
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__OVL_1D_INTEGRALS
C                OED__OVL_INT1D_TO_E0
C  DESCRIPTION : This operation calculates a batch of unnormed 3-center
C                overlap integrals between primitive cartesian gaussians
C                for the shell range:
C
C                              [F00]     , F = A to Q = A + B + C
C                                   ij,k
C
C                for the block of ij exponent pairs and the block
C                of k exponents. The total number of 3-center overlap
C                integrals generated here is thus given by the total
C                number of cartesian monomials NXYZFT times the total
C                number of exponents MIJK in the present block.
C
C                On exit, the batch elements will be stored as:
C
C                             batch (ij,k,nxyzft)
C
C                Note especially that the labels corresponding to the
C                cartesian monomials part refer to the A,B,C shells,
C                while the the labels corresponding to the primitives
C                and contraction indices refer to the I,J,K shells.
C
C
C                  Input:
C
C                    NBATCH       =  size of the primitive cartesian
C                                    integral batch
C                    NINT1D       =  space needed for each of the 1D
C                                    X,Y,Z integral arrays
C                    ATOMIC       =  indicates, if purely atomic
C                                    integrals will be evaluated
C                    MIJ          =  current # of ij primitive index
C                                    pairs corresponding to the
C                                    contracted shell pair I,J
C                    MK           =  current # of k exponent indices
C                                    corresponding to the contracted
C                                    shell K
C                    MIJK         =  current # of ij primitive index
C                                    pairs times current # of k
C                                    exponent indices
C                    NIJ          =  total # of ij primitive index
C                                    pairs for the contracted shell
C                                    pair I,J
C                    NIJBEG(END)  =  first(last) ij primitive index
C                                    defining the ij block
C                    NK           =  total # of k exponent indices
C                                    corresponding to the contracted
C                                    shell K
C                    NKBEG(END)   =  first(last) k primitive index
C                                    defining the k block
C                    NPGTOx       =  # of primitives per contraction
C                                    for contraction shells x = I,J,K
C                    NXYZFT       =  sum of # of cartesian monomials
C                                    for all shells in the range
C                                    F = A,...,Q=A+B+C
C                    NXYZQ        =  # of cartesian monomials for
C                                    the Q=A+B+C shell
C                    SHELLx       =  the shell type for contraction
C                                    shells x = A and Q=A+B+C
C                    Xx,Yx,Zx     =  the x,y,z-coordinates for centers
C                                    x = A,I,J,K
C                    ALPHAx       =  the primitive exponents for
C                                    contraction shells x = I,J,K
C                    PRIMx        =  i,j,k labels of primitives for
C                                    the respective contraction shells
C                                    x = I,J,K
C                    NORMx        =  the normalization factors due to
C                                    the primitive exponents for the
C                                    contraction shells x = I,J,K
C                    RHOIJK       =  the complete set of NIJ * NK
C                                    exponential prefactors between
C                                    contraction shells I,J and K
C                    QAx          =  will hold current MIJK coordinate
C                                    x=X,Y,Z differences Q-A between
C                                    centers Q and A
C                    QINVHF       =  will hold current MIJK values of
C                                    1/(2*Q), where Q are the exponent
C                                    sums for contraction shells I,J
C                                    and K
C                    SCALE        =  will hold current MIJK values of
C                                    scaling factors
C                    INT1Dx       =  will hold all current 1D integrals
C                                    for each cartesian component
C                                    (x = X,Y,Z)
C
C                  Output:
C
C                    BATCH        =  current batch of primitive
C                                    cartesian 3-center overlap [F00]
C                                    integrals
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     ATOMIC

         INTEGER     I,J,K,M,P
         INTEGER     IJ
         INTEGER     MIJ,MK,MIJK
         INTEGER     NBATCH,NINT1D
         INTEGER     NIJ,NIJBEG,NIJEND,NK,NKBEG,NKEND
         INTEGER     NPGTOI,NPGTOJ,NPGTOK
         INTEGER     NXYZFT,NXYZQ
         INTEGER     SHELLA,SHELLQ

         INTEGER     PRIMI (1:NIJ)
         INTEGER     PRIMJ (1:NIJ)
         INTEGER     PRIMK (1:NK)

         DOUBLE PRECISION  EXPI,EXPJ,EXPK
         DOUBLE PRECISION  KNORM
         DOUBLE PRECISION  QINV
         DOUBLE PRECISION  QXI,QXJ,QXK,QYI,QYJ,QYK,QZI,QZJ,QZK
         DOUBLE PRECISION  XA,YA,ZA,XI,YI,ZI,XJ,YJ,ZJ,XK,YK,ZK
         DOUBLE PRECISION  ZERO,HALF,ONE,ONEP5

         DOUBLE PRECISION  ALPHAI  (1:NPGTOI)
         DOUBLE PRECISION  ALPHAJ  (1:NPGTOJ)
         DOUBLE PRECISION  ALPHAK  (1:NPGTOK)

         DOUBLE PRECISION  NORMI   (1:NPGTOI)
         DOUBLE PRECISION  NORMJ   (1:NPGTOJ)
         DOUBLE PRECISION  NORMK   (1:NPGTOK)

         DOUBLE PRECISION  RHOIJK  (1:NIJ,1:NK)

         DOUBLE PRECISION  BATCH   (1:NBATCH)

         DOUBLE PRECISION  QAX     (1:MIJK)
         DOUBLE PRECISION  QAY     (1:MIJK)
         DOUBLE PRECISION  QAZ     (1:MIJK)
         DOUBLE PRECISION  QINVHF  (1:MIJK)
         DOUBLE PRECISION  SCALE   (1:MIJK)

         DOUBLE PRECISION  INT1DX  (1:NINT1D)
         DOUBLE PRECISION  INT1DY  (1:NINT1D)
         DOUBLE PRECISION  INT1DZ  (1:NINT1D)

         PARAMETER  (ZERO  = 0.D0)
         PARAMETER  (HALF  = 0.5D0)
         PARAMETER  (ONE   = 1.D0)
         PARAMETER  (ONEP5 = 1.5D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...calculate the quantities needed to establish the
C                1D overlap integrals.
C
C
         IF (ATOMIC) THEN
             M = 0
             DO 100 P = NKBEG,NKEND
                K = PRIMK (P)
                EXPK = ALPHAK (K)
                KNORM = NORMK (K)
                DO 110 IJ = NIJBEG,NIJEND
                   M = M + 1
                   I = PRIMI (IJ)
                   J = PRIMJ (IJ)
                   QINV = ONE / (ALPHAI (I) + ALPHAJ (J) + EXPK)
                   QINVHF (M) = HALF * QINV
                   SCALE  (M) = NORMI (I) * NORMJ (J) * KNORM
     +                           * (QINV ** ONEP5)
  110           CONTINUE
  100        CONTINUE
         ELSE
             M = 0
             DO 200 P = NKBEG,NKEND
                K = PRIMK (P)
                EXPK = ALPHAK (K)
                KNORM = NORMK (K)
                QXK = EXPK * XK
                QYK = EXPK * YK
                QZK = EXPK * ZK
                DO 220 IJ = NIJBEG,NIJEND
                   M = M + 1
                   I = PRIMI (IJ)
                   J = PRIMJ (IJ)
                   EXPI = ALPHAI (I)
                   EXPJ = ALPHAJ (J)
                   QXI = EXPI * XI
                   QYI = EXPI * YI
                   QZI = EXPI * ZI
                   QXJ = EXPJ * XJ
                   QYJ = EXPJ * YJ
                   QZJ = EXPJ * ZJ
                   QINV = ONE / (EXPI + EXPJ + EXPK)
                   QAX (M) = QINV * (QXI + QXJ + QXK) - XA
                   QAY (M) = QINV * (QYI + QYJ + QYK) - YA
                   QAZ (M) = QINV * (QZI + QZJ + QZK) - ZA
                   QINVHF (M) = HALF * QINV
                   SCALE  (M) = NORMI (I) * NORMJ (J) * KNORM
     +                           * RHOIJK (IJ,K) * (QINV ** ONEP5)
  220           CONTINUE
  200        CONTINUE
         END IF
C
C
C             ...perform the following steps:
C
C                1) construct all 1D x,y,z overlap integrals for
C                   all present ij,k triples.
C
C                2) assemble the complete [F00] 3-center overlap
C                   batch for all ij,k triples using the 1D integrals.
C                   Arrays QAX and QAY are passed as scratch arrays.
C                   Note, that for that purpose we can use the routine
C                   which generates the [E|0] batches for the normal
C                   2-center overlap integrals.
C
C
         CALL    OED__OVL_1D_INTEGRALS
     +
     +                ( SHELLQ,0,
     +                  ATOMIC,
     +                  MIJK,
     +                  QAX,QAY,QAZ,
     +                  QINVHF,
     +                  ZERO,ZERO,ZERO,
     +
     +                            INT1DX,
     +                            INT1DY,
     +                            INT1DZ )
     +
     +
         CALL    OED__OVL_INT1D_TO_E0
     +
     +                ( SHELLA,SHELLQ,
     +                  ATOMIC,
     +                  MIJK,
     +                  NXYZFT,NXYZQ,
     +                  INT1DX,INT1DY,INT1DZ,
     +                  QAX,QAY,
     +                  SCALE,
     +
     +                            BATCH )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
