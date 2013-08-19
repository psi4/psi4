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
         SUBROUTINE  ERD__PREPARE_CTR
     +
     +                    ( NCSIZE,
     +                      NIJ,NKL,
     +                      NPGTOA,NPGTOB,NPGTOC,NPGTOD,
     +                      SHELLA,SHELLB,SHELLC,SHELLD,
     +                      ALPHAA,ALPHAB,ALPHAC,ALPHAD,
     +                      PREFACT,SPNORM,
     +                      EQUALAB,EQUALCD,
     +                      BLOCKED,
     +                      RHO,
     +
     +                              NORMA,NORMB,NORMC,NORMD,
     +                              RHOAB,RHOCD,
     +                              CBATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__PREPARE_CTR
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine prepares the for contracting the
C                primitive batches. Everything that needs still to be
C                done at the stage when calling this routine should be
C                placed here. At the moment, we have the following:
C
C                       i) copy the exponential prefactors
C                      ii) generate all the A,B,C,D norms
C                     iii) add the contribution due to the
C                          overall integral prefactor and the
C                          s-/p-shell type norm into one of
C                          the A,B,C,D norms, depending on
C                          which has the least elements
C                      iv) initialize the contraction batch
C
C
C                  Input:
C
C                    NCSIZE       =  size of the contraction batch
C                    NIJ(KL)      =  total # of ij(kl) primitive index
C                                    pairs for the contracted shell
C                                    pair A,B(C,D)
C                    NPGTOx       =  # of primitives per contraction
C                                    for contraction shells x = A,B,C,D
C                    SHELLx       =  the shell types for contraction
C                                    shells x = A,B,C,D
C                    ALPHAx       =  the primitive exponents for
C                                    contraction shells x = A,B,C,D
C                    PREFACT      =  overall prefactor for all integrals
C                    SPNORM       =  normalization factor due to
C                                    presence of s- and p-type shells.
C                                    For each s-type shell there is a
C                                    factor of 1 and for each p-type
C                                    shell a factor of 2
C                    EQUALxy      =  indicates, if csh x and csh y are
C                                    considered to be equal for the
C                                    pairs xy = AB and CD
C                    BLOCKED      =  if false, there will be no need
C                                    to block the contraction step over
C                                    the set of primitives and thus as
C                                    a consequence there is no need to
C                                    initialize the contraction batch
C                    RHO          =  NIJ exponential prefactors rho(a,b)
C                                    + NKL exponential prefactors
C                                    rho(c,d), in that order
C
C                  Output:
C
C                    NORMx        =  the normalization factors due to
C                                    the primitive exponents for the
C                                    contraction shells x = A,B,C,D
C                    RHOAB(CD)    =  the complete set of NIJ (NKL)
C                                    exponential prefactors between
C                                    contraction shells A and B
C                                    (C and D)
C                    CBATCH       =  contraction batch initialized
C                                    to zero (if needed)
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT   NONE

         LOGICAL    BLOCKED
         LOGICAL    EQUALAB,EQUALCD

         INTEGER    N
         INTEGER    NCSIZE
         INTEGER    NIJ,NKL
         INTEGER    NPGTOA,NPGTOB,NPGTOC,NPGTOD
         INTEGER    NPMIN
         INTEGER    SHELLA,SHELLB,SHELLC,SHELLD

         DOUBLE PRECISION   FACTOR
         DOUBLE PRECISION   POWER
         DOUBLE PRECISION   PREFACT,SPNORM
         DOUBLE PRECISION   ZERO,HALF,ZP75

         DOUBLE PRECISION   CBATCH (1:NCSIZE)
         DOUBLE PRECISION   ALPHAA (1:NPGTOA)
         DOUBLE PRECISION   ALPHAB (1:NPGTOB)
         DOUBLE PRECISION   ALPHAC (1:NPGTOC)
         DOUBLE PRECISION   ALPHAD (1:NPGTOD)
         DOUBLE PRECISION   NORMA  (1:NPGTOA)
         DOUBLE PRECISION   NORMB  (1:NPGTOB)
         DOUBLE PRECISION   NORMC  (1:NPGTOC)
         DOUBLE PRECISION   NORMD  (1:NPGTOD)
         DOUBLE PRECISION   RHO    (1:NIJ+NKL)
         DOUBLE PRECISION   RHOAB  (1:NIJ)
         DOUBLE PRECISION   RHOCD  (1:NIJ)
         double precision   rho_temp(1:nij+nkl)
         double precision   rho_temp2(1:nij+nkl)

         PARAMETER   (ZERO = 0.D0)
         PARAMETER   (HALF = 0.5D0)
         PARAMETER   (ZP75 = 0.75D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...copy the exponential prefactors. This has to be
C                done with care, as the memory location for the
C                final RHOAB and RHOCD arrays might overlap with
C                the input RHO array. The following diagram shows
C                how such overlap might happen:
C
C                       RHO  ->  |  NIJ  |  NKL  |
C                                 \       \       \
C                                  \       \       \
C                                   | RHOAB | RHOCD |
C
C                We are always safe, if we start copying from the
C                last element of RHO downwards.
C
C
c         DO N = NKL,1,-1
c            RHOCD (N) = RHO (NIJ+N)
c         END DO

         do n = 1, nkl
            rho_temp(n) = rho(nij+n)
         enddo

 
c         DO N = NIJ,1,-1
c            RHOAB (N) = RHO (N)
c         END DO

         do n = 1, nij
            rho_temp2(n) = rho(n)
         enddo

#ifndef ERDFIX
         do n = 1, nkl
            rhocd(n) = rho_temp(n)
         enddo

         do n = 1, nij
            rhoab(n) = rho_temp2(n)
         enddo
#endif
C
C
C             ...calculate the A,B,C,D norms.
C
C
         POWER = DFLOAT (SHELLA) * HALF + ZP75
         IF (EQUALAB) THEN
             DO N = 1,NPGTOA
                NORMA (N) = ALPHAA (N) ** POWER
                NORMB (N) = NORMA (N)
             END DO
         ELSE
             DO N = 1,NPGTOA
                NORMA (N) = ALPHAA (N) ** POWER
             END DO
             POWER = DFLOAT (SHELLB) * HALF + ZP75
             DO N = 1,NPGTOB
                NORMB (N) = ALPHAB (N) ** POWER
             END DO
         END IF

         POWER = DFLOAT (SHELLC) * HALF + ZP75
         IF (EQUALCD) THEN
             DO N = 1,NPGTOC
                NORMC (N) = ALPHAC (N) ** POWER
                NORMD (N) = NORMC (N)
             END DO
         ELSE
             DO N = 1,NPGTOC
                NORMC (N) = ALPHAC (N) ** POWER
             END DO
             POWER = DFLOAT (SHELLD) * HALF + ZP75
             DO N = 1,NPGTOD
                NORMD (N) = ALPHAD (N) ** POWER
             END DO
         END IF
C
C
C             ...rescale one of the A,B,C,D norms, which has the
C                least number of elements.
C
C
         FACTOR = PREFACT * SPNORM

         NPMIN = MIN (NPGTOA,NPGTOB,NPGTOC,NPGTOD)

         IF (NPGTOA.EQ.NPMIN) THEN
             DO N = 1,NPGTOA
                NORMA (N) = FACTOR * NORMA (N)
             END DO
         ELSE IF (NPGTOB.EQ.NPMIN) THEN
             DO N = 1,NPGTOB
                NORMB (N) = FACTOR * NORMB (N)
             END DO
         ELSE IF (NPGTOC.EQ.NPMIN) THEN
             DO N = 1,NPGTOC
                NORMC (N) = FACTOR * NORMC (N)
             END DO
         ELSE
             DO N = 1,NPGTOD
                NORMD (N) = FACTOR * NORMD (N)
             END DO
         END IF
C
C
C             ...initialize contraction batch (if necessary).
C
C
         IF (BLOCKED) THEN
             DO N = 1,NCSIZE
                CBATCH (N) = ZERO
             END DO
         END IF
C
C
C             ...ready!
C
C

#ifdef ERDFIX
         do n = 1, nkl
            rhocd(n) = rho_temp(n)
         enddo

         do n = 1, nij
            rhoab(n) = rho_temp2(n)
         enddo
#endif

         RETURN
         END
