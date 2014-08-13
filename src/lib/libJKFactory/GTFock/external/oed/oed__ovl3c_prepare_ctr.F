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
         SUBROUTINE  OED__OVL3C_PREPARE_CTR
     +
     +                    ( NCSIZE,
     +                      NIJK,
     +                      NPGTOI,NPGTOJ,NPGTOK,
     +                      SHELLI,SHELLJ,SHELLK,
     +                      ALPHAI,ALPHAJ,ALPHAK,
     +                      PREFACT,SPNORM,
     +                      EQUALIJ,EQUALIK,EQUALJK,
     +                      BLOCKED,
     +                      RHO,
     +
     +                              NORMI,NORMJ,NORMK,
     +                              RHOIJK,
     +                              CBATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__OVL3C_PREPARE_CTR
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine prepares for contracting the 3-center
C                overlap primitive batches. Everything that needs still
C                to be done at the stage when calling this routine
C                should be placed here. At the moment, we have the
C                following:
C
C                       i) copy the IJK exponential prefactors
C                      ii) generate all the I,J and K norms
C                     iii) add the contribution due to the
C                          overall integral prefactor and the
C                          s-/p-shell type norm into one of
C                          the I,J or K norms, depending on
C                          which has the least elements
C                      iv) initialize the 3-center overlap
C                          contraction batch
C
C
C                  Input:
C
C                    NCSIZE       =  size of the contraction batch
C                    NIJK         =  # of ij primitive index pairs
C                                    times # of k primitive indices
C                    NPGTOx       =  # of primitives per contraction
C                                    for contraction shells x = I,J,K
C                    SHELLx       =  the shell types for contraction
C                                    shells x = I,J,K
C                    ALPHAx       =  the primitive exponents for
C                                    contraction shells x = I,J,K
C                    PREFACT      =  overall prefactor for all integrals
C                    SPNORM       =  normalization factor due to
C                                    presence of s- and p-type shells.
C                                    For each s-type shell there is a
C                                    factor of 1 and for each p-type
C                                    shell a factor of 2
C                    EQUALxy      =  indicates, if csh x and csh y are
C                                    considered to be equal for the
C                                    csh pairs xy = IJ,IK,JK
C                    BLOCKED      =  if false, there will be no need
C                                    to block the contraction step over
C                                    the set of primitives and thus as
C                                    a consequence there is no need to
C                                    initialize the contraction batch
C                    RHO          =  the NIJK exponential prefactors
C                                    rho (i,j,k)
C
C                  Output:
C
C                    NORMx        =  the normalization factors due to
C                                    the primitive exponents for the
C                                    contraction shells x = I,J,K
C                    RHOIJK       =  the copied set of NIJK exponential
C                                    prefactors between contraction
C                                    shells I,J and K
C                    CBATCH       =  contraction batch initialized
C                                    to zero (if needed)
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT   NONE

         LOGICAL    BLOCKED
         LOGICAL    EQUALIJ,EQUALIK,EQUALJK

         INTEGER    N
         INTEGER    NCSIZE
         INTEGER    NIJK
         INTEGER    NPGTOI,NPGTOJ,NPGTOK
         INTEGER    NPMIN
         INTEGER    SHELLI,SHELLJ,SHELLK

         DOUBLE PRECISION   FACTOR
         DOUBLE PRECISION   POWER
         DOUBLE PRECISION   PREFACT,SPNORM
         DOUBLE PRECISION   ZERO,HALF,ZP75

         DOUBLE PRECISION   CBATCH (1:NCSIZE)
         DOUBLE PRECISION   ALPHAI (1:NPGTOI)
         DOUBLE PRECISION   ALPHAJ (1:NPGTOJ)
         DOUBLE PRECISION   ALPHAK (1:NPGTOK)
         DOUBLE PRECISION   NORMI  (1:NPGTOI)
         DOUBLE PRECISION   NORMJ  (1:NPGTOJ)
         DOUBLE PRECISION   NORMK  (1:NPGTOK)
         DOUBLE PRECISION   RHO    (1:NIJK)
         DOUBLE PRECISION   RHOIJK (1:NIJK)
         double precision   rho_temp(1:nijk)

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
C                final RHOIJK array might overlap with the input RHO
C                array. The following diagram shows how such overlap
C                might happen:
C
C                       RHO  ->  |  NIJK  |
C                                 \        \
C                                  \        \
C                                   | RHOIJK |
C
C                We are always safe, if we start copying from the
C                last element of RHO downwards.
C
C
c         DO N = NIJK,1,-1
c            RHOIJK (N) = RHO (N)
c         END DO
    
          do n = 1, nijk
             rho_temp(n) = rho(n) 
          enddo 
#ifndef ERDFIX
          do n = 1, nijk
             rhoijk(n) = rho_temp(n)
          enddo
#endif

C
C
C             ...calculate the I,J and K norms.
C
C
         IF (EQUALIJ .AND. EQUALIK) THEN

             POWER = DFLOAT (SHELLI) * HALF + ZP75
             DO N = 1,NPGTOI
                NORMI (N) = ALPHAI (N) ** POWER
                NORMJ (N) = NORMI (N)
                NORMK (N) = NORMK (N)
             END DO

         ELSE IF (EQUALIJ) THEN

             POWER = DFLOAT (SHELLI) * HALF + ZP75
             DO N = 1,NPGTOI
                NORMI (N) = ALPHAI (N) ** POWER
                NORMJ (N) = NORMI (N)
             END DO
             POWER = DFLOAT (SHELLK) * HALF + ZP75
             DO N = 1,NPGTOK
                NORMK (N) = ALPHAK (N) ** POWER
             END DO

         ELSE IF (EQUALIK) THEN

             POWER = DFLOAT (SHELLI) * HALF + ZP75
             DO N = 1,NPGTOI
                NORMI (N) = ALPHAI (N) ** POWER
                NORMK (N) = NORMI (N)
             END DO
             POWER = DFLOAT (SHELLJ) * HALF + ZP75
             DO N = 1,NPGTOJ
                NORMJ (N) = ALPHAJ (N) ** POWER
             END DO

         ELSE IF (EQUALJK) THEN

             POWER = DFLOAT (SHELLJ) * HALF + ZP75
             DO N = 1,NPGTOJ
                NORMJ (N) = ALPHAJ (N) ** POWER
                NORMK (N) = NORMJ (N)
             END DO
             POWER = DFLOAT (SHELLI) * HALF + ZP75
             DO N = 1,NPGTOI
                NORMI (N) = ALPHAI (N) ** POWER
             END DO

         ELSE

             POWER = DFLOAT (SHELLI) * HALF + ZP75
             DO N = 1,NPGTOI
                NORMI (N) = ALPHAI (N) ** POWER
             END DO
             POWER = DFLOAT (SHELLJ) * HALF + ZP75
             DO N = 1,NPGTOJ
                NORMJ (N) = ALPHAJ (N) ** POWER
             END DO
             POWER = DFLOAT (SHELLK) * HALF + ZP75
             DO N = 1,NPGTOK
                NORMK (N) = ALPHAK (N) ** POWER
             END DO

         END IF
C
C
C             ...rescale one of the I,J or K norms, which has the
C                least number of elements.
C
C
         FACTOR = PREFACT * SPNORM

         NPMIN = MIN (NPGTOI,NPGTOJ,NPGTOK)

         IF (NPGTOI.EQ.NPMIN) THEN
             DO N = 1,NPGTOI
                NORMI (N) = FACTOR * NORMI (N)
             END DO
         ELSE IF (NPGTOJ.EQ.NPMIN) THEN
             DO N = 1,NPGTOJ
                NORMJ (N) = FACTOR * NORMJ (N)
             END DO
         ELSE
             DO N = 1,NPGTOK
                NORMK (N) = FACTOR * NORMK (N)
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
          do n = 1, nijk
             rhoijk(n) = rho_temp(n)
          enddo
#endif
         RETURN
         END
