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
         SUBROUTINE  OED__NAI_PREPARE_CTR
     +
     +                    ( NCSIZE,
     +                      NIJ,
     +                      NPGTOA,NPGTOB,
     +                      SHELLA,SHELLB,
     +                      ALPHAA,ALPHAB,
     +                      PREFACT,SPNORM,
     +                      EQUALAB,
     +                      BLOCKED,
     +                      RHO,
     +
     +                              NORMA,NORMB,
     +                              RHOAB,
     +                              CBATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_PREPARE_CTR
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine prepares for contracting the nuclear
C                attraction primitive batches. Everything that needs
C                still to be done at the stage when calling this
C                routine should be placed here. At the moment, we have
C                the following:
C
C                       i) copy the exponential prefactors
C                      ii) generate all the A and B norms
C                     iii) add the contribution due to the
C                          overall integral prefactor and the
C                          s-/p-shell type norm into one of
C                          the A or B norms, depending on
C                          which has the least elements
C                      iv) initialize the nuclear attraction
C                          contraction batch
C
C
C                  Input:
C
C                    NCSIZE       =  size of the contraction batch
C                    NIJ          =  total # of ij primitive index
C                                    pairs for the contracted shell
C                                    pair A,B
C                    NCEN         =  # of ij nuclear attraction centers
C                    NPGTOx       =  # of primitives per contraction
C                                    for contraction shells x = A,B
C                    SHELLx       =  the shell types for contraction
C                                    shells x = A,B
C                    ALPHAx       =  the primitive exponents for
C                                    contraction shells x = A,B
C                    PREFACT      =  overall prefactor for all integrals
C                    SPNORM       =  normalization factor due to
C                                    presence of s- and p-type shells.
C                                    For each s-type shell there is a
C                                    factor of 1 and for each p-type
C                                    shell a factor of 2
C                    EQUALAB      =  indicates, if csh A and csh A are
C                                    considered to be equal
C                    BLOCKED      =  if false, there will be no need
C                                    to block the contraction step over
C                                    the set of primitives and thus as
C                                    a consequence there is no need to
C                                    initialize the contraction batch
C                    RHO          =  NIJ exponential prefactors rho(a,b)
C
C                  Output:
C
C                    NORMx        =  the normalization factors due to
C                                    the primitive exponents for the
C                                    contraction shells x = A,B
C                    RHOAB        =  the complete set of NIJ exponential
C                                    prefactors between contraction
C                                    shells A and B
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
         LOGICAL    EQUALAB

         INTEGER    N
         INTEGER    NCSIZE
         INTEGER    NIJ
         INTEGER    NPGTOA,NPGTOB
         INTEGER    SHELLA,SHELLB

         DOUBLE PRECISION   FACTOR
         DOUBLE PRECISION   POWER
         DOUBLE PRECISION   PREFACT,SPNORM
         DOUBLE PRECISION   ZERO,HALF,ZP75

         DOUBLE PRECISION   CBATCH (1:NCSIZE)
         DOUBLE PRECISION   ALPHAA (1:NPGTOA)
         DOUBLE PRECISION   ALPHAB (1:NPGTOB)
         DOUBLE PRECISION   NORMA  (1:NPGTOA)
         DOUBLE PRECISION   NORMB  (1:NPGTOB)
         DOUBLE PRECISION   RHO    (1:NIJ)
         DOUBLE PRECISION   RHOAB  (1:NIJ)
         double precision   rho_temp(1:nij)

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
C                final RHOAB array might overlap with the input RHO
C                array. The following diagram shows how such overlap
C                might happen:
C
C                       RHO  ->  |  NIJ  |
C                                 \       \
C                                  \       \
C                                   | RHOAB |
C
C                We are always safe, if we start copying from the
C                last element of RHO downwards.
C
C
c         DO N = NIJ,1,-1
c            RHOAB (N) = RHO (N)
c         END DO

          do n = 1, nij
             rho_temp(n) = rho(n)
          enddo
#ifndef ERDFIX
          do n = 1, nij
             rhoab(n) = rho_temp(n)
          enddo
#endif

C
C
C             ...calculate the A and B norms.
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
C
C
C             ...rescale one of the A or B norms, which has the
C                least number of elements.
C
C
         FACTOR = PREFACT * SPNORM

         IF (NPGTOA.LE.NPGTOB) THEN
             DO N = 1,NPGTOA
                NORMA (N) = FACTOR * NORMA (N)
             END DO
         ELSE
             DO N = 1,NPGTOB
                NORMB (N) = FACTOR * NORMB (N)
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
          do n = 1, nij
             rhoab(n) = rho_temp(n)
          enddo
#endif
         RETURN
         END
