C------------------------------------------------------------------------
C  INCLUDE FILE: ERD__TUNING
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  DESCRIPTION : Include file used to set all tuning parameters.
C                For optimum performance of the present two electron
C                integral package, there are 3 parameters that are
C                set and which allow for tuning of performance:
C
C
C                  L1CACHE   =   Size of level 1 cache in units of
C                                8 Byte
C
C                  TILE      =   Number of rows and columns in units
C                                of 8 Byte of level 1 cache square
C                                tile array used for performing optimum
C                                out-of-place matrix transposition
C                                operations
C                                
C                  NCTROW    =   Minimum number of rows that are
C                                still accepted for block contractions.
C                                This quantity controls extremly low
C                                # of row values due to small level
C                                1 cache size.
C
C                The value for TILE can be estimated as follows. Assume
C                the level 1 cache is 3/4 free. We have to fit 2 square
C                matrices (out-of-place transpositions!) of tile size
C                into it, giving:
C
C                       TILE = Int ( SQRT ((3/4)*L1CACHE/2)
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C                      >>>  COPYRIGHT NOTICE  <<<
C                             Program ERD
C
C    This routine is part of the integral program package ERD and as
C    such is copyrighted by N. Flocke, 2004; all rights reserved.
C    The full copyright notice is printed in the main routine of the
C    ERD package and all terms and conditions given therein are
C    incorporated here.
C
C------------------------------------------------------------------------
C
C
C             ...declare and set the tuning parameters.
C
C
         INTEGER   L1CACHE
         INTEGER   NCTROW
         INTEGER   TILE

         PARAMETER  (L1CACHE = 16384)
         PARAMETER  (TILE    = 64)
         PARAMETER  (NCTROW  = 50)
