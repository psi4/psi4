#if HAVE_CONFIG_H
#   include "config.fh"
#endif
#include "galinalg.fh"
      SUBROUTINE TEST_DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,
     $                      LDB,BETA,C,LDC)
C     .. Parameters ..
      DOUBLE PRECISION ONE, ZERO
      PARAMETER        (ONE=1.0D+0,ZERO=0.0D+0)
      INTEGER          MB, NB, KB
      PARAMETER        (MB=64,NB=MB,KB=64)
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA, BETA
      INTEGER          K, LDA, LDB, LDC, M, N
      CHARACTER        TRANSA, TRANSB
C     .. Array Arguments ..
C
C     Purpose
C     =======
C
C     DGEMM  performs one of the matrix-matrix operations
C
C     C := alpha*op( A )*op( B ) + beta*C,
C
C     where  op( X ) is one of
C
C     op( X ) = X   or   op( X ) = X',
C
C     alpha and beta are scalars, and A, B and C are matrices, with op( A )
C     an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
C
C     Parameters
C     ==========
C
C     TRANSA - CHARACTER*1.
C     On entry, TRANSA specifies the form of op( A ) to be used in
C     the matrix multiplication as follows:
C
C     TRANSA = 'N' or 'n',  op( A ) = A.
C
C     TRANSA = 'T' or 't',  op( A ) = A'.
C
C     TRANSA = 'C' or 'c',  op( A ) = A'.
C
C     Unchanged on exit.
C
C     TRANSB - CHARACTER*1.
C     On entry, TRANSB specifies the form of op( B ) to be used in
C     the matrix multiplication as follows:
C
C     TRANSB = 'N' or 'n',  op( B ) = B.
C
C     TRANSB = 'T' or 't',  op( B ) = B'.
C
C     TRANSB = 'C' or 'c',  op( B ) = B'.
C
C     Unchanged on exit.
C
C     M      - INTEGER.
C     On entry,  M  specifies  the number  of rows  of the  matrix
C     op( A )  and of the  matrix  C.  M  must  be at least  zero.
C     Unchanged on exit.
C
C     N      - INTEGER.
C     On entry,  N  specifies the number  of columns of the matrix
C     op( B ) and the number of columns of the matrix C. N must be
C     at least zero.
C     Unchanged on exit.
C
C     K      - INTEGER.
C     On entry,  K  specifies  the number of columns of the matrix
C     op( A ) and the number of rows of the matrix op( B ). K must
C     be at least  zero.
C     Unchanged on exit.
C
C     ALPHA  - DOUBLE PRECISION.
C     On entry, ALPHA specifies the scalar alpha.
C     Unchanged on exit.
C
C     A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
C     k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
C     Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
C     part of the array  A  must contain the matrix  A,  otherwise
C     the leading  k by m  part of the array  A  must contain  the
C     matrix A.
C     Unchanged on exit.
C
C     LDA    - INTEGER.
C     On entry, LDA specifies the first dimension of A as declared
C     in the calling (sub) program. When  TRANSA = 'N' or 'n' then
C     LDA must be at least  max( 1, m ), otherwise  LDA must be at
C     least  max( 1, k ).
C     Unchanged on exit.
C
C     B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
C     n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
C     Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
C     part of the array  B  must contain the matrix  B,  otherwise
C     the leading  n by k  part of the array  B  must contain  the
C     matrix B.
C     Unchanged on exit.
C
C     LDB    - INTEGER.
C     On entry, LDB specifies the first dimension of B as declared
C     in the calling (sub) program. When  TRANSB = 'N' or 'n' then
C     LDB must be at least  max( 1, k ), otherwise  LDB must be at
C     least  max( 1, n ).
C     Unchanged on exit.
C
C     BETA   - DOUBLE PRECISION.
C     On entry,  BETA  specifies the scalar  beta.  When  BETA  is
C     supplied as zero then C need not be set on input.
C     Unchanged on exit.
C
C     C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
C     Before entry, the leading  m by n  part of the array  C must
C     contain the matrix  C,  except when  beta  is zero, in which
C     case C need not be set on entry.
C     On exit, the array  C  is overwritten by the  m by n  matrix
C     ( alpha*op( A )*op( B ) + beta*C ).
C
C     LDC    - INTEGER.
C     On entry, LDC specifies the first dimension of C as declared
C     in  the  calling  (sub)  program.   LDC  must  be  at  least
C     max( 1, m ).
C     Unchanged on exit.
C
C
C     Level 3 Blas routine.
C
C     -- Written on 8-February-1989.
C     Jack Dongarra, Argonne National Laboratory.
C     Iain Duff, AERE Harwell.
C     Jeremy Du Croz, Numerical Algorithms Group Ltd.
C     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*     This code comes from a report entitled:
*     The IBM RISC System/6000 and Linear Algebra Operations, by
*     Jack J. Dongarra, Peter Mayes, and Giuseppe Radicati di Brozolo,
*     University of Tennessee Computer Science Tech Report: CS - 90 - 122.
C
C
      DOUBLE PRECISION A(LDA,*), B(LDB,*), C(LDC,*)
C     .. External Functions ..
      LOGICAL          LSAME
      EXTERNAL         LSAME
C     .. External Subroutines ..
      EXTERNAL         XERBLA
C     .. Intrinsic Functions ..
      INTRINSIC        MAX, MIN
C     .. Local Scalars ..
      DOUBLE PRECISION T11, T12, T21, T22
      INTEGER          I, IDEPTH, II, ILEN, INFO, ISPAN, J, JDEPTH, JJ,
     *                 JLEN, JSPAN, L, LL, LSPAN, NCOLA, NROWA, NROWB
      LOGICAL          NOTA, NOTB
C     .. Local Arrays ..
      DOUBLE PRECISION CH(KB,MB), CH1(KB), CH2(KB)
C     .. Executable Statements ..
C
C     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
C     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
C     and  columns of  A  and the  number of  rows  of  B  respectively.
C
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      IF (NOTA) THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF (NOTB) THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
C
C     Test the input parameters.
C
      INFO = 0
      IF (( .NOT. NOTA) .AND. ( .NOT. LSAME(TRANSA,'C'))
     *    .AND. ( .NOT. LSAME(TRANSA,'T'))) THEN
         INFO = 1
      ELSE IF (( .NOT. NOTB) .AND. ( .NOT. LSAME(TRANSB,'C'))
     *         .AND. ( .NOT. LSAME(TRANSB,'T'))) THEN
         INFO = 2
      ELSE IF (M.LT.0) THEN
         INFO = 3
      ELSE IF (N.LT.0) THEN
         INFO = 4
      ELSE IF (K.LT.0) THEN
         INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
         INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
         INFO = 13
      END IF
      IF (INFO.NE.0) THEN
         CALL XERBLA('DGEMM ',INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (((ALPHA.EQ.ZERO) .OR. (K.EQ.0))
     *     .AND. (BETA.EQ.ONE))) RETURN
      IF (BETA.EQ.ZERO) THEN
         DO 40 J = 1, N
            DO 20 I = 1, M
               C(I,J) = ZERO
   20       CONTINUE
   40    CONTINUE
      ELSE
         DO 80 J = 1, N
            DO 60 I = 1, M
               C(I,J) = BETA*C(I,J)
   60       CONTINUE
   80    CONTINUE
      END IF
C
C     And if  alpha.eq.zero.
C
      IF (ALPHA.EQ.ZERO) RETURN
C
C     Start the operations.
C
      IF (NOTB) THEN
         IF (NOTA) THEN
C
C           Form  C := C + alpha*A*B.
C
            DO 380 L = 1, K, KB
               LSPAN = MIN(KB,K-L+1)
               DO 360 I = 1, M, MB
                  IDEPTH = 2
                  ISPAN = MIN(MB,M-I+1)
                  ILEN = IDEPTH*(ISPAN/IDEPTH)
                  DO 120 II = I, I + ISPAN - 1
                     DO 100 LL = L, L + LSPAN - 1
                        CH(LL-L+1,II-I+1) = ALPHA*A(II,LL)
  100                CONTINUE
  120             CONTINUE
                  DO 340 J = 1, N, NB
                     JDEPTH = 2
                     JSPAN = MIN(NB,N-J+1)
                     JLEN = JDEPTH*(JSPAN/JDEPTH)
                     DO 220 JJ = J, J + JLEN - 1, JDEPTH
                        DO 160 II = I, I + ILEN - 1, IDEPTH
                           T11 = ZERO
                           T21 = ZERO
                           T12 = ZERO
                           T22 = ZERO
                           DO 140 LL = L, L + LSPAN - 1
                              T11 = T11 + CH(LL-L+1,II-I+1)*B(LL,JJ)
                              T21 = T21 + CH(LL-L+1,II-I+2)*B(LL,JJ)
                              T12 = T12 + CH(LL-L+1,II-I+1)*B(LL,JJ+1)
                              T22 = T22 + CH(LL-L+1,II-I+2)*B(LL,JJ+1)
  140                      CONTINUE
                           C(II,JJ) = C(II,JJ) + T11
                           C(II+1,JJ) = C(II+1,JJ) + T21
                           C(II,JJ+1) = C(II,JJ+1) + T12
                           C(II+1,JJ+1) = C(II+1,JJ+1) + T22
  160                   CONTINUE
                        IF (ILEN.LT.ISPAN) THEN
                           DO 200 II = I + ILEN, I + ISPAN - 1
                              T11 = ZERO
                              T12 = ZERO
                              DO 180 LL = L, L + LSPAN - 1
                                 T11 = T11 + CH(LL-L+1,II-I+1)*B(LL,JJ)
                                 T12 = T12 + CH(LL-L+1,II-I+1)*B(LL,
     *                                 JJ+1)
  180                         CONTINUE
                              C(II,JJ) = C(II,JJ) + T11
                              C(II,JJ+1) = C(II,JJ+1) + T12
  200                      CONTINUE
                        END IF
  220                CONTINUE
                     IF (JLEN.LT.JSPAN) THEN
                        DO 320 JJ = J + JLEN, J + JSPAN - 1
                           DO 260 II = I, I + ILEN - 1, IDEPTH
                              T11 = ZERO
                              T21 = ZERO
                              DO 240 LL = L, L + LSPAN - 1
                                 T11 = T11 + CH(LL-L+1,II-I+1)*B(LL,JJ)
                                 T21 = T21 + CH(LL-L+1,II-I+2)*B(LL,JJ)
  240                         CONTINUE
                              C(II,JJ) = C(II,JJ) + T11
                              C(II+1,JJ) = C(II+1,JJ) + T21
  260                      CONTINUE
                           IF (ILEN.LT.ISPAN) THEN
                              DO 300 II = I + ILEN, I + ISPAN - 1
                                 T11 = ZERO
                                 DO 280 LL = L, L + LSPAN - 1
                                    T11 = T11 + CH(LL-L+1,II-I+1)*B(LL,
     *                                    JJ)
  280                            CONTINUE
                                 C(II,JJ) = C(II,JJ) + T11
  300                         CONTINUE
                           END IF
  320                   CONTINUE
                     END IF
  340             CONTINUE
  360          CONTINUE
  380       CONTINUE
         ELSE
C
C           Form  C := C + alpha*A'*B
C
            DO 680 I = 1, M, MB
               IDEPTH = 2
               ISPAN = MIN(MB,M-I+1)
               ILEN = IDEPTH*(ISPAN/IDEPTH)
               DO 660 L = 1, K, KB
                  LSPAN = MIN(KB,K-L+1)
                  DO 420 II = I, I + ISPAN - 1
                     DO 400 LL = L, L + LSPAN - 1
                        CH(LL-L+1,II-I+1) = ALPHA*A(LL,II)
  400                CONTINUE
  420             CONTINUE
                  DO 640 J = 1, N, NB
                     JDEPTH = 2
                     JSPAN = MIN(NB,N-J+1)
                     JLEN = JDEPTH*(JSPAN/JDEPTH)
                     DO 520 JJ = J, J + JLEN - 1, JDEPTH
                        DO 460 II = I, I + ILEN - 1, IDEPTH
                           T11 = ZERO
                           T21 = ZERO
                           T12 = ZERO
                           T22 = ZERO
                           DO 440 LL = L, L + LSPAN - 1
                              T11 = T11 + CH(LL-L+1,II-I+1)*B(LL,JJ)
                              T21 = T21 + CH(LL-L+1,II-I+2)*B(LL,JJ)
                              T12 = T12 + CH(LL-L+1,II-I+1)*B(LL,JJ+1)
                              T22 = T22 + CH(LL-L+1,II-I+2)*B(LL,JJ+1)
  440                      CONTINUE
                           C(II,JJ) = C(II,JJ) + T11
                           C(II+1,JJ) = C(II+1,JJ) + T21
                           C(II,JJ+1) = C(II,JJ+1) + T12
                           C(II+1,JJ+1) = C(II+1,JJ+1) + T22
  460                   CONTINUE
                        IF (ILEN.LT.ISPAN) THEN
                           DO 500 II = I + ILEN, I + ISPAN - 1
                              T11 = ZERO
                              T12 = ZERO
                              DO 480 LL = L, L + LSPAN - 1
                                 T11 = T11 + CH(LL-L+1,II-I+1)*B(LL,JJ)
                                 T12 = T12 + CH(LL-L+1,II-I+1)*B(LL,
     *                                 JJ+1)
  480                         CONTINUE
                              C(II,JJ) = C(II,JJ) + T11
                              C(II,JJ+1) = C(II,JJ+1) + T12
  500                      CONTINUE
                        END IF
  520                CONTINUE
                     IF (JLEN.LT.JSPAN) THEN
                        DO 620 JJ = J + JLEN, J + JSPAN - 1
                           DO 560 II = I, I + ILEN - 1, IDEPTH
                              T11 = ZERO
                              T21 = ZERO
                              DO 540 LL = L, L + LSPAN - 1
                                 T11 = T11 + CH(LL-L+1,II-I+1)*B(LL,JJ)
                                 T21 = T21 + CH(LL-L+1,II-I+2)*B(LL,JJ)
  540                         CONTINUE
                              C(II,JJ) = C(II,JJ) + T11
                              C(II+1,JJ) = C(II+1,JJ) + T21
  560                      CONTINUE
                           IF (ILEN.LT.ISPAN) THEN
                              DO 600 II = I + ILEN, I + ISPAN - 1
                                 T11 = ZERO
                                 DO 580 LL = L, L + LSPAN - 1
                                    T11 = T11 + CH(LL-L+1,II-I+1)*B(LL,
     *                                    JJ)
  580                            CONTINUE
                                 C(II,JJ) = C(II,JJ) + T11
  600                         CONTINUE
                           END IF
  620                   CONTINUE
                     END IF
  640             CONTINUE
  660          CONTINUE
  680       CONTINUE
         END IF
      ELSE
         IF (NOTA) THEN
C
C           Form  C := C + alpha*A*B'
C
            DO 1000 J = 1, N, NB
               JDEPTH = 2
               JSPAN = MIN(NB,N-J+1)
               JLEN = JDEPTH*(JSPAN/JDEPTH)
               DO 980 L = 1, K, KB
                  LSPAN = MIN(KB,K-L+1)
                  DO 720 JJ = J, J + JSPAN - 1
                     DO 700 LL = L, L + LSPAN - 1
                        CH(LL-L+1,JJ-J+1) = ALPHA*B(JJ,LL)
  700                CONTINUE
  720             CONTINUE
                  DO 960 I = 1, M, MB
                     IDEPTH = 2
                     ISPAN = MIN(MB,M-I+1)
                     ILEN = IDEPTH*(ISPAN/IDEPTH)
                     DO 840 II = I, I + ILEN - 1, IDEPTH
                        DO 740 LL = L, L + LSPAN - 1
                           CH1(LL-L+1) = A(II,LL)
                           CH2(LL-L+1) = A(II+1,LL)
  740                   CONTINUE
                        DO 780 JJ = J, J + JLEN - 1, JDEPTH
                           T11 = ZERO
                           T21 = ZERO
                           T12 = ZERO
                           T22 = ZERO
                           DO 760 LL = L, L + LSPAN - 1
                              T11 = T11 + CH1(LL-L+1)*CH(LL-L+1,JJ-J+1)
                              T21 = T21 + CH2(LL-L+1)*CH(LL-L+1,JJ-J+1)
                              T12 = T12 + CH1(LL-L+1)*CH(LL-L+1,JJ-J+2)
                              T22 = T22 + CH2(LL-L+1)*CH(LL-L+1,JJ-J+2)
  760                      CONTINUE
                           C(II,JJ) = C(II,JJ) + T11
                           C(II+1,JJ) = C(II+1,JJ) + T21
                           C(II,JJ+1) = C(II,JJ+1) + T12
                           C(II+1,JJ+1) = C(II+1,JJ+1) + T22
  780                   CONTINUE
                        IF (JLEN.LT.JSPAN) THEN
                           DO 820 JJ = J + JLEN, J + JSPAN - 1
                              T11 = ZERO
                              T21 = ZERO
                              DO 800 LL = L, L + LSPAN - 1
                                 T11 = T11 + A(II,LL)*CH(LL-L+1,JJ-J+1)
                                 T21 = T21 + A(II+1,LL)*CH(LL-L+1,
     *                                 JJ-J+1)
  800                         CONTINUE
                              C(II,JJ) = C(II,JJ) + T11
                              C(II+1,JJ) = C(II+1,JJ) + T21
  820                      CONTINUE
                        END IF
  840                CONTINUE
                     IF (ILEN.LT.ISPAN) THEN
                        DO 940 II = I + ILEN, I + ISPAN - 1
                           DO 880 JJ = J, J + JLEN - 1, JDEPTH
                              T11 = ZERO
                              T12 = ZERO
                              DO 860 LL = L, L + LSPAN - 1
                                 T11 = T11 + A(II,LL)*CH(LL-L+1,JJ-J+1)
                                 T12 = T12 + A(II,LL)*CH(LL-L+1,JJ-J+2)
  860                         CONTINUE
                              C(II,JJ) = C(II,JJ) + T11
                              C(II,JJ+1) = C(II,JJ+1) + T12
  880                      CONTINUE
                           IF (JLEN.LT.JSPAN) THEN
                              DO 920 JJ = J + JLEN, J + JSPAN - 1
                                 T11 = ZERO
                                 DO 900 LL = L, L + LSPAN - 1
                                    T11 = T11 + A(II,LL)*CH(LL-L+1,
     *                                    JJ-J+1)
  900                            CONTINUE
                                 C(II,JJ) = C(II,JJ) + T11
  920                         CONTINUE
                           END IF
  940                   CONTINUE
                     END IF
  960             CONTINUE
  980          CONTINUE
 1000       CONTINUE
         ELSE
C
C           Form  C := C + alpha*A'*B'
C
            DO 1300 J = 1, N, NB
               JDEPTH = 2
               JSPAN = MIN(NB,N-J+1)
               JLEN = JDEPTH*(JSPAN/JDEPTH)
               DO 1280 L = 1, K, KB
                  LSPAN = MIN(KB,K-L+1)
                  DO 1040 JJ = J, J + JSPAN - 1
                     DO 1020 LL = L, L + LSPAN - 1
                        CH(LL-L+1,JJ-J+1) = ALPHA*B(JJ,LL)
 1020                CONTINUE
 1040             CONTINUE
                  DO 1260 I = 1, M, MB
                     IDEPTH = 2
                     ISPAN = MIN(MB,M-I+1)
                     ILEN = IDEPTH*(ISPAN/IDEPTH)
                     DO 1140 II = I, I + ILEN - 1, IDEPTH
                        DO 1080 JJ = J, J + JLEN - 1, JDEPTH
                           T11 = ZERO
                           T21 = ZERO
                           T12 = ZERO
                           T22 = ZERO
                           DO 1060 LL = L, L + LSPAN - 1
                              T11 = T11 + A(LL,II)*CH(LL-L+1,JJ-J+1)
                              T21 = T21 + A(LL,II+1)*CH(LL-L+1,JJ-J+1)
                              T12 = T12 + A(LL,II)*CH(LL-L+1,JJ-J+2)
                              T22 = T22 + A(LL,II+1)*CH(LL-L+1,JJ-J+2)
 1060                      CONTINUE
                           C(II,JJ) = C(II,JJ) + T11
                           C(II+1,JJ) = C(II+1,JJ) + T21
                           C(II,JJ+1) = C(II,JJ+1) + T12
                           C(II+1,JJ+1) = C(II+1,JJ+1) + T22
 1080                   CONTINUE
                        IF (JLEN.LT.JSPAN) THEN
                           DO 1120 JJ = J + JLEN, J + JSPAN - 1
                              T11 = ZERO
                              T21 = ZERO
                              DO 1100 LL = L, L + LSPAN - 1
                                 T11 = T11 + A(LL,II)*CH(LL-L+1,JJ-J+1)
                                 T21 = T21 + A(LL,II+1)*CH(LL-L+1,
     *                                 JJ-J+1)
 1100                         CONTINUE
                              C(II,JJ) = C(II,JJ) + T11
                              C(II+1,JJ) = C(II+1,JJ) + T21
 1120                      CONTINUE
                        END IF
 1140                CONTINUE
                     IF (ILEN.LT.ISPAN) THEN
                        DO 1240 II = I + ILEN, I + ISPAN - 1
                           DO 1180 JJ = J, J + JLEN - 1, JDEPTH
                              T11 = ZERO
                              T12 = ZERO
                              DO 1160 LL = L, L + LSPAN - 1
                                 T11 = T11 + A(LL,II)*CH(LL-L+1,JJ-J+1)
                                 T12 = T12 + A(LL,II)*CH(LL-L+1,JJ-J+2)
 1160                         CONTINUE
                              C(II,JJ) = C(II,JJ) + T11
                              C(II,JJ+1) = C(II,JJ+1) + T12
 1180                      CONTINUE
                           IF (JLEN.LT.JSPAN) THEN
                              DO 1220 JJ = J + JLEN, J + JSPAN - 1
                                 T11 = ZERO
                                 DO 1200 LL = L, L + LSPAN - 1
                                    T11 = T11 + A(LL,II)*CH(LL-L+1,
     *                                    JJ-J+1)
 1200                            CONTINUE
                                 C(II,JJ) = C(II,JJ) + T11
 1220                         CONTINUE
                           END IF
 1240                   CONTINUE
                     END IF
 1260             CONTINUE
 1280          CONTINUE
 1300       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of DGEMM .
C
      END

      SUBROUTINE TEST_ZGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A,
     $                        LDA, B, LDB, BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      COMPLEX*16         ALPHA, BETA
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - COMPLEX*16       array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - COMPLEX*16      .
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
*     .. Local Scalars ..
      LOGICAL            CONJA, CONJB, NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      COMPLEX*16         TEMP
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Executable Statements ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
*     B  respectively are to be  transposed but  not conjugated  and set
*     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
*     and the number of rows of  B  respectively.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      CONJA = LSAME( TRANSA, 'C' )
      CONJB = LSAME( TRANSB, 'C' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.CONJA                ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.CONJB                ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZGEMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE IF( CONJA )THEN
*
*           Form  C := alpha*conjg( A' )*B + beta*C.
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
            DO 150, J = 1, N
               DO 140, I = 1, M
                  TEMP = ZERO
                  DO 130, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  130             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  140          CONTINUE
  150       CONTINUE
         END IF
      ELSE IF( NOTA )THEN
         IF( CONJB )THEN
*
*           Form  C := alpha*A*conjg( B' ) + beta*C.
*
            DO 200, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 160, I = 1, M
                     C( I, J ) = ZERO
  160             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 170, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  170             CONTINUE
               END IF
               DO 190, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*CONJG( B( J, L ) )
                     DO 180, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  180                CONTINUE
                  END IF
  190          CONTINUE
  200       CONTINUE
         ELSE
*
*           Form  C := alpha*A*B'          + beta*C
*
            DO 250, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 210, I = 1, M
                     C( I, J ) = ZERO
  210             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 220, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  220             CONTINUE
               END IF
               DO 240, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 230, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  230                CONTINUE
                  END IF
  240          CONTINUE
  250       CONTINUE
         END IF
      ELSE IF( CONJA )THEN
         IF( CONJB )THEN
*
*           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
*
            DO 280, J = 1, N
               DO 270, I = 1, M
                  TEMP = ZERO
                  DO 260, L = 1, K
                     TEMP = TEMP +
     $                      CONJG( A( L, I ) )*CONJG( B( J, L ) )
  260             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  270          CONTINUE
  280       CONTINUE
         ELSE
*
*           Form  C := alpha*conjg( A' )*B' + beta*C
*
            DO 310, J = 1, N
               DO 300, I = 1, M
                  TEMP = ZERO
                  DO 290, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*B( J, L )
  290             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  300          CONTINUE
  310       CONTINUE
         END IF
      ELSE
         IF( CONJB )THEN
*
*           Form  C := alpha*A'*conjg( B' ) + beta*C
*
            DO 340, J = 1, N
               DO 330, I = 1, M
                  TEMP = ZERO
                  DO 320, L = 1, K
                     TEMP = TEMP + A( L, I )*CONJG( B( J, L ) )
  320             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  330          CONTINUE
  340       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
            DO 370, J = 1, N
               DO 360, I = 1, M
                  TEMP = ZERO
                  DO 350, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  350             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  360          CONTINUE
  370       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZGEMM .
*
      END
