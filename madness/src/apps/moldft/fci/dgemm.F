*
* $Id: dgemm.F,v 1.2 1999/07/28 00:23:38 d3e129 Exp $
*
      SUBROUTINE DGEMM (TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
     $     BETA,C,LDC)
*     .. Scalar Arguments ..
C     $Id: dgemm.F,v 1.2 1999/07/28 00:23:38 d3e129 Exp $
      CHARACTER*1        TRANSA,TRANSB
      INTEGER            M,N,K,LDA,LDB,LDC
      DOUBLE PRECISION   ALPHA,BETA
*     .. Array Arguments ..
      DOUBLE PRECISION   A(LDA,*),B(LDB,*),C(LDC,*)
*     
*     !!!!!!!!!!!!!!!!!!!!!!!
*     This version of DGEMM has been hacked by RJH to use blocking
*     (loop unrolling removed)
*     The original version is in dgemm.orig.f
*     !!!!!!!!!!!!!!!!!!!!!!!
*     ..
*     
*     Purpose
*     =======
*     
*     DGEMM  performs one of the matrix-matrix operations
*     
*     C := alpha*op(A)*op(B) + beta*C,
*     
*     where  op(X) is one of
*     
*     op(X) = X   or   op(X) = X',
*     
*     alpha and beta are scalars,and A,B and C are matrices,with op(A)
*     an m by k matrix, op(B)  a  k by n matrix and  C an m by n matrix.
*     
*     Parameters
*     ==========
*     
*     TRANSA - CHARACTER*1.
*     On entry,TRANSA specifies the form of op(A) to be used in
*     the matrix multiplication as follows:
*     
*     TRANSA = 'N' or 'n', op(A) = A.
*     
*     TRANSA = 'T' or 't', op(A) = A'.
*     
*     TRANSA = 'C' or 'c', op(A) = A'.
*     
*     Unchanged on exit.
*     
*     TRANSB - CHARACTER*1.
*     On entry,TRANSB specifies the form of op(B) to be used in
*     the matrix multiplication as follows:
*     
*     TRANSB = 'N' or 'n', op(B) = B.
*     
*     TRANSB = 'T' or 't', op(B) = B'.
*     
*     TRANSB = 'C' or 'c', op(B) = B'.
*     
*     Unchanged on exit.
*     
*     M      - INTEGER.
*     On entry, M  specifies  the number  of rows  of the  matrix
*     op(A)  and of the  matrix  C.  M  must  be at least  zero.
*     Unchanged on exit.
*     
*     N      - INTEGER.
*     On entry, N  specifies the number  of columns of the matrix
*     op(B) and the number of columns of the matrix C. N must be
*     at least zero.
*     Unchanged on exit.
*     
*     K      - INTEGER.
*     On entry, K  specifies  the number of columns of the matrix
*     op(A) and the number of rows of the matrix op(B). K must
*     be at least  zero.
*     Unchanged on exit.
*     
*     ALPHA  - DOUBLE PRECISION.
*     On entry,ALPHA specifies the scalar alpha.
*     Unchanged on exit.
*     
*     A      - DOUBLE PRECISION array of DIMENSION (LDA,ka),where ka is
*     k  when  TRANSA = 'N' or 'n', and is  m  otherwise.
*     Before entry with  TRANSA = 'N' or 'n', the leading  m by k
*     part of the array  A  must contain the matrix  A, otherwise
*     the leading  k by m  part of the array  A  must contain  the
*     matrix A.
*     Unchanged on exit.
*     
*     LDA    - INTEGER.
*     On entry,LDA specifies the first dimension of A as declared
*     in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*     LDA must be at least  max(1,m),otherwise  LDA must be at
*     least  max(1,k).
*     Unchanged on exit.
*     
*     B      - DOUBLE PRECISION array of DIMENSION (LDB,kb),where kb is
*     n  when  TRANSB = 'N' or 'n', and is  k  otherwise.
*     Before entry with  TRANSB = 'N' or 'n', the leading  k by n
*     part of the array  B  must contain the matrix  B, otherwise
*     the leading  n by k  part of the array  B  must contain  the
*     matrix B.
*     Unchanged on exit.
*     
*     LDB    - INTEGER.
*     On entry,LDB specifies the first dimension of B as declared
*     in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*     LDB must be at least  max(1,k),otherwise  LDB must be at
*     least  max(1,n).
*     Unchanged on exit.
*     
*     BETA   - DOUBLE PRECISION.
*     On entry, BETA  specifies the scalar  beta.  When  BETA  is
*     supplied as zero then C need not be set on input.
*     Unchanged on exit.
*     
*     C      - DOUBLE PRECISION array of DIMENSION (LDC,n).
*     Before entry,the leading  m by n  part of the array  C must
*     contain the matrix  C, except when  beta  is zero,in which
*     case C need not be set on entry.
*     On exit,the array  C  is overwritten by the  m by n  matrix
*     (alpha*op(A)*op(B) + beta*C).
*     
*     LDC    - INTEGER.
*     On entry,LDC specifies the first dimension of C as declared
*     in  the  calling  (sub)  program.   LDC  must  be  at  least
*     max(1,m).
*     Unchanged on exit.
*     
*     
*     Level 3 Blas routine.
*     
*     -- Written on 8-February-1989.
*     Jack Dongarra,Argonne National Laboratory.
*     Iain Duff,AERE Harwell.
*     Jeremy Du Croz,Numerical Algorithms Group Ltd.
*     Sven Hammarling,Numerical Algorithms Group Ltd.
*     
*     
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
*      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            NOTA,NOTB
      INTEGER            I,INFO,J,L,NCOLA,NROWA,NROWB
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         ,ZERO
      PARAMETER        (ONE = 1.0D+0,ZERO = 0.0D+0)
c
      double precision work(0:1023) ! RJH To handle TT case
      double precision ddot
      external ddot
*     ..
*     .. Executable Statements ..
*     
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA,NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*     
      NOTA  = LSAME(TRANSA,'N')
      NOTB  = LSAME(TRANSB,'N')
      IF(NOTA)THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF(NOTB)THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*     
*     Test the input parameters.
*     
      INFO = 0
      IF(     (.NOT.NOTA                ).AND.
     $     (.NOT.LSAME(TRANSA,'C')).AND.
     $     (.NOT.LSAME(TRANSA,'T'))     )THEN
         INFO = 1
      ELSE IF((.NOT.NOTB                ).AND.
     $        (.NOT.LSAME(TRANSB,'C')).AND.
     $        (.NOT.LSAME(TRANSB,'T'))     )THEN
         INFO = 2
      ELSE IF(M  .LT.0              )THEN
         INFO = 3
      ELSE IF(N  .LT.0              )THEN
         INFO = 4
      ELSE IF(K  .LT.0              )THEN
         INFO = 5
      ELSE IF(LDA.LT.MAX(1,NROWA))THEN
         INFO = 8
      ELSE IF(LDB.LT.MAX(1,NROWB))THEN
         INFO = 10
      ELSE IF(LDC.LT.MAX(1,M    ))THEN
         INFO = 13
      END IF
      IF(INFO.NE.0)THEN
*         CALL XERBLA('DGEMM ',INFO)
         call fci_err('dgemm: ', info)
         RETURN
      END IF
*     
*     Quick return if possible.
*     
      IF((M.EQ.0).OR.(N.EQ.0).OR.
     $     (((ALPHA.EQ.ZERO).OR.(K.EQ.0)).AND.(BETA.EQ.ONE)))
     $     RETURN
*     
*     And if  alpha.eq.zero.
*     
      IF(ALPHA.EQ.ZERO)THEN
         IF(BETA.EQ.ZERO)THEN
            DO 20,J = 1,N
               DO 10,I = 1,M
                  C(I,J) = ZERO
 10            CONTINUE
 20         CONTINUE
         ELSE
            DO 40,J = 1,N
               DO 30,I = 1,M
                  C(I,J) = BETA*C(I,J)
 30            CONTINUE
 40         CONTINUE
         END IF
         RETURN
      END IF
*     
*     Start the operations.
*     
c     
c     rjh ... determine blocking parameters from the cachesize
c     .       and allowing inner loop to be 3x longer for 
c     .       better vectorization
c
*      icachesize = 131072*0.5   ! 1/2 MB cache in doubles
      icachesize = 131072*0.25  ! 1/4 MB cache in doubles
*      icachesize = 2048         ! 16KB cache in doubles
c
      jchunk = sqrt(icachesize/7.0d0)
      if (nota) then
         lchunk = jchunk
         ichunk = jchunk*3      ! For max vector length
      else
         if (notb) then
            lchunk = jchunk*3   ! For max vectors length
            ichunk = jchunk
         else
            jchunk = min(1024,jchunk)
            lchunk = jchunk
            ichunk = 3*jchunk   ! To minimize cache thrash loading B
         endif
      endif
c
*      write(6,111) ichunk, jchunk, lchunk
* 111  format(' ijlchunk ', 3i6)
c

c
      do jlo = 1,n,jchunk         ! blocking by rjh
         jhi = min(jlo+jchunk-1,n)
         do ilo = 1,m,ichunk
            ihi = min(ilo+ichunk-1,m)
            do j = jlo,jhi
               IF(BETA.EQ.ZERO)THEN
                  DO I = ilo,ihi
                     C(I,J) = ZERO
                  enddo
               ELSE IF(BETA.NE.ONE)THEN
                  DO I = ilo,ihi
                     C(I,J) = BETA*C(I,J)
                  enddo
               END IF
            enddo
            do llo = 1,k,lchunk
               lhi = min(llo+lchunk-1,k)
c     
               IF(NOTB)THEN
                  IF(NOTA)THEN
*     
*     Form  C := alpha*A*B + beta*C.
*     
                     DO 90,J = jlo,jhi
                        DO 80,L = llo,lhi
                           IF(B(L,J).NE.ZERO)THEN
                              TEMP = ALPHA*B(L,J)
#ifdef LINUX
                              call daxpy(ihi-ilo+1,temp,
     $                             a(ilo,l),1,c(ilo,j),1)
#else
                              do i = ilo, ihi
                                 C(I  ,J) = C(I  ,J) + TEMP*A(I  ,L)
                              enddo
#endif
                           END IF
 80                     CONTINUE
 90                  CONTINUE
                  ELSE
*     
*     Form  C := alpha*A'*B + beta*C
*     
                     DO 120,J = jlo,jhi
                        DO 110,I = ilo,ihi
#ifdef LINUX
                           temp = ddot(lhi-llo+1,a(llo,i),1,b(llo,j),1)
#else
                           TEMP = ZERO
                           DO L = llo,lhi
                              TEMP = TEMP + A(L,I)*B(L,J)
                           enddo
#endif
                           C(I,J) = ALPHA*TEMP + C(I,J)
 110                    CONTINUE
 120                 CONTINUE
                  END IF
               ELSE
                  IF(NOTA)THEN
*     
*     Form  C := alpha*A*B' + beta*C
*     
                     DO 170,J = jlo,jhi
                        DO 160,L = llo,lhi
                           IF(B(J,L).NE.ZERO)THEN
                              TEMP = ALPHA*B(J,L)
                              DO I = ilo,ihi
                                 C(I,J) = C(I,J) + TEMP*A(I,L)
                              enddo
                           END IF
 160                    CONTINUE
 170                 CONTINUE
                  ELSE
*     
*     Form  C := alpha*A'*B' + beta*C
*     
                     DO 200,J = jlo,jhi
                        do l = llo,lhi
                           work(l-llo) = b(j,l)
                        enddo
                        DO 190,I = ilo,ihi
                           TEMP = ZERO
                           DO L = llo,lhi
*                              TEMP = TEMP + A(L,I)*B(J,L)
                              TEMP = TEMP + A(L,I)*work(l-llo)
                           enddo
                           C(I,J) = ALPHA*TEMP + C(I,J)
 190                    CONTINUE
 200                 CONTINUE
                  END IF
               END IF
            enddo
         enddo
      enddo
*     
      RETURN
*     
*     End of DGEMM .
*     
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
* $Id: dgemm.F,v 1.2 1999/07/28 00:23:38 d3e129 Exp $
*
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     January 31, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
