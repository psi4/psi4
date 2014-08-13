      subroutine mxv(a,ncol,b,nrow,c)
C$Id: mxv_daxpy1.f,v 1.2 1995-02-02 23:24:19 d3g681 Exp $
      implicit double precision (a-h, o-z)
      dimension a(ncol, nrow), b(nrow), c(ncol)
      parameter (nchunk = 800)
c
c     matrix vector product stripmined to optimize cache usage
c     when inner loop is replaced with a daxpy that uses pipelined
c     loads for a to avoid writing over c in the cache.
c
      do 10 ilo = 1, ncol, nchunk
         ihi = min(ncol, ilo+nchunk-1)
     ndo = ihi - ilo + 1
         do 20 i = ilo, ihi
            c(i) = 0.0d0
 20      continue
         do 30 j = 1, nrow
c           do 40 i = ilo, ihi
c              c(i) = c(i) + a(i,j)*b(j)
c40         continue
            call daxpy1(ndo, b(j), a(ilo, j), c(ilo))
 30      continue
 10   continue
c     
      end
