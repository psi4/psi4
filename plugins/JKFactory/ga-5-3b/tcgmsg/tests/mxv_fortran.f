      subroutine mxv(a,ncol,b,nrow,c)
      implicit double precision (a-h, o-z)
      dimension a(ncol, nrow), b(nrow), c(ncol)
      parameter (nchunk = 127)
c
c     matrix vector product stripmined to optimize cache usage
c     when inner loop is replaced with a daxpy that uses pipelined
c     loads for a to avoid writing over c in the cache.
c
      do 10 ilo = 1, ncol, nchunk
         ihi = min(ncol, ilo+nchunk-1)
         do 20 i = ilo, ihi
            c(i) = 0.0d0
 20      continue
         do 30 j = 1, nrow
            do 40 i = ilo, ihi
               c(i) = c(i) + a(i,j)*b(j)
 40         continue
 30      continue
 10   continue
c     
      end
