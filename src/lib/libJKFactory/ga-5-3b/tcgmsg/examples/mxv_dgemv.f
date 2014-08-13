      subroutine mxv(a,ncol,b,nrow,c)
C$Id: mxv_dgemv.f,v 1.2 1995-02-02 23:24:20 d3g681 Exp $
      implicit double precision (a-h, o-z)
      double precision a(ncol, nrow), b(nrow), c(ncol)
      parameter (ilen=500, jlen=60)
c
      call dgemv('n', ncol, nrow, 1.0d0, a, ncol, b, 1, 0.0d0, c, 1)
c$$$      do 10 i = 1, ncol
c$$$         c(i) = 0.0d0
c$$$ 10   continue
c$$$c     
c$$$      do 40 jlo = 1, nrow, jlen
c$$$         jhi = min(jlo+jlen-1, nrow)
c$$$         do 30 ilo = 1, ncol, ilen
c$$$            ihi = min(ilo+ilen-1, ncol)
c$$$            ndo = ihi - ilo + 1
c$$$            do 20 j = jlo, jhi
c$$$               call daxpy2(ndo, b(j), a(ilo,j), c(ilo))
c$$$ 20         continue
c$$$ 30      continue
c$$$ 40   continue
c
      end
