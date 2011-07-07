*
* $Id: timer.f,v 1.2 1999/07/28 00:23:43 d3e129 Exp $
*
      implicit none
      integer maxn
      parameter (maxn = 100000)
      double precision a(maxn), b(maxn), c(maxn), s
c
      integer i, j, n, nloop, loop
      double precision linux_cputime, start, usedf, usedb, ratef, rateb
      external linux_cputime
c
      do i = 1, n
         a(i) = 2.0d0
         b(i) = 1.0d0
         c(i) = 3.0d0
      enddo
c
      n = 1
      s = 4.5d0
c
 10   nloop = 10000000/n + 1
      start = linux_cputime()
      do loop = 1, nloop
         do i = 1, n
            a(i) = a(i) + s*b(i)
         enddo
      enddo
      usedf = linux_cputime() - start
      ratef = 1d-6 * n * nloop / usedf
      start = linux_cputime()
      do loop = 1, nloop
         call daxpy(n, s, b, 1, a, 1)
      enddo
      usedb = linux_cputime() - start
      rateb = 1d-6 * n * nloop / usedb
      write(6,1) n, usedf, ratef, usedb, rateb
 1    format(1x,i8,2f10.3,2x,2f10.3)
      n = n * 2
      if (n .lt. maxn) goto 10
c
      end


