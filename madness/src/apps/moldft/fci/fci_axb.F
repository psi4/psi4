*
* $Id: fci_axb.F,v 1.2 1999/07/28 00:23:39 d3e129 Exp $
*
      subroutine fci_axb(ni, nj, nk, a, b, c)
      implicit none
c
      integer ni, nj, nk
      double precision c(ni, nj), a(ni, nk), b(nk, nj)
c
c     c(i,j) = sum(k) a(i,k)*b(k,j)
c
      integer ilo, ihi, jlo, jhi, klo, khi
      integer i, j, k, ndoi
      integer ichunk, jchunk, kchunk
      parameter (ichunk=192, jchunk=32, kchunk=32) ! IBM 755CX
*      parameter (ichunk=100, jchunk=8, kchunk=8)
c
      do jlo = 1,nj,jchunk
         jhi = min(nj,jlo+jchunk-1)
         do ilo = 1,ni,ichunk
            ihi = min(ni,ilo+ichunk-1)
            ndoi = ihi - ilo + 1
            do j = jlo, jhi
               do i = ilo, ihi
                  c(i,j) = 0.0d0
               enddo
            enddo
            do klo = 1,nk,kchunk
               khi = min(nk,klo+kchunk-1)
c
               do j = jlo, jhi
                  do k = klo, khi
                     if (b(k,j) .ne. 0.0d0) then
#ifdef LINUX
                        call daxpy(ndoi,b(k,j),a(ilo,k),1,c(ilo,j),1)
#else
                        do i = ilo, ihi
                           c(i,j) = c(i,j) + a(i,k)*b(k,j)
                        enddo
#endif
                     endif
                  enddo
               enddo
c$$$c OK
c$$$               do k = klo, khi
c$$$                  do j = jlo, jhi
c$$$                     if (b(k,j) .ne. 0.0d0) then
c$$$                        do i = ilo, ihi
c$$$                           c(i,j) = c(i,j) + a(i,k)*b(k,j)
c$$$                        enddo
c$$$                     endif
c$$$                  enddo
c$$$               enddo
c
            enddo
         enddo
      enddo
c
      end

         
