#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine factor(p,idx,idy,idz)
      implicit none
      integer i,j,p,idx,idy,idz,it
      integer ip,ifac,pmax,prime(1000)
      integer fac(1000)
c
      i = 1
      ip = p
c
c    factor p completely
c    first, find all prime numbers less than or equal to p
c
      pmax = 0
      do i = 2, p
        do j = 1, pmax
          if (mod(i,prime(j)).eq.0) go to 100
        end do
        pmax = pmax + 1
        prime(pmax) = i
  100   continue
      end do
c
c    find all prime factors of p
c
      ifac = 0
      do i = 1, pmax
  200   if (mod(ip,prime(i)).eq.0) then
          ifac = ifac + 1
          fac(ifac) = prime(i)
          ip = ip/prime(i)
          go to 200
        endif
      end do
c
c    determine three factors of p of approximately the
c    same size
c
      idx = 1
      idy = 1
      idz = 1
      do i = ifac, 1, -1
        if (idx.le.idy.and.idx.le.idz) then
          idx = fac(i)*idx
        elseif (idy.le.idx.and.idy.le.idz) then
          idy = fac(i)*idy
        elseif (idz.le.idx.and.idz.le.idy) then
          idz = fac(i)*idz
        endif
      end do
c      it = idy
c      idy = idx
c      idx = it
c
c      it = idx
c      idx = idz
c      idz = it
c
c      it = idy
c      idy = idz
c      idz = it
      return
      end
c
      subroutine i_proc_to_xyz(p,ix,iy,iz,idx,idy,idz)
      implicit none
      integer p,ix,iy,iz,ip,it
      integer idx,idy,idz
c
      ip = p
c
      it = mod(ip,idx)
      ix = it
      ip = (ip - it)/idx
      it = mod(ip,idy)
      iy = it
      ip = (ip - it)/idy
      it = mod(ip,idz)
      iz = it
c
      return
      end
c
      subroutine i_xyz_to_proc(p,ix,iy,iz,idx,idy,idz)
      implicit none
      integer p,ix,iy,iz,idx,idy,idz
c
      p = ix + idx*iy + idx*idy*iz
c
      return
      end
