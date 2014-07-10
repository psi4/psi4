#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine newcfg(natom1,natom2,itask,rcluster,rcut)
      implicit none 
      integer natom1, natom2, itask
      integer     n, nc, nx, ny, nz, ntot
 
      parameter ( nc = 50, n = 2*nc ** 3)
 
      double precision        lj(3,2)
      double precision        cell,cell2,mult,xadd,yadd,zadd
      double precision        multx,multy,multz,cmx,cmy,cmz
      double precision        x,y,z,dx,dy,dz,rcluster,r,rx,ry,rz,rmax
      double precision        volume, radius, pi, rcut
      double precision        coords(n, 3)
      integer     i, j, ix, iy, iz, iref, m, mp, mp1, mp4
      integer     one,two,nside, irad, ixcent, iycent, izcent
      integer     ixmin,ixmax,iymin,iymax,izmin,izmax
      integer     fill(nc,nc,nc,2), icnt
      character*32 filename
c
c   *******************************************************
c
      mult = 1.5d00
      lj(1,1) = 0.0d00
      lj(2,1) = 0.0d00
      lj(3,1) = 0.0d00
      lj(1,2) = mult/2.0d00
      lj(2,2) = mult/2.0d00
      lj(3,2) = mult/2.0d00
c
      one = 1
      two = 2
c
      if (itask.lt.10) then
        write(filename,100) itask
      else if (itask.ge.10.and.itask.lt.100) then
        write(filename,101) itask
      else if (itask.ge.100.and.itask.lt.1000) then
        write(filename,102) itask
      else if (itask.ge.1000.and.itask.lt.10000) then
        write(filename,103) itask
      endif
  100 format('md.cfg',i1)
  101 format('md.cfg',i2)
  102 format('md.cfg',i3)
  103 format('md.cfg',i4)

      open (unit=2,file=filename,status='unknown',form='formatted')
c
      nside = nint((0.5d00*dble(natom1+natom2))**(1.0d00/3.0d00))
      if (2*nside**3.lt.natom1+natom2) nside = nside+1
      pi = 4.0d00*atan(1.0d00)
      volume = 0.5d00*dble(natom2)*mult**3
      radius = 3.0d00*volume/(4.0d00*pi)
      radius = radius**(1.0d00/3.0d00)
      irad = int(radius/mult) + 1
      if (2*irad.ge.nside) nside = 2*irad+2
      do ix = 1, nside
        do iy = 1, nside
          do iz = 1, nside
            fill(ix,iy,iz,1) = 0
            fill(ix,iy,iz,2) = 0
          end do
        end do
      end do
      multx = mult * float(nside)
      multy = mult * float(nside)
      multz = mult * float(nside)
      cell = mult
      cell2 = 0.5d00 * cell
      dx = multx/float(nside)
      dy = multy/float(nside)
      dz = multz/float(nside)

      volume = 0.5d00*dble(natom2)*mult**3
      radius = 3.0d00*volume/(4.0d00*pi)
      radius = radius**(1.0d00/3.0d00)
      irad = int(radius/mult) + 1
      ixcent = nside/2
      iycent = nside/2
      izcent = nside/2
      ixmin = ixcent - irad
      ixmax = ixcent + irad
 300  icnt = 0
      do ix = ixmin, ixmax
        iy = int(sqrt(dble(irad)**2-dble(ix-ixcent)**2))
        iymin = iycent - iy
        iymax = iycent + iy
        do iy = iymin, iymax
          iz = int(sqrt(dble(irad)**2-dble(ix-ixcent)**2
     +                               -dble(iy-iycent)**2))
          izmin = izcent - iz
          izmax = izcent + iz
          do iz = izmin, izmax
            xadd = dx*float(ix-1) - multx/2.0d00
            yadd = dy*float(iy-1) - multy/2.0d00
            zadd = dz*float(iz-1) - multz/2.0d00
            if (fill(ix,iy,iz,1).eq.0.and.icnt.lt.natom2) then
              icnt = icnt + 1
              coords(natom1+icnt,1) = xadd + lj(1,1)
              coords(natom1+icnt,2) = yadd + lj(2,1)
              coords(natom1+icnt,3) = zadd + lj(3,1)
              fill(ix,iy,iz,1) = 1
            endif
            if (fill(ix,iy,iz,2).eq.0.and.icnt.lt.natom2) then
              icnt = icnt + 1
              coords(natom1+icnt,1) = xadd + lj(1,2)
              coords(natom1+icnt,2) = yadd + lj(2,2)
              coords(natom1+icnt,3) = zadd + lj(3,2)
              fill(ix,iy,iz,2) = 1
            endif
          end do
        end do
      end do
      if (icnt.lt.natom2) then
        irad = irad+1
        go to 300
      endif
c
c  Find center of mass of cluster
c
      rx = 0.0d00
      ry = 0.0d00
      rz = 0.0d00
      do icnt = 1, natom2
        rx = rx + coords(natom1+icnt,1)
        ry = ry + coords(natom1+icnt,2)
        rz = rz + coords(natom1+icnt,3)
      end do
      cmx = rx/dble(natom2)
      cmy = ry/dble(natom2)
      cmz = rz/dble(natom2)
      rmax = rmax + sqrt(dx**2+dy**2+dz**2)
c
c  Find radius of cluster
c
      rmax = 0.0d00
      do icnt = 1, natom2
        rx = coords(natom1+icnt,1) - cmx
        ry = coords(natom1+icnt,2) - cmy
        rz = coords(natom1+icnt,3) - cmz
        r = sqrt(rx**2 + ry**2 + rz**2)
        if (r.gt.rmax) rmax = r
      end do
      rmax = rmax + 0.1
c
      icnt = 0
      do ix = 1, nside
        do iy = 1, nside
          do iz = 1, nside
            xadd = dx*float(ix-1) - multx/2.0d00
            yadd = dy*float(iy-1) - multy/2.0d00
            zadd = dz*float(iz-1) - multz/2.0d00
            if (fill(ix,iy,iz,1).eq.0.and.icnt.lt.natom1) then
              icnt = icnt + 1
              coords(icnt,1) = xadd + lj(1,1)
              coords(icnt,2) = yadd + lj(2,1)
              coords(icnt,3) = zadd + lj(3,1)
              fill(ix,iy,iz,1) = 1
            endif
            if (fill(ix,iy,iz,2).eq.0.and.icnt.lt.natom1) then
              icnt = icnt + 1
              coords(icnt,1) = xadd + lj(1,2)
              coords(icnt,2) = yadd + lj(2,2)
              coords(icnt,3) = zadd + lj(3,2)
              fill(ix,iy,iz,2) = 1
            endif
          end do
        end do
      end do
c
      rcluster = rmax
c
c  If there are no atoms of type 1, then assume cluster is
c  isolated and set size to 100 LJ units
c
      if (natom1.eq.0) then
        multx = 25.0
        multy = 25.0
        multz = 25.0
      endif
      if (rcluster.lt.rcut) then
        rcluster = rcut
      endif
      write(2,2100) natom1+natom2
      write(2,2000) multx,multy,multz,rcluster
      do i = 1, natom1
        coords(i,1) = coords(i,1) - cmx
        coords(i,2) = coords(i,2) - cmy
        coords(i,3) = coords(i,3) - cmz
        write(2,200) one,(coords(i,j),j=1,3)
      end do
      do i = natom1+1, natom1+natom2
        coords(i,1) = coords(i,1) - cmx
        coords(i,2) = coords(i,2) - cmy
        coords(i,3) = coords(i,3) - cmz
        write(2,200) two,(coords(i,j),j=1,3)
      end do
c
 200  format (i5,3f12.5)
2000  format (4f16.8)
2100  format (i8)
      close(2)
 
      return
      end
