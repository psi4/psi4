#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine input
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "mpif.h"
c.......................................................................
c Input configuration from an XYZ format file call be.inpt and set up
c initial data structures. Atomic numbers from XYZ file are ignored and
c an atomic number of 4 (Beryllium) is used instead.
c.......................................................................
c
c   initialize variables
c
      natom = 0
      do i = 1, maxatom
        ax(i) = 0.0d00
        ay(i) = 0.0d00
        az(i) = 0.0d00
      end do
c
      if (ga_nodeid().eq.0) then
        open(2,file='be.inpt',status='old')
        read(2,*) natom
        read(2,*)
c
c Read in coordinates
c
        do i = 1, natom
          read(2,*) j,ax(i),ay(i),az(i)
        end do
        close(2)
      endif
      call ga_igop(1,natom,1,'+')
      call ga_dgop(2,ax,natom,'+')
      call ga_dgop(3,ay,natom,'+')
      call ga_dgop(4,az,natom,'+')
c
c Set up s-function centers and nuclear charges
c
      ifcnt = 1
      do i = 1, natom
        q(i) = 4.0d00
c
        expnt(ifcnt) = 1741.0d00
        expnt(ifcnt+1) = 262.1d00
        expnt(ifcnt+2) = 60.33d00
        expnt(ifcnt+3) = 17.62d00
        expnt(ifcnt+4) = 5.933d00
        expnt(ifcnt+5) = 2.185d00
        expnt(ifcnt+6) = 0.859d00
        expnt(ifcnt+7) = 0.1806d00
        expnt(ifcnt+8) = 0.05835d00
        expnt(ifcnt+9) = 0.3d00
        expnt(ifcnt+10) = 0.3d00
        expnt(ifcnt+11) = 0.3d00
        expnt(ifcnt+12) = 0.3d00
        expnt(ifcnt+13) = 0.3d00
        expnt(ifcnt+14) = 0.3d00
c
        do j = 1, 15
          x(ifcnt) = ax(i)
          y(ifcnt) = ay(i)
          z(ifcnt) = az(i)
          if (j.eq.10) then
            x(ifcnt) = x(ifcnt) + 1.6d00
          endif
          if (j.eq.11) then
            x(ifcnt) = x(ifcnt) - 1.6d00
          endif
          if (j.eq.12) then
            y(ifcnt) = y(ifcnt) + 1.6d00
          endif
          if (j.eq.13) then
            y(ifcnt) = y(ifcnt) - 1.6d00
          endif
          if (j.eq.14) then
            z(ifcnt) = z(ifcnt) + 1.6d00
          endif
          if (j.eq.15) then
            z(ifcnt) = z(ifcnt) - 1.6d00
          endif
          ifcnt = ifcnt + 1
        end do
      end do
c
c  evaluate repulsion energy
c
      enrep = 0.0d00
      do i = 1, natom
        do j = i+1, natom
          r = sqrt((ax(i)-ax(j))**2 + (ay(i)-ay(j))**2
     +      + (az(i)-az(j))**2)
          enrep = enrep + q(i)*q(j)/r
        end do
      end do
      nocc = 2*natom
      nbfn =  15*natom
      nnbfn = nbfn*(nbfn+1)/2
      return
      end
