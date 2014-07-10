#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c
c                                   NOTICE
c
c   This software is being made available for internal testing and
c   evaluation purposes only. This software is a pre-release test version
c   which has not yet been authenticated and cleared for publication. Adherence
c   to this notice may be necessary for the author, Battelle Memorial
c   Institute, to successfully assert copyright in and commercialize this
c   software. This software is not intended for duplication or distribution
c   to third parties without the permission of the Manager of Software
c   Products at Pacific Northwest Laboratory, Richland, Washington,  99352.
c
      subroutine kin
#include "common.fh"
c
      integer i,j
c
c   This subroutine calculates the kinetic energy of the current
c   configuration. It assumes that the center of mass motion in
c   the xy plane has already been subtracted out.
c
      nrg(4) = 0.0d00
      nrg(8) = 0.0d00
      do 200 j = 1, 3
        do 100 i = 1, antot
          nrg(4) = nrg(4) + 0.5d00 * ra(i,j,5)**2 / mass(i)
c          nrg(4) = nrg(4) + 0.5d00 * mass(i) * ra(i,j,2)**2
  100   continue
  200 continue
c
      return
      end
c
      subroutine com
#include "common.fh"
c
      double precision px,py,pz,vx,vy,vz,tcmass
      double precision comm(4)
      integer i
c
c   This subroutine calculates the contributions to the center of mass
c   and the center of mass momentum from each processor
c
c   Calculate center of mass motion momentum
c
      px = 0.0d00
      py = 0.0d00
      pz = 0.0d00
      tcmass = 0.0d00
      do 100 i = 1, antot
        px = px + ra(i,1,5)
        py = py + ra(i,2,5)
        pz = pz + ra(i,3,5)
        tcmass = tcmass + mass(i)
  100 continue
c
      comm(1) = px
      comm(2) = py
      comm(3) = pz
      comm(4) = tcmass
      call ga_dgop(7,comm,4,'+')
      px = comm(1)
      py = comm(2)
      pz = comm(3)
      tcmass = comm(4)
c
      nrg(22) = px
      nrg(23) = py
      nrg(24) = pz
c
      nrg(25) = tcmass
c
      return
      end
c
      subroutine com_rmv
#include "common.fh"
c
      double precision px,py,pz,vx,vy,vz,tcmass
      integer i
c
c   This subroutine subtracts out the center of mass motion
c   Calculate center of mass motion momentum
c
       px = nrg(22)
       py = nrg(23)
       pz = nrg(24)
       tcmass = nrg(25)
c
c   Calculate center of mass velocity
c
       vx = px / tcmass
       vy = py / tcmass
       vz = pz / tcmass
c
c   Subtract out center of mass velocity
c
       do 200 i = 1, antot
         ra(i,1,2) = ra(i,1,2) - vx
         ra(i,2,2) = ra(i,2,2) - vy
         ra(i,3,2) = ra(i,3,2) - vz
         ra(i,1,5) = ra(i,1,5) - mass(i) * vx
         ra(i,2,5) = ra(i,2,5) - mass(i) * vy
         ra(i,3,5) = ra(i,3,5) - mass(i) * vz
  200 continue
c
      return
      end
c
      subroutine fixper
#include "common.fh"
c
      integer i
c
c   This subroutine fixes up the periodic boundary conditions
c   in the xy plane after the particle positions have been updated
c
      do 100 i = 1, antot
        ra(i,1,1) = ra(i,1,6) - xbox * anint(ra(i,1,6)/xbox)
        if (ra(i,1,1).eq.xbox2) ra(i,1,1) = -xbox2
        ra(i,2,1) = ra(i,2,6) - ybox * anint(ra(i,2,6)/ybox)
        if (ra(i,2,1).eq.ybox2) ra(i,2,1) = -ybox2
        ra(i,3,1) = ra(i,3,6) - zbox * anint(ra(i,3,6)/zbox)
        if (ra(i,3,1).eq.zbox2) ra(i,3,1) = -zbox2
  100 continue
      return
      end
