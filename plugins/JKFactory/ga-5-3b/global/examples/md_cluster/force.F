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
      subroutine force
#include "common.fh"
c
      double precision tbeg,wraptime
      integer i,j
      integer inode
c
      tbeg = wraptime()
c
c   clear forces
c
      do 100 j = 1, 3
        do 200 i = 1, antot
          ra(i,j,4) = 0.0d00
  200   continue
  100 continue
c
      do 300 i = 1, antot
        asev(i) = 0.0d00
  300 continue
c
      esvr = 0.0d00
      esvd = 0.0d00
      esvoh = 0.0d00
      esva = 0.0d00
      presf = 0.0d00
      presfx = 0.0d00
      presfy = 0.0d00
      presfz = 0.0d00
c
      call gather
c
c   clear buffers
c
      do i = 1, btot
        xfrc(i) = 0.0d00
        yfrc(i) = 0.0d00
        zfrc(i) = 0.0d00
      end do
c
c   assign atoms to individual molecules
c
      call pairs
      call scatter
c
      nrg(6) =esvr + esvd + esvoh + esva
      nrg(13) = esvr
      nrg(14) = esvd
      nrg(17) = esvoh
      nrg(21) = esva
      nrg(10) = presf
      nrg(18) = presfx
      nrg(19) = presfy
      nrg(20) = presfz
      presf = presf / (xbox * ybox * zbox)
      nrg(15) = presf/ 3.0d00
c
      tmstat(2) = tmstat(2) + wraptime() - tbeg
c
      return
      end
