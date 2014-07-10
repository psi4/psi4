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
      function wraptime()
      double precision wraptime
#if HAVE_ETIME      
      real*4 etime,tarray(2)
      external etime
c
c   This function returns the UNIX clock time in seconds.
c   It is designed to minimize the effort required to 
c   convert the timing statistics in going from one
c   machine to another.
c
c      SGI
c
      wraptime = dble(etime(tarray))
#else
c
c      NO TIME STATISTICS
c
      wraptime = 0.0d00
#endif
c
      return
      end
c
      subroutine tmout
#include "common.fh"
      double precision tmax(MAXTIM),tmin(MAXTIM)
      integer i,pnum
c
c  This routine writes out the timekeeping statistics to
c  the output file 'md.out'
c
      do i = 1, MAXTIM
        tmin(i) = tmstat(i)
        tmax(i) = tmstat(i)
      end do
      call ga_dgop(5,tmstat,MAXTIM,'+')
      call ga_dgop(5,tmax,MAXTIM,'max')
      call ga_dgop(5,tmin,MAXTIM,'min')
      pnum = ga_nnodes()
      do i = 1, MAXTIM
        tmstat(i) = tmstat(i)/dble(pnum)
      end do
      if (ga_nodeid().ne.0) return
c
      if (l_stdio) then
        write(6,100) tmstat(1),tmin(1),tmax(1)
        write(6,200) tmstat(2),tmin(2),tmax(2)
        write(6,300) tmstat(3),tmin(3),tmax(3)
        write(6,400) tmstat(4),tmin(4),tmax(4)
        write(6,500) tmstat(5),tmin(5),tmax(5)
        write(6,2300) tmstat(23),tmin(23),tmax(23)
        write(6,550) tmstat(6)+tmstat(7)+tmstat(8)
        write(6,600) tmstat(6),tmin(6),tmax(6)
        write(6,700) tmstat(7),tmin(7),tmax(7)
        write(6,750)
        write(6,800) tmstat(8),tmin(8),tmax(8)
        write(6,2800) tmstat(28),tmin(28),tmax(28)
        write(6,1400) tmstat(14),tmin(14),tmax(14)
        write(6,850) tmstat(9)+tmstat(10)
        write(6,900) tmstat(9),tmin(9),tmax(9)
        write(6,1000) tmstat(10),tmin(10),tmax(10)
        write(6,1050) tmstat(11)+tmstat(12)+tmstat(13)
        write(6,1100) tmstat(11),tmin(11),tmax(11)
        write(6,1200) tmstat(12),tmin(12),tmax(12)
        write(6,1300) tmstat(13),tmin(13),tmax(13)
        write(6,1500) tmstat(15),tmin(15),tmax(15)
        write(6,1600) tmstat(16),tmin(16),tmax(16)
      endif
c
  100 format('Total elapsed time for calculation     ',3f8.1,' (sec)')
  200 format('Total time in force calculation        ',3f8.1,' (sec)')
  300 format('          Pairwise interactions        ',3f8.1,' (sec)')
  400 format('          Bending interactions         ',3f8.1,' (sec)')
  500 format('          Stretching interactions      ',3f8.1,' (sec)')
  550 format('Total time exchanging particle lists   ',f8.1,' (sec)')
  600 format('          Gathering list               ',3f8.1,' (sec)')
  700 format('          Scattering list              ',3f8.1,' (sec)')
  750 format('          Updating lists               ',3f8.1,' (sec)')
  800 format('             Coordinates               ',3f8.1,' (sec)')
  850 format('Total time sorting lists               ',f8.1,' (sec)')
  900 format('          Sorting local list           ',3f8.1,' (sec)')
 1000 format('          Sorting buffer list          ',3f8.1,' (sec)')
 1050 format('Total time sending messages            ',f8.1,' (sec)')
 1100 format('          Gather                       ',3f8.1,' (sec)')
 1200 format('          Scatter                      ',3f8.1,' (sec)')
 1300 format('          Update                       ',3f8.1,' (sec)')
 1400 format('             Angles                    ',3f8.1,' (sec)')
 1500 format('Total time calculating stress tensor   ',3f8.1,' (sec)')
 1600 format('Total time in scatter synchronization  ',3f8.1,' (sec)')
 2300 format('Total time constructing neighbor list  ',3f8.1,' (sec)')
 2800 format('             Bonds                     ',3f8.1,' (sec)')
c
      return
      end
