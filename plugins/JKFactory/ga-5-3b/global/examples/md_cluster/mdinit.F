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
      subroutine mdinit
#include "common.fh"
c
      double precision six,twelve,ac
      integer i,j
c
c  turn of all cluster related functions
c
      nocluster = .false.
c
c    This subroutine sets up the potentials, assigns atomic masses
c    and finishes initializing the calculation
c
      six = 6.0d00
      twelve = 12.0d00
c
c      initialize timing and energy statistics arrays
c
      do i = 1, MAXNRG
        nrg(i) = 0.0d00
        anrg(i) = 0.0d00
        anrg2(i) = 0.0d00
      end do
      nrgcnt = 0
      nrg2ct = 0
      do i = 1, MAXTIM
        tmstat(i) = 0.0d00
      end do
      failcount = 0
c
c     assign atom masses
c
      do i = 1, antot
        mass(i) = amass(at(i))
      end do
c
c   Initialize cutoffs for  potential arrays.
c
      rcmax = 0.0d00
      do i = 1, atnum
        do j = 1, atnum
          ac = acut(i,j)
          if (ac.gt.rcmax) rcmax = ac
          rcorr(i,j) = e12(i,j) / ac**12
          dcorr(i,j) = -e6(i,j) / ac**6
          if (icut.eq.1) then
            frcorr(i,j) = - twelve * e12(i,j) / ac**13
            fdcorr(i,j) = six * e6(i,j) / ac**7
          else
            frcorr(i,j) = 0.0d00
            fdcorr(i,j) = 0.0d00
          endif
        end do
      end do
      rcmax = rcmax + 1.0d00
c
c   calculate total number of degrees of freedom
c
      nrg(1) = dble(3*atot-3)
      nrg(2) = nrg(1) / 2.0
c
c   Convert algorithm parameters into system units
c
      istep = 0
      tavg = 0.0d00
      itavg = 0
      tmstrt = 0
      t_done = 0.0d00
      t_rmndr = tau
c
c   Initialize parameters for mass and temperature degrees of freedom
c
      vol1 = xbox * ybox * zbox
      alen1(1) = xbox
      alen1(2) = ybox
      alen1(3) = zbox
      vol2 = 0.0d00
      vol3 = 0.0d00
      do j = 1, 3
        alen2(j) = 0.0d00
        alen3(j) = 0.0d00
      end do
      scal1 = 1.0d00
      scal2 = 0.0d00
      scal3 = 0.0d00
c
c   Initialize bins accumulating r_cluster values
c
      do j = 1, 10
        do i = 1, mcbins
          r_distr(i,j) = 0
        end do
      end do
c
c   Initialize absolute coordinates
c
      call cluster_center
      call update
      call force
      mbflg = .true.
      call cluster_com
c
      tavg = 0.0d00
      do j = 1, 3
        do i = 1, antot
          ra(i,j,7) = ra(i,j,8)
          ra(i,j,3) = ra(i,j,4)/mass(i)
        end do
      end do
c
      return
      end
