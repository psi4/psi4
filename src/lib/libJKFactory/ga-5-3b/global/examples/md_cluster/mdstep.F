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
      subroutine mdstep
#include "common.fh"
c
      double precision vbox, rmax
      integer i,j,me
      logical newcfg
      double precision cluster_check_radius
      logical debug
      if (istep.gt.3930438) then
        debug = .false.
      else
        debug = .false.
      endif
c
c   This routine guides the MD steps.
c   Begin the main loop through the MD steps
c
      me = ga_nodeid()
      r_confine = 0.0d00
      do 5000 istep = 1, nstep
c
        newcfg = .false.
        mbflg = .false.
        cmflg = .false.
        ipmode = 0
        t_rmndr = tau
        t_done = 0.0d00
        cllsn_cnt = 0
c
c  Check to see if there are any special instructions
c
        do 100 i = 1, nsc
c
c  is end >= istep >= beg
c
          if ((istep.ge.isc(i,1)).and.(istep.le.isc(i,2))) then
c
c   is mod(istep-beg,inc) = 0
c
            if (mod(istep,isc(i,3)).eq.0) then
c
c   get next configuration using the appropriate algorithm
c
              if (isc(i,4).eq.1) then
                call estep
                newcfg = .true.
              elseif (isc(i,4).eq.2) then
                prssr = rsc(i,2) 
                pmass = rsc(i,4)
                call pstep
                newcfg = .true.
              elseif (isc(i,4).eq.3) then
                tmprtr = rsc(i,1)
                prssr = rsc(i,2)
                pmass = rsc(i,4)
                call sstep
                newcfg = .true.
              elseif (isc(i,4).eq.4) then
                tmprtr = rsc(i,1)
                tmass = rsc(i,3)
                call tstep
                newcfg = .true.
              elseif (isc(i,4).eq.5) then
                tmprtr = rsc(i,1)
                prssr = rsc(i,2) 
                tmass = rsc(i,3)
                pmass = rsc(i,4)
                call ptstep
                newcfg = .true.
              elseif (isc(i,4).eq.6) then
                tmprtr = rsc(i,1)
                call mbstep
                newcfg = .true.
              elseif (isc(i,4).eq.7) then
                itarg = isc(i,2)
                tmprtr = rsc(i,1)
                tvol = rsc(i,2)
                tmass = rsc(i,3)
                call vlstep
                newcfg = .true.
              elseif (isc(i,4).eq.8) then
                tmprtr = rsc(i,1)
                call kstep
                newcfg = .true.
              elseif (isc(i,4).eq.9) then
                tmprtr = rsc(i,1)
                prssr = rsc(i,2) 
                tmass = rsc(i,3)
                pmass = rsc(i,4)
                ipmode = 1
                call ptstep
                newcfg = .true.
              elseif (isc(i,4).eq.10) then
                tmprtr = rsc(i,1)
                prssr = rsc(i,2) 
                tmass = rsc(i,3)
                pmass = rsc(i,4)
                ipmode = 2
                call ptstep
                newcfg = .true.
              endif
            endif
          endif
  100   continue
c
c  get next configuration if no special step is taken
c
        if (.not.newcfg) then
          if (dflalg.eq.1) then
            call estep
          elseif (dflalg.eq.2) then
            prssr = dfprs
            pmass = dfpm
            call pstep
          elseif (dflalg.eq.3) then
            tmprtr = dftmp
            prssr = dfprs
            pmass = dfpm
            call sstep
          elseif (dflalg.eq.4) then
            tmprtr = dftmp
            tmass = dftm
            call tstep
          elseif (dflalg.eq.5) then
            tmprtr = dftmp
            prssr = dfprs
            tmass = dftm
            pmass = dfpm
            call ptstep
          elseif (dflalg.eq.6) then
            tmprtr = dftmp
            call mbstep
          elseif (dflalg.eq.9) then
            tmprtr = dftmp
            prssr = dfprs
            tmass = dftm
            pmass = dfpm
            ipmode = 1
            call ptstep
          elseif (dflalg.eq.10) then
            tmprtr = dftmp
            prssr = dfprs
            tmass = dftm
            pmass = dfpm
            ipmode = 2
            call ptstep
          endif
        endif
            if (debug) then
              write(6,*) ga_nodeid(),' Got to 1 at step ',istep
            endif
c
c  Update remaining energy quantities
c
        nrg(3) = nrg(4) + nrg(6)
        vbox = xbox*ybox*zbox
        nrg(7) = nrg(5) * dble(atot-1) / vbox + nrg(15)

        if (istep.eq.equil_1) then
          call fixper
          do i = 1, antot
            do j = 1, 3
              ra(i,j,6) = ra(i,j,1)
            end do
          end do
          call cluster_com
          call cluster_center
          rmax = cluster_check_radius()
          if (rmax.gt.r_cluster) r_cluster = rmax + 0.01
        endif
        if (mod(istep,mcfreq).eq.0.and.istep.gt.equil_1) then
          call cluster_mc
          if (me.eq.0.and.l_rad) write(7,7100) dble(istep)*tau,r_cluster
        endif
            if (debug) then
              write(6,*) ga_nodeid(),' Got to 2 at step ',istep
            endif
        if (istep.gt.equil_2.and.r_cluster.le.cl_upper)
     +    call cluster_binr
            if (debug) then
              write(6,*) ga_nodeid(),' Got to 3 at step ',istep
            endif
        if (istep.eq.window_1) call cluster_reset_binr(1)
        if (istep.eq.window_2) call cluster_reset_binr(2)
c
c   Perform all statistical operations
c   on the new configuration.
c
c   print pressure
c
         if (mod(istep,istat).eq.0.and.l_stdio) then
           call header(istep)
           if (me.eq.0) write(6,6000) nrg(7) 
           if (me.eq.0) write(6,6300) nrg(3)
           if (me.eq.0) write(6,6700) nrg(6),nrg(4),nrg(5)
           if (me.eq.0) write(6,6800) xbox,ybox,zbox,scal1
           if (me.eq.0) write(6,6100) nrg(13),nrg(14),
     +       nrg(17),nrg(21)
           if (me.eq.0) write(6,6900) nrg(9)
         endif
c
c   accumulate energy statistics
c
         if (istep.gt.equil_2) call estat
c
         if (me.eq.0.and.l_step.and.mod(istep,1000).eq.0) then
           open(unit=2,file='step.cnt',status='unknown')
           write(2,*) 'proc : ',ga_pgroup_nodeid(ga_pgroup_get_world())
           write(2,*) 'istep : ',istep
           write(2,6900) nrg(9)
           write(2,6100) nrg(13),nrg(14),nrg(17),nrg(21)
           write(2,6000) nrg(7) 
           write(2,6300) nrg(3)
           write(2,6700) nrg(6),nrg(4),nrg(5)
           write(2,6800) xbox,ybox,zbox,scal1
           write(2,7000) r_cluster
           write(2,7200) cl_lower
           write(2,7300) cl_upper
           close(2)
         endif
 5000 continue
      return
 6000 format(1x,'The instantaneous pressure is ',f12.4)
 6100 format(1x,'Current energy statistics'/
     +          '       repulsion: ',f16.4,/
     +          '      dispersion: ',f16.4,/
     +          '           bonds: ',f16.4,/
     +          '          angles: ',f16.4)
 6200 format(1x,'Statistics at time ',i6,' ps')
 6300 format(1x,'The total energy is ',f12.4)
 6700 format('     potential      kinetic'/
     +       '     energy         energy   temperature'/,
     +          1x,3f13.3)
 6800 format(1x,'The current simulation cell dimensions:'/
     +         '              x:   ',f16.4,/
     +         '              y:   ',f16.4,/
     +         '              z:   ',f16.4,/
     +         '              s:   ',f16.4)
 6820 format(1x,'The current simulation cell dimensions:'/
     +         '              x:   ',f16.4,/
     +         '              y:   ',f16.4,/
     +         '              s:   ',f16.4)
 6900 format(1x,'The instantaneous value of the Hamiltonian is ',f12.4)
 7000 format(1x,'The current value of confining sphere is      ',f12.4)
 7100 format(2f16.8)
 7200 format(1x,'Lower bound of confining sphere               ',f12.4)
 7300 format(1x,'Upper bound of confining sphere               ',f12.4)
      end
