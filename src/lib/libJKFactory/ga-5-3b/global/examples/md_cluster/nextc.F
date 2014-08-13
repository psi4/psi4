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
      subroutine estep
#include "common.fh"
      integer iter
      iter = 0
c
      prsflg = .false.
      tmpflg = .false.
      ptflg = .false.

      do while (t_rmndr.gt.0.0d00)
        iter = iter + 1
        call pred
        call force
        call corr
        if (iter.gt.100) then
          write(6,*) 'istep: ',istep
          write(6,*) 't_rmndr: ',t_rmndr
        endif
        if (iter.gt.1000) then
          call ga_error('Failed in estep',istep)
        endif
      end do
c
      return
      end
c
      subroutine pstep
#include "common.fh"
      integer iter
      iter = 0
c
      prsflg = .true.
      tmpflg = .false.
      ptflg = .false.
      do while (t_rmndr.gt.0.0d00)
        iter = iter + 1
        call pred
        call force
        call corr
        if (iter.gt.100) then
          write(6,*) 'istep: ',istep
          write(6,*) 't_rmndr: ',t_rmndr
        endif
        if (iter.gt.1000) then
          call ga_error('Failed in pstep',istep)
        endif
      end do
c
      return
      end
c
      subroutine tstep
#include "common.fh"
      integer iter
      iter = 0
c
      prsflg = .false.
      tmpflg = .true.
      ptflg = .false.
      do while (t_rmndr.gt.0.0d00)
        iter = iter + 1
        call pred
        call force
        call corr
        if (iter.gt.100) then
          write(6,*) 'istep: ',istep
          write(6,*) 't_rmndr: ',t_rmndr
        endif
        if (iter.gt.1000) then
          call ga_error('Failed in tstep',istep)
        endif
      end do
c
      return
      end
c
      subroutine sstep
#include "common.fh"
c
      double precision scale
      integer i,j
      integer iter
      iter = 0
c
c   This subroutine scales the atomic velocities and momenta
c   so that the system temperature corresponds to the desired
c   temperature.
c
      if (nrg(5).gt.0.0d00) then
        scale = sqrt(tmprtr / nrg(5))
      else
        scale = 1.0d00
      endif
      do 200 j = 1, 3
        do 100 i = 1, antot
          ra(i,j,2) = scale * ra(i,j,2)
          ra(i,j,5) = scale * ra(i,j,5)
  100   continue
  200 continue
      prsflg = .false.
      tmpflg = .false.
      ptflg = .false.
      mbflg = .true.
      do while (t_rmndr.gt.0.0d00)
        iter = iter + 1
        call pred
        call force
        call corr
        if (iter.gt.100) then
          write(6,*) 'istep: ',istep
          write(6,*) 't_rmndr: ',t_rmndr
        endif
        if (iter.gt.1000) then
          call ga_error('Failed in sstep',istep)
        endif
      end do
c
      return
      end
c
      subroutine kstep
#include "common.fh"
c
      double precision scale,ke
      integer i,j
      integer iter
      iter = 0
c
c   This subroutine adjusts the kinetic energy so that the
c   average temperature corresponds to the target temperature
c
      if (itavg.gt.0) then
        tavg = tavg / dble(itavg)
        ke = nrg(2) * (tmprtr - tavg)
        scale = (nrg(4) + ke)/nrg(4)
      else
        scale = 0.0d00
      endif
      if (scale.gt.0.0d00) then
        scale = sqrt(scale)
      else
        scale = 1.0d00
      endif
      do 200 j = 1, 3
        do 100 i = 1, antot
          ra(i,j,2) = scale * ra(i,j,2)
          ra(i,j,5) = scale * ra(i,j,5)
  100   continue
  200 continue
      tavg = 0.0d00
      itavg = 0
      prsflg = .false.
      tmpflg = .false.
      ptflg = .false.
      mbflg = .true.
      do while (t_rmndr.gt.0.0d00)
        iter = iter + 1
        call pred
        call force
        call corr
        if (iter.gt.100) then
          write(6,*) 'istep: ',istep
          write(6,*) 't_rmndr: ',t_rmndr
        endif
        if (iter.gt.1000) then
          call ga_error('Failed in kstep',istep)
        endif
      end do
c
      return
      end
c
      subroutine ptstep
#include "common.fh"
      integer iter
      iter = 0
c
      prsflg = .false.
      tmpflg = .false.
      ptflg = .true.
      do while (t_rmndr.gt.0.0d00)
        iter = iter + 1
        call pred
        call force
        call corr
        if (iter.gt.100) then
          write(6,*) 'istep: ',istep
          write(6,*) 't_rmndr: ',t_rmndr
        endif
        if (iter.gt.1000) then
          call ga_error('Failed in ptstep',istep)
        endif
      end do
      return
      end
c
      subroutine mbstep
#include "common.fh"
c
      double precision scl,gasdev
      integer i,j
      integer iter
      iter = 0
c
c   This subroutine updates configuration if only a single 
c   configuration is available. It starts by assigning a
c   Boltzmann distribution of velocities and momenta to all 
c   the atoms and then calculating the next configuration
c   using a regular Gear step.
c
c    generate atomic velocities
c
      do 200 j = 1, 3
        do 100 i = 1, antot
          scl = sqrt(tmprtr / mass(i))
          ra(i,j,2) = scl * gasdev(0)
          ra(i,j,5) = mass(i) * ra(i,j,2)
          ra(i,j,3) = ra(i,j,4)/mass(i)
  100   continue
  200 continue
c
c   Re-initialize absolute coordinates
c
      call com
      call com_rmv
      call cluster_com
      tmstrt = istep
      scal2 = 0.0d00
      scal3 = 0.0d00
      vol2 = 0.0d00
      vol3 = 0.0d00
      prsflg = .false.
      tmpflg = .false.
      ptflg = .false.
      mbflg = .true.
      cmflg = .true.
      do while (t_rmndr.gt.0.0d00)
        iter = iter + 1
        call pred
        call force
        call corr
        if (iter.gt.100) then
          write(6,*) 'istep: ',istep
          write(6,*) 't_rmndr: ',t_rmndr
        endif
        if (iter.gt.1000) then
          call ga_error('Failed in mbstep',istep)
        endif
      end do
#if 0
      if (mod(istep,ilist).eq.0) then
        if (istep.gt.equil_1) then
          call cluster_com
          call cluster_center
        endif
        call fixper
        do i = 1, antot
          ra(i,1,6) = ra(i,1,1)
          ra(i,2,6) = ra(i,2,1)
          ra(i,3,6) = ra(i,3,1)
        end do
      endif
#endif
      return
      end
c
      subroutine vlstep
#include "common.fh"
c
      double precision scale,sclinc
      integer i,j
      integer iter
      iter = 0
c
      prsflg = .false.
      tmpflg = .false.
      ptflg = .false.
      if (itarg.ge.istep) then
        scale = xbox * ybox * zbox
        sclinc = (tvol - scale ) / float(itarg - istep + 1)
        scale = (scale + sclinc) / (xbox * ybox * zbox)
        scale = exp(log(scale)/3.0d00)
        xbox = scale * xbox
        ybox = scale * ybox
        zbox = scale * zbox
        xbox2 = scale * xbox2
        ybox2 = scale * ybox2
        zbox2 = scale * zbox2
        do 300 j = 1, 3
          do 100 i = 1, antot
            ra(i,j,1) = scale * ra(i,j,1)
            ra(i,j,6) = scale * ra(i,j,6)
  100     continue
  300   continue
      endif
      mbflg = .true.
      do while (t_rmndr.gt.0.0d00)
        iter = iter + 1
        call pred
        call force
        call corr
        if (iter.gt.100) then
          write(6,*) 'istep: ',istep
          write(6,*) 't_rmndr: ',t_rmndr
        endif
        if (iter.gt.1000) then
          call ga_error('Failed in vlstep',istep)
        endif
      end do
      if (nrg(5).gt.0.0d00) then
        scale = sqrt(tmprtr / nrg(5))
      else
        scale = 1.0d00
      endif
      do 600 j = 1, 3
        do 400 i = 1, antot
          ra(i,j,2) = scale * ra(i,j,2)
          ra(i,j,5) = scale * ra(i,j,5)
  400   continue
  600 continue
c
      return
      end
c
      subroutine pred
#include "common.fh"
c
      double precision scale,c1,c2
      double precision rx,ry,rz,dt
      integer i,j
      double precision cluster_check_radius
      logical debug
c
c   This subroutine implements the predictor for a 3 point Gear
c   predictor-corrector algorithm
c
      if (istep.gt.3930438) then
        debug = .false.
      else
        debug = .false.
      endif
      call cluster_com
      call cluster_old_at
      call cluster_therm
      if (debug) then
        write(6,*) ga_nodeid(),' (pred) Got to 1 at step ',istep
      endif
c      rx = cluster_check_radius()
c      if (rx.gt.r_cluster.and.istep.gt.equil_3) then
c        write(6,*) ga_nodeid(),' Step ',istep
c        write(6,*) ga_nodeid(),' t_done               ',t_done
c        write(6,*) ga_nodeid(),' Radius out of bounds ',rx
c        write(6,*) ga_nodeid(),' r_cluster            ',r_cluster
c      endif
      if (t_done.eq.tau) then
        t_done = 0.0d00
        dt = tau
        t_rmndr = tau
      else
        dt = t_rmndr
      endif
c
      c1 = dt
      c2 = dt * c1 / 2.0d00
c
c   predictor step for atoms
c
      do 200 j = 1, 3
        do 100 i = 1, antot
          ra(i,j,6) = ra(i,j,6) + c1 * ra(i,j,2) + c2 * ra(i,j,3)
          ra(i,j,2) = ra(i,j,2) + c1 * ra(i,j,3) 
  100   continue
  200 continue
c
c   predictor step for volume
c
      if ((prsflg.or.ptflg).and.ipmode.eq.0) then
        vol1 = vol1 + c1 * vol2 + c2 * vol3
        vol2 = vol2 + c1 * vol3
        scale = vol1 / (xbox * ybox * zbox)
        scale = exp(log(scale)/3.0d00)
        xbox = scale * xbox
        ybox = scale * ybox
        zbox = scale * zbox
        xbox2 = scale * xbox2
        ybox2 = scale * ybox2
        zbox2 = scale * zbox2
      endif
c
      if ((prsflg.or.ptflg).and.ipmode.eq.1) then
        do 500 j = 1, 3
          alen1(j) = alen1(j) + c1 * alen2(j) + c2 * alen3(j)
          alen2(j) = alen2(j) + c1 * alen3(j)
  500   continue
        xbox = alen1(1)
        ybox = alen1(2)
        zbox = alen1(3)
        xbox2 = xbox/2.0d00
        ybox2 = ybox/2.0d00
        zbox2 = zbox/2.0d00
      endif
c
      if ((prsflg.or.ptflg).and.ipmode.eq.2) then
        do 600 j = 1, 2
          alen1(j) = alen1(j) + c1 * alen2(j) + c2 * alen3(j)
          alen2(j) = alen2(j) + c1 * alen3(j)
  600   continue
        xbox = alen1(1)
        ybox = alen1(2)
        xbox2 = xbox/2.0d00
        ybox2 = ybox/2.0d00
      endif
c
c   predictor step for time scale
c
      if (tmpflg.or.ptflg) then
        scal1 = scal1 + c1 * scal2 + c2 * scal3
        scal2 = scal2 + c1 * scal3 
      endif
c
      if (debug) then
        write(6,*) ga_nodeid(),' (pred) Got to 2 at step ',istep
      endif
      call fixper
      if (istep.gt.equil_1) then
        call cluster_com
      endif
      if (debug) then
        write(6,*) ga_nodeid(),' (pred) Got to 3 at step ',istep
      endif
      if (istep.gt.equil_1) then
        call cluster_check_cllsn
      else
        t_rmndr = 0.0d00
        t_done = tau
      endif
      if (debug) then
        write(6,*) ga_nodeid(),' (pred) Got to 4 at step ',istep
      endif
      if (.not.(ptflg.or.prsflg)) then
        if (t_rmndr.eq.0.0d00) call cluster_center
      endif
      if (debug) then
        write(6,*) ga_nodeid(),' (pred) Got to 5 at step ',istep
      endif
      if (mod(istep,ilist).eq.0.and.t_rmndr.eq.0.0d00) then
        call update
      endif
      if (cllsn_cnt.gt.100) then
        call ga_error("Too many collisions ",istep)
      endif
      if (debug) then
        write(6,*) ga_nodeid(),' (pred) Got to 6 at step ',istep
      endif
c
      return
      end
c
      subroutine corr
#include "common.fh"
c
      double precision c1,c2,coef1,xcorr,vcorr,scorr,dt
      double precision v1,s1,pmassi,tmassi,dlnvdt(3),itmp
      double precision mv1,cmr,kvt,kvtl(3),l1(3)
      double precision rbuf(8),virial,virc(3)
      double precision alfx,alfy,alfpx,alfpy,alfpz
      integer i,j,iat,nreal,nbyte
c
c   This subroutine implements the corrector for a 3 point Gear
c   predictor-corrector algorithm
c
      dt = t_done
      c1 = dt
      c2 = dt * c1 / 2.0d00
      coef1 = c2 / c1
c
      if (prsflg.or.ptflg) then
        pmassi = 1.0d00 / pmass
      else
        pmassi = 0.0d00
        vol1 = xbox * ybox * zbox
        vol2 = 0.0
        vol3 = 0.0
c
        alen1(1) = xbox
        alen1(2) = ybox
        alen1(3) = zbox
        do  20 j = 1, 3
          alen2(j) = 0.0d00
          alen3(j) = 0.0d00
   20   continue
      endif
      if (tmpflg.or.ptflg) then
        tmassi = 1.0d00 / tmass
      else
        tmassi = 0.0d00
        scal1 = 1.0d00
        scal2 = 0.0d00
        scal3 = 0.0d00
      endif
c
      if (ipmode.eq.0) then
        dlnvdt(1) = vol2 / (3.0d00 * vol1)
        dlnvdt(2) = vol2 / (3.0d00 * vol1)
        dlnvdt(3) = vol2 / (3.0d00 * vol1)
      elseif (ipmode.eq.1) then
        dlnvdt(1) = alen2(1) / alen1(1)
        dlnvdt(2) = alen2(2) / alen1(2)
        dlnvdt(3) = alen2(3) / alen1(3)
      elseif (ipmode.eq.2) then
        dlnvdt(1) = alen2(1) / alen1(1)
        dlnvdt(2) = alen2(2) / alen1(2)
        dlnvdt(3) = 0.0d00
      endif
      kvt = 0.0d00
      kvtl(1) = 0.0d00
      kvtl(2) = 0.0d00
      kvtl(3) = 0.0d00
      do 200 j = 1, 3
        do 150 iat = 1, antot
          mv1 = mass(iat) * (ra(iat,j,2) - ra(iat,j,6) * dlnvdt(j))
          kvt = kvt + mv1**2/mass(iat)
          kvtl(j) = kvtl(j) + mv1**2/mass(iat)
  150   continue
  200 continue
c
c   do global sum on quantities needed for corrector, if necessary
c
      if (prsflg.or.ptflg.or.tmpflg) then
        rbuf(1) = kvt
        rbuf(2) = kvtl(1)
        rbuf(3) = kvtl(2)
        rbuf(4) = kvtl(3)
        rbuf(5) = nrg(10)
        rbuf(6) = nrg(18)
        rbuf(7) = nrg(19)
        rbuf(8) = nrg(20)
c
        call ga_dgop(3,rbuf,8,'+')
c
        kvt = rbuf(1)
        kvtl(1) = rbuf(2)
        kvtl(2) = rbuf(3)
        kvtl(3) = rbuf(4)
        virial = rbuf(5)
        virc(1) = rbuf(6)
        virc(2) = rbuf(7)
        virc(3) = rbuf(8)
      endif
c
      itmp = (kvt + virial)/(3.0d00 * vol1)
      if (prsflg.or.ptflg) then
        if (ipmode.eq.0) then
          v1 = scal1**2 * pmassi * (itmp - prssr) + vol2 * scal2 / scal1
        else
          l1(1) = scal1**2*pmassi*((kvtl(1)+virc(1))/alen1(1)
     +          - prssr*alen1(2)*alen1(3)) + alen2(1)*scal2/scal1
          l1(2) = scal1**2*pmassi*((kvtl(2)+virc(2))/alen1(2)
     +          - prssr*alen1(1)*alen1(3)) + alen2(2)*scal2/scal1
          if (ipmode.eq.1) then
            l1(3) = scal1**2*pmassi*((kvtl(3)+virc(3))/alen1(3)
     +            - prssr*alen1(1)*alen1(2)) + alen2(3)*scal2/scal1
          else
            l1(3) = 0.0d00
          endif
          alfx = l1(1)
          alfy = l1(2)
        endif
      else
        v1 = 0.0d00
      endif
c
      do 400 j = 1, 3
        if (ipmode.eq.0) then
          do 300 i = 1, antot
            cmr = ra(i,j,6)
            ra(i,j,8) = ra(i,j,4) / mass(i) 
     +         - (dlnvdt(j) + scal2 / scal1)
     +         * (ra(i,j,2) - dlnvdt(j) * cmr)
     +         - (vol2**2 * cmr / (3.0d00*vol1**2) 
     +         - v1 * cmr / (3.0d00 * vol1)
     +         - dlnvdt(j) * ra(i,j,2)) 
  300     continue
        elseif (ipmode.ge.1) then
          do 340 i = 1, antot
            cmr = ra(i,j,6)
            ra(i,j,8) = ra(i,j,4) / mass(i)
     +        - scal2 * ra(i,j,2) / scal1
     +        + cmr * scal2 * dlnvdt(j) / scal1
     +        + cmr * l1(j) / alen1(j)
  340     continue
        endif
  400 continue
c
      do 425 j = 1, 3
        do 450 i = 1, antot
          xcorr = ra(i,j,8) - ra(i,j,3)
          ra(i,j,2) = ra(i,j,2) + coef1 * xcorr
          ra(i,j,3) = ra(i,j,8)
          ra(i,j,7) = ra(i,j,8)
  450   continue
  425 continue
c
      if (tmpflg.or.ptflg) then
        s1 = scal1 * tmassi * (kvt - nrg(1) * tmprtr)
     +     + scal2**2 / scal1
      endif
      if ((prsflg.or.ptflg).and.ipmode.eq.0) then
        vcorr = v1 - vol3
        vol2 = vol2 + coef1 * vcorr
        vol3 = v1
      endif
      if ((prsflg.or.ptflg).and.ipmode.ge.1) then
        l1(1) = l1(1) + alfpx
        l1(2) = l1(2) + alfpy
        l1(3) = l1(3) + alfpz
        if (ipmode.eq.1) then
          do 475 j = 1, 3
            vcorr = l1(j) - alen3(j)
            alen2(j) = alen2(j) + coef1 * vcorr
            alen3(j) = l1(j)
  475     continue
        elseif (ipmode.eq.2) then
          do 485 j = 1, 2
            vcorr = l1(j) - alen3(j)
            alen2(j) = alen2(j) + coef1 * vcorr
            alen3(j) = l1(j)
  485     continue
        endif
        alfpx = l1(1)
        alfpy = l1(2)
        alfpz = l1(3)
      endif
      if (tmpflg.or.ptflg) then
        scorr = s1 - scal3
        scal2 = scal2 + coef1 * scorr
        scal3 = s1
      endif
c
c   calculate value of the Hamiltonian
c
      kvt = 0.0d00
      if (ipmode.eq.0) then
        dlnvdt(1) = vol2 / (3.0d00 * vol1)
        dlnvdt(2) = vol2 / (3.0d00 * vol1)
        dlnvdt(3) = vol2 / (3.0d00 * vol1)
      elseif (ipmode.eq.1) then
        dlnvdt(1) = alen2(1) / alen1(1)
        dlnvdt(2) = alen2(2) / alen1(2)
        dlnvdt(3) = alen2(3) / alen1(3)
        vol1 = alen1(1)*alen1(2)*alen1(3)
      elseif (ipmode.eq.2) then
        dlnvdt(1) = alen2(1) / alen1(1)
        dlnvdt(2) = alen2(2) / alen1(2)
        dlnvdt(3) = 0.0d00
        vol1 = alen1(1)*alen1(2)*zbox
      endif
      do 600 j = 1, 3
        do 550 iat = 1, antot
          ra(iat,j,5) = mass(iat) * (ra(iat,j,2)-ra(iat,j,6)*dlnvdt(j))
          kvt = kvt + ra(iat,j,5)**2 / mass(iat)
  550   continue
  600 continue
c
      call com
      call kin
      call nrgsum
      call com_rmv
c
      kvt = 2.0 * nrg(4)
      nrg(9) = 0.5d00 * (kvt + tmass * scal2**2 / scal1**2) + nrg(6)
      if (ipmode.eq.0) nrg(9) = nrg(9)+0.5d00*pmass*vol2**2/scal1**2
      if (ipmode.eq.1) nrg(9) = nrg(9)+0.5d00*pmass*alen2(1)**2/scal1**2
     +                                +0.5d00*pmass*alen2(2)**2/scal1**2
     +                                +0.5d00*pmass*alen2(3)**2/scal1**2
      if (ipmode.eq.2) nrg(9) = nrg(9)+0.5d00*pmass*alen2(1)**2/scal1**2
     +                                +0.5d00*pmass*alen2(2)**2/scal1**2
      itmp = (kvt + nrg(10))/(3.0d00 * vol1)
      if (prsflg.or.ptflg) then
        nrg(11) = 0.5d00 * kvt + nrg(6) + vol1 * prssr
      else
        nrg(11) = 0.5d00 * kvt + nrg(6) + vol1 * itmp
      endif
      if (prsflg.or.ptflg) nrg(9) = nrg(9) + prssr * vol1
      if (tmpflg.or.ptflg) then
        nrg(9) = nrg(9) + nrg(1) * tmprtr * log(scal1)
      endif
c
      nrg(5) = nrg(4) / nrg(2)
      tavg = tavg + nrg(5)
      itavg = itavg + 1
c
c    make sure volume and temperature scale factor is the same on
c    all processors
c
      if ((mod(istep,100).eq.0).and.(prsflg.or.ptflg.or.tmpflg)) then
        if (ga_nodeid().eq.0) then
          rbuf(1) = scal1
          rbuf(2) = scal2
          rbuf(3) = xbox
          rbuf(4) = ybox
          rbuf(5) = zbox
          rbuf(6) = vol2
        else
          rbuf(1) = 0.0d00
          rbuf(2) = 0.0d00
          rbuf(3) = 0.0d00
          rbuf(4) = 0.0d00
          rbuf(5) = 0.0d00
          rbuf(6) = 0.0d00
        endif
c
        nreal = 6
        call ga_dgop(1,rbuf,6,'+')
c
        scal1 = rbuf(1)
        scal2 = rbuf(2)
        xbox = rbuf(3)
        ybox = rbuf(4)
        zbox = rbuf(5)
        vol2 = rbuf(6)
      endif
      if (l_cllsn) then 
        call cluster_do_cllsn
      endif
      call cluster_com
c
      return
      end
