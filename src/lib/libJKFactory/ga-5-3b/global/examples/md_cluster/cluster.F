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
      subroutine cluster_com
#include "common.fh"
c
      integer i
      double precision masstot, dx, dy, dz
      double precision r,rdot,com(10)
c
c   This subroutine calculates the center of mass (COM) of the cluster of
c   particles of type 2, excluding the lone particle.
c
      if (nocluster) return
      cl_cmx = 0.0d00
      cl_cmy = 0.0d00
      cl_cmz = 0.0d00
      cl_acmx = 0.0d00
      cl_acmy = 0.0d00
      cl_acmz = 0.0d00
      cl_vcmx = 0.0d00
      cl_vcmy = 0.0d00
      cl_vcmz = 0.0d00
      masstot = 0.0d00
c
      do i = 1, antot
c      if (at(i).eq.2.and.aidx(i).ne.cl_lone_particle) then
        if (at(i).eq.2) then
          cl_cmx = cl_cmx + mass(i)*ra(i,1,6)
          cl_cmy = cl_cmy + mass(i)*ra(i,2,6)
          cl_cmz = cl_cmz + mass(i)*ra(i,3,6)
          cl_acmx = cl_acmx + ra(i,1,4)
          cl_acmy = cl_acmy + ra(i,2,4)
          cl_acmz = cl_acmz + ra(i,3,4)
          cl_vcmx = cl_vcmx + mass(i)*ra(i,1,2)
          cl_vcmy = cl_vcmy + mass(i)*ra(i,2,2)
          cl_vcmz = cl_vcmz + mass(i)*ra(i,3,2)
          masstot = masstot + mass(i)
        endif
      end do
      com(1) = cl_cmx
      com(2) = cl_cmy
      com(3) = cl_cmz
      com(4) = cl_vcmx
      com(5) = cl_vcmy
      com(6) = cl_vcmz
      com(7) = cl_acmx
      com(8) = cl_acmy
      com(9) = cl_acmz
      com(10) = masstot
      call ga_dgop(3,com,10,'+')
      cl_cmx  = com(1)
      cl_cmy  = com(2)
      cl_cmz  = com(3)
      cl_acmx  = com(7)
      cl_acmy  = com(8)
      cl_acmz  = com(9)
      cl_vcmx  = com(4)
      cl_vcmy  = com(5)
      cl_vcmz  = com(6)
      masstot = com(10)
      if (masstot.gt.0.0d00) then
        cl_cmx = cl_cmx / masstot
        cl_cmy = cl_cmy / masstot
        cl_cmz = cl_cmz / masstot
        cl_acmx = cl_acmx / masstot
        cl_acmy = cl_acmy / masstot
        cl_acmz = cl_acmz / masstot
        cl_vcmx = cl_vcmx / masstot
        cl_vcmy = cl_vcmy / masstot
        cl_vcmz = cl_vcmz / masstot
      endif
      cl_mass = masstot
c
      return
      end
c
      subroutine cluster_therm
#include "common.fh"
c
      integer i
      double precision dx, dy, dz
      double precision r,rdot
c
c  remove component of velocity that lies along the vector from COM
c  to the particle, if particle is outside cluster radius. This gradually
c  forces particles to form a drop.
c
      if (nocluster) return
      if (istep.ge.equil_1) return
      call cluster_com
      do i = 1, antot
c        if (at(i).eq.2.and.aidx(i).ne.cl_lone_particle) then
        if (at(i).eq.2) then
          dx = ra(i,1,1) - cl_cmx
          dy = ra(i,2,1) - cl_cmy
          dz = ra(i,3,1) - cl_cmz
          r = sqrt(dx**2 + dy**2 + dz**2)
          dx = dx/r
          dy = dy/r
          dz = dz/r
          if (r.gt.r_cluster) then
            rdot = ra(i,1,2)*dx+ra(i,2,2)*dy+ra(i,3,2)*dz
            if (rdot.gt.0.0d00) then
              ra(i,1,2) = ra(i,1,2) - 2.0d00*rdot*dx
              ra(i,2,2) = ra(i,2,2) - 2.0d00*rdot*dy
              ra(i,3,2) = ra(i,3,2) - 2.0d00*rdot*dz
            endif
          endif
        endif
      end do
c
      return
      end
c
      subroutine cluster_center
#include "common.fh"
c
c  move center of mass to origin (only perform this if it is immediately
c  followed by a call to update subroutine)
c
      integer i,me
      logical debug
      if (nocluster) return
      if (istep.eq.0.or.(mod(istep,ilist).eq.0.and.
     +    t_rmndr.eq.0.0d00)) then
        me = ga_nodeid()
        if (istep.ge.6932366) then
          debug = .false.
        else
          debug = .false.
        endif
        if (debug) then
          write(6,*) me,' mod(istep,ilist) = ',mod(istep,ilist)
          write(6,*) me,' ilist = ',ilist
        endif
        call cluster_com
        if (debug) write(6,*) me,'cl_cmx (a) = ',cl_cmx,istep
        if (debug) write(6,*) me,'cl_cmy (a) = ',cl_cmy,istep
        if (debug) write(6,*) me,'cl_cmz (a) = ',cl_cmz,istep
        do i = 1, antot
          ra(i,1,6) = ra(i,1,6) - cl_cmx
          ra(i,2,6) = ra(i,2,6) - cl_cmy
          ra(i,3,6) = ra(i,3,6) - cl_cmz
        end do
        call fixper
        call cluster_com
        if (debug) write(6,*) me,'cl_cmx (b) = ',cl_cmx,istep
        if (debug) write(6,*) me,'cl_cmy (b) = ',cl_cmy,istep
        if (debug) write(6,*) me,'cl_cmz (b) = ',cl_cmz,istep
      endif
      return
      end
c
      subroutine cluster_old_at
#include "common.fh"
c
c   store original coordinates of cluster atoms
c
      integer i, icnt, jcnt
      if (nocluster) return
      icnt = 0
      jcnt = 0
      do i = 1, antot
        if (at(i).eq.2) then
          icnt = icnt + 1
          cl_at(icnt) = i
          cl_old(icnt,1,1) = ra(i,1,1)
          cl_old(icnt,2,1) = ra(i,2,1)
          cl_old(icnt,3,1) = ra(i,3,1)
          cl_old(icnt,1,2) = ra(i,1,6)
          cl_old(icnt,2,2) = ra(i,2,6)
          cl_old(icnt,3,2) = ra(i,3,6)
          cl_old(icnt,1,3) = ra(i,1,2)
          cl_old(icnt,2,3) = ra(i,2,2)
          cl_old(icnt,3,3) = ra(i,3,2)
        endif
        if (at(i).eq.1) then
          jcnt = jcnt + 1
          sl_at(jcnt) = i
          sl_old(jcnt,1,1) = ra(i,1,1)
          sl_old(jcnt,2,1) = ra(i,2,1)
          sl_old(jcnt,3,1) = ra(i,3,1)
          sl_old(jcnt,1,2) = ra(i,1,6)
          sl_old(jcnt,2,2) = ra(i,2,6)
          sl_old(jcnt,3,2) = ra(i,3,6)
          sl_old(jcnt,1,3) = ra(i,1,2)
          sl_old(jcnt,2,3) = ra(i,2,2)
          sl_old(jcnt,3,3) = ra(i,3,2)
        endif
      end do
c
      cl_cm_old(1) = cl_cmx
      cl_cm_old(2) = cl_cmy
      cl_cm_old(3) = cl_cmz
      cl_vcm_old(1) = cl_vcmx
      cl_vcm_old(2) = cl_vcmy
      cl_vcm_old(3) = cl_vcmz
c
      do i = 1, 3
        cl_alen1_old(i) = alen1(i)
        cl_alen2_old(i) = alen2(i)
      end do
      cl_vol1_old = vol1
      cl_vol2_old = vol2
      cl_box_old(1) = xbox
      cl_box_old(2) = ybox
      cl_box_old(3) = zbox
      cl_scal1_old = scal1
      cl_scal2_old = scal2
c
      cl_tot = icnt
      sl_tot = jcnt
      return
      end
c
      subroutine cluster_check_cllsn
#include "common.fh"
c
c   This subroutine checks to see if any of the cluster particles have moved
c   beyond the cutoff radius from the center of mass of the cluster. If they
c   have, then the initial guess for the particle coordinates and velocities
c   is modified to reflect the collision with the restraining sphere.
c
      integer icnt, iat, i, j, ii, jj, iloc, imin, scndat, reset
      double precision r, r2, rx, ry, rz, dtmax, tsav
      double precision rrx, rry, rrz, vvx, vvy, vvz
      double precision vn, vx, vy, vz, fx, fy, fz, mu, v2, f2
      double precision vrdot, frdot, vfdot, a, b, c
      double precision t1, t2, tcut, comm(4), tmax, htausq
      double precision rnx, rny, rnz, vpx, vpy, vpz, vllx, vlly, vllz
      double precision t3, t4, ax, ay, az, tchk, tmin, t_est
      double precision r1x, r1y, r1z, r2x, r2y, r2z, rmax, rskin, rn
      double precision cluster_find_tau
      integer failsafe, ibuf(MD_MAXPROC), me, nproc, twohit
      logical debug
      integer iter
c
      iter = 0
      l_cllsn = .false.
      if (nocluster) then
        t_rmndr = 0.0d00
        t_done = tau
        return
      endif
      me = ga_nodeid()
      nproc = ga_nnodes()
      if (istep.ge.68365721.and.istep.le.68365723) then
        debug = .false.
      else
        debug = .false.
      endif
      rskin = 0.02d00
      tmax = t_rmndr
      if (t_rmndr.eq.tau) cllsn_isav = 0
  100 dtmax = 2.0d00*t_rmndr
      twohit = 0
      failsafe = 0
      iat = 0
      iloc = 0
      tsav = 0.0d00
      rmax = r_cluster
      if (debug) write(6,*) me,' atot = ',atot,istep
      if (debug) write(6,*) me,' antot = ',antot,istep
      if (debug) write(6,*) me, ' r_cluster = ',r_cluster,istep
      if (debug) write(6,*) me,' t_rmndr = ',tmax,istep
      if (debug) write(6,*) me,' cl_cmx(1) = ',cl_cmx,istep
      if (debug) write(6,*) me,' cl_cmy(1) = ',cl_cmy,istep
      if (debug) write(6,*) me,' cl_cmz(1) = ',cl_cmz,istep
      do i = 1, cl_tot
        ii = cl_at(i)
        rx = ra(ii,1,6) - cl_cmx
        ry = ra(ii,2,6) - cl_cmy
        rz = ra(ii,3,6) - cl_cmz
        r2 = rx**2 + ry**2 + rz**2
        r = sqrt(r2)
c        if (debug) write(6,*) me, 'ra(ii,1,6) ',ra(ii,1,6),istep
c        if (debug) write(6,*) me, 'ra(ii,2,6) ',ra(ii,2,6),istep
c        if (debug) write(6,*) me, 'ra(ii,3,6) ',ra(ii,3,6),istep
        if (r.gt.r_cluster) then
          if (r.gt.rmax) rmax = r
          if (debug) write(6,*) me,' tmax = ',tmax,istep
          if (debug) write(6,*) me,' dtmax = ',dtmax,istep
          if (debug) write(6,*) me,' ii = ',ii,istep
          if (debug) write(6,*) me,' antot = ',antot,istep
          if (debug) write(6,*) me,' r-r_cluster ',r-r_cluster,istep
          if (debug) write(6,*) me,' r_cluster = ',r_cluster,istep
          if (debug) write(6,*) me,' sqrt(r2) = ',sqrt(r2),istep
          if (debug) write(6,*) me,' ra(1) = ',ra(ii,1,6),istep
          if (debug) write(6,*) me,' ra(2) = ',ra(ii,2,6),istep
          if (debug) write(6,*) me,' ra(3) = ',ra(ii,3,6),istep
          if (debug) write(6,*) me,' cl_cmx = ',cl_cmx,istep
          if (debug) write(6,*) me,' cl_cmy = ',cl_cmy,istep
          if (debug) write(6,*) me,' cl_cmz = ',cl_cmz,istep
          if (debug) write(6,*) me,' r_cluster_old = ',r_cluster_old,
     +      istep
          if (debug) write(6,*) me,' cl_old(1) = ',cl_old(i,1,2),istep
          if (debug) write(6,*) me,' cl_old(2) = ',cl_old(i,2,2),istep
          if (debug) write(6,*) me,' cl_old(3) = ',cl_old(i,3,2),istep
          if (debug) write(6,*) me,' cl_cm_old(1) = ',cl_cm_old(1),istep
          if (debug) write(6,*) me,' cl_cm_old(2) = ',cl_cm_old(2),istep
          if (debug) write(6,*) me,' cl_cm_old(3) = ',cl_cm_old(3),istep
          if (debug) write(6,*) me,' r_old = ', sqrt(
     +      + (cl_old(i,1,2)-cl_cm_old(1))**2
     +      + (cl_old(i,2,2)-cl_cm_old(2))**2
     +      + (cl_old(i,3,2)-cl_cm_old(3))**2)
          rn = sqrt(rx**2+ry**2+rz**2)
          vn = ra(ii,1,2)*rx+ra(ii,2,2)*ry+ra(ii,3,2)*rz
          vn = vn/rn
          t_est = vn*(rn-r_cluster)
          tcut = cluster_find_tau(ii,i,cllsn_isav,tmax,.false.)
          if (debug) write(6,*) me,' tcut = ',tcut,istep
          if (debug) write(6,*) me,' t_est = ',t_est,istep
          if (debug) write(6,155) me,1,tcut,istep
  155     format(i3,' tcut at ',i1,': ',f12.4,' step: ',i8)
c
c   Check to see if collision is earlier than previously found collisions
c   and that it is also greater than zero.
c
          if (tcut.lt.dtmax.and.tcut.gt.0.0d00) then
c
c   Try to protect against numerical roundoff by checking that different
c   particle collides with sphere or that collision time is significantly
c   greater than zero if it is the same particle
c
            if (.not.(ii.eq.cllsn_isav.and.tcut.lt.1.0d-03)) then
              dtmax = tcut
              tsav = tcut
              iat = ii
              iloc = i
            else
              write(6,101) ii, tcut, istep
  101         format('Rejected collision of atom ',i8,' at time ',
     +               f16.8,' at step ',i8)
              failsafe = 1
            endif
          else if (tcut.lt.0.0d00) then
            write(6,*) ga_nodeid(),' vn ',vn,istep
            write(6,*) ga_nodeid(),' t_est ',t_est,istep
            write(6,*) ga_nodeid(),' r-r_cluster ',r-r_cluster,istep
            write(6,*) ga_nodeid(),' Returned negative value at 1',istep
            failsafe = 1
          endif
#if 0
c
c   Check for trajectories that may have gone out and then back into
c   confining sphere
c
        else if (r_cluster-r.lt.rskin.and.r_cluster-r.gt.0.0d00) then
          vn = rx*ra(ii,1,2)+ry*ra(ii,2,2)+rz*ra(ii,3,2)
          if (vn.le.0.0d00) then
            tcut = cluster_find_tau(ii,i,cllsn_isav,tmax,.true.)
            if (tcut.lt.dtmax.and.tcut.gt.0.0d00) then
              write(6,*) 'Found looping trajectory',istep
              write(6,*) 'dtmax = ',dtmax,istep
              write(6,*) 'tcut = ',tcut,istep
              write(6,*) 'ii = ',ii,istep
              write(6,*) 'cllsn_isav = ',cllsn_isav,istep
              write(6,*) 'r = ',r,istep
              write(6,*) 'r_cluster = ',r_cluster,istep
              if (.not.(ii.eq.cllsn_isav.and.tcut.lt.1.0d-10)) then
                write(6,*) 'Resetting trajectory',istep
                dtmax = tcut
                tsav = tcut
                iat = ii
                iloc = i
              else
                write(6,*)'Not resetting trajectory (loop is too short)'
              endif
            endif
          endif
#endif
        endif
      end do 
c
c   Check to see if a collision occured on another processor
c
      ibuf(1) = iat
      ibuf(2) = failsafe
      if (debug) write(6,*) me,' iat = ',iat,istep
      if (debug) write(6,*) me,' failsafe = ',failsafe,istep
      if (debug) write(6,*) me,' ibuf(1) = ',ibuf(1),istep
      if (debug) write(6,*) me,' ibuf(2) = ',ibuf(2),istep
      call ga_igop(5,ibuf,2,'+')
      if (debug) write(6,*) me,' ibuf(1) = ',ibuf(1),istep
      if (debug) write(6,*) me,' ibuf(2) = ',ibuf(2),istep
c
c   No collisions detected. Step is complete so just return.
c
      if (ibuf(1).eq.0.and.ibuf(2).eq.0) then
        t_rmndr = 0.0d00
        t_done = tau
        if (debug) write(6,*) me,' Returning with no hits ',istep
        return
      endif
c
c   If no solution found for atom outside confining sphere, then bail
c   completely and increase diameter of confining sphere so that all
c   atoms are enclosed.
c
      if (ibuf(2).gt.0) then
        call ga_dgop(4,rmax,1,'max')
        if (debug) write(6,*) me,' r_cluster (a) = ',r_cluster
        r_cluster = r_cluster + (rmax - r_cluster)*1.1d00
        if (debug) write(6,*) me,' r_cluster (b) = ',r_cluster
        if (ga_nodeid().eq.0) then
          write(6,*) ga_nodeid(),' Returned negative value at 2',istep
        endif
        failcount = failcount + 1
        t_rmndr = 0.0d00
        t_done = tau
        return
      endif
c
c   At this point, a collision has been detected on at least on
c   processor. Check to find minimum time if collisions occur on more
c   than one processor.
c
      if (iat.eq.0) then
        tsav = tau
      endif
      a = tsav
      call ga_dgop(6,a,1,'min')
c
c   Find mass of particle that has first collision with confining sphere
c
      if (a.eq.tsav) then
        b = mass(iat)
      else
        b = 0.0d00
      endif
      call ga_dgop(7,b,1,'max')
 1000 mmass = cl_mass - b
      tcut = a
      if (debug) write(6,155) me,2,tcut,istep
      if (tsav.ne.a) then
        iat = 0
        iloc = 0
        cllsn_isav = 0
        if (debug) write(6,*) me,' cllsn_isav(1): ',cllsn_isav
      else
        cllsn_isav = iat
        if (debug) write(6,*) me,' cllsn_isav(2): ',cllsn_isav
      endif
c
c  Handle degenerate case of only two atoms in cluster, which will
c  generally have two simultaneous collisions
c
      if (ctot.eq.2) then
        do i = 1, nproc
          if (i-1.eq.me) then
            ibuf(i) = iat
          else
            ibuf(i) = 0
          endif
        end do
        call ga_igop(3,ibuf,nproc,'+')
        j = 0
        iat = 0
        cllsn_isav = 0
        do i = 1, nproc
          if (ibuf(i).gt.0.and.i-1.eq.me.and.j.eq.0) then
            iat = ibuf(i)
            cllsn_isav = ibuf(i)
            j = 1
          else if (ibuf(i).gt.0) then
            j = 1
          endif
        end do
c
c  Make sure that correct mass is being used for collision
c
        if (iat.ne.0) then
          b = mass(iat)
        else
          b = 0.0d00
        endif
        call ga_dgop(6,b,1,'max')
        mmass = cl_mass - b
      endif
c
c   Now recalculate coordinates of all particles in the cluster
c
      call cluster_reset(tcut)
      call cluster_com !CHECK
c
c   Check to see if there were any trajectories that pass through confining
c   sphere twice (i.e. goes out and then comes back in). Don't bother with
c   check if there are only 2 particles in system.
c
      if (ctot.eq.2) then
        go to 2000
      endif
      reset = 0
      tchk = 0.0d00
      scndat = 0
      tmin = tcut
      tsav = tcut
      failsafe = 0
      do i = 1, cl_tot
        ii = cl_at(i)
        rx = ra(ii,1,6) - cl_cmx
        ry = ra(ii,2,6) - cl_cmy
        rz = ra(ii,3,6) - cl_cmz
        r = sqrt(rx**2+ry**2+rz**2)
        if (r.ge.r_cluster.and.ii.ne.cllsn_isav) then
          tchk = cluster_find_tau(ii,i,cllsn_isav,tcut,.false.)
          if (debug) write(6,155) me,7,tchk,istep
          if (debug) write(6,*) me, 'ii: ',ii
          if (debug) write(6,*) me, 'cllsn_isav: ',cllsn_isav
          reset = 1
          if (tchk.lt.tmin.and.tchk.gt.0.0d00) then
            if (debug) write(6,*) me,' Shorter collision found'
            tmin = tchk
            imin = ii
            if (debug) write(6,155) me,4,tchk,istep
            if (debug) write(6,155) me,5,tmin,istep
          else if (tchk.lt.0.0d00) then
            write(6,*) ga_nodeid(),' Returned negative value at 3',istep
            failsafe = 1
          else
c
c  This case should theoretically be impossible, but it may occur because of
c  roundoff
c
            write(6,*) me,' tchk greater than or equal to tmin',istep
            write(6,*) me,' tchk = ',tchk,istep
            write(6,*) me,' tmin = ',tmin,istep
            failsafe = 1
          endif
        endif
      end do
      ibuf(1) = reset
      ibuf(2) = failsafe
      call ga_igop(2,ibuf,2,'+')
c
c   No solution for tau found, so increase confining sphere and then bail
c
      if (ibuf(2).gt.0) then
        if (debug) write(6,*) me,' rmax (c) = ',rmax
        call ga_dgop(4,rmax,1,'max')
        if (debug) write(6,*) me,' rmax (d) = ',rmax
        if (debug) write(6,*) me,' r_cluster (c) = ',r_cluster
        r_cluster = r_cluster + (rmax - r_cluster)*1.1d00
        if (debug) write(6,*) me,' r_cluster (d) = ',r_cluster
        if (ga_nodeid().eq.0) then
          write(6,*) 'Bogus value of collision time at 2'
        endif
        call cluster_reset(t_rmndr)
        t_rmndr = 0.0d00
        t_done = tau
        return
      endif
      reset = ibuf(1)
c
c   Shorter collision found so recalculate positions
c
      if (reset.gt.0) then
        if (tmin.le.tsav) then
          if (debug) write(6,*) 'tmin = ',tmin,istep
          if (debug) write(6,*) 'tsav = ',tsav,istep
          a = tmin
        else
          a = tau
        endif
        call ga_dgop(2,a,1,'min')
c
        if (a.le.tsav) then
          if (a.ne.tmin) then
            iat = 0
            iloc = 0
            cllsn_isav = 0
            if (debug) write(6,*) me,'cllsn_isav(3): ',cllsn_isav
          else
            cllsn_isav = imin
            if (debug) write(6,*) me,'cllsn_isav(4): ',cllsn_isav
          endif
          tcut = a
          tsav = a   !BJP is this correct?
          if (debug) write(6,155) me,3,tcut,istep
          if (iat.gt.0) then
            b = mass(iat)
          else
            b = 0.0d00
          endif
          call ga_dgop(7,b,1,'+')
          call cluster_reset(0.0d00)
          iter = iter + 1
          if (iter.gt.100) debug = .true.
          if (iter.gt.1000) call ga_error('Too many iterations',0)
          go to 1000
        endif
      endif
c
 2000 l_cllsn = .true.
      cllsn_idx = cllsn_isav
c 
      t_rmndr = tmax - tcut
      t_done = tau - t_rmndr
      if (debug) write(6,*) me, 't_rmndr = ',t_rmndr
      if (debug) write(6,*) me, 't_done = ',t_done
c
      return
      end
c
      double precision function cluster_find_tau(iat,iloc,isav,tmax,
     +                                           forward)
#include "common.fh"
c
c   This function calculates the time at which a particle outside the confining
c   sphere would have collided with the sphere.
c
      integer iat, iloc, isav
      double precision r, r2, rx, ry, rz, dr
      double precision rrx, rry, rrz, vvx, vvy, vvz
      double precision vx, vy, vz, fx, fy, fz, v2, f2
      double precision clcut, frdot, vrdot, vfdot, a, b, c
      double precision t, t1, t2, tcut, tmax, thi, tlo, tmid
      double precision ttmin,ttmax
      integer i,j,me,iter
      logical forward, debug
c
c  A particle lies outside the cutoff. Find particle that crosses
c  cutoff first.
c
      if (istep.ge.68365721.and.istep.le.68365723) then
        debug = .false.
      else
        debug = .false.
      endif
      me = ga_nodeid()
      mmass = cl_mass - mass(iat)
c
c  Calculate center of mass and velocity of center of mass of remaining
c  particles in cluster (these are the rr and vv vectors)
c
      rrx = (cl_mass*cl_cm_old(1) - mass(iat)*cl_old(iloc,1,2))/mmass
      rry = (cl_mass*cl_cm_old(2) - mass(iat)*cl_old(iloc,2,2))/mmass
      rrz = (cl_mass*cl_cm_old(3) - mass(iat)*cl_old(iloc,3,2))/mmass
      vvx = (cl_mass*cl_vcm_old(1) - mass(iat)*cl_old(iloc,1,3))/mmass
      vvy = (cl_mass*cl_vcm_old(2) - mass(iat)*cl_old(iloc,2,3))/mmass
      vvz = (cl_mass*cl_vcm_old(3) - mass(iat)*cl_old(iloc,3,3))/mmass
c
c  The r and v vectors are the relative coordinates and velocities of
c  particle iat and the center of mass vectors. The vector f is the
c  force along this relative coordinate.
c
      rx = cl_old(iloc,1,2) - cl_cm_old(1)
      ry = cl_old(iloc,2,2) - cl_cm_old(2)
      rz = cl_old(iloc,3,2) - cl_cm_old(3)
      vx = cl_old(iloc,1,3) - cl_vcm_old(1)
      vy = cl_old(iloc,2,3) - cl_vcm_old(2)
      vz = cl_old(iloc,3,3) - cl_vcm_old(3)
      fx = ra(iat,1,3) - cl_acmx
      fy = ra(iat,2,3) - cl_acmy
      fz = ra(iat,3,3) - cl_acmz
c
      r = sqrt(rx**2 + ry**2 + rz**2)
c
c  Note that we are using acceleration instead of force,
c  so factors of reduced mass disappear
c
      vrdot = vx*rx + vy*ry + vz*rz
      frdot = fx*rx + fy*ry + fz*rz
      vfdot = vx*fx + vy*fy + vz*fz
      clcut = r_cluster
      v2 = vx**2 + vy**2 + vz**2
      a = v2+frdot
      b = 2.0d00*vrdot
c
c      c = rx**2 + ry**2 + rz**2 - clcut**2
c
c  Calculate c to minimize roundoff error from subtracting large numbers
c
      dr = r - clcut
      c = 2.0d00*r*dr+dr**2
      f2 = fx**2 + fy**2 + fz**2
c
c    Scan forwards from 0 to find a final value of t2. If no final
c    value found, then return with value of -1.
c
              if (debug) write(6,*) me, 'iat: ',iat
              if (debug) write(6,*) me, 'a: ',a
              if (debug) write(6,*) me, 'b: ',b
              if (debug) write(6,*) me, 'c: ',c
              if (debug) write(6,*) me, 'f2: ',f2
              if (debug) write(6,*) me, 'vfdot: ',vfdot
      if (forward) then
        t = tmax
        t1 = 0.0d00
        ttmin = t1
        ttmax = t
c        tlo = c + t1*(b + t1*(a+t1*(vfdot+t1*0.25d00*f2)))
        tlo = (-c-t1**2*(a+t1*(vfdot+t1*0.25d00*f2)))/b - t1
        do i = 1, 1000
          if (t.eq.tmax) then
            t2= tmax*dble(i)/1000.0d00
c            thi = c + t2*(b + t2*(a+t2*(vfdot+t2*0.25d00*f2)))
            thi = (-c-t2**2*(a+t2*(vfdot+t2*0.25d00*f2)))/b - t2
            if ((thi.ge.0.0d00.and.tlo.le.0.0d00).or.
     +          (thi.le.0.0d00.and.tlo.ge.0.0d00)) then
              t = t2
              write(6,*) ga_nodeid(),' tmin = ',ttmin,istep
              write(6,*) ga_nodeid(),' tmax = ',ttmax,istep
              write(6,*) ga_nodeid(),' tlo = ',tlo,istep
              write(6,*) ga_nodeid(),' thi = ',thi,istep
              write(6,*) ga_nodeid(),' t = ',t,istep
              go to 500
            endif
          endif
        end do
        cluster_find_tau = -1.0d00
        return
c
c    scan backwards from tmax to find an initial value of t1
c
      else
        t = 0.0d00
        t2 = tmax
c        thi = c + t2*(b + t2*(a+t2*(vfdot+t2*0.25d00*f2)))
        thi = (-c-t2**2*(a+t2*(vfdot+t2*0.25d00*f2)))/b - t2
        do i = 1, 1000
          if (t.eq.0.0d00) then
            t1= tmax - tmax*dble(i)/1000.0d00
c            tlo = c + t1*(b + t1*(a+t1*(vfdot+t1*0.25d00*f2)))
            tlo = (-c-t1**2*(a+t1*(vfdot+t1*0.25d00*f2)))/b - t1
            if ((thi.ge.0.0d00.and.tlo.le.0.0d00).or.
     +          (thi.le.0.0d00.and.tlo.ge.0.0d00)) then
              t = t1
              if (debug) write(6,*) me, 'tlo: ',tlo
              if (debug) write(6,*) me, 'thi: ',thi
              if (debug) write(6,*) me, 't: ',t
              go to 500
            endif
          endif
        end do
      endif
c
c  Use bisections to find accurate solution
c
  500 if (forward) then
        t1 = 0.0d00
        t2 = t
      else
        t1 = t
        t2 = tmax
        if (debug) write(6,*) me, 't1: ',t1
        if (debug) write(6,*) me, 't2: ',t2
      endif
      iter = 0
c  600 tlo = c + t1*(b + t1*(a+t1*(vfdot+t1*0.25d00*f2)))
c      thi = c + t2*(b + t2*(a+t2*(vfdot+t2*0.25d00*f2)))
  600 tlo = (-c-t1**2*(a+t1*(vfdot+t1*0.25d00*f2)))/b - t1
      thi = (-c-t2**2*(a+t2*(vfdot+t2*0.25d00*f2)))/b - t2
      if ((thi.ge.0.0d00.and.tlo.le.0.0d00).or.
     +    (thi.le.0.0d00.and.tlo.ge.0.0d00)) then
        t = 0.5d00*(t1+t2)
c        tmid = c + t*(b + t*(a+t*(vfdot+t*0.25d00*f2)))
        tmid = (-c-t**2*(a+t*(vfdot+t*0.25d00*f2)))/b - t
        if ((tmid.ge.0.0d00.and.thi.le.0.0d00).or.
     +      (tmid.le.0.0d00.and.thi.ge.0.0d00)) then
            t1 = t
        else if ((tmid.ge.0.0d00.and.tlo.le.0.0d00).or.
     +           (tmid.le.0.0d00.and.tlo.ge.0.0d00)) then
            t2 = t
        endif
        tcut = 0.5d00*(t1+t2)
        iter = iter + 1
        if (abs(t1-t2)/tau.gt.1.0d-14.and.iter.lt.1000) go to 600
c
c  Protect against numerical roundoff errors
c
        if (iter.ge.1000) then
          write(6,*) ga_nodeid(),' Maxed out on iterations'
        endif
        if (forward.and.(tcut.lt.1.0d-10.or.
     +      abs(tcut-tmax).lt.1.0d-10)) then
          cluster_find_tau = -1.0d00
          return
        endif
      else
        open(unit=2,file='tau.dat',status='unknown')
        do j = 0, 1000
          t = tmax*dble(j)/1000.0d00
          t1 = (-c-t**2*(a+t*(vfdot+t*0.25d00*f2)))/b - t
          tmid = c + t*(b + t*(a+t*(vfdot+t*0.25d00*f2)))
          write(2,300) t,t1, tmid
  300     format(3('    ',f16.8))
        end do
        close(2)
        if (iat.ne.isav) then
          if (abs(thi).lt.0.000001d00) then
            tcut = tmax
          else if (abs(tlo).lt.0.000001d00) then
            if (vrdot.gt.0.0d00) then
              tcut = 0.0d00
            else
              tcut = tmax
            endif
          else
            if (forward) write(6,*) 'Searching forward'
            write(6,*) 'Collision time does not converge at step ',istep
            tcut = -1.0d00
          endif
        else
          tcut = 4.0d00*tau
        endif
      endif
c
      if (forward) write(6,*) ga_nodeid(),' tcut = ',tcut
      cluster_find_tau = tcut 
      return
      end
c
      double precision function cluster_check_radius()
#include "common.fh"
      double precision rx, ry, rz, r2, cl2
      integer i
c
c   check to make sure that all cluster particles are within cluster
c   radius
c
      cl2 = 0.0d00
      do i = 1, antot
        if (at(i).eq.2) then
          rx = ra(i,1,6) - cl_cmx
          ry = ra(i,2,6) - cl_cmy
          rz = ra(i,3,6) - cl_cmz
          r2 = rx**2 + ry**2 + rz**2
          if (r2.gt.cl2) cl2 = r2
        endif
      end do
      call ga_dgop(4,cl2,1,'max')
c
      cluster_check_radius = sqrt(cl2)
      return
      end
c
      subroutine cluster_mc
#include "common.fh"
c
c  Perform Monte Carlo adjustment of volume on confining volume
c
      double precision r_old, delta_r, r2, cl2, rx, ry, rz, pi
      double precision vol_new, vol_old, x, ran1, ratio
      integer i, icnt, iat, me
      logical force_move,debug
c
      if (nocluster) return
      if (istep.ge.6932366) then
        debug = .false.
      else
        debug = .false.
      endif
      call cluster_com
      me = ga_nodeid()
      r_cluster_old = r_cluster
      r_old = r_cluster
      force_move = .false.
      delta_r = 0.2
      if (debug) write(6,*)me,' (1) r_cluster ',r_cluster,istep
      if (me.eq.0) then
        r_cluster = r_cluster + delta_r*(ran1(0)-0.5d00)
      else
        r_cluster = 0.0d00
      endif
      if (debug) write(6,*)me,' (2) r_cluster ',r_cluster,istep
      call ga_dgop(1,r_cluster,1,'+')
      if (debug) write(6,*)me,' (3) r_cluster ',r_cluster,istep
c
c  If new value of r_cluster falls outside of allowed interval then
c  it is rejected, unless old value was aready outside interval.
c
      if (r_cluster.ge.cl_upper.or.r_cluster.le.cl_lower) then
c
c  If new values move cluster radius closer to cutoffs from the outside,
c  then accept them, reject them otherwise. Before accepting move, make
c  sure that there are no illegal overlaps.
c
        force_move = .true.
        if (r_cluster.gt.r_old.and.r_cluster.gt.cl_upper) then
          r_cluster = r_old
        else if (r_cluster.lt.r_old.and.r_cluster.lt.cl_lower) then
          r_cluster = r_old
        endif
        if (r_cluster.eq.r_old) return
      endif
      if (debug) write(6,*)me,' (4) r_cluster ',r_cluster,istep
c
c  Check to see if any particles are outside new value of r_cluster.
c  If so, then new value is rejected.
c
      if (r_cluster.lt.r_old) then
        icnt = 0
        cl2 = r_cluster**2
        do i = 1, antot
          if (at(i).eq.2) then
            rx = ra(i,1,6) - cl_cmx
            ry = ra(i,2,6) - cl_cmy
            rz = ra(i,3,6) - cl_cmz
            r2 = rx**2 + ry**2 + rz**2
            if (r2.gt.cl2) then
              icnt = icnt + 1
            endif
          endif
        end do
        call ga_igop(9,icnt,1,'+')
        if (icnt.gt.0) then
          r_cluster = r_old
          return
        endif
      endif
      if (debug) write(6,*)me,' (5) r_cluster ',r_cluster,istep
      if (force_move) return
c
c  Accept with Monte Carlo probability.
c
      pi = 4.0d00*atan(1.0d00)
      vol_new = 4.0d00*pi*r_cluster**3/3.0d00
      vol_old = 4.0d00*pi*r_old**3/3.0d00
      ratio = (r_cluster/r_old)**2
c
      x = ratio*exp(-cl_prssr*(vol_new-vol_old)/mc_tmprtr)
      if (me.eq.0) then
        if (x.ge.1.0d00) then
          icnt = 1
        else 
          if (ran1(0).lt.x) then
            icnt = 1
          else
            icnt = 0
          endif
        endif
      else
        icnt = 0
      endif
      call ga_igop(1,icnt,1,'+')
      if (icnt.eq.0) then
        r_cluster = r_old
      endif
      if (debug) write(6,*)me,' (6) r_cluster ',r_cluster,istep
      return
      end
c
      subroutine cluster_binr
#include "common.fh"
c
c  Bin the current value of r_cluster
c
      double precision dr
      integer ir, isample
      if (nocluster) return
c
c  Find bin
c
      dr = r_cluster - cl_lower
      ir = int(dr/mc_step) + 1
      if (ir.le.0) ir = 1
      if (ir.gt.mcbins) ir = mcbins
      mc_cnt = mc_cnt + 1
c
c  Find slice for statistics
c
      if (istep.gt.mc_start) then
        isample = int(10.0d00*(dble(istep-mc_start)-0.5d00)
     +          / dble(nstep-mc_start)) + 1
        if (isample.gt.10) isample = 10
      else
        isample = 1
      endif
      r_cnt(isample) = r_cnt(isample) + 1
      r_distr(ir,isample) = r_distr(ir,isample) + 1
      return
      end
c
      subroutine cluster_print_binr
#include "common.fh"
c
c  Print out the distribution of r_cluster
c
      double precision rdist(MAXBINS), norm, r, normi
      double precision rdist2(MAXBINS), sum
      double precision sigr(MAXBINS),sigr2(MAXBINS)
      double precision sig2(MAXBINS),sig2_2(MAXBINS)
      double precision rdisti(MAXBINS,10),rdisti2(MAXBINS,10)
      integer itot, i, j, isample
      character*32 filename
c
      if (nocluster) return
      if (task_id.lt.10) then
        write(filename,100) task_id
      else if (task_id.ge.10.and.task_id.lt.100) then
        write(filename,101) task_id
      else if (task_id.ge.100.and.task_id.lt.1000) then
        write(filename,102) task_id
      else if (task_id.ge.1000.and.task_id.lt.10000) then
        write(filename,103) task_id
      endif
  100 format('cl.distr.',i1)
  101 format('cl.distr.',i2)
  102 format('cl.distr.',i3)
  103 format('cl.distr.',i4)
c
c  Evaluate normalized distribution and uncertainty for P(r)
c
      if (mc_cnt.lt.0) return
      norm = dble(mc_cnt)
      do i = 1, mcbins
        rdist(i) = 0.0d00
        rdist2(i) = 0.0d00
        sig2(i) = 0.0d00
        sig2_2(i) = 0.0d00
        sigr(i) = 0.0d00
        sigr2(i) = 0.0d00
        do j = 1, 10
          rdisti(i,j) = 0.0d00
          rdisti2(i,j) = 0.0d00
        end do
      end do
      do j = 1, 10
        normi = dble(r_cnt(j))
        do i = 1, mcbins
          rdist(i) = rdist(i) + dble(r_distr(i,j))/norm
          rdisti(i,j) = dble(r_distr(i,j))/normi
        end do
      end do
      do j = 1, 10
        do i = 1, mcbins
          sig2(i) = sig2(i) + rdisti(i,j)**2
        end do
      end do
      do i = 1, mcbins
        sigr(i) = (sig2(i)-10.0d00*rdist(i)**2)/9.0d00
      end do
c
c  Uncertainty at the 95% confidence level
c
      do i = 1, mcbins
        sigr(i) = 2.23d00*sqrt(sigr(i)/10.d00)
      end do
      sum = 0.0d00
      do i = 1, mcbins
        sum = sum + rdist(i)
      end do
      do i = 1, mcbins
        rdist(i) = rdist(i)/(sum*mc_step)
        sigr(i) = sigr(i)/(sum*mc_step)
      end do
c
c  Evaluate normalized distribution and uncertainty for P(V)
c
      do j = 1, 10
        normi = dble(r_cnt(j))
        do i = 1, mcbins
          r = cl_lower + (dble(i)-0.5d00)*mc_step
          rdist2(i) = rdist2(i) + dble(r_distr(i,j))/(r**2*norm)
          rdisti2(i,j) = dble(r_distr(i,j))/(r**2*normi)
        end do
      end do
      do j = 1, 10
        do i = 1, mcbins
          sig2_2(i) = sig2_2(i) + rdisti2(i,j)**2
        end do
      end do
      do i = 1, mcbins
        sigr2(i) = (sig2_2(i)-10.0d00*rdist2(i)**2)/9.0d00
      end do
c
c  Uncertainty at the 95% confidence level
c
      do i = 1, mcbins
        sigr2(i) = 2.23d00*sqrt(sigr2(i)/10.d00)
      end do
      sum = 0.0d00
      do i = 1, mcbins
        sum = sum + rdist2(i)
      end do
      do i = 1, mcbins
        rdist2(i) = rdist2(i)/(sum*mc_step)
        sigr2(i) = sigr2(i)/(sum*mc_step)
      end do
c
      if (ga_nodeid().eq.0) then
        open(unit=2,file=filename,status='unknown')
        do i = 1, mcbins
          r = cl_lower + (dble(i)-0.5d00)*mc_step
          write(2,200) r, rdist(i), rdist2(i),sigr(i),sigr2(i)
  200     format(f12.4,'    ',f16.8,'    ',f16.8,'       ',
     +           f16.8,'       'f16.8)
        end do
        close(2)
      endif
c
      return
      end
c
      subroutine cluster_clear_binr
#include "common.fh"
c
c  clear the distribution of r_cluster
c
      integer itot, i, j
c
      if (nocluster) return
      do j = 1, 10
        do i = 1, MAXBINS
          r_distr(i,j) = 0
        end do
        r_cnt(j) = 0
      end do
      mc_cnt = 0
      mc_start = istep
c
      return
      end
c
      subroutine cluster_reset_binr(iflg)
#include "common.fh"
c
c   recalculate bounds for confining sphere
c
      integer iflg,ilow,ihi,imax,i
      if (nocluster) return
      if (iflg.eq.1) then
c
c   find minimum and maximum values of distribution and reset
c   windows to just bound these values
c
        if (ga_nodeid().eq.0.and.l_rad) then
          open(unit=3,file="win1.dat",status='unknown')
          do i = 1, mcbins
            write(3,*) i, r_distr(i,1)
          end do
        endif
        ilow = 0
        ihi = 0
        do i = 1, mcbins
          if (r_distr(i,1).gt.0.and.ilow.eq.0) then
            ilow = i-1
            if (ilow.lt.1) ilow = 1
          endif
          if (ilow.gt.0.and.ihi.eq.0.and.r_distr(i,1).eq.0) then
            ihi = i
          endif
        end do
        if (ihi.eq.0) ihi = mcbins
        if (ihi.lt.mcbins) then
          cl_upper = cl_lower + mc_step*dble(ihi)
        endif
        if (ihi.gt.1) then
          cl_lower = cl_lower + mc_step*dble(ilow-1)
        endif
      else if (iflg.eq.2) then
c
c   reset bounds so that maximum is set at radius value that
c   is the maximum value of distribution
c
        ilow = 1
        ihi = 0
        imax = 0
        do i = 1, mcbins
c          if (r_distr(i,1).gt.0.and.ilow.eq.0) then
c            ilow = i-1
c            if (ilow.lt.1) ilow = 1
c          endif
          if (r_distr(i,1).gt.imax) then
            ihi = i
            imax = r_distr(i,1)
          endif
        end do
        if (ihi.eq.0) ihi = mcbins
        if (ihi.lt.mcbins) then
          cl_upper = cl_lower + mc_step*dble(ihi)
        endif
        if (ilow.gt.1) then
          cl_lower = cl_lower + mc_step*dble(ilow-1)
        endif
      endif
      cl_upper = 3.5d00
      cl_lower = 0.5d00
      cl_upper = 4.5d00
      cl_lower = 0.5d00
      mc_step = (cl_upper-cl_lower)/dble(mcbins)
      call cluster_clear_binr
      return
      end
c
      subroutine cluster_do_cllsn
#include "common.fh"
      integer iat, i, j, ii, jj, iloc, isav, scndat
      double precision r, r2, rc2, rx, ry, rz
      double precision rrx, rry, rrz, vvx, vvy, vvz
      double precision vx, vy, vz, fx, fy, fz, mu, v2, f2
      double precision comm(4), tmax, tcut, vrdot
      double precision rnx, rny, rnz, vpx, vpy, vpz, vllx, vlly, vllz
c
c  Re-evaluate all particle velocities. First, recalculate center of
c  mass and center of mass velocity of cluster at new coordinates
c  and velocities.
c
      if (nocluster) return
      cllsn_cnt = cllsn_cnt + 1
      iat = cllsn_idx
      call cluster_com
      if (iat.gt.0) then
        rrx = (cl_mass*cl_cmx - mass(iat)*ra(iat,1,6))/mmass
        rry = (cl_mass*cl_cmy - mass(iat)*ra(iat,2,6))/mmass
        rrz = (cl_mass*cl_cmz - mass(iat)*ra(iat,3,6))/mmass
        vvx = (cl_mass*cl_vcmx - mass(iat)*ra(iat,1,2))/mmass
        vvy = (cl_mass*cl_vcmy - mass(iat)*ra(iat,2,2))/mmass
        vvz = (cl_mass*cl_vcmz - mass(iat)*ra(iat,3,2))/mmass
c
c  The r and v vectors are the relative coordinates and velocities of
c  particle iat and the center of mass vectors. The vector f is the
c  force along this relative coordinate.
c
        rx = ra(iat,1,6) - rrx
        ry = ra(iat,2,6) - rry
        rz = ra(iat,3,6) - rrz
        vx = ra(iat,1,2) - vvx
        vy = ra(iat,2,2) - vvy
        vz = ra(iat,3,2) - vvz
c
c  Decompose velocity v into components that are parallel and perpendicular
c  to r. Save this information for computing post-collision velocities.
c
        r = sqrt(rx**2 + ry**2 + rz**2)
        rnx = rx/r
        rny = ry/r
        rnz = rz/r
        vrdot = vx*rnx + vy*rny + vz*rnz
        vllx = vrdot * rnx
        vlly = vrdot * rny
        vllz = vrdot * rnz
        vpx = vx - vllx
        vpy = vy - vlly
        vpz = vz - vllz
        comm(1) = (vpx-vllx-vx)/cl_mass
        comm(2) = (vpy-vlly-vy)/cl_mass
        comm(3) = (vpz-vllz-vz)/cl_mass
        comm(4) = mmass
      else
        comm(1) = 0.0d00
        comm(2) = 0.0d00
        comm(3) = 0.0d00
        comm(4) = 0.0d00
      endif
      call ga_dgop(8,comm,4,'+')
      vvx = comm(1)
      vvy = comm(2)
      vvz = comm(3)
      mmass = comm(4)
c
      do i = 1, cl_tot
        ii = cl_at(i)
        if (ii.eq.iat) then
          ra(ii,1,2) = ra(ii,1,2) + mmass*vvx
          ra(ii,2,2) = ra(ii,2,2) + mmass*vvy
          ra(ii,3,2) = ra(ii,3,2) + mmass*vvz
        else
          ra(ii,1,2) = ra(ii,1,2) - mass(ii)*vvx
          ra(ii,2,2) = ra(ii,2,2) - mass(ii)*vvy
          ra(ii,3,2) = ra(ii,3,2) - mass(ii)*vvz
        endif
      end do
      call cluster_old_at
c 
      return
      end
c
      subroutine cluster_reset(dt)
#include "common.fh"
      double precision dt, htausq, scale
      integer i, ii
c
c  Reset coordinates and velocities using old values from beginning of
c  step.
c
      htausq = 0.5d00*dt**2
      do i = 1, cl_tot
        ii = cl_at(i)
        ra(ii,1,1) = cl_old(i,1,1)+dt*cl_old(i,1,3)+htausq*ra(ii,1,3)
        ra(ii,2,1) = cl_old(i,2,1)+dt*cl_old(i,2,3)+htausq*ra(ii,2,3)
        ra(ii,3,1) = cl_old(i,3,1)+dt*cl_old(i,3,3)+htausq*ra(ii,3,3)
        ra(ii,1,6) = cl_old(i,1,2)+dt*cl_old(i,1,3)+htausq*ra(ii,1,3)
        ra(ii,2,6) = cl_old(i,2,2)+dt*cl_old(i,2,3)+htausq*ra(ii,2,3)
        ra(ii,3,6) = cl_old(i,3,2)+dt*cl_old(i,3,3)+htausq*ra(ii,3,3)
        ra(ii,1,2) = cl_old(i,1,3)+dt*ra(ii,1,3)
        ra(ii,2,2) = cl_old(i,2,3)+dt*ra(ii,2,3)
        ra(ii,3,2) = cl_old(i,3,3)+dt*ra(ii,3,3)
      end do
#if 1
      do i = 1, sl_tot
        ii = sl_at(i)
        ra(ii,1,1) = sl_old(i,1,1)+dt*sl_old(i,1,3)+htausq*ra(ii,1,3)
        ra(ii,2,1) = sl_old(i,2,1)+dt*sl_old(i,2,3)+htausq*ra(ii,2,3)
        ra(ii,3,1) = sl_old(i,3,1)+dt*sl_old(i,3,3)+htausq*ra(ii,3,3)
        ra(ii,1,6) = sl_old(i,1,2)+dt*sl_old(i,1,3)+htausq*ra(ii,1,3)
        ra(ii,2,6) = sl_old(i,2,2)+dt*sl_old(i,2,3)+htausq*ra(ii,2,3)
        ra(ii,3,6) = sl_old(i,3,2)+dt*sl_old(i,3,3)+htausq*ra(ii,3,3)
        ra(ii,1,2) = sl_old(i,1,3)+dt*ra(ii,1,3)
        ra(ii,2,2) = sl_old(i,2,3)+dt*ra(ii,2,3)
        ra(ii,3,2) = sl_old(i,3,3)+dt*ra(ii,3,3)
      end do
#endif
#if 1
c
c   resete extended Hamiltonian parameters
c
      if ((prsflg.or.ptflg).and.ipmode.eq.0) then
        vol1 = cl_vol1_old + dt * cl_vol2_old + htausq * vol3
        vol2 = cl_vol2_old + dt * vol3
        scale = vol1 / (cl_box_old(1) * cl_box_old(2) * cl_box_old(3))
        scale = exp(log(scale)/3.0d00)
        xbox = scale * cl_box_old(1)
        ybox = scale * cl_box_old(2)
        zbox = scale * cl_box_old(3)
        xbox2 = 0.5d00 * xbox
        ybox2 = 0.5d00 * ybox
        zbox2 = 0.5d00 * zbox
      endif
c
      if ((prsflg.or.ptflg).and.ipmode.eq.1) then
        do 500 i = 1, 3
          alen1(i) = cl_alen1_old(i) + dt * cl_alen2_old(i)
     +             + htausq * alen3(i)
          alen2(i) = cl_alen2_old(i) + dt * alen3(i)
  500   continue
        xbox = alen1(1)
        ybox = alen1(2)
        zbox = alen1(3)
        xbox2 = 0.5d00 * xbox
        ybox2 = 0.5d00 * ybox
        zbox2 = 0.5d00 * zbox
      endif
c
      if ((prsflg.or.ptflg).and.ipmode.eq.2) then
        do 600 i = 1, 2
          alen1(i) = cl_alen1_old(i) + dt * cl_alen2_old(i)
     +             + htausq * alen3(i)
          alen2(i) = cl_alen2_old(i) + dt * alen3(i)
  600   continue
        xbox = alen1(1)
        ybox = alen1(2)
        xbox2 = 0.5d00 * xbox
        ybox2 = 0.5d00 * ybox
      endif
c
c   predictor step for time scale
c
      if (tmpflg.or.ptflg) then
        scal1 = cl_scal1_old + dt * cl_scal2_old + htausq * scal3
        scal2 = cl_scal2_old + dt * scal3 
      endif
c
c   center of mass parameters
c
      cl_cmx = cl_cm_old(1)
      cl_cmy = cl_cm_old(2)
      cl_cmz = cl_cm_old(3)
      cl_vcmx = cl_vcm_old(1)
      cl_vcmy = cl_vcm_old(2)
      cl_vcmz = cl_vcm_old(3)
#endif
      return
      end
