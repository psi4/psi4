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
      subroutine rdpar
#include "common.fh"
      double precision x,ran1
      integer i,j
      integer me
c
c    This subroutine reads in the molecular simulation parameters
c
      me = ga_nodeid()
      if (me.eq.0) open (unit=5,file='md_lj.in',status='old')
c
c     Line 1
c
c       nstep: total number of steps in simulation
c
c       tau: time step increment
c
c       nsc: number of special control steps
c
c       iseed: random number seed (this should be a NEGATIVE
c              integer)
c
      nstep = 0
      tau = 0.0d00
      nsc = 0
      iseed = 0
      if (me.eq.0) then
        read(5,*) nstep,tau,nsc,iseed
      endif
      call ga_igop(2,nstep,1,'+')
      call ga_dgop(3,tau,1,'+')
      call ga_igop(4,nsc,1,'+')
      call ga_igop(5,iseed,1,'+')
c
c     Line 2
c
c       rcut: cutoff for use in neighbor list calculation
c
c       ilist: number of steps taken before updating
c              neighbor list
c
c       icut: cutoff parameter for potential
c
c           (0) only potential set to zero at cutoff
c
c           (1) potential and forces set to zero at cutoff
c
      rcut = 0.0d00
      ilist = 0
      icut = 0
      if (me.eq.0) then
        read(5,*) rcut,ilist,icut
      endif
      call ga_dgop(6,rcut,1,'+')
      call ga_igop(7,ilist,1,'+')
      call ga_igop(8,icut,1,'+')
c
c     Line 3
c
c       dflalg: default algorithm
c
c            1: constant energy 3 point Gear algorithm
c            2: constant pressure algorithm
c            3: constant temperature algorithm (velocity scaling)
c            4: constant temperature algorithm (Nose' dynmamics)
c            5: constant pressure and temperature algorithm
c            6: Boltzmann velocity step
c            7: adjust volume to specified target volume
c               (available only as a special step)
c            8: adjust kinetic energy based on average temperature
c               since last call to this step so that temperature
c               matches target temperature.
c               (available only as a special step)
c            9: constant pressure and temperature algorithm with
c               independent box adjustments
c            10: constant pressure and temperature algorithm with
c               independent box adjustments in only the x and y
c               dimensions
c
c       dftmp: default temperature
c
c       dfprs: default pressure
c
c       dftm: default mass for temperature
c
c       dfpm: default mass for pressure
c
      dflalg = 0
      dftmp = 0.0d00
      dfprs = 0.0d00
      dftm = 0.0d00
      dfpm = 0.0d00
      if (me.eq.0) then
        read(5,*) dflalg,dftmp,dfprs,dftm,dfpm
      endif
      call ga_igop(9,dflalg,1,'+')
      call ga_dgop(1,dftmp,1,'+')
      call ga_dgop(2,dfprs,1,'+')
      call ga_dgop(3,dftm,1,'+')
      call ga_dgop(4,dfpm,1,'+')
c
c     Line 4.1 - 4.nsc
c
c     special step control instructions
c
c       isc(i,1): begining step for special instruction
c
c       isc(i,2): final step for special instruction
c
c       isc(i,3): increment for special instruction
c
c       isc(i,4): stepping algorithm for special instruction
c                 (see documentation on dflalg)
c
c       rsc(i,1): temperature (in K) if temperature scaling or
c                 Boltzmann algorithm is specified
c
c       rsc(i,2): pressure if constant pressure
c                 algorithms are specified. If volume adjustment
c                 algorithm is specified then this is the target
c                 volume in A**3.
c
c       rsc(i,3): mass (in g A**2 /mole) for Nose temperature
c                 algorithm
c
c       rsc(i,4): mass (in g/(mole A**4)) for constant pressure
c                 algorithm
c
      do 100 i=1,nsc
        do j=1,4 
          isc(i,j) = 0.0d00
          rsc(i,j) = 0.0d00
        end do
        if (me.eq.0) then
          read(5,*) (isc(i,j),j=1,4),(rsc(i,j),j=1,4)
        endif
        do j=1,4 
          call ga_igop(4,isc(i,j),1,'+')
          call ga_dgop(5,rsc(i,j),1,'+')
        end do
  100 continue
c
c     Line 5
c
c       istart: startup format
c               1: start from single configuration
c               2: start from configuration plus velocity
c
c       istop: write-out format
c               1: single configuration
c               2: configuration plus velocity
c
      istart = 0
      istop = 0
      if (me.eq.0) then
        read(5,*) istart,istop
      endif
      call ga_igop(6,istart,1,'+')
      call ga_igop(7,istop,1,'+')
c
c     Line 6
c
c     istat: frequency to print out simulation information 
c
      istat = 0
      if (me.eq.0) then
        read(5,*) istat
      endif
      call ga_igop(8,istat,1,'+')
c
c     Line 7
c
c     equil_1: number of steps in first equilibration regime
c
c     equil_2: number of steps in second equilibration regime
c
c     equil_3: number of steps in third equilibration regime
c
      equil_1 = 0
      equil_2 = 0
      equil_3 = 0
      if (me.eq.0) then
        read(5,*) equil_1, equil_2, equil_3
      endif
      call ga_igop(9,equil_1,1,'+')
      call ga_igop(1,equil_2,1,'+')
      call ga_igop(2,equil_3,1,'+')
      window_1 = int(0.2d00*dble(equil_3-equil_2))
      window_1 = equil_2 + window_1
      window_2 = int(0.4d00*dble(equil_3-equil_2))
      window_2 = window_1 + window_2
c
c     Line 7
c
c     cl_prssr: pressure that is applied to cluster
c
c     mc_tmprtr: temperature used in Monte Carlo steps
c
      cl_prssr = 0.0d00
      mc_tmprtr = 0.0d00
      if (me.eq.0) then
        read(5,*) cl_prssr, mc_tmprtr
      endif
      call ga_dgop(3,cl_prssr,1,'+')
      call ga_dgop(4,mc_tmprtr,1,'+')
c
c     Line 8
c
c     mcfreq: frequency to apply Monte Carlo step on radius
c
c     mcbins: number of bins to accumulate radius in
c
      mcfreq = 0
      mcbins = 0
      if (me.eq.0) then
        read(5,*) mcfreq, mcbins
      endif
      call ga_igop(7,mcfreq,1,'+')
      call ga_igop(8,mcbins,1,'+')
c
      if (me.eq.0) close(5)
c
c    Write out simulation information to output file
c
      if (me.eq.0.and.l_stdio) then
        write(6,1300) nstep
        write(6,1400) 1000.0 * tau
        write(6,4325) rcut
        write(6,4350) ilist
        if (icut.eq.0) then
          write(6,4400)
        else
          write(6,4420)
        endif
        write(6,1500) dflalg
        write(6,1600)
        write(6,1700) dftmp
        write(6,1800) dfprs
        write(6,1900) dftm
        write(6,2000) dfpm
        write(6,2100)
        do 300 i = 1, nsc
          if (i.eq.1) then
             write(6,2200)
          else
             if (isc(i-1,4).eq.7) write(6,2200)
          endif
          if (isc(i,4).ne.7) then
            write(6,2300) (isc(i,j),j=1,4),(rsc(i,j),j=1,4)
          else
            write(6,2400) 
            write(6,2500) (isc(i,j),j=1,4)
            write(6,2600) rsc(i,1)
            write(6,2700) rsc(i,2)
          endif
  300   continue
        write(6,3400) istat
        write(6,3800) mcfreq
        write(6,3500) equil_1
        write(6,3600) equil_2
        write(6,3700) equil_3
        if (istart.eq.1) then
          write(6,2800)
        else
          write(6,2900)
        endif
        if (istop.eq.1) then
          write(6,3000)
        else
          write(6,3100)
        endif
        write(6,3200) iseed
        write(6,4500) ga_pgroup_nodeid(ga_pgroup_get_world())
      endif
      iseed = iseed - 10*me
      x = ran1(iseed)
      r_cluster = 0.0d00
      mc_cnt = 0
c
      return
 1300 format('Total number of steps in simulation                :',i9)
 1400 format('Time step interval                                 :'
     +       ,f16.6)
 1500 format('Default algorithm                                  :',i9)
 1600 format('     Default algorithm parameters ')
 1700 format('        Default temperature                        :'
     +       ,f16.6)
 1800 format('        Default pressure                           :'
     +       ,f16.6)
 1900 format('        Default temperature mass                   :'
     +       ,f16.6)
 2000 format('        Default pressure mass                      :'
     +       ,f16.6)
 2100 format('Special step instructions:')
 2200 format('   begin     end  incrmt algrthm    T       P      ',
     +'T mass  P mass  ')
 2300 format(4i8,2f8.3,1pe9.2,1pe9.2)
 2400 format('   begin     end  incrmt algrthm')
 2500 format(4i8)
 2600 format('     Temperature for volume adjustment             :'
     +       ,f16.6)
 2700 format('     Target volume for volume adjustment           :'
     +       ,f16.6)
 2800 format('Initial configuration is coordinates only')
 2900 format('Initial configuration is coordinates and velocities')
 3000 format('Final configuration is coordinates only')
 3100 format('Final configuration is coordinates and velocities')
 3200 format('Random number seed                                 :',i9)
 3400 format('Frequency to print out simulation information      :',i9)
 3500 format('Number of steps in first equilibration regime      :',i9)
 3600 format('Number of steps in second equilibration regime     :',i9)
 3700 format('Number of steps in third equilibration regime      :',i9)
 3800 format('Frequency to take Monte Carlo step                 :',i9)
 4325 format('Cutoff distance                                    :',
     +       f16.6)
 4350 format('Neighbor list update frequency                     :',i9)
 4400 format('Potential set to zero at cutoff distance')
 4420 format('Potential and forces set to zero at cutoff distance')
 4500 format('Writing from node                                  :',i9)
      end
c
      subroutine rdcfg
#include "common.fh"
c
      double precision scl,gasdev
      integer inode,icnt,i,j,k,pnum,ilast,ifirst,me
      character*32 filename
c
c   This subroutine reads in the initial coordinates of all atoms in
c   the system. First determine which coordinates each processor
c   should read.
c
      me = ga_nodeid()
      pnum = ga_nnodes()
c
c  clean up everything
c
      atot = 0
      ctot = 0
      antot = 0
      do k = 1, 8
        do j = 1, 3
          do i = 1, MAXAT
            ra(i,j,k) = 0.0d00
          end do
        end do
      end do
      do i = 1, MAXAT
        mass(i) = 0.0d00
        at(i) = 0
        aidx(i) = 0
c
        xcrd(i) = 0.0d00
        ycrd(i) = 0.0d00
        zcrd(i) = 0.0d00
        xfrc(i) = 0.0d00
        yfrc(i) = 0.0d00
        zfrc(i) = 0.0d00
        xacc(i) = 0.0d00
        yacc(i) = 0.0d00
        zacc(i) = 0.0d00
        mbuf(i) = 0.0d00
        bat(i) = 0
        bidx(i) = 0
      end do
      btot = 0
c
      if (task_id.lt.10) then
        write(filename,100) task_id
      else if (task_id.ge.10.and.task_id.lt.100) then
        write(filename,101) task_id
      else if (task_id.ge.100.and.task_id.lt.1000) then
        write(filename,102) task_id
      else if (task_id.ge.1000.and.task_id.lt.10000) then
        write(filename,103) task_id
      endif
  100 format('md.cfg',i1)
  101 format('md.cfg',i2)
  102 format('md.cfg',i3)
  103 format('md.cfg',i4)
      atot = 0
      if (me.eq.0) then
        open (unit=2,file=filename,status='old')
        read(2,*) atot
        close(2)
      endif
      call ga_igop(5,atot,1,'+')
      ilast = nint(dble((me+1)*atot)/dble(pnum))
      ifirst = nint(dble(me*atot)/dble(pnum))
      ifirst = ifirst + 1
c
c   read in portions of MD configuration on each node
c
      do inode = 0, pnum - 1
        call ga_sync()
        if (me.eq.inode) then
          open (unit=2,file=filename,status='old')
c
          read(2,*) atot
          read(2,*) xbox,ybox,zbox
          xbox2 = xbox/2.0d00
          ybox2 = ybox/2.0d00
          zbox2 = zbox/2.0d00
c
          cl_lower = 0.0d00
          if (xbox.lt.ybox.and.xbox.lt.zbox) then
            cl_upper = xbox2
          else if (ybox.lt.xbox.and.ybox.lt.zbox) then
            cl_upper = ybox2
          else
            cl_upper = zbox2
          endif
          mc_step = (cl_upper-cl_lower)/dble(mcbins)
c
          icnt = 0
          if (istart.eq.1.or.(.not.l_oldcfg)) then
            do i = 1, ilast
              if (i.lt.ifirst) then
                read(2,*)
              else
                icnt = icnt + 1
                read(2,*) at(icnt),(ra(icnt,j,1),j=1,3)
                if (at(icnt).eq.2) ctot = ctot + 1
c
c    generate atomic velocities from a Maxwell-Boltzmann distribution
c
                do j = 1, 3
                  scl = sqrt(dftmp / amass(at(icnt)))
                  ra(icnt,j,2) = scl * gasdev(0)
                end do
                aidx(icnt) = i
              endif
            end do
          else
            do i = 1, ilast
              if (i.lt.ifirst) then
                read(2,*)
              else
                icnt = icnt + 1
                read(2,*) at(icnt),(ra(icnt,j,1),j=1,3),
     +                    (ra(icnt,j,2),j=1,3)
                if (at(icnt).eq.2) ctot = ctot + 1
                aidx(icnt) = i
              endif
            end do
          endif
          close(2)
        endif
      end do
      antot = icnt
      call ga_igop(9,ctot,1,'+')
c
      if (me.eq.0.and.l_stdio) then
        write(6,1100) atot
        write(6,1000) xbox,ybox,zbox
      endif
c
c   initialize absolute coordinates
c
      do j = 1, 3
        do i = 1, antot
          ra(i,j,6) = ra(i,j,1)
        end do
      end do
c
      call fixper
c
      return
 1000 format('The initial box size parameters are '/
     +       '      x dimension = ',f9.4,' A'/
     +       '      y dimension = ',f9.4,' A'/
     +       '      z dimension = ',f9.4,' A') 
 1100 format('The total number of atoms in simulation            :',i9)
      end
c
      subroutine atomin
#include "common.fh"
c
      double precision epsln,sigma
      integer i,j,ipairs,ip,ktrnc
      integer me
      logical ichk(50,50)
c
c    This subroutine reads in atomic parameters from a parameter
c    file. These parameters will be used to construct the
c    system potential function.
c
      me = ga_nodeid()
      do i = 1, 50
        do j = 1, 50
          ichk(i,j) = .false.
        end do
      end do
c
      if (me.eq.0) then
        open(unit=2,file='atom.inp',status='old')
      endif
c
      if (me.eq.0) then
        read (2,*) amass(1)
        read (2,*) amass(2)
        read (2,*) amass(3)
        read (2,*) amass(4)
        atnum = 4
      else
        amass(1) = 0.0d00
        amass(2) = 0.0d00
        amass(3) = 0.0d00
        amass(4) = 0.0d00
        atnum = 0
      endif
      call ga_dgop(1,amass(1),4,'+')
      call ga_igop(2,atnum,1,'+')
c
c   read in atomic parameters
c
      if (me.eq.0) then
        read(2,*) ipairs
      else
        ipairs = 0
      endif
      call ga_igop(3,ipairs,1,'+')
c
      do 100 ip = 1, ipairs
        i = 0
        j = 0
        ktrnc = 0
        epsln = 0.0d00
        sigma = 0.0d00
        if (me.eq.0) then
          read (2,*) i,j,ktrnc,epsln,sigma
        endif
        call ga_igop(1,i,1,'+')
        call ga_igop(2,j,1,'+')
        call ga_igop(3,ktrnc,1,'+')
        call ga_dgop(4,epsln,1,'+')
        call ga_dgop(5,sigma,1,'+')
        ichk(i,j) = .true.
        if (i.ne.j) ichk(j,i) = .true.
        e12(i,j) = 4.0d00*epsln*sigma**12
        e6(i,j) = 4.0d00*epsln*sigma**6
        if (i.ne.j) e12(j,i) = e12(i,j)
        if (i.ne.j) e6(j,i) = e6(i,j)
        if (ktrnc.eq.1) then
          acut(i,j) = exp(log(2.0d00)/6.0d00)*sigma
          acut(j,i) = exp(log(2.0d00)/6.0d00)*sigma
        else
          acut(i,j) = rcut * sigma
          acut(j,i) = rcut * sigma
        endif
        acut2(i,j) = acut(i,j)**2
        acut2(j,i) = acut(j,i)**2
  100 continue
c
      if (me.eq.0) close(2)
      return
      end
