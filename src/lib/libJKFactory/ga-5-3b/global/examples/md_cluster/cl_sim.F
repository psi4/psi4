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
      subroutine cl_sim(natom1,natom2,task)
#include "common.fh"
c
      double precision tbeg,wraptime
      integer i,j,nprocs,ndim,me,task
      integer heap, stack, natom1, natom2
      integer blocks(2),map0(MD_MAXPROC)
      integer map1(MD_MAXPROC),map2(MD_MAXPROC+1)
      character*32 filename
      logical status
c
c   This is the main calling program for the Molecular Dynamics
c   calculation.
c
c   set io control flags
c
      l_stdio = .true.
      l_rad = .false.
      l_rst = .false.
      l_step = .true.
      l_oldcfg = .false.
c
      do i = 1, MAXTIM
        tmstat(i) = 0.0d00
      end do
      tbeg = wraptime()
c
c   create communication arrays
c
      nprocs = ga_nnodes()
      me = ga_nodeid()
c
c   create map and block arrays to use irregular distributions
c
      map2(1) = 1
      blocks(1) = 1
      blocks(2) = nprocs
      do i = 1, nprocs
        map0(i) = i
        map1(i) = (i-1)*MAXAT + 1
        map2(i+1) = (i-1)*MAXAT + 1
      end do
c
      ndim = 1
      g_size = ga_create_handle()
      call ga_set_data(g_size,ndim,nprocs,MT_INT)
      call ga_set_irreg_distr(g_size,map0,blocks(2))
c      call ga_set_chunk(g_size,ndim)
      status = ga_allocate(g_size)
c
      g_coords = ga_create_handle()
      ndim = 2
      gcoords_hi(1) = 3
      gcoords_hi(2) = MAXAT*nprocs
      call ga_set_data(g_coords,ndim,gcoords_hi,MT_DBL)
      call ga_set_irreg_distr(g_coords,map2,blocks)
c      gcoords_hi(2) = MAXAT
c      call ga_set_chunk(g_coords,gcoords_hi)
      status = ga_allocate(g_coords)
c
      g_frc = ga_create_handle()
      ndim = 2
      gfrc_hi(1) = 3
      gfrc_hi(2) = MAXAT*nprocs
      call ga_set_data(g_frc,ndim,gfrc_hi,MT_DBL)
      call ga_set_irreg_distr(g_frc,map2,blocks)
c      gfrc_hi(2) = MAXAT
c      call ga_set_chunk(g_frc,gfrc_hi)
      status = ga_allocate(g_frc)
c
      g_acc = ga_create_handle()
      ndim = 2
      gacc_hi(1) = 3
      gacc_hi(2) = MAXAT*nprocs
      call ga_set_data(g_acc,ndim,gacc_hi,MT_DBL)
      call ga_set_irreg_distr(g_acc,map2,blocks)
c      gacc_hi(2) = MAXAT
c      call ga_set_chunk(g_acc,gacc_hi)
      status = ga_allocate(g_acc)
c
      g_index = ga_create_handle()
      ndim = 1
      gindex_hi = MAXAT*nprocs
      call ga_set_data(g_index,ndim,gindex_hi,MT_INT)
      call ga_set_irreg_distr(g_index,map1,blocks(2))
c      gindex_hi = MAXAT
c      call ga_set_chunk(g_index,gindex_hi)
      status = ga_allocate(g_index)
c
      g_rvec = ga_create_handle()
      ndim = 1
      grvec_hi = MAXAT*nprocs
      call ga_set_data(g_rvec,ndim,grvec_hi,MT_DBL)
      call ga_set_irreg_distr(g_rvec,map1,blocks(2))
c      grvec_hi = MAXAT
c      call ga_set_chunk(g_rvec,grvec_hi)
      status = ga_allocate(g_rvec)
c
      g_iat = ga_create_handle()
      ndim = 2
      giat_hi(1) = 2
      giat_hi(2) = MAXAT*nprocs
      call ga_set_data(g_iat,ndim,giat_hi,MT_INT)
      call ga_set_irreg_distr(g_iat,map2,blocks)
c      giat_hi(2) = MAXAT
c      call ga_set_chunk(g_iat,giat_hi)
      status = ga_allocate(g_iat)
c
c   set up arrays for describing data held by each processor
c
      do i = 0, nprocs-1
        gsize_lo(i) = 0
        gindex_lo(i) = 0
        grvec_lo(i) = 0
        do j = 1, 2
          gcoords_lo(j,i) = 0
          gfrc_lo(j,i) = 0
          gacc_lo(j,i) = 0
          giat_lo(j,i) = 0
        end do
      end do
      call nga_distribution(g_size,me,gsize_lo(me),gsize_hi)
      call ga_igop(1,gsize_lo,nprocs,'+')
      call nga_distribution(g_coords,me,gcoords_lo(1,me),gcoords_hi)
      call ga_igop(2,gcoords_lo,2*nprocs,'+')
      call nga_distribution(g_frc,me,gfrc_lo(1,me),gfrc_hi)
      call ga_igop(3,gfrc_lo,2*nprocs,'+')
      call nga_distribution(g_acc,me,gacc_lo(1,me),gacc_hi)
      call ga_igop(4,gacc_lo,2*nprocs,'+')
      call nga_distribution(g_index,me,gindex_lo(me),gindex_hi)
      call ga_igop(5,gindex_lo,nprocs,'+')
      call nga_distribution(g_rvec,me,grvec_lo(me),grvec_hi)
      call ga_igop(6,grvec_lo,nprocs,'+')
      call nga_distribution(g_iat,me,giat_lo(1,me),giat_hi)
      call ga_igop(7,giat_lo,2*nprocs,'+')
c
c   Read in variables and initial configuration
c
      task_id = task
      if (me.eq.0.and.l_stdio) then
        if (task.lt.10) then
          write(filename,100) task
        else if (task.ge.10.and.task.lt.100) then
          write(filename,101) task
        else if (task.ge.100.and.task.lt.1000) then
          write(filename,102) task
        else if (task.ge.1000.and.task.lt.10000) then
          write(filename,103) task
        endif
  100   format('md.out',i1)
  101   format('md.out',i2)
  102   format('md.out',i3)
  103   format('md.out',i4)
        open(unit=6,file=filename,status='unknown')
      endif
      if (me.eq.0.and.l_rad)
     +  open(unit=7,file='rad.dat',status='unknown')
      call rdpar
      call atomin
      r_cluster = 0.0d00
      if (.not.l_oldcfg) then
        if (me.eq.0) call newcfg(natom1,natom2,task,r_cluster,rcut)
      endif
      call ga_dgop(3,r_cluster,1,'+')
      call rdcfg
      call shuffle
c
c   Initialize MD calculations
c
      call mdinit
c
c   Begin MD calculation
c
      call mdstep
      tmstat(1) = tmstat(1) + wraptime() - tbeg
      call mdout
c
c   Close out calculation
c
      if (l_stdio) close(6)
      if (l_rad) close(7)
      status = ga_destroy(g_size)
      status = ga_destroy(g_coords)
      status = ga_destroy(g_frc)
      status = ga_destroy(g_acc)
      status = ga_destroy(g_index)
      status = ga_destroy(g_rvec)
      status = ga_destroy(g_iat)
c
      return
      end
