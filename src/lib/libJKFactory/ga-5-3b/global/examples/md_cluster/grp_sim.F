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
#define USE_SUBGROUPS 1
      program grp_sim
#include "common.fh"
c
      integer MAXTASKS
      parameter (MAXTASKS=2000)
      integer i,j,me,icnt,ndim,one
      integer heap, stack, group_size, nprocs, ngroups, my_grp
      integer max_task, itask, jtask, natom1, natom2
      integer group_list(2000), proc_list(2000)
      double precision tbeg, wraptime, elapsed
      double precision task_time(MAXTASKS),task_beg
      double precision proc_time(MD_MAXPROC), delta_t
      double precision maxtime, mintime,mingtime,maxgtime
      logical status
c
c   This is the main calling program for the Molecular Dynamics
c   calculation.
c
c
c   Initialize message passing
c
#ifdef MPI
      integer ierr
      call mpi_init(ierr)
#else
      call pbeginf
#endif
      tbeg = wraptime()
c
      call ga_initialize()
c
c   Initialize global arrays
c
      heap = 2000000
      stack = 2000000
c      heap = 2000
c      stack = 2000
      if (.not.ma_init(MT_DBL, stack, heap))
     +  call ga_error("ma_init failed",-1)
c
c   Create process groups
c
#if USE_SUBGROUPS
      group_size = 1
#else
      group_size = ga_nnodes()
#endif
      max_task = 8
      one = 1
      nprocs = ga_nnodes() 
      me = ga_nodeid()
      my_grp = (me-mod(me,group_size))/group_size
      do i = 1, min(max_task, MAXTASKS)
        task_time(i) = 0.0d00
      end do

#if USE_SUBGROUPS
c      write(6,101) my_grp,me
  101 format('My group is ',i2,' on proc ',i3)
      ngroups = nprocs/group_size
      do i = 1, min(nprocs, MD_MAXPROC)
        proc_time(i) = 0.0d00
      end do
c      write(6,102) ngroups,me
  102 format('Ngroups is  ',i2,' on proc ',i3)
      icnt = 0
      do i = 1, ngroups
        do j = 1, group_size
          proc_list(j) = icnt
          icnt = icnt + 1
        end do
        group_list(i) = ga_pgroup_create(proc_list,group_size)
      end do
#endif
c
c  Create global array to use for master-worker algorithm
c
      g_counter = ga_create_handle()
      ndim = 1
      call ga_set_data(g_counter,ndim,ndim,MT_INT)
      status = ga_allocate(g_counter)
      call ga_zero(g_counter)
c      write(6,103) me
  103 format('Created counter array on   ',i3)
#if USE_SUBGROUPS
c      write(6,104) me,group_list(my_grp+1)
  104 format('Default group on ',i3,' is ',i2)
      call ga_pgroup_set_default(group_list(my_grp+1))
c      write(6,105) me
  105 format('Set default group on   ',i3)
#endif
  100 if (ga_nodeid().eq.0) then
        itask = nga_read_inc(g_counter,one,one)
      else
        itask = 0
      endif
      call ga_igop(1,itask,1,'+')
      if (itask.lt.max_task) then
c        write(6,106) itask,me
  106   format('Executing task ',i3,' on proc ',i3)
        natom1 = 400
        natom1 = 0
        jtask = max_task - 1 - itask
        natom2 = jtask*5 + 25
        natom2 = jtask + 2
        task_beg = wraptime()
        call cl_sim(natom1,natom2,jtask)
        delta_t = wraptime() - task_beg
        if (itask.lt.MAXTASKS) task_time(itask+1) = delta_t
#if USE_SUBGROUPS
        i = ga_pgroup_nodeid(ga_pgroup_get_world())
        i = i/group_size+1
        if (i.le.MD_MAXPROC) proc_time(i) = delta_t
#endif
        go to 100
      endif
#if USE_SUBGROUPS
      call ga_pgroup_set_default(ga_pgroup_get_world())
#endif
      call ga_dgop(3,task_time,MAXTASKS,'+')
      task_time(1) = task_time(1)/dble(group_size)
      mintime = task_time(1)
      maxtime = task_time(1)
      do i = 2, min(max_task, MAXTASKS)
        task_time(i) = task_time(i)/dble(group_size)
        if (task_time(i).gt.maxtime) maxtime=task_time(i)
        if (task_time(i).lt.mintime) mintime=task_time(i)
      end do
      call ga_dgop(4,proc_time,MD_MAXPROC,'+')
      proc_time(1) = proc_time(1)/dble(group_size)
      mingtime = proc_time(1)
      maxgtime = proc_time(1)
      do i = 2, min(ngroups, MD_MAXPROC)
        proc_time(i) = proc_time(i)/dble(group_size)
        if (proc_time(i).gt.maxgtime) maxgtime=proc_time(i)
        if (proc_time(i).lt.mingtime) mingtime=proc_time(i)
      end do
      elapsed = wraptime()-tbeg
      call ga_dgop(2,elapsed,1,'+')
      elapsed = elapsed/dble(nprocs)
      if (me.eq.0) then
        do i = 1, min(max_task, MAXTASKS)
          write(6,300) i,task_time(i)
        end do
        do i = 1, min(ngroups, MD_MAXPROC)
          write(6,301) i,proc_time(i)
        end do
        write(6,201) max_task
        write(6,202) group_size
        write(6,203) mintime
        write(6,204) maxtime
        write(6,205) mingtime
        write(6,206) maxgtime
        write(6,200) elapsed
      endif
  200 format('Elapsed time for simulation  : ',f16.4)
  201 format('Number of tasks              : ',i12)
  202 format('Number of processors in group: ',i12)
  203 format('Minimum time for task        : ',f16.4)
  204 format('Maximum time for task        : ',f16.4)
  205 format('Minimum time for group       : ',f16.4)
  206 format('Maximum time for group       : ',f16.4)
  300 format('Time for task[',i3,']           : ',f16.4)
  301 format('Time for group[',i3,']          : ',f16.4)
      call ga_terminate()
#ifdef MPI
c      write(6,*) 'Calling mpi_finalize'
      call mpi_finalize()
#else
      call pend
#endif
c      write(6,*) 'Called mpi_finalize'
c
c      close(6)
      stop
      end
