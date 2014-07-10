#if HAVE_CONFIG_H
#   include "config.fh"
#endif
C
C     Test the minval, minloc, maxval, maxloc, and enum functions in GA.
C
      program main
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer heap, stack, fudge, ma_heap, me, nproc
      logical status
      parameter (heap=200*200*4, fudge=100, stack=200*200)
c     
c***  Intitialize a message passing library
c
#include "mp3.fh"
c
c***  Initialize GA
c
c     There are 2 choices: ga_initialize or ga_initialize_ltd.
c     In the first case, there is no explicit limit on memory usage.
c     In the second, user can set limit (per processor) in bytes.
c
      call ga_initialize()
      nproc = ga_nnodes()
      me = ga_nodeid()
c     we can also use GA_set_memory_limit BEFORE first ga_create call
c
      if(ga_nodeid().eq.0)then
         print *,' GA initialized '
         call ffflush(6)
      endif
c
c***  Initialize the MA package
c     MA must be initialized before any global array is allocated
c
      ma_heap = heap/nproc + fudge
      status = ma_init(MT_DCPL, stack, ma_heap)
      if (.not. status) call ga_error('ma_init failed',-1) 
c
      if(me.eq.0)then
        print *, 'using ', nproc, ' process(es)'
        call ffflush(6)
      endif
c
      call test_nga_pack()                               ! Test PACK/UNPACK
      call test_nga_scan()                               ! Test SCAN_COPY/SCAN_ADD

c     if(me.eq.0) call ga_print_stats()
c
c***  Tidy up the GA package
c
      call ga_terminate()
c
c***  Tidy up after message-passing library
c
      call MP_FINALIZE()
c
      end

      subroutine test_nga_pack()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
C
      integer g_ielm1, g_result1                           ! handles to global arrays
      integer g_ielm2, g_sbit2, g_src2, g_sink2, g_result2 ! handles to global arrays
      integer g_irowmat, g_icolmat                         ! handles to global arrays
      integer g_velm4, g_velm5                             ! handles to global arrays
C
      integer ilo, ihi, num
      integer ilo4, ihi4, num4
      integer i, it, count
      integer me, nproc                       ! my processor & number of procs 
      integer ndim,dims(1),chunk(1)
C
      integer nelements_max, nen_max
      parameter (nelements_max=100000, nen_max=4)
      integer ipack(nelements_max), 
     *        jpack(nelements_max), 
     *        kpack(nelements_max), 
     *        ielm1(nelements_max), 
     *        ielm4(nen_max*nelements_max)
      double precision velm4(nen_max*nelements_max), 
     *                 velm5(nen_max*nelements_max)
c
      integer ilocmax, ilocmin
      integer ilocmax_ga, ilocmin_ga
      integer imax_ga, imin_ga
C
      integer nelements, nen
      integer isum, icount, itoff
      double precision result
      integer i1, i2
C
c
c***  check parallel environment
      me = ga_nodeid()
      nproc = ga_nnodes()
      if(me.eq.0)then
        print *, 'testing pack/unpack'
        call ffflush(6)
      endif
c
      nelements=20000
      nen=4
c
c***  create a global 1-D array
      ndim=1
      dims(1)=nen*nelements
      chunk(1)=min(1,dims(1)/nproc)
      if (.not. nga_create(MT_INT, ndim, dims, 'ielm1', chunk, 
     $                     g_ielm1))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'result1', chunk, 
     $                     g_result1))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'sbit2', chunk, 
     $                     g_sbit2))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'src2', chunk, 
     $                     g_src2))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'sink2', chunk, 
     $                     g_sink2))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'result2', chunk, 
     $                     g_result2))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'irowmat', chunk, 
     $                     g_irowmat))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'icolmat', chunk, 
     $                     g_icolmat))
     $     call ga_error(' ga_create failed ',0)
C
      if (.not. nga_create(MT_DBL, ndim, dims, 'velm4', chunk, 
     $                     g_velm4))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_DBL, ndim, dims, 'velm5', chunk, 
     $                     g_velm5))
     $     call ga_error(' ga_create failed ',0)
c
c     Enumerate the sbit arrary to get a sequential vector.
      call ga_fill(g_ielm1,0)
      call ga_patch_enum(g_ielm1,1,nelements,1,1)
c     *** ielm1: 1 2 3 4 5 ... nelements ... 0 0 0 0 0
c
c***  compute local ilo, ihi, num for each processor
      call nga_distribution(g_ielm1,me,ilo,ihi,chunk)
      do i=0,nproc-1
         if(me.eq.i) then
            write(6,'(a,i4,a,i8,a,i8)') 'me: ',me,
     +            ' lo: ', ilo,' hi: ',ihi
            call ffflush(6)
         endif
         call ga_sync()
      enddo
      num=ihi-ilo+1
      if(ihi.le.0) num=0
      if(ihi.gt.nelements) then
         ihi=nelements
         num=max(0,ihi-ilo+1)
      endif
      if(num.gt.nelements_max) then
         print *,"Too many elements "
         goto 9999
      endif
cbjp      print *, 'me=',me, num,ilo,ihi
      do i=0,nproc-1
         if(me.eq.i) then
            write(6,'(a,i4,a,i8,a,i8,a,i8)') 'me: ',me,
     +            ' num: ',num, ' lo: ', ilo,' hi: ',ihi
            call ffflush(6)
         endif
         call ga_sync()
      enddo
c
c***  scatter some values into the global array
      call ga_fill(g_sbit2,0)
      call ga_fill(g_src2,0)
      if(num.gt.0) call nga_get(g_ielm1,ilo,ihi,ielm1,1)
      do i=ilo,ihi
         ipack(i-ilo+1)=1+nen*(ielm1(i-ilo+1)-1)
         jpack(i-ilo+1)=1
         kpack(i-ilo+1)=ielm1(i-ilo+1)
      enddo
      if(num.gt.0) call NGA_scatter(g_sbit2,jpack,ipack,num)
      if(num.gt.0) call NGA_scatter(g_src2 ,kpack,ipack,num)
c     *** sbit2: 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 .....
c     ***  src2: 1 0 0 0 2 0 0 0 3 0 0 0 4 0 0 0 5 0 0 0 6 .....
c
c     ..................................................................
c     TEST THE "PACK" FUNCTION.
c
      itoff=nen*nelements
      call ga_fill(g_sink2,0)
      call ga_pack(g_src2, g_sink2,g_sbit2,1,itoff,icount)
c     *** sbit2: 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 .....
c     ***  src2: 1 0 0 0 2 0 0 0 3 0 0 0 4 0 0 0 5 0 0 0 6 .....
c     *** sink2: 1 2 3 4 5 ... icount ... 0 0 0 0 0 0 0 0 0 ..... 0
c         [NOTE: icount should equal nelements.]
      if(icount.ne.nelements) then
         print *,"Wrong number of elements from PACK operation: ",
     *           "count = ",icount," should be = ",nelements
      endif
      call ga_add(1,g_ielm1,-1,g_sink2,g_result1)
c     ***   ielm1: 1 2 3 4 5 ... nelements ... 0 0 0 0 0 0 0 0
c     ***   sink2: 1 2 3 4 5 ... nelements ... 0 0 0 0 0 0 0 0
c     *** result1: 0 0 0 0 0 0 0 ... 0  (result1=ielm1-sink2)
c***  Find the maximum value and the index of the maximum value
      call nga_select_elem(g_result1,'max',imax_ga,ilocmax_ga)
c***  Find the minimum value and the index of the minimum value
      call nga_select_elem(g_result1,'min',imin_ga,ilocmin_ga)
      if(imax_ga.eq.0.and.imin_ga.eq.0) then
         if (me.eq.0) then
           print *,"Pack successful"
         endif
      else
         call ga_error("Pack unsuccessful",0)
      endif
      call ga_sync()
C
C     ..................................................................
C     TEST THE "UNPACK" FUNCTION.
C
      itoff=nen*nelements
      call ga_fill(g_src2,0)
      call ga_unpack(g_sink2,g_src2,g_sbit2,1,itoff,count)
C     *** sink2: 1 2 3 4 5 ... icount ... 0 0 0 0 0 0 0 0 0 ..... 0
C     *** sbit2: 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 .....
C     ***  src2: 1 0 0 0 2 0 0 0 3 0 0 0 4 0 0 0 5 0 0 0 6 .....
      call ga_fill(g_sink2,0)
      call ga_pack(g_src2, g_sink2,g_sbit2,1,itoff,icount)
C     ***  src2: 1 0 0 0 2 0 0 0 3 0 0 0 4 0 0 0 5 0 0 0 6 .....
C     *** sbit2: 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 .....
C     *** sink2: 1 2 3 4 5 ... icount ... 0 0 0 0 0 0 0 0 0 ..... 0
C         [NOTE: icount should equal nelements.]
      if(icount.ne.nelements) then
         print *,"Wrong number of elements from PACK operation: ",
     *           "count = ",icount," should be = ",nelements
      endif

      call ga_add(1,g_ielm1,-1,g_sink2,g_result1)
C     ***   ielm1: 1 2 3 4 5 ... nelements ... 0 0 0 0 0 0 0 0
C     ***   sink2: 1 2 3 4 5 ... nelements ... 0 0 0 0 0 0 0 0
C     *** result1: 0 0 0 0 0 0 0 ... 0  (result1=ielm1-sink2)
c***  Find the maximum value and the index of the maximum value
      call nga_select_elem(g_result1,'max',imax_ga,ilocmax_ga)
c***  Find the minimum value and the index of the minimum value
      call nga_select_elem(g_result1,'min',imin_ga,ilocmin_ga)
      if(imax_ga.eq.0.and.imin_ga.eq.0) then
         if (me.eq.0) then
           print *,"Unpack successful"
         endif
      else
         print *,"Unpack unsuccessful", imin_ga,ilocmin_ga
         print *,"Unpack unsuccessful", imax_ga,ilocmax_ga
         call ga_error("Unpack unsuccessful", -1)
      endif
      call ga_sync()
C
      goto 9999
 9999 continue
C
C     ..................................................................
C
      return
      end
      subroutine test_nga_scan()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
C
      integer g_ielm1, g_result1                           ! handles to global arrays
      integer g_ielm2, g_sbit2, g_src2, g_sink2, g_result2 ! handles to global arrays
      integer g_irowmat, g_icolmat                         ! handles to global arrays
      integer g_velm4, g_velm5                             ! handles to global arrays
C
      integer ilo, ihi, num
      integer ilo4, ihi4, num4
      integer i, j, it, count
      integer me, nproc                       ! my processor & number of procs 
      integer ndim,dims(1),chunk(1)
C
      integer nelements_max, nen_max
      parameter (nelements_max=10000, nen_max=4)
      integer ipack(nelements_max), 
     *        jpack(nelements_max), 
     *        kpack(nelements_max), 
     *        ielm1(nelements_max), 
     *        ielm2(nen_max*nelements_max),
     *        ielm3(nen_max*nelements_max),
     *        ielm4(nen_max*nelements_max),
     *        ielm5(nen_max*nelements_max)
      integer itest(nen_max*nelements_max), one
      double precision velm4(nen_max*nelements_max), 
     *                 velm5(nen_max*nelements_max)
c
      integer ilocmax, ilocmin
      integer ilocmax_ga, ilocmin_ga
      integer imax_ga, imin_ga
C
      integer nelements, nen
      integer isum, icount, itoff
      double precision result
      integer i1, i2, ichk
C
c
c***  check parallel environment
      me = ga_nodeid()
      nproc = ga_nnodes()
c
      nelements=20000
      nelements=20
      nen=4
c
c***  create a global 1-D array
      ndim=1
      dims(1)=nen*nelements
      chunk(1)=min(1,dims(1)/nproc)
      if (.not. nga_create(MT_INT, ndim, dims, 'ielm1', chunk, 
     $                     g_ielm1))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'result1', chunk, 
     $                     g_result1))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'sbit2', chunk, 
     $                     g_sbit2))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'src2', chunk, 
     $                     g_src2))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'sink2', chunk, 
     $                     g_sink2))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'result2', chunk, 
     $                     g_result2))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'irowmat', chunk, 
     $                     g_irowmat))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'icolmat', chunk, 
     $                     g_icolmat))
     $     call ga_error(' ga_create failed ',0)
C
      if (.not. nga_create(MT_DBL, ndim, dims, 'velm4', chunk, 
     $                     g_velm4))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_DBL, ndim, dims, 'velm5', chunk, 
     $                     g_velm5))
     $     call ga_error(' ga_create failed ',0)
c
c     Enumerate the sbit arrary to get a sequential vector.
      call ga_fill(g_ielm1,0)
      call ga_patch_enum(g_ielm1,1,nelements,1,1)
c     *** ielm1: 1 2 3 4 5 ... nelements ... 0 0 0 0 0
c
c***  compute local ilo, ihi, num for each processor
      call nga_distribution(g_ielm1,me,ilo,ihi,chunk)
      do i=0,nproc-1
         if(me.eq.i) then
            write(6,'(a,i4,a,i8,a,i8)') 'me: ',me,
     +            ' lo: ', ilo,' hi: ',ihi
            call ffflush(6)
         endif
         call ga_sync()
      enddo
      num=ihi-ilo+1
      if(ihi.le.0) num=0
      if(ihi.gt.nelements) then
         ihi=nelements
         num=max(0,ihi-ilo+1)
      endif
      if(num.gt.nelements_max) then
         print *,"Too many elements "
         goto 9999
      endif
cbjp      print *, 'me=',me, num,ilo,ihi
      do i=0,nproc-1
         if(me.eq.i) then
            write(6,'(a,i4,a,i8,a,i8,a,i8)') 'me: ',me,
     +            ' num: ',num, ' lo: ', ilo,' hi: ',ihi
            call ffflush(6)
         endif
         call ga_sync()
      enddo
c
c***  scatter some values into the global array
      call ga_fill(g_sbit2,0)
      call ga_fill(g_src2,0)
      if(num.gt.0) call nga_get(g_ielm1,ilo,ihi,ielm1,1)
      do i=ilo,ihi
         ipack(i-ilo+1)=1+nen*(ielm1(i-ilo+1)-1)
         jpack(i-ilo+1)=1
         kpack(i-ilo+1)=ielm1(i-ilo+1)
      enddo
      if(num.gt.0) call NGA_scatter(g_sbit2,jpack,ipack,num)
      if(num.gt.0) call NGA_scatter(g_src2 ,kpack,ipack,num)
c     *** sbit2: 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 .....
c     ***  src2: 1 0 0 0 2 0 0 0 3 0 0 0 4 0 0 0 5 0 0 0 6 .....
c
c     ..................................................................
c     TEST THE "PACK" FUNCTION.
c
      itoff=nen*nelements
      call ga_fill(g_sink2,0)
      call ga_pack(g_src2, g_sink2,g_sbit2,1,itoff,icount)
c     *** sbit2: 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 .....
c     ***  src2: 1 0 0 0 2 0 0 0 3 0 0 0 4 0 0 0 5 0 0 0 6 .....
c     *** sink2: 1 2 3 4 5 ... icount ... 0 0 0 0 0 0 0 0 0 ..... 0
c         [NOTE: icount should equal nelements.]
      if(icount.ne.nelements) then
         print *,"Wrong number of elements from PACK operation: ",
     *           "count = ",icount," should be = ",nelements
      endif
      call ga_add(1,g_ielm1,-1,g_sink2,g_result1)
c     ***   ielm1: 1 2 3 4 5 ... nelements ... 0 0 0 0 0 0 0 0
c     ***   sink2: 1 2 3 4 5 ... nelements ... 0 0 0 0 0 0 0 0
c     *** result1: 0 0 0 0 0 0 0 ... 0  (result1=ielm1-sink2)
c***  Find the maximum value and the index of the maximum value
      call nga_select_elem(g_result1,'max',imax_ga,ilocmax_ga)
c***  Find the minimum value and the index of the minimum value
      call nga_select_elem(g_result1,'min',imin_ga,ilocmin_ga)
      if(imax_ga.eq.0.and.imin_ga.eq.0) then
         if (me.eq.0) then
           print *,"Pack successful"
         endif
      else
         call ga_error("Pack unsuccessful",-1)
      endif
      call ga_sync()
c
c     ..................................................................
c     TEST THE "UNPACK" FUNCTION.
c
      itoff=nen*nelements
      call ga_fill(g_src2,0)
      call ga_unpack(g_sink2,g_src2,g_sbit2,1,itoff,count)
*      call ga_unpack(g_src2,g_sink2,g_sbit2,1,nelements,count)
C     *** sink2: 1 2 3 4 5 ... icount ... 0 0 0 0 0 0 0 0 0 ..... 0
C     *** sbit2: 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 .....
C     ***  src2: 1 0 0 0 2 0 0 0 3 0 0 0 4 0 0 0 5 0 0 0 6 .....
      call ga_fill(g_sink2,0)
      call ga_pack(g_src2, g_sink2,g_sbit2,1,itoff,icount)
*     call ga_pack(g_sink2,g_src2,g_sbit2,1,itoff,icount)
C     *** sbit2: 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 .....
C     ***  src2: 1 0 0 0 2 0 0 0 3 0 0 0 4 0 0 0 5 0 0 0 6 .....
C     *** sink2: 1 2 3 4 5 ... icount ... 0 0 0 0 0 0 0 0 0 ..... 0
C         [NOTE: icount should equal nelements.]
      if(icount.ne.nelements) then
         print *,"Wrong number of elements from PACK operation: ",
     *           "count = ",icount," should be = ",nelements
      endif
      call ga_add(1,g_ielm1,-1,g_sink2,g_result1)
C     ***   ielm1: 1 2 3 4 5 ... nelements ... 0 0 0 0 0 0 0 0
C     ***   sink2: 1 2 3 4 5 ... nelements ... 0 0 0 0 0 0 0 0
C     *** result1: 0 0 0 0 0 0 0 ... 0  (result1=ielm1-sink2)
c***  Find the maximum value and the index of the maximum value
      call nga_select_elem(g_result1,'max',imax_ga,ilocmax_ga)
c***  Find the minimum value and the index of the minimum value
      call nga_select_elem(g_result1,'min',imin_ga,ilocmin_ga)
      if(imax_ga.eq.0.and.imin_ga.eq.0) then
         if (me.eq.0) then
           print *,"Unpack successful"
         endif
      else
         call ga_error("Unpack unsuccessful",-1)
      endif
      call ga_sync()
*     call ga_print(g_src2)
c
c     ..................................................................
c     TEST THE SCAN COPY FUNCTION.
c
c***  compute local ilo, ihi, num for each processor
      call nga_distribution(g_irowmat,me,ilo4,ihi4,chunk)
      num4=ihi4-ilo4+1
      if(ihi.le.0) num4=0
      if(num4.gt.nen*nelements_max) then
         print *,"Too many elements "
         goto 9999
      endif
cbjp      print *, 'me=',me, num4,ilo4,ihi4
      do i=0,nproc-1
         if(me.eq.i) then
            write(6,'(a,i4,a,i8,a,i8,a,i8)') 'me: ',me,
     +            ' num4: ',num4, ' lo4: ', ilo4,' h4i: ',ihi4
            call ffflush(6)
         endif
         call ga_sync()
      enddo
      call ga_fill(g_irowmat,-1)
      call GA_scan_copy(g_src2, g_irowmat, g_sbit2,1,itoff)
c     ***   sbit2: 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 .....
c     ***    src2: 1 0 0 0 2 0 0 0 3 0 0 0 4 0 0 0 5 0 0 0 6 .....
c     *** irowmat: 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 .....
c
      call nga_get(g_irowmat,ilo4,ihi4,ielm4,1)
      do i=ilo4,ihi4
         velm4(i-ilo4+1)=1.0d+00
      enddo
      call ga_fill(g_velm4,0.0d+00)
cbjp      print *, num4, (ielm4(i),i=1,min(4,num4))
      do j=0,nproc-1
         if(me.eq.j) then
            write(6,'(a,i4,a,i8,a,i8,a,i8)') 'me: ',me, ' num4: ',num4,
     +            ' ielm4: ', (ielm4(i),i=1,min(4,num4))
            call ffflush(6)
         endif
         call ga_sync()
      enddo
      call nga_scatter_acc(g_velm4,velm4,ielm4,num4,1.0d+00)
C     *** velm4: 4 4 4 4 4 4 4 4 4 4 ... nelements
      result=ga_ddot(g_velm4,g_velm4)
      if(result.eq.nen*nen*nelements) then
         if (me.eq.0) then
           print *,"Scan_copy successful."
         endif
      else
         call ga_error("Scan_copy failed.",-1)
      endif
      call ga_sync()
c
c     ..................................................................
c     TEST THE SCAN ADD FUNCTION.
c
      call ga_fill(g_src2,1)
      call ga_fill(g_icolmat,0)
*     call ga_print(g_icolmat)
*     call ga_print(g_src2)
*     call ga_print(g_sbit2)
c      if (me.eq.0) then
c        one = 1
c        call nga_get(g_src2,one,dims,itest,one)
c        write(6,'(25i4)') (itest(i),i=1,dims(1))
c        write(6,*)
c        call nga_get(g_sbit2,one,dims,itest,one)
c        write(6,'(25i4)') (itest(i),i=1,dims(1))
c      endif
      call ga_sync
      call GA_scan_add(g_src2, g_icolmat,g_sbit2,1,itoff,0)
c      if (me.eq.0) then
c        one = 1
c        write(6,*)
c        call nga_get(g_icolmat,one,dims,itest,one)
c        write(6,'(25i4)') (itest(i),i=1,dims(1))
c      endif
c     ***   sbit2: 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 .....
c     ***    src2: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
c     *** icolmat: 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 .....
      call ga_fill(g_velm4,1.0d+00)
      call nga_get(g_icolmat,ilo4,ihi4,ielm4,1)
      do i=ilo4,ihi4
         velm4(i-ilo4+1)=ielm4(i-ilo4+1)
         velm5(i-ilo4+1)=1.0d+00
      enddo
      call nga_put(g_velm4,ilo4,ihi4,velm4,1)
      call nga_put(g_velm5,ilo4,ihi4,velm5,1)
      result=ga_ddot(g_velm4,g_velm5)
c     *** velm4: 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 .....
c     *** velm5: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 .....
c     *** result = sum(velm4 * velm5) = (1+2+3+4+...+nen)*nelements
      isum=0
      do i=1,nen
         isum=isum+i
      enddo
      if(result.eq.isum*nelements) then
         if (me.eq.0) then
           print *,"Scan_add successful."
         endif
      else
         call ga_error("Scan_add failed.",-1)
      endif
      call ga_sync()
      call ga_sync
      call GA_scan_add(g_src2, g_icolmat,g_sbit2,1,itoff,1)
c     ***   sbit2: 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 .....
c     ***    src2: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
c     *** icolmat: 0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3 0 .....
c      if (me.eq.0) then
c        one = 1
c        write(6,*)
c        call nga_get(g_icolmat,one,dims,itest,one)
c        write(6,'(25i4)') (itest(i),i=1,dims(1))
c      endif
      call nga_get(g_src2, 1, itoff, ielm2, one)
      call nga_get(g_sbit2, 1, itoff, ielm3, one)
      call nga_get(g_icolmat, 1, itoff, ielm4, one)
      call nga_get(g_icolmat, 1, itoff, ielm5, one)
      ichk = 1
      do i = 1, itoff
        if (ielm3(i).eq.1.or.i.eq.1) then
          ielm4(i) = 0
        else
          ielm4(i) = ielm4(i-1) + ielm2(i-1)
        endif
        if (ielm4(i).ne.ielm5(i)) then
          ichk = 0
        endif
      end do
cc      ichk=ga_idot(g_src2,g_icolmat)
      isum = 0
      do i=1,nen
         isum=isum+i-1
      enddo
cc      if (ichk.eq.nelements*isum) then
      if (ichk.eq.1) then
        if (me.eq.0) then
          print *,"Exclusive scan_add successful."
        endif
      else
       call ga_error("Exclusive scan_add failed.",-1)
      endif
c
      goto 9999
 9999 continue
c
c     ..................................................................
c
      return
      end
