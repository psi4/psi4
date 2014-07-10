#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      program ga_test
c $Id: pgtestmatmult.F,v 1.2.6.1 2007-04-25 21:49:59 manoj Exp $
c
c     test ga_dgemm
c     compile with make FLD_REN="optimized blas" testmatmult.x
c     Note: - change nummax for large arrays
c           - turn off "verify" for large arrays due to memory 
c             limitations, as verify=TRUE for large arrays produces 
c             segfault, dumps core,or any crap.
c     
c
      implicit none
      integer num_m,num_n,num_k,nummax,howmany,ntrans
      parameter (nummax=1024,howmany=2,ntrans=4)
      integer num1
      parameter(num1=nummax)
      integer i,ii
      real*8 h0(num1*num1)
      integer g_c,g_b,g_a
      real*8 a
      real*8 t1,mf,avg_t(ntrans),avg_mf(ntrans)
      integer itime,ntimes,me
      integer nums_m(howmany),nums_n(howmany),nums_k(howmany)
      character*1 transa(ntrans),transb(ntrans),ta,tb
      real *8 tmpa(nummax,nummax), tmpb(nummax,nummax)
      real *8 tmpc(nummax,nummax)
      logical verify,status
      integer grp_me, nproc, nprocs, proc_group(0:100), proclist(100)
      data transa/'n','t','n','t'/
      data transb/'n','n','t','t'/

      data nums_m/512,1024/
      data nums_n/512,1024/
      data nums_k/512,1024/

      

#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
#include "mp3.fh"
      if (.not.ma_init(MT_DBL,1,20000000)) then
         call ga_error('failed: ma_init(MT_DBL,1,20000000)',10)
      endif
      status = ma_set_auto_verify(.true.)
!      status = ma_set_hard_fail(.true.)
!      status = ma_set_error_print(.true.)
      call ga_initialize()
      me = ga_nodeid()
      nproc = ga_nnodes()
      if(nproc .ne. 8) then
         write(*,*) 'nproc =', nproc
         call ga_error('failed: need exactly 8 processors',10)
      endif
      call ga_sync()
      proclist(1) = 0
      proclist(2) = 1
      proclist(3) = 2
      proclist(4) = 3
      nprocs = 4
      proc_group(0) = ga_pgroup_create(proclist, nprocs)
      proclist(1) = 4
      proclist(2) = 5
      proclist(3) = 6
      proclist(4) = 7
      nprocs = 4
      proc_group(1) = ga_pgroup_create(proclist, nprocs)

      if (me.lt.4) then
         call ga_pgroup_set_default(proc_group(0))
      endif
      if (me.ge.4) then
         call ga_pgroup_set_default(proc_group(1))
      endif
      
      
c     get new proc ids for your default group 
      me = ga_nodeid()
      nproc = ga_nnodes()

      ii = 0
      do i = 1, num1*num1
         ii = ii + 1
         if(ii.gt.num1) then
            ii = 0
         endif
         h0(i) = ii
      enddo
      
c     
c     Compute times assuming 500 mflops and 5 second target time
c     
c     ntimes = max(3.0d0,5.0d0/(4.0d-9*num**3))
      ntimes = 5

      verify = .TRUE.           ! to verify ga_dgemm
c      verify = .FALSE.           ! to verify ga_dgemm

#if 1
      do ii=1,howmany
         num_m = nums_m(ii)
         num_n = nums_n(ii)
         num_k = nums_k(ii)
         a = 0.5d0/(num_m*num_n)
         if(num_m.gt.nummax .OR. num_n.gt.nummax .OR.
     &        num_k.gt.nummax) then
            call ga_error('Insufficient memory: check nummax', 0)
         endif
c     
c     
         if (.not.ga_create(MT_DBL,num_m,num_n,'g_c',
     &        0,0,g_c)) then
            call ga_error('failed: create g_c',20)
         endif
         if (.not.ga_create(MT_DBL,num_k,num_n,'g_b',
     &        0,0,g_b)) then
            call ga_error('failed: create g_b',30)
         endif
         if (.not.ga_create(MT_DBL,num_m,num_k,'g_a',
     &        0,0,g_a)) then
            call ga_error('failed: create g_a',40)
         endif
c     
c     Initialize matrices A and B
c     call load_ga_from_square(g_a,num_m,h0,num1)
c     call load_ga_from_square(g_b,num_m,h0,num1)
c     
         if(me.eq.0) then 
            call load_ga(g_a, h0, num_m, num_k)
            call load_ga(g_b, h0, num_k, num_n)
         endif
         call ga_zero(g_c)
         call ga_sync()

         if (ga_nodeid().eq.0) then
            write(*,*) ! for new line
            write(*,*) 'Matrix Multiplication C = A[', num_m, ',', 
     .           num_k, '] x B[', num_k, ',', num_n, ']'
            write(*,*) ! for new line
            call ffflush(6) 
         endif

         do i = 1, ntrans
            avg_t(i)  = 0.0d0
            avg_mf(i) = 0.0d0
         enddo

         do itime = 1, ntimes
            do i = 1, ntrans
               call ga_sync()
               ta = transa(i)
               tb = transb(i)
               t1 = util_timer()
               call ga_dgemm(ta,tb,num_m,num_n,num_k,
     &              1.0d0, g_a, g_b, 0.0d0, g_c) 
               t1 = util_timer() - t1
               if (me.eq.0) then
                  mf = 2d0*num_m*num_n*num_k/t1*1d-6/ga_nnodes()
                  avg_t(i)  = avg_t(i)+t1
                  avg_mf(i) = avg_mf(i) + mf
                  write(6,200) ' Run#', itime, t1, mf, ta, tb
                  call ffflush(6) 
                  if (verify .AND. itime.eq.1) then
                     call verify_ga_dgemm(ta, tb, num_m, num_n, num_k, 
     &                    1.0d0, g_a, g_b, 0.0d0, g_c,tmpa,tmpb,tmpc)
                  endif
               endif
            enddo
         enddo

         if (ga_nodeid().eq.0) then
            write(*,*)          ! for new line
            do i = 1, ntrans
               write(6,200) 'Average:',0,avg_t(i)/ntimes,
     .              avg_mf(i)/ntimes,transa(i),transb(i)
            enddo
            if(verify) write(*,*) 'All ga_dgemms are verified...O.K.'
            write(*,*)          ! for new line
            call ffflush(6) 
         endif

c
c         call ga_print(g_a)
c         call ga_print(g_b)
c         call ga_print(g_c)
c
         if (.not.ga_destroy(g_c)) then
            call ga_error('failed: destroy g_c',20)
         endif
         if (.not.ga_destroy(g_b)) then
            call ga_error('failed: destroy g_b',30)
         endif
         if (.not.ga_destroy(g_a)) then
            call ga_error('failed: destroy g_a',40)
         endif
      enddo
#endif     

 200  format(a15, i2, ': ', e12.4, ' seconds ',f12.1, 
     .     ' mflops/proc ', 3a2)
      
      call ga_pgroup_set_default(ga_pgroup_get_world())
      call ga_sync() 
      call ga_terminate
      call MP_FINALIZE()
c     
      end
c
c-----------------------------------------------------------------------
c     Verify for correctness. Process 0 computes BLAS dgemm 
c     locally. For larger arrays, disbale this test as memory
c     might not be sufficient
c
      subroutine verify_ga_dgemm(xt1, xt2, num_m, num_n, num_k, 
     $     alpha, g_a, g_b, beta, g_c, tmpa, tmpb, tmpc)
      implicit none
      character *1 xt1, xt2
      integer num_m, num_n, num_k, g_a, g_b, g_c
      double precision alpha, beta
      real *8 tmpa(num_m,num_k), tmpb(num_k,num_n), tmpc(num_m,num_n)

      integer i,j,type,dim1, dim2
      real *8 abs_value
#include "mafdecls.fh"
#include "testutil.fh"
#include "galinalg.fh"

      do i = 1,num_n
         do j = 1, num_m
            tmpc(j,i) = -1.0
            tmpa(j,i) = -2.0
         enddo
      enddo
      
      call ga_inquire(g_a, type, dim1, dim2)
      call ga_get(g_a, 1, dim1, 1, dim2, tmpa, dim1)
      
      call ga_inquire(g_b, type, dim1, dim2)
      call ga_get(g_b, 1, dim1, 1, dim2, tmpb, dim1)

c     compute dgemm sequentially
      call dgemm(xt1, xt2, num_m, num_n, num_k, alpha,
     &     tmpa, num_m, tmpb, num_k, beta, tmpc, num_m)
      
c     after computing c locally, verify it with the values in g_c
      call ga_inquire(g_c, type, dim1, dim2)
      call ga_get(g_c, 1, dim1, 1, dim2, tmpa, dim1)
      
      do i = 1,num_n
         do j = 1, num_m
            abs_value = abs(tmpc(j,i)-tmpa(j,i))
            if(abs_value .gt. 1.0 .OR. abs_value .lt. -1.0) then
               write(*,*) 'Values are = ', tmpc(j,i), tmpa(j,i)
               write(*,*) 'Values are = ', abs(tmpc(j,i)-tmpa(j,i)),
     .              abs_value
               call ffflush(6)
               call ga_error('verify ga_dgemm failed', 0)
            endif
         enddo
      enddo
      
      
      return 
      end
c
c-----------------------------------------------------------------------
c called by process '0' (or your master process )
c
      subroutine load_ga(handle,f, dim1,dim2)
      implicit none
      integer handle
      integer dim1,dim2,i
      real*8 f(dim1,dim2)
#include "mafdecls.fh"
      
      if(dim1.le.0)return
      if(dim2.le.0)return

      call ga_put(handle, 1, dim1, 1, dim2, f, dim1)

      return
      end
c
c-----------------------------------------------------------------------
c     must be called by all processors, if you need to fillup the 
c     entire array
c
      subroutine load_ga_from_square(handle,num,f,ndim)
      implicit none
      integer handle, memhandle
      integer num,ndim
      real*8 f(ndim,ndim)
      integer ilo, ihi, jlo, jhi, nx, ny, ibuff
      integer ga_nodeid, i1, i2, i, j, ix, jx
#include "mafdecls.fh"

      call ga_distribution(handle, ga_nodeid(), ilo, ihi, jlo, jhi)

      if(ihi.le.0)return
      if(jhi.le.0)return

c     nx = ihi - ilo + 1
c     ny = jhi - jlo + 1

      do i = ilo,ihi,ndim
         do j = jlo,jhi,ndim
            call ga_put(handle,i,min(ihi,i+ndim),j,min(jhi,j+ndim),
     &                    f,ndim)
         enddo
      enddo

      return
      end
c
c-----------------------------------------------------------------------
c     must be called by all processors, if you need to fillup the 
c     entire array
c
      subroutine load_ga_from_triangle(handle,f,ndim)
      implicit none
      integer handle, memhandle
      real*8 f(*)
      integer ndim
      integer ilo, ihi, jlo, jhi, nx, ny, ibuff
      integer ga_nodeid, i1, i2, i, j, ix, jx
#include "mafdecls.fh"

      call ga_distribution(handle, ga_nodeid(), ilo, ihi, jlo, jhi)

      if(ihi.le.0)return
      if(jhi.le.0)return

      nx = ihi - ilo + 1
      ny = jhi - jlo + 1

      if (.not.ma_alloc_get(MT_DBL,nx*ny,'flap',memhandle,ibuff)) then
         call ga_error('failed: allocate triangle',100)
      endif

      do i = 1,nx
         do j = 1,ny
            ix = i + ilo - 1
            jx = j + jlo - 1
            i1 = min(ix,jx)
            i2 = max(ix,jx)
            dbl_mb(ibuff + nx*(j-1) + (i-1) ) = f(i2*(i2-1)/2 + i1)
         enddo
      enddo

      call ga_put(handle,ilo,ihi,jlo,jhi,
     &              dbl_mb(ibuff),nx)

      if (.not.ma_free_heap(memhandle)) then
         call ga_error('failed: free triangle',100)
      endif

      return
      end
