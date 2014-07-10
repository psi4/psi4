#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c vector boxes lack arithmetic precision 
#ifdef CRAY_YMP
# define THRESH 1d-10
#elif defined(FUJITSU)
# define THRESH 1d-12
# define THRESHF 1e-5
#else
# define THRESH 1d-13
# define THRESHF 1e-5
#endif
#define MISMATCH(x,y) abs(x-y)/max(1d0,abs(x)).gt.THRESH
#define MISMATCHF(x,y) abs(x-y)/max(1.0,abs(x)).gt.THRESHF

#define USE_SCALAPACK 0
#define SMALL_TEST 1

      program main
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
      integer heap, stack, fudge, ma_heap, me
      integer nmax, DIM
#if SMALL_TEST
      parameter (nmax = 8, DIM = 2)
#else
      parameter (nmax = 2375, DIM = 2)
#endif
      integer ndim, nprocs, type, length
      integer g_a, g_b, g_c, g_d, g_e, g_f, g_h, g_i, g_j, inode
      MA_ACCESS_INDEX_TYPE index
      integer i, j, k, nb, dims(7)
      integer lo(7), hi(7), tlo(7), thi(7), t2lo(7), t2hi(7)
      integer block_list(10000), block_map(10000), nblock
      integer chunk(7), ld(7), block(7), proc_grid(7)
      integer a(nmax, nmax),b(nmax,nmax),e(nmax,nmax),f(nmax,nmax)
      integer skip(7), i_one, ialpha, ibeta, check_int
      double precision c(nmax,nmax),d(nmax,nmax), one, ddot
      double precision alpha, beta
      logical status
      parameter (heap=400*400*4, fudge=100, stack=400*400*4)
c
c***  Intitialize a message passing library
c
#include "mp3.fh"
c
c***  Initialize GA
c
      call ga_initialize()
      me = ga_nodeid()
c
      ma_heap = heap + fudge 
c
      if(me.eq.0)then
        write(6,100)
        call ffflush(6)
  100   format(' GA initialized')
      endif
c
c***  Initialize the MA package
c     MA must be initialized before any global array is allocated
c
      status = ma_init(MT_DBL, stack, ma_heap)
      if (.not. status) call ga_error('ma_init failed',-1) 
c
      if(me.eq.0)then
        write(6,101) ga_nnodes()
        call ffflush(6)
  101   format(' Using ',i3,' process(es)')
      endif
c
c  Create GA
c
      ndim = DIM
      do i = 1, ndim
        dims(i) = nmax
        ld(i) = nmax
        lo(i) = 1
        hi(i) = nmax
        chunk(i) = -1
#if SMALL_TEST
        block(i) = 2 
#else
        block(i) = 87
#endif
      end do
      proc_grid(1) = 2
      proc_grid(2) = ga_nnodes()/2
c
      g_a = ga_create_handle() 
      call ga_set_data(g_a, ndim, dims, MT_F_INT)
      call ga_set_chunk(g_a, chunk)
#if USE_SCALAPACK
      call ga_set_block_cyclic_proc_grid(g_a, block, proc_grid)
#else
      call ga_set_block_cyclic(g_a, block)
#endif
      status = ga_allocate(g_a)
      call ga_zero(g_a)
c
      g_b = ga_create_handle() 
      call ga_set_data(g_b, ndim, dims, MT_F_DBL)
      call ga_set_chunk(g_b, chunk)
#if USE_SCALAPACK
      call ga_set_block_cyclic_proc_grid(g_b, block, proc_grid)
#else
      call ga_set_block_cyclic(g_b, block)
#endif
      status = ga_allocate(g_b)
      call ga_zero(g_b)
c
      g_c = ga_create_handle() 
      call ga_set_data(g_c, ndim, dims, MT_F_INT)
      call ga_set_chunk(g_c, chunk)
      status = ga_allocate(g_c)
      call ga_zero(g_c)
c
      g_d = ga_create_handle() 
      call ga_set_data(g_d, ndim, dims, MT_F_INT)
      call ga_set_chunk(g_d, chunk)
#if USE_SCALAPACK
      call ga_set_block_cyclic_proc_grid(g_d, block, proc_grid)
#else
      call ga_set_block_cyclic(g_d, block)
#endif
      status = ga_allocate(g_d)
      call ga_zero(g_d)
c
      g_e = ga_create_handle() 
      call ga_set_data(g_e, ndim, dims, MT_F_DBL)
      call ga_set_chunk(g_e, chunk)
#if USE_SCALAPACK
      call ga_set_block_cyclic_proc_grid(g_e, block, proc_grid)
#else
      call ga_set_block_cyclic(g_e, block)
#endif
      status = ga_allocate(g_e)
      call ga_zero(g_e)
c
      g_f = ga_create_handle() 
      call ga_set_data(g_f, ndim, dims, MT_F_DBL)
      call ga_set_chunk(g_f, chunk)
#if USE_SCALAPACK
      call ga_set_block_cyclic_proc_grid(g_f, block, proc_grid)
#else
      call ga_set_block_cyclic(g_f, block)
#endif
      status = ga_allocate(g_f)
      call ga_zero(g_f)
c
      g_h = ga_create_handle() 
      call ga_set_data(g_h, ndim, dims, MT_F_DBL)
      call ga_set_chunk(g_h, chunk)
      status = ga_allocate(g_h)
      call ga_zero(g_h)
c
      g_i = ga_create_handle() 
      call ga_set_data(g_i, ndim, dims, MT_F_DBL)
      call ga_set_chunk(g_i, chunk)
      status = ga_allocate(g_i)
      call ga_zero(g_i)
c
      g_j = ga_create_handle() 
      call ga_set_data(g_j, ndim, dims, MT_F_DBL)
      call ga_set_chunk(g_j, chunk)
      status = ga_allocate(g_j)
      call ga_zero(g_j)
c
      if(me.eq.0)then
        write(6,102)
        call ffflush(6)
  102   format(' Completed allocation of GAs')
      endif
c
c  Initialize local arrays
c
#if SMALL_TEST
      tlo(1) = 3
      thi(1) = 4
      tlo(2) = 4
      thi(2) = 7
      t2lo(1) = 4
      t2hi(1) = 5
      t2lo(2) = 5
      t2hi(2) = 8
#else
      tlo(1) = 3
      thi(1) = 2277
      tlo(2) = 4
      thi(2) = 1832
      t2lo(1) = 4
      t2hi(1) = 2278
      t2lo(2) = 5
      t2hi(2) = 1833
#endif
c
      one = 1.0d00
      i_one = 1
      k = 1
      do i = 1, nmax
        do j = 1, nmax
          a(i,j) = k
          b(i,j) = 0
          c(i,j) = dble(k)
          d(i,j) = 0.0d00
          e(i,j) = 0
          f(i,j) = 0
          k = k + 1
        end do
      end do
c
c  Test PUT and GET
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Put and GA_Get...'
        write(6,*)
      endif
c
c  Copy data from local array to global array
c
      if (me.eq.0) then
        call nga_put(g_a,lo,hi,a,ld)
      endif
c
c  Copy data from global array back to local array
c
      call ga_sync
      call nga_get(g_a,lo,hi,e,ld)
      call ga_sync
c
      do j = 1, nmax
        do i = 1, nmax
          if (a(i,j).ne.e(i,j)) then
            write(6,103) me,a(i,j),e(i,j)
            call ga_error('GA PUT and GET failed',me)
          endif
          e(i,j) = 0
        end do
      end do
  103 format('proc: ',i4,' a(i,j): ',i8,' e(i,j): ',i8)
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'GA_Put and GA_Get are OK'
        write(6,*)
      endif
c
c  Test ga_add
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Add...'
        write(6,*)
      endif
      if (me.eq.0) then
        call nga_put(g_c,lo,hi,a,ld)
      endif
      ialpha = 1
      ibeta = -1
      call ga_add(ialpha, g_a, ibeta, g_c, g_c)
      if (ga_idot(g_c,g_c).eq.0) then
        if(me.eq.0)then
          write(6,*)
          write(6,*) 'GA_Add operation for regular and block-cyclic',
     +               ' arrays is OK'
          write(6,*)
          call ffflush(6)
        endif
      else
        if(me.eq.0)then
          write(6,*)
          write(6,*) 'GA_Add operation for regular and block-cyclic',
     +               ' arrays is not OK'
          write(6,*)
          call ga_error('exiting', 1)
          call ffflush(6)
        endif
      endif
c
c   Check add between two block-cyclic arrays
c
      if (me.eq.0) then
        call nga_put(g_d,lo,hi,a,ld)
      endif
      call ga_add(ialpha, g_a, ialpha, g_d, g_d)
c
c   Copy data to local buffer and then back to regular array
c
      call nga_get(g_d,lo,hi,b,ld)
      if (me.eq.0) then
        call nga_get(g_c,lo,hi,b,ld)
      endif
      call ga_sync
      if (ga_idot(g_c,g_c).eq.0) then
        if(me.eq.0)then
          write(6,*)
          write(6,*) 'GA_Add operation for two block-cyclic',
     +               ' arrays is OK'
          write(6,*)
          call ffflush(6)
        endif
      else
        if(me.eq.0)then
          write(6,*)
          write(6,*) 'GA_Add operation for two block-cyclic arrays',
     +               ' is not OK'
          write(6,*)
          call ga_error('exiting', 1)
          call ffflush(6)
        endif
      endif
c
c  Check to find which blocks correspond to local patch
c
#if SMALL_TEST
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Check nga_locate_region and nga_locate_num_blocks'
        write(6,*) 'functions'
        write(6,*)
        write(6,*) 'Printing original matrix'
        do i = 1, min(nmax,8)
          write(6,200) (a(i,j),j=1,min(nmax,8))
        end do
      endif
      status=nga_locate_region(g_a,tlo,thi,block_map,block_list,nblock)
      nb = nga_locate_num_blocks(g_a,tlo,thi)
      if (me.eq.0) then
        write(6,*)
        write(6,111) nb
        do i = 1, nblock
          write(6,110) i,block_list(i),block_map(4*(i-1)+1),
     +      block_map(4*(i-1)+3),block_map(4*(i-1)+2),
     +      block_map(4*(i-1)+4)
        end do
      endif
  110 format(i3,' block(',i3,') tlo(1): ',i3,' thi(1): ',i3,
     +                        ' tlo(2): ',i3,' thi(2): ',i3)
  111 format(' Number of blocks: ',i3)
      if (me.eq.0) then
        call nga_access_block_segment(g_a,me,index,length)
        call print_block(int_mb(index),length)
      endif
#endif
c
c  Check onesided accumulate
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Acc...'
        write(6,*)
      endif
      if (me.eq.0) then
        call nga_acc(g_b,lo,hi,c,ld,one)
      endif
      call ga_sync
      if (me.eq.0) then
        call nga_acc(g_b,lo,hi,c,ld,one)
      endif
      call ga_sync
      if (me.eq.0) then
        call nga_put(g_h,lo,hi,c,ld)
      endif
      alpha = 1.0
      beta = -2.0
      call ga_add(alpha,g_b,beta,g_h,g_h)
      if (ga_ddot(g_h,g_h).eq.0.0d00) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Acc is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Acc is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif

c
c  Check copy
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Copy...'
        write(6,*)
      endif
      call ga_copy(g_a, g_c)
      call ga_copy(g_c, g_d)
      call ga_add(ialpha,g_a,ibeta,g_d,g_c)
      if (ga_idot(g_c,g_c).eq.0) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Copy is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Copy is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
c
c  Check zero
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Zero...'
        write(6,*)
        call nga_get(g_a,lo,hi,a,ld)
      endif
      call ga_sync
      call ga_zero(g_a)
      call ga_copy(g_a,g_c)
      if (ga_idot(g_c,g_c).eq.0) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Zero is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Zero is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
c
c  Check scale
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Scale...'
        write(6,*)
      endif
      if (me.eq.0) then
        call nga_put(g_b,lo,hi,c,ld)
        call nga_put(g_h,lo,hi,c,ld)
      endif
      call ga_scale(g_b,2.0d00)
      call ga_scale(g_h,2.0d00)
      alpha = 1.0d00
      beta = -1.0d00
      call ga_add(alpha,g_b,beta,g_h,g_h)
      if (ga_ddot(g_h,g_h).eq.0.0d00) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Scale is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Scale is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
c
c  Check fill
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Fill...'
        write(6,*)
      endif
      call ga_fill(g_a,1)
      call ga_fill(g_c,1)
      call ga_add(ialpha,g_a,ibeta,g_c,g_c)
      if (ga_idot(g_c,g_c).eq.0) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Fill is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Fill is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
c
c  Check zero_patch
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Zero_patch...'
        write(6,*)
      endif
      call ga_copy(g_a,g_c)
      call nga_zero_patch(g_a,tlo,thi)
      call nga_zero_patch(g_c,tlo,thi)
      call ga_add(ialpha,g_a,ibeta,g_c,g_c)
      if (ga_idot(g_c,g_c).eq.0) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Zero_patch is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Zero_patch is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
c
c  Check fill_patch
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Fill_patch...'
        write(6,*)
      endif
      call ga_copy(g_a,g_c)
      call nga_fill_patch(g_a,tlo,thi,2)
      call nga_fill_patch(g_c,tlo,thi,2)
      call ga_add(ialpha,g_a,ibeta,g_c,g_c)
      if (ga_idot(g_c,g_c).eq.0) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Fill_patch is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Fill_patch is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
c
c  Check scale_patch
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Scale_patch...'
        write(6,*)
      endif
      call ga_copy(g_a,g_c)
      call nga_scale_patch(g_a,tlo,thi,2)
      call nga_scale_patch(g_c,tlo,thi,2)
      call ga_add(ialpha,g_a,ibeta,g_c,g_c)
      if (ga_idot(g_c,g_c).eq.0) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Scale_patch is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Scale_patch is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
c
c  Check copy_patch
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Copy_patch...'
        write(6,*)
      endif
      call ga_fill(g_a,1)
      call ga_fill(g_c,1)
      call nga_fill_patch(g_a,tlo,thi,2)
      call nga_copy_patch('n',g_a,tlo,thi,g_c,t2lo,t2hi)
      call ga_fill(g_a,1)
      call nga_fill_patch(g_a,t2lo,t2hi,2)
      call ga_add(ialpha,g_a,ibeta,g_c,g_c)
      if (ga_idot(g_c,g_c).eq.0) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Copy_patch from block-cyclic to regular is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Copy_patch from block-cyclic to regular is',
     +               ' not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
      call ga_fill(g_a,1)
      call ga_fill(g_c,1)
      call nga_fill_patch(g_c,tlo,thi,2)
      call nga_copy_patch('n',g_c,tlo,thi,g_a,t2lo,t2hi)
      call ga_fill(g_c,1)
      call nga_fill_patch(g_c,t2lo,t2hi,2)
      call ga_add(ialpha,g_a,ibeta,g_c,g_c)
      if (ga_idot(g_c,g_c).eq.0) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Copy_patch from regular to block-cyclic is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Copy_patch from regular to block-cyclic is',
     +               ' not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
c
c   test ga_add_patch
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Add_patch...'
        write(6,*)
      endif
      call ga_zero(g_b)
      call ga_zero(g_e)
      call ga_zero(g_f)
      call ga_zero(g_h)
      call ga_zero(g_i)
      call ga_zero(g_j)
      if (me.eq.0) then
        call nga_put(g_e,t2lo,t2hi,c,ld)
        call nga_put(g_f,t2lo,t2hi,c,ld)
        call nga_put(g_i,t2lo,t2hi,c,ld)
        call nga_put(g_j,t2lo,t2hi,c,ld)
      endif
      call ga_sync
      call nga_add_patch(one,g_e,t2lo,t2hi,one,g_f,tlo,thi,g_b,tlo,thi)
      call nga_add_patch(one,g_i,t2lo,t2hi,one,g_j,tlo,thi,g_h,tlo,thi)
      call ga_add(alpha,g_b,beta,g_h,g_i)
      if (ga_ddot(g_i,g_i).eq.0.0d00) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Add_patch is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Add_patch is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
c
c   test ga_ddot
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Ddot...'
        write(6,*)
      endif
      call ga_fill(g_b, 2.0d00)
      call ga_fill(g_b, 2.0d00)
      ddot = ga_ddot(g_b, g_b)
      if (me.eq.0) then
        write(6,112) ddot,dble(nmax*nmax*4)
      endif
  112 format(' Value of DDOT: ',f12.2,' Expected value: ',f12.2)
      if (me.eq.0) write(6,*)
c
c   test ga_ddot_patch
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Ddot_patch...'
        write(6,*)
      endif
      ddot = nga_ddot_patch(g_b, 'n', tlo, thi, g_b, 'n', tlo, thi)
      if (me.eq.0) then
        write(6,113) ddot,dble((thi(1)-tlo(1)+1)*(thi(2)-tlo(2)+1)*4)
      endif
  113 format(' Value of DDOT_PATCH: ',f12.2,' Expected value: ',f12.2)
c
c   test ga_abs_value_patch
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Abs_patch...'
        write(6,*)
      endif
      call ga_fill(g_a,-1)
      call ga_abs_value_patch(g_a, tlo, thi)
      call ga_fill(g_c,-1)
      call ga_abs_value_patch(g_c, tlo, thi)
      call ga_add(ialpha,g_a,ibeta,g_c,g_c)
      if (ga_idot(g_c,g_c).eq.0) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Abs_value_patch is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Abs_value_patch is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
c
c   test ga_elem_multiply
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Elem_multiply...'
        write(6,*)
      endif
      call ga_fill(g_b, 2.0d00)
      call ga_fill(g_e, 3.0d00)
      call ga_zero(g_f)
      call ga_elem_multiply(g_b, g_e, g_f)
      call ga_fill(g_h, 2.0d00)
      call ga_fill(g_i, 3.0d00)
      call ga_zero(g_j)
      call ga_elem_multiply(g_h, g_i, g_j)
      call ga_add(alpha,g_f,beta,g_j,g_h)
      if (ga_ddot(g_h,g_h).eq.0.0d00) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Elem_multiply is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Elem_multiply is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
c
c   test ga_elem_divide_patch
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Elem_divide_patch...'
        write(6,*)
      endif
      call ga_fill(g_b, 3.0d00)
      call ga_elem_divide_patch(g_f,t2lo,t2hi,g_e,tlo,thi,g_b,tlo,thi)
      call ga_fill(g_h, 3.0d00)
      call ga_elem_divide_patch(g_j,t2lo,t2hi,g_i,tlo,thi,g_h,tlo,thi)
      call ga_add(alpha,g_b,beta,g_h,g_h)
      if (ga_ddot(g_h,g_h).eq.0.0d00) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Elem_divide_patch is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Elem_divide_patch is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
c
c    test strided get
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Strided_get...'
        write(6,*)
      endif
      skip(1) = 2
      skip(2) = 2
      call ga_zero(g_a)
      call ga_zero(g_c)
      call nga_get(g_a,lo,hi,a,ld)
      call nga_get(g_c,lo,hi,b,ld)
      call ga_fill(g_a, 1)
      call ga_fill(g_c, 1)
      if (me.eq.0) then
        call nga_strided_get(g_a,lo,hi,skip,a,ld)
        call nga_strided_get(g_c,lo,hi,skip,b,ld)
      endif
      call ga_sync
      nb = 0
      if (me.eq.0) then
        do i = 1, nmax
          do j = 1, nmax
            if (a(i,j).ne.b(i,j)) nb = nb + 1
          end do
        end do
        if (nb.eq.0) then
          if (me.eq.0) then
            write(6,*)
            write(6,*) 'GA_Strided_get is OK'
            write(6,*)
          endif
        else
          if (me.eq.0) then
            write(6,*)
            write(6,*) 'GA_Strided_get is not OK'
            write(6,*)
            call ga_error('exiting', 1)
          endif
        endif
      endif
c
c    test strided put
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Strided_put...'
        write(6,*)
      endif
      call ga_fill(g_a, 1)
      if (me.eq.0) call nga_get(g_a,lo,hi,a,ld)
      call ga_zero(g_a)
      call ga_fill(g_c, 1)
      if (me.eq.0) call nga_get(g_c,lo,hi,b,ld)
      call ga_zero(g_c)
      if (me.eq.0) then
        call nga_strided_put(g_a,lo,hi,skip,a,ld)
        call nga_strided_put(g_c,lo,hi,skip,b,ld)
      endif
      call ga_sync
      call ga_add(ialpha,g_a,ibeta,g_c,g_c)
      if (ga_idot(g_c,g_c).eq.0) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Strided_put is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Strided_put is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
c
c    test strided accumulate
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Strided_acc...'
        write(6,*)
      endif
      call ga_fill(g_b, 1.0d00)
      call nga_get(g_b,lo,hi,c,ld)
      call ga_fill(g_h, 1.0d00)
      call nga_get(g_h,lo,hi,d,ld)
      call ga_sync
      if (me.eq.0) then
        call nga_strided_acc(g_b,lo,hi,skip,c,ld,one)
        call nga_strided_acc(g_h,lo,hi,skip,d,ld,one)
      endif
      call ga_sync
      call ga_add(alpha,g_b,beta,g_h,g_h)
      if (ga_ddot(g_h,g_h).eq.0.0d00) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Strided_acc is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Strided_acc is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
c
c    test transpose
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_transpose...'
        write(6,*)
      endif
      k = 0
      do i = 1, nmax
        do j = 1, nmax
          k = k+1
          c(i,j) = dble(k)
        end do
      end do
      if (me.eq.0) then
        call nga_put(g_b,lo,hi,c,ld)
        call nga_put(g_h,lo,hi,c,ld)
      endif
      call ga_sync
      call ga_transpose(g_b,g_i)
      call ga_transpose(g_h,g_j)
      call ga_add(alpha,g_i,beta,g_j,g_h)
      if (ga_ddot(g_h,g_h).eq.0.0d00) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Transpose is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Transpose is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
c
c    test symmetrize
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_symmetrize...'
        write(6,*)
      endif
      if (me.eq.0) then
        call nga_put(g_b,lo,hi,c,ld)
        call nga_put(g_h,lo,hi,c,ld)
      endif
      call ga_sync
      call ga_symmetrize(g_b)
      call ga_symmetrize(g_h)
      call ga_add(alpha,g_b,beta,g_h,g_h)
      if (ga_ddot(g_h,g_h).eq.0.0d00) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Symmetrize is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Symmetrize is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
c
c    test periodic get
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Periodic_get...'
        write(6,*)
      endif
      tlo(1) = lo(1) + 2
      thi(1) = hi(1) + 2
      tlo(2) = lo(2) + 2
      thi(2) = hi(2) + 2
      k = 0
      do i = 1, nmax
        do j = 1, nmax
          k = k + 1
          a(i,j) = k
        end do
      end do
      if (me.eq.0) then
        call nga_put(g_a,lo,hi,a,ld)
        call nga_put(g_c,lo,hi,a,ld)
      endif
      call ga_sync
      if (me.eq.0) then
        call nga_periodic_get(g_a,tlo,thi,a,ld)
        call nga_periodic_get(g_c,tlo,thi,b,ld)
      endif
      call ga_sync
      nb = 0
      if (me.eq.0) then
        do i = 1, nmax
          do j = 1, nmax
            if (a(i,j).ne.b(i,j)) nb = nb + 1
          end do
        end do
        if (nb.eq.0) then
          if (me.eq.0) then
            write(6,*)
            write(6,*) 'GA_Periodic_get is OK'
            write(6,*)
          endif
        else
          if (me.eq.0) then
            write(6,*)
            write(6,*) 'GA_Periodic_get is not OK'
            write(6,*)
            call ga_error('exiting', 1)
          endif
        endif
      endif
c
c    test periodic put
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Periodic_put...'
        write(6,*)
      endif
      k = 0
      do i = 1, nmax
        do j = 1, nmax
          k = k + 1
          a(i,j) = k
        end do
      end do
      if (me.eq.0) then
        call nga_periodic_put(g_a,tlo,thi,a,ld)
        call nga_periodic_put(g_c,tlo,thi,a,ld)
      endif
      call ga_sync
      call ga_add(ialpha,g_a,ibeta,g_c,g_c)
      if (ga_idot(g_c,g_c).eq.0) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Periodic_put is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Periodic_put is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
c
c    test periodic accumulate
c
      if (me.eq.0) then
        write(6,*)
        write(6,*) 'Testing GA_Periodic_acc...'
        write(6,*)
      endif
      call ga_fill(g_b, 1.0d00)
      call ga_fill(g_h, 1.0d00)
      call ga_sync
      if (me.eq.0) then
        call nga_periodic_acc(g_b,tlo,thi,c,ld,one)
        call nga_periodic_acc(g_h,tlo,thi,c,ld,one)
      endif
      call ga_sync
      call ga_add(alpha,g_b,beta,g_h,g_h)
      if (ga_ddot(g_h,g_h).eq.0.0d00) then
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Periodic_acc is OK'
          write(6,*)
        endif
      else
        if (me.eq.0) then
          write(6,*)
          write(6,*) 'GA_Periodic_acc is not OK'
          write(6,*)
          call ga_error('exiting', 1)
        endif
      endif
      
  200 format(8i8)
c 300 format(8f8.1)
      if (me.eq.0) then
        write(6,302)
      endif
  302 format(' Finished tests: success')
      status = ga_destroy(g_a)
      status = ga_destroy(g_b)
      status = ga_destroy(g_c)
      status = ga_destroy(g_d)
      status = ga_destroy(g_e)
      status = ga_destroy(g_f)
      call ga_terminate()
c
c***  Tidy up after message-passing library
c
      call MP_FINALIZE()
c
      stop
      end
c
      subroutine print_block(a,length)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer length, i, j
      integer a(length)
      write(6,*)
      write(6,200)
      do i=1, 4
        write(6,100) (a(4*(j-1)+i), j=1,4)
  100   format(8i8)
      end do
  200 format(' Write contents of processor')
      return
      end
