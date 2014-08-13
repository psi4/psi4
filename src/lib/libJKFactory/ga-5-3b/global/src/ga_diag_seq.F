#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine gai_diag_seq(g_a, g_s, g_v, evals)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer g_a               ! Matrix to diagonalize
      integer g_s               ! Metric
      integer g_v               ! Global matrix to return evecs
      double precision evals(*) ! Local array to return evals
c
      integer n, ierr
      integer l_fv1, l_fv2
      MA_ACCESS_INDEX_TYPE k_fv1, k_fv2
      integer l_a, l_s, l_v
      MA_ACCESS_INDEX_TYPE k_a, k_s, k_v
      integer dim1, dim2, type, me
      logical status
c
c
c     Solve the generalized eigen-value problem returning
c     all eigen-vectors and values in ascending order
c
c     The input matrices may be destroyed
c
      call ga_check_handle(g_a, 'ga_diag_seq a')
      call ga_check_handle(g_s, 'ga_diag_seq s')
      call ga_check_handle(g_v, 'ga_diag_seq v')
      call ga_sync()
c
c     Only process 0 does the diag
c
      call ga_inquire(g_a, type, dim1, dim2)
      if(dim1.ne.dim2)
     $  call ga_error('ga_diag_seq: nonsquare matrix ',0)
 
      n = dim1
      me = ga_nodeid()
      if (me .eq. 0) then
c
c     allocate scratch space
c     
         status=ma_push_get(MT_DBL, n*n,
     $        'diag_seq:a', l_a, k_a)
         status=status.and.ma_push_get(MT_DBL, n*n,
     $        'diag_seq:s', l_s, k_s)
         status=status.and.ma_push_get(MT_DBL, n*n,
     $        'diag_seq:v', l_v, k_v)
         status=status.and.ma_push_get(MT_DBL, n,
     $        'diag_seq:fv1', l_fv1, k_fv1)
         status=status.and.ma_push_get(MT_DBL, n,
     $        'diag_seq:fv2', l_fv2, k_fv2)
         if (.not. status)
     $        call ga_error('ga_diag_seq: insufficient memory', n)
c     
c     Fill local arrays from global arrays
c     
         call ga_get(g_a, 1, n, 1, n, dbl_mb(k_a), n)
         call ga_get(g_s, 1, n, 1, n, dbl_mb(k_s), n)
c
*         write(6,*) ' fock matrix '
*         call output(dbl_mb(k_a), 1, n, 1, n, n, n, 1)
*         write(6,*) ' overlap matrix '
*         call output(dbl_mb(k_s), 1, n, 1, n, n, n, 1)
c     
         call rsg(n, n, dbl_mb(k_a), dbl_mb(k_s), evals, 1,
     $        dbl_mb(k_v), dbl_mb(k_fv1), dbl_mb(k_fv2), ierr)
c
*         write(6,*) 'eigen vectors '
*         call output(dbl_mb(k_v), 1, n, 1, n, n, n, 1)
*         write(6,*) 'eigen values '
*         call output(evals, 1, n, 1, 1, n, 1, 1)
*         call forflush(6)
c     
c     Shove eigen-vectors back into global array
c     
         call ga_put(g_v, 1, n, 1, n, dbl_mb(k_v), n)
c     
c     Free scratch space
c     
         status = ma_pop_stack(l_fv2)
         status = ma_pop_stack(l_fv1)
         status = ma_pop_stack(l_v)
         status = ma_pop_stack(l_s)
         status = ma_pop_stack(l_a)
      endif
c     
c     Broadcast the eigenvalues to all processes
c
      call ga_brdcst(32500, evals, 
     $               ma_sizeof(MT_DBL,n,MT_BYTE), 0)
      call ga_sync()
c
      end
c
      subroutine gai_diag_std_seq(g_a,  g_v, evals)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer g_a               ! Matrix to diagonalize
      integer g_v               ! Global matrix to return evecs
      double precision evals(*) ! Local array to return evals
c
      integer n, ierr
      integer l_fv1, l_fv2
      MA_ACCESS_INDEX_TYPE k_fv1, k_fv2
      integer l_a, l_v
      MA_ACCESS_INDEX_TYPE k_a, k_v
      integer dim1, dim2, type, me
      logical status
c
c
c     Solve the standard eigen-value problem returning
c     all eigen-vectors and values in ascending order
c
c     The input matrices may be destroyed
c
      call ga_check_handle(g_a, 'ga_diag_std a')
      call ga_check_handle(g_v, 'ga_diag_std v')
      call ga_sync()
c
c     Only process 0 does the diag
c
      call ga_inquire(g_a, type, dim1, dim2)
      if(dim1.ne.dim2)
     $  call ga_error('ga_diag_std_seq: nonsquare matrix ',0)
 
      n = dim1
      me = ga_nodeid()
      if (me .eq. 0) then
c
c     allocate scratch space
c     
         status=ma_push_get(MT_DBL, n*n,
     $        'diag_std_seq:a', l_a, k_a)
         status=status.and.ma_push_get(MT_DBL, n*n,
     $        'diag_std_seq:v', l_v, k_v)
         status=status.and.ma_push_get(MT_DBL, n,
     $        'diag_std_seq:fv1', l_fv1, k_fv1)
         status=status.and.ma_push_get(MT_DBL, n,
     $        'diag_std_seq:fv2', l_fv2, k_fv2)
         if (.not. status)
     $      call ga_error('ga_diag_std_seq: insufficient memory', n)
c     
c     Fill local arrays from global arrays
c     
         call ga_get(g_a, 1, n, 1, n, dbl_mb(k_a), n)
c     
         call rs(n, n, dbl_mb(k_a),  evals, 1,
     $        dbl_mb(k_v), dbl_mb(k_fv1), dbl_mb(k_fv2), ierr)
c
c     
c     Shove eigen-vectors back into global array
c     
         call ga_put(g_v, 1, n, 1, n, dbl_mb(k_v), n)
c     
c     Free scratch space
c     
         status = ma_pop_stack(l_fv2)
         status = ma_pop_stack(l_fv1)
         status = ma_pop_stack(l_v)
         status = ma_pop_stack(l_a)
      endif
c     
c     Broadcast the eigenvalues to all processes
c
      call ga_brdcst(32500, evals, 
     $               ma_sizeof(MT_DBL,n,MT_BYTE), 0)
      call ga_sync()
c
      end
