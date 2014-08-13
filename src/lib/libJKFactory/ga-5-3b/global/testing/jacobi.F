#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c $Id: jacobi.F,v 1.9 2003-06-24 19:01:59 d3h325 Exp $
      program jacobi 
*
*     Jacobi iterative method for solving Uxx + Uyy = 1 on [0,1]x[0,1] area
*     Global Arrays are used to store and access U and exact solution vectors 
*     J.N., 06.01.94, works with GA 1.2
*
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer heap, stack
c
c***  Intitialize a message passing library
c
#include "mp3.fh"
c
      heap = 80000
      stack = heap
      call ga_initialize()                                   ! initialize GAs 
      if (.not. ma_init(MT_DBL, heap, stack))                ! initialize MA
     $    call ga_error("ma init failed",heap+stack)
#ifdef GA_TRACE
      call trace_init(10000)                     ! initialize trace
#endif
      call iterate()                                         ! do the work      
#ifdef GA_TRACE
      call trace_end(ga_nodeid())                 ! end trace
#endif
c
      if(ga_nodeid().eq.0) print *,'All tests successful '
c
      call ga_terminate()                                    ! terminate GAs 
c
      call MP_FINALIZE()
      end


      subroutine iterate()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer m                               ! grid size
      integer g_u, g_ux, g_diff               ! handles to global arrays
      parameter (m = 10)
      integer ilo, ihi, jlo, jhi, ld,  ni, nj
      integer i, k
      MA_ACCESS_INDEX_TYPE indexl, index, index_first, indexd    ! MA addressing indices
      integer h_l, h_d                        ! handles to MA objects
      integer me, nproc                       ! my processor & number of procs 
      double precision eps                    ! stopping number
      double precision S, W, C, N, E          ! coef. of G matrix and d
      logical have_data                          ! do I "own" any data ?
      double precision norm, product, uxnorm
      data eps /1d-5/
c
c***  check parallel environment
      me = ga_nodeid()
      nproc = ga_nnodes()
c
c***  create global arrays: g_u - approx. solution,  g_ux - exact solution
c                           g_diff will store g_u - g_ux 
      if (.not. ga_create(MT_DBL, m, m, 'u', 1, 1, g_u))
     $     call ga_error(' ga_create failed ',0)
      if (.not. ga_create(MT_DBL, m, m, 'ux', 1, 1, g_ux))
     $     call ga_error(' ga_create failed ',0)
      if (.not. ga_create(MT_DBL, m, m, 'diff', 1, 1, g_diff))
     $     call ga_error(' ga_create failed ',0)
c
c***  check which piece, if any, of g_ux (and g_u) I own
      call ga_distribution(g_ux, me, ilo, ihi, jlo, jhi)
      ni = ihi - ilo + 1                 ! number of 'local' rows 
      nj = jhi - jlo + 1                 ! number of 'local' columns 
      have_data = ni .gt. 0 .and. nj .gt. 0
c
      if(have_data)then
c
c***     allocate memory for D 
         if(.not. ma_push_get(MT_DBL,ni*nj,'d',h_d, indexd))
     &            call ga_error('memory allocation failed',0)
c
c***     access local data
         call ga_access(g_ux, ilo, ihi, jlo, jhi, index, ld)
c
c***     create coefficients of the jacobi iteration matrix G and exact solution
         call jacinit(m, S, W, C, N, E, DBL_MB(indexd), ni, 
     &                DBL_MB(index), ilo, ihi, jlo, jhi, ld)
c
c***     release access to the local data (the data were updated !)
         call ga_release_update(g_ux, ilo, ihi, jlo, jhi)
c
c***     allocate memory for a copy of my piece of g_u and neighb. grid points
         if(.not. ma_push_get(MT_DBL,(ni+2)*(nj+2),'loc',h_l, indexl))
     &            call ga_error('memory allocation failed',0)
c
c***     zero the allocated memory
         do i = 1, (ni+2)*(nj+2)
            DBL_MB(indexl - 1 + i) = 0d0        ! indexl points to first element
         enddo
      endif
c
c***  initial guess for u -- zero
      call ga_zero(g_u)
c
c***  the stopping test is ||u-ux||/||ux||, ||.|| - second norm
      uxnorm = sqrt(ga_ddot(g_ux, g_ux))                   ! ||g_ux||
c
c***  Now, synchronize and then iterate 
      k = 0
10    call ga_sync()
           k = k + 1
c
           if(have_data)then
c
c...          first, get access to the local piece of g_u
              call ga_access(g_u, ilo, ihi, jlo, jhi, index, ld)
c
c...          compute next approximation of u
              call nextu(ni, nj, dbl_mb(index), ld, dbl_mb(indexl),
     &                   ni+1, S, W, C, N, E, dbl_mb(indexd) ) 
c
c...          release access to the local data (the data were updated !)
              call ga_release_update(g_u, ilo, ihi, jlo, jhi)
           endif
c
c...       compute the stopping number -- second norm
           call ga_add(1d0, g_u, -1d0, g_ux, g_diff)    ! g_diff = g_u - g_ux 
           product = ga_ddot(g_diff, g_diff)             ! <g_diff, g_diff>
           norm    = sqrt(product)/uxnorm
           if(me.eq.0 .and. Mod(k,10).eq.1) print *,k,' error= ',norm
c
      if(norm .gt. eps) then
         if(.not. have_data) goto 10
c
c...     not done yet, get a copy of my piece of g_u and neighboring grid points
c        ... determine where from we should copy data -- consider topology 
         index_first = indexl
         if(jlo.eq.1) index_first = index_first + ni+2 
         if(ilo.eq.1) index_first = index_first + 1 
c
         call ga_get(g_u, max(1,ilo-1), min(m,ihi+1),
     &                    max(1,jlo-1), min(m,jhi+1),
     &                    dbl_mb(index_first), ni+2 ) 
         goto 10
      endif
c
      if(me.eq.0) then
         print *,' converged in ',k, ' iterations, error= ',norm
      endif 
c
c***  deallocate MA memory and destroy global arrays
      if(have_data)then
        if(.not. ma_pop_stack(h_l))call ga_error('invalid handle ?',0)
        if(.not. ma_pop_stack(h_d))call ga_error('invalid handle ?',0)
      endif
      if(.not. ga_destroy(g_u)) call ga_error('invalid handle ?',0)
      if(.not. ga_destroy(g_ux)) call ga_error('invalid handle ?',0)
      if(.not. ga_destroy(g_diff)) call ga_error('invalid handle ?',0)
c
      call ga_sync()
      end


      subroutine nextu(ni, nj, u, ld, u0, ld0, S, W, C, N, E, D ) 
      implicit none
      integer ni, nj, ld, ld0
      double precision u(1:ld,1:nj), u0(0:ld0,0:nj), D(1:ni,1:nj)
      double precision S, W, C, N, E
      integer i, j      
c
      do j = 1, nj
         do i = 1, ni
            u(i,j) = S*u0(i+1,j) + W*u0(i,j-1) +  
     &               N*u0(i-1,j) + E*u0(i,j+1) + D(i,j)
         enddo
      enddo
      end


      subroutine jacinit(m, S, W, C, N, E, D, ldd, ux,
     &                   ilo, ihi, jlo, jhi, ldu )
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      double precision S, W, C, N, E, temp
      integer ilo, ihi, jlo, jhi, ldu, ldd, i, j
      integer m, ni, nj, jj, ii 
      double precision h, ux(ldu,m), D(ldd,m)
      if(m.le.1)call ga_error('jacinit: wrong value of m',0) 
      S = .25
      W = .25
      N = .25
      E = .25
      h = 1d0/(m+1)
      ni = ihi - ilo + 1                 ! number of 'local' rows 
      nj = jhi - jlo + 1                 ! number of 'local' columns 
      do j = 1, nj 
         jj = jlo + j -1
         do i = 1, ni 
            temp = 0d0
            ii = ilo + i -1
            if(jj .eq. 1) temp = temp + 1 
            if(ii .eq. 1) temp = temp + 1 
            if(ii .eq. m) temp = temp + jj*h + 1
            if(jj .eq. m) temp = temp + ii*h + 1
            D(i,j) = temp/4d0
            ux(i,j)= 1d0 + h*h*ii*jj
         enddo
      enddo
      end
