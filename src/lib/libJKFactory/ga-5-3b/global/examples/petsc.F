#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      program petsc
      implicit none
#include "mafdecls.fh"
#include "global.fh"

      integer heap, stack
      data heap,stack /2*40000/
c
c***  Intitialize a message passing library
c
#ifdef MPI
      call mpi_init
#else
      call pbeginf
#endif
c
      call ga_initialize()                               ! initialize GAs 
      if (.not. ma_init(MT_DBL, heap, stack))            ! initialize MA
     $    call ga_error("ma init failed",heap+stack)
#ifdef ENABLE_TRACE
      call trace_init(10000)                 ! initialize trace
#endif
      call iterate()                                     ! do the work
#ifdef ENABLE_TRACE
      call trace_end(ga_nodeid())             ! end trace
#endif
      call ga_terminate()                                ! terminate GAs 
c
#ifdef MPI
      call mpi_finalize()
#else
      call pend()
#endif
      end

      subroutine iterate()
      implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/vec.h"
#include "include/finclude/mat.h"
#include "include/finclude/pc.h"
#include "include/finclude/ksp.h"
#include "include/finclude/sles.h"
#include "include/finclude/sys.h"

#include "mafdecls.fh"
#include "global.fh"

      double precision  norm
      integer     i, j, II, JJ, ierr, m, n
      parameter (m = 3)
      parameter (n = 3)
      integer     its, Istart, Iend, flg
      integer     me, nproc
      Scalar      v, one, neg_one
      Vec         x, b, u
      Mat         A 
      SLES        sles
      KSP         ksp
      PetscRandom rctx

      integer g_x
      integer ld
      double precision buf_v(1)
      PetscOffset idx

c$$$      PC          pc
c$$$      PCType      ptype 
c$$$      double precision tol

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                 Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c***  check parallel environment
      me = ga_nodeid()
      nproc = ga_nnodes()

      one  = 1.0
      neg_one = -1.0
      ld = m*n

c***  Global Array
c***  create global arrays: g_x - approx. solution
      if (.not. ga_create(MT_DBL, m*n, 1, 'x', 1, 1, g_x))
     $     call ga_error(' ga_create failed ',0)
c
c***  initial guess for x -- zero
      call ga_zero(g_x)
c
c$$$      do i=1,m*n
c$$$         buf(i) = i
c$$$      enddo
c$$$      call ga_put(g_x,1,m*n,1,1,buf,ld)
c
c***  PETSC      
      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!      Compute the matrix and right-hand-side vector that define
!      the linear system, Ax = b.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

c***  Matrix A
      call MatCreate(PETSC_COMM_WORLD,m*n,m*n,A,ierr)
      call MatGetOwnershipRange(A,Istart,Iend,ierr)

      do II=Istart,Iend-1
        v = -1.0
        i = II/n
        j = II - i*n  
        if ( i.gt.0 ) then
          JJ = II - n
          call MatSetValues(A,1,II,1,JJ,v,ADD_VALUES,ierr)
        endif
        if ( i.lt.m-1 ) then
          JJ = II + n
          call MatSetValues(A,1,II,1,JJ,v,ADD_VALUES,ierr)
        endif
        if ( j.gt.0 ) then
          JJ = II - 1
          call MatSetValues(A,1,II,1,JJ,v,ADD_VALUES,ierr)
        endif
        if ( j.lt.n-1 ) then
          JJ = II + 1
          call MatSetValues(A,1,II,1,JJ,v,ADD_VALUES,ierr)
        endif
        v = 4.0
        call  MatSetValues(A,1,II,1,II,v,ADD_VALUES,ierr)
      enddo

      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

c***  Vector b
      call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,m*n,u,ierr)
      call VecSetFromOptions(u,ierr)
      call VecDuplicate(u,b,ierr)
      call VecDuplicate(b,x,ierr)

c***  u is the exact solution
      call VecSet(one,u,ierr)
      if (me .eq. 0) print *, 'Exact solution:'
      call VecView(u, VIEWER_STDOUT_WORLD, ierr)

c***  b is the right hand side
      call MatMult(A,u,b,ierr)

c***  Manage to make connection of ga to petsc: g_x -> x
      call VecGetOwnershipRange(x,Istart,Iend,ierr)
      call VecGetArray(x,buf_v,idx,ierr)

      call ga_get(g_x,Istart+1,Iend,1,1,buf_v(idx+1),ld)

      call VecRestoreArray(x,buf_v,idx,ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!         Create the linear solver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      call SLESCreate(PETSC_COMM_WORLD,sles,ierr)

      call SLESSetOperators(sles,A,A,DIFFERENT_NONZERO_PATTERN,
     &                      ierr)

c$$$      call SLESGetPC(sles,pc,ierr)
c$$$      ptype = PCJACOBI
c$$$      call PCSetType(pc,ptype,ierr)
c$$$      call SLESGetKSP(sles,ksp,ierr)
c$$$      tol = 1.e-7
c$$$      call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_DOUBLE_PRECISION,
c$$$     &     PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,ierr)

      call SLESSetFromOptions(sles,ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                      Solve the linear system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      call SLESSolve(sles,b,x,its,ierr)

      call ga_sync()
c***  write the approx solution back to ga
      call VecGetArray(x,buf_v,idx,ierr)

      call ga_put(g_x,Istart+1,Iend,1,1,buf_v(idx+1),ld)

      call VecRestoreArray(x,buf_v,idx,ierr)
      
      if (me .eq. 0) print *, 'Approx solution:'
      call ga_print(g_x)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!                     Check solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!  Check the error

      call VecAXPY(neg_one,u,x,ierr)
      call VecNorm(x,NORM_2,norm,ierr)
      if (me .eq. 0) then
        if (norm .gt. 1.e-12) then
           write(6,100) norm, its
        else
           write(6,110) its
        endif
      endif
  100 format('Norm of error ',e10.4,' iterations ',i5)
  110 format('Norm of error < 1.e-12, iterations ',i5)

c***  clean up
      call SLESDestroy(sles,ierr)
      call VecDestroy(u,ierr)
      call VecDestroy(x,ierr)
      call VecDestroy(b,ierr)
      call MatDestroy(A,ierr)

      call PetscFinalize(ierr)

      if(.not. ga_destroy(g_x)) call ga_error('invalid handle ?',0)
c
      call ga_sync()
      end

