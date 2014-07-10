#if HAVE_CONFIG_H
#   include "config.fh"
#endif
# define THRESH  1.0d-10
#define MISMATCH(x,y) abs(x-y)/max(1,abs(x)).gt.THRESH
c
c#define BLOCK_CYCLIC
c#define USE_SCALAPACK_DISTR
#ifdef USE_SCALAPACK_DISTR
#define BLOCK_CYCLIC
#endif
c
      program test
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer TESTDIM
      parameter(TESTDIM = 256)
      logical status
      integer g_a, g_b, g_c
      double precision alpha, beta
      integer ndim, adims(2), bdims(2), cdims(2), tlo(2), thi(2)
      integer alo(2), ahi(2), blo(2), bhi(2), clo(2), chi(2)
      integer ald, bld, cld, i_inc, j_inc
      double precision val
      integer me, nproc, i, j, ii, jj
      GA_ACCESS_INDEX_TYPE idx, inc
c
c***  Initialize a message passing library
c
#include "mp3.fh"
c
      call nga_initialize()
      if(ga_nodeid().eq.0)then
         write(6,*)
         write(6,'(a)') ' GA initialized'
         write(6,*)
         call ffflush(6)
      endif
c
      status = ma_init(MT_DBL, 500000, 900000/ga_nnodes())
      if (.not. status)call ga_error( 'ma_init failed', -1)
c
      me = ga_nodeid()
c
c   create test arrays
c
      g_a = nga_create_handle() 
      ndim = 2
      adims(1) = TESTDIM
      adims(2) = TESTDIM
      call nga_set_data(g_a,ndim,adims,MT_DBL)
      status = ga_allocate(g_a)
c
      g_b = nga_create_handle() 
      ndim = 2
      bdims(1) = TESTDIM + 1
      bdims(2) = TESTDIM + 1
      call nga_set_data(g_b,ndim,bdims,MT_DBL)
      status = ga_allocate(g_b)
c
      g_c = nga_create_handle() 
      ndim = 2
      cdims(1) = TESTDIM + 2
      cdims(2) = TESTDIM + 2
      call nga_set_data(g_c,ndim,cdims,MT_DBL)
      status = ga_allocate(g_c)
c
c  initialize a and b
c
      i_inc = TESTDIM/2
      j_inc = TESTDIM/2
c
      call nga_distribution(g_a,me,alo,ahi)
      call nga_access(g_a,alo,ahi,idx,ald)
      do j = alo(2), ahi(2)
        do i = alo(1), ahi(1)
          dbl_mb(idx) = dble((j-1)*adims(1) + i-1)
          idx = idx + 1
        end do
      end do
      call nga_release(g_a, alo, ahi)
c
      call nga_distribution(g_b,me,blo,bhi)
      call nga_access(g_b,blo,bhi,idx,bld)
      do j = blo(2), bhi(2)
        do i = blo(1), bhi(1)
          dbl_mb(idx) = dble((j-1)*bdims(1) + i-1)
          idx = idx + 1
        end do
      end do
      call nga_release(g_b, blo, bhi)
c
      alo(1) = TESTDIM/4
      alo(2) = TESTDIM/4
      ahi(1) = alo(1) + i_inc
      ahi(2) = alo(2) + j_inc
      blo(1) = TESTDIM/4 + 1
      blo(2) = TESTDIM/4 + 1
      bhi(1) = blo(1) + i_inc
      bhi(2) = blo(2) + j_inc
      clo(1) = TESTDIM/4 + 2
      clo(2) = TESTDIM/4 + 2
      chi(1) = clo(1) + i_inc
      chi(2) = clo(2) + j_inc
      alpha = 1.0d00
      beta = 1.0d00
c
      call nga_add_patch(alpha, g_a, alo, ahi, beta, g_b, blo, bhi,
     +                   g_c, clo, chi)
c      call ga_print(g_a)
c      call ga_print(g_b)
c      call ga_print(g_c)
c
c    check C for answer
c
      call nga_distribution(g_c,me,tlo,thi)
      if (tlo(1).lt.clo(1)) tlo(1) = clo(1)
      if (tlo(2).lt.clo(2)) tlo(2) = clo(2)
      if (thi(1).gt.chi(1)) thi(1) = chi(1)
      if (thi(2).gt.chi(2)) thi(2) = chi(2)
c
      if (tlo(1).le.thi(1).and.tlo(2).le.thi(2)) then
        call nga_access(g_c,tlo,thi,idx,cld)
        do j = tlo(2), thi(2)
          jj = j - tlo(2)
          do i = tlo(1), thi(1)
            ii = i - tlo(1)
            val = alpha*dble((j-3)*adims(1)+i-3)
     +          + beta*dble((j-2)*bdims(1)+i-2)
            if (dbl_mb(idx+jj*cld+ii).ne.val) then
              write(6,'(i4,a,2i8,2f8.0)') me,' Mismatch for values: ',
     +            i,j,dbl_mb(idx+jj*cld+ii),val
            endif
          end do
        end do
        call nga_release(g_c, tlo, thi)
      endif
c
      if (me.eq.0) then
        write(6,'(a)') 'Successfully completed test of nga_add_patch'
      endif
c
      status = nga_destroy(g_a)
      status = nga_destroy(g_b)
      status = nga_destroy(g_c)
      call nga_terminate()
c
      call MP_FINALIZE()
      end
