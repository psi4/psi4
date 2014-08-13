#if HAVE_CONFIG_H
#   include "config.fh"
#endif

#define VERBOSE 0
#define NELEM 20

      program main
      implicit none
#include "global.fh"
      integer nproc, me
#include "mp3.fh"
      call ga_initialize()
      nproc = ga_nnodes()
      me = ga_nodeid()
#if VERBOSE
      if (0.eq.me) then
          write(6,*) 'TESTING INTEGER'
      endif
#endif
      call test_int()
#if VERBOSE
      if (0.eq.me) then
          write(6,*) 'TESTING REAL'
      endif
#endif
      call test_real()
#if VERBOSE
      if (0.eq.me) then
          write(6,*) 'TESTING DOUBLE PRECISION'
      endif
#endif
      call test_dbl()
#if VERBOSE
      if (0.eq.me) then
          write(6,*) 'TESTING COMPLEX'
      endif
#endif
      call test_scpl()
#if VERBOSE
      if (0.eq.me) then
          write(6,*) 'TESTING DOUBLE COMPLEX'
      endif
#endif
      call test_dcpl()
      call ga_terminate()
      call MP_FINALIZE()
      end

      subroutine test_int()
#include "global.fh"
#include "mafdecls.fh"
      integer g_a
      integer lo(1), dims(1), chunk(1), ld(1)
      logical rc
      integer buf(NELEM)
      integer get(NELEM)
      integer start
      integer inc
      start = 1
      inc = 1
      lo(1) = 1
      dims(1) = NELEM
      chunk(1) = 0
      rc = nga_create(MT_INT, 1, dims, "g_a", chunk, g_a)
      if (.not. rc) then
          write(6,*) 'nga_create failed'
      endif
      call ga_patch_enum(g_a, 1, NELEM, start, inc)
      call nga_get(g_a, lo, dims, get, ld)
      do i = 1, NELEM
          buf(i) = start + (i-1)*inc
      end do
      do i = 1, NELEM
          if (buf(i).ne.get(i)) then
              write(6,*) 'at',i,'value mismatch',buf(i),'.ne.',get(i)
              call ga_error('error', i)
          endif
      end do
      rc = ga_destroy(g_a)
      if (.not. rc) then
          write(6,*) 'ga_destroy failed'
      endif
      end

      subroutine test_real()
#include "global.fh"
#include "mafdecls.fh"
      integer g_a
      integer lo(1), dims(1), chunk(1), ld(1)
      logical rc
      real buf(NELEM)
      real get(NELEM)
      real start
      real inc
      start = 1
      inc = 1
      lo(1) = 1
      dims(1) = NELEM
      chunk(1) = 0
      rc = nga_create(MT_REAL, 1, dims, "g_a", chunk, g_a)
      if (.not. rc) then
          write(6,*) 'nga_create failed'
      endif
      call ga_patch_enum(g_a, 1, NELEM, start, inc)
      call nga_get(g_a, lo, dims, get, ld)
      do i = 1, NELEM
          buf(i) = start + (i-1)*inc
      end do
      do i = 1, NELEM
          if (buf(i).ne.get(i)) then
              write(6,*) 'at',i,'value mismatch',buf(i),'.ne.',get(i)
              call ga_error('error', i)
          endif
      end do
      rc = ga_destroy(g_a)
      if (.not. rc) then
          write(6,*) 'ga_destroy failed'
      endif
      end

      subroutine test_dbl()
#include "global.fh"
#include "mafdecls.fh"
      integer g_a
      integer lo(1), dims(1), chunk(1), ld(1)
      logical rc
      double precision buf(NELEM)
      double precision get(NELEM)
      double precision start
      double precision inc
      start = 1
      inc = 1
      lo(1) = 1
      dims(1) = NELEM
      chunk(1) = 0
      rc = nga_create(MT_DBL, 1, dims, "g_a", chunk, g_a)
      if (.not. rc) then
          write(6,*) 'nga_create failed'
      endif
      call ga_patch_enum(g_a, 1, NELEM, start, inc)
      call nga_get(g_a, lo, dims, get, ld)
      do i = 1, NELEM
          buf(i) = start + (i-1)*inc
      end do
      do i = 1, NELEM
          if (buf(i).ne.get(i)) then
              write(6,*) 'at',i,'value mismatch',buf(i),'.ne.',get(i)
              call ga_error('error', i)
          endif
      end do
      rc = ga_destroy(g_a)
      if (.not. rc) then
          write(6,*) 'ga_destroy failed'
      endif
      end

      subroutine test_scpl()
#include "global.fh"
#include "mafdecls.fh"
      integer g_a
      integer lo(1), dims(1), chunk(1), ld(1)
      logical rc
      complex buf(NELEM)
      complex get(NELEM)
      complex start
      complex inc
      start = cmplx(1,1)
      inc = cmplx(1,1)
      lo(1) = 1
      dims(1) = NELEM
      chunk(1) = 0
      rc = nga_create(MT_SCPL, 1, dims, "g_a", chunk, g_a)
      if (.not. rc) then
          write(6,*) 'nga_create failed'
      endif
      call ga_patch_enum(g_a, 1, NELEM, start, inc)
      call nga_get(g_a, lo, dims, get, ld)
      do i = 1, NELEM
          buf(i) = start + (i-1)*inc
      end do
      do i = 1, NELEM
          if (buf(i).ne.get(i)) then
              write(6,*) 'at',i,'value mismatch',buf(i),'.ne.',get(i)
              call ga_error('error', i)
          endif
      end do
      rc = ga_destroy(g_a)
      if (.not. rc) then
          write(6,*) 'ga_destroy failed'
      endif
      end

      subroutine test_dcpl()
#include "global.fh"
#include "mafdecls.fh"
      integer g_a
      integer lo(1), dims(1), chunk(1), ld(1)
      logical rc
      double complex buf(NELEM)
      double complex get(NELEM)
      double complex start
      double complex inc
      start = dcmplx(1,1)
      inc = dcmplx(1,1)
      lo(1) = 1
      dims(1) = NELEM
      chunk(1) = 0
      rc = nga_create(MT_DCPL, 1, dims, "g_a", chunk, g_a)
      if (.not. rc) then
          write(6,*) 'nga_create failed'
      endif
      call ga_patch_enum(g_a, 1, NELEM, start, inc)
      call nga_get(g_a, lo, dims, get, ld)
      do i = 1, NELEM
          buf(i) = start + (i-1)*inc
      end do
      do i = 1, NELEM
          if (buf(i).ne.get(i)) then
              write(6,*) 'at',i,'value mismatch',buf(i),'.ne.',get(i)
              call ga_error('error', i)
          endif
      end do
      rc = ga_destroy(g_a)
      if (.not. rc) then
          write(6,*) 'ga_destroy failed'
      endif
      end
