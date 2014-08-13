#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c
c     $Id: maf.F,v 1.14 2003-07-10 19:19:28 d3h325 Exp $
c

c
c     FORTRAN routines for a portable dynamic memory allocator.
c

#define MAF_INTERNAL

cc
cc    private routines
cc

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      integer function ma_set_sizes ()

      implicit none

#include "maf2c.fh"
#include "mafdecls.fh"

      ma_set_sizes = 0
      if (f2c_inform_base(MT_BYTE, byte_mb(1), byte_mb(2)) .eq.
     $    MA_FALSE) return
      if (f2c_inform_base(MT_INT, int_mb(1), int_mb(2)) .eq.
     $    MA_FALSE) return
      if (f2c_inform_base(MT_LOG, log_mb(1), log_mb(2)) .eq.
     $    MA_FALSE) return
      if (f2c_inform_base(MT_REAL, real_mb(1), real_mb(2)) .eq.
     $    MA_FALSE) return
      if (f2c_inform_base(MT_DBL, dbl_mb(1), dbl_mb(2)) .eq.
     $    MA_FALSE) return
      if (f2c_inform_base(MT_SCPL, scpl_mb(1), scpl_mb(2)) .eq.
     $    MA_FALSE) return
      if (f2c_inform_base(MT_DCPL, dcpl_mb(1), dcpl_mb(2)) .eq.
     $    MA_FALSE) return

      ma_set_sizes = 1

      return
      end

cc
cc    public routines
cc

c     In general, each routine simply calls its corresponding f2c_ C
c     wrapper routine, which performs any necessary argument munging
c     and then calls the corresponding C routine.

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_alloc_get (datatype, nelem, name, memhandle,
     $    index)

      implicit none
#include "mafdecls.fh"

      integer datatype
      integer nelem
      character*(*) name
      integer memhandle
      MA_ACCESS_INDEX_TYPE index

#include "maf2c.fh"

      if (f2c_alloc_get(datatype, nelem, name, memhandle, index) .eq.
     $    MA_TRUE) then
          MA_alloc_get = .true.
      else
          MA_alloc_get = .false.
      endif

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_allocate_heap (datatype, nelem, name,
     $    memhandle)

      implicit none

      integer datatype
      integer nelem
      character*(*) name
      integer memhandle

#include "maf2c.fh"

      if (f2c_allocate_heap(datatype, nelem, name, memhandle) .eq.
     $    MA_TRUE) then
          MA_allocate_heap = .true.
      else
          MA_allocate_heap = .false.
      endif

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_chop_stack (memhandle)

      implicit none

      integer memhandle

#include "maf2c.fh"

      if (f2c_chop_stack(memhandle) .eq. MA_TRUE) then
          MA_chop_stack = .true.
      else
          MA_chop_stack = .false.
      endif

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_free_heap (memhandle)

      implicit none

      integer memhandle

#include "maf2c.fh"

      if (f2c_free_heap(memhandle) .eq. MA_TRUE) then
          MA_free_heap = .true.
      else
          MA_free_heap = .false.
      endif

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_free_heap_piece (memhandle, nelem)

      implicit none

      integer memhandle
      integer nelem

#include "maf2c.fh"

      if (f2c_free_heap_piece(memhandle, nelem) .eq. MA_TRUE) then
          MA_free_heap_piece = .true.
      else
          MA_free_heap_piece = .false.
      endif

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_get_index (memhandle, index)

      implicit none
#include "mafdecls.fh"

      integer memhandle
      MA_ACCESS_INDEX_TYPE index

#include "maf2c.fh"

      if (f2c_get_index(memhandle, index) .eq. MA_TRUE) then
          MA_get_index = .true.
      else
          MA_get_index = .false.
      endif

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_get_next_memhandle (ithandle, memhandle)

      implicit none

      integer ithandle
      integer memhandle

#include "maf2c.fh"

      if (f2c_get_next_memhandle(ithandle, memhandle) .eq. MA_TRUE)
     $    then
          MA_get_next_memhandle = .true.
      else
          MA_get_next_memhandle = .false.
      endif

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_get_numalign (numalign)

      implicit none

      integer numalign

#include "maf2c.fh"

      MA_get_numalign = (f2c_get_numalign(numalign) .eq. MA_TRUE)

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_init (datatype, stack, heap)

      implicit none

      integer datatype
      integer stack
      integer heap

#include "maf2c.fh"

      if (f2c_init(datatype, stack, heap) .eq. MA_TRUE) then
          MA_init = .true.
      else
          MA_init = .false.
      endif

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_initialized ()

      implicit none

#include "maf2c.fh"

      if (f2c_initialized() .eq. MA_TRUE) then
          MA_initialized = .true.
      else
          MA_initialized = .false.
      endif

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_init_memhandle_iterator (ithandle)

      implicit none

      integer ithandle

#include "maf2c.fh"

      if (f2c_init_memhandle_iterator(ithandle) .eq. MA_TRUE)
     $    then
          MA_init_memhandle_iterator = .true.
      else
          MA_init_memhandle_iterator = .false.
      endif

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      integer function MA_inquire_avail (datatype)

      implicit none

      integer datatype

#include "maf2c.fh"

      MA_inquire_avail = f2c_inquire_avail(datatype)

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      integer function MA_inquire_heap (datatype)

      implicit none

      integer datatype

#include "maf2c.fh"

      MA_inquire_heap = f2c_inquire_heap(datatype)

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      integer function MA_inquire_heap_check_stack (datatype)

      implicit none

      integer datatype

#include "maf2c.fh"

      MA_inquire_heap_check_stack =
     $    f2c_inquire_heap_check_stack(datatype)

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      integer function MA_inquire_heap_no_partition (datatype)

      implicit none

      integer datatype

#include "maf2c.fh"

      MA_inquire_heap_no_partition =
     $    f2c_inquire_heap_no_partition(datatype)

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      integer function MA_inquire_stack (datatype)

      implicit none

      integer datatype

#include "maf2c.fh"

      MA_inquire_stack = f2c_inquire_stack(datatype)

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      integer function MA_inquire_stack_check_heap (datatype)

      implicit none

      integer datatype

#include "maf2c.fh"

      MA_inquire_stack_check_heap =
     $    f2c_inquire_stack_check_heap(datatype)

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      integer function MA_inquire_stack_no_partition (datatype)

      implicit none

      integer datatype

#include "maf2c.fh"

      MA_inquire_stack_no_partition =
     $    f2c_inquire_stack_no_partition(datatype)

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_pop_stack (memhandle)

      implicit none

      integer memhandle

#include "maf2c.fh"

      if (f2c_pop_stack(memhandle) .eq. MA_TRUE) then
          MA_pop_stack = .true.
      else
          MA_pop_stack = .false.
      endif

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      subroutine MA_print_stats(printroutines)

      implicit none

      logical printroutines

#include "maf2c.fh"

      if (printroutines) then
          call f2c_print_stats(1)
      else
          call f2c_print_stats(0)
      endif

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_push_get (datatype, nelem, name, memhandle,
     $    index)

      implicit none
#include "mafdecls.fh"

      integer datatype
      integer nelem
      character*(*) name
      integer memhandle
      MA_ACCESS_INDEX_TYPE index

#include "maf2c.fh"

      if (f2c_push_get(datatype, nelem, name, memhandle, index) .eq.
     $    MA_TRUE) then
          MA_push_get = .true.
      else
          MA_push_get = .false.
      endif

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_push_stack (datatype, nelem, name,
     $    memhandle)

      implicit none

      integer datatype
      integer nelem
      character*(*) name
      integer memhandle

#include "maf2c.fh"

      if (f2c_push_stack(datatype, nelem, name, memhandle) .eq.
     $    MA_TRUE) then
          MA_push_stack = .true.
      else
          MA_push_stack = .false.
      endif

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_set_auto_verify (value)

      implicit none

      logical value
      integer ivalue

#include "maf2c.fh"

      if (value) then
          ivalue = MA_TRUE
      else
          ivalue = MA_FALSE
      endif

      if (f2c_set_auto_verify(ivalue) .eq. MA_TRUE) then
          MA_set_auto_verify = .true.
      else
          MA_set_auto_verify = .false.
      endif

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_set_error_print (value)

      implicit none

      logical value
      integer ivalue

#include "maf2c.fh"

      if (value) then
          ivalue = MA_TRUE
      else
          ivalue = MA_FALSE
      endif

      if (f2c_set_error_print(ivalue) .eq. MA_TRUE) then
          MA_set_error_print = .true.
      else
          MA_set_error_print = .false.
      endif

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_set_hard_fail (value)

      implicit none

      logical value
      integer ivalue

#include "maf2c.fh"

      if (value) then
          ivalue = MA_TRUE
      else
          ivalue = MA_FALSE
      endif

      if (f2c_set_hard_fail(ivalue) .eq. MA_TRUE) then
          MA_set_hard_fail = .true.
      else
          MA_set_hard_fail = .false.
      endif

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_set_numalign (numalign)

      implicit none

      integer numalign

#include "maf2c.fh"

      MA_set_numalign = (f2c_set_numalign(numalign) .eq. MA_TRUE)

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      integer function MA_sizeof (datatype1, nelem1, datatype2)

      implicit none

      integer datatype1
      integer nelem1
      integer datatype2

#include "maf2c.fh"

      MA_sizeof = f2c_sizeof(datatype1, nelem1, datatype2)

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      integer function MA_sizeof_overhead (datatype)

      implicit none

      integer datatype

#include "maf2c.fh"

      MA_sizeof_overhead = f2c_sizeof_overhead(datatype)

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      subroutine MA_summarize_allocated_blocks

      implicit none

#include "maf2c.fh"

      call f2c_summarize_allocated_blocks()

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      subroutine ma_trace(value)

      implicit none

      integer value

#include "maf2c.fh"

      call f2c_trace(value)

      return
      end

c     --------------------------------------------------------------- c
c     --------------------------------------------------------------- c

      logical function MA_verify_allocator_stuff ()

      implicit none

#include "maf2c.fh"

      if (f2c_verify_allocator_stuff() .eq. MA_TRUE) then
          MA_verify_allocator_stuff = .true.
      else
          MA_verify_allocator_stuff = .false.
      endif

      return
      end

#undef MAF_INTERNAL
