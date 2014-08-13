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


      program main
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c#include "tcgmsg.fh"
      integer heap, stack, fudge, ma_heap, me
      integer nmax, DIM, nwidth, MAXPROC, nloop, nskip
      parameter (nmax = 174, DIM = 2, nwidth = 0, MAXPROC = 2000)
      parameter (nloop = 1, nskip = 2)
      integer ndim, nproc, pdims(7), type, g_a, maxval
      integer i, j, ichk, ii, jj, ival, dims(7), width(7)
      integer lo(7), hi(7), ld(7)
      integer chunk(7), dim_strided(7), lop(7), hip(7)
      integer a(nmax, nmax), skip(7)
      parameter (heap=500*500*4, fudge=100, stack=500*500)
      logical status
      integer one
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
      print*
      call ga_initialize()
      nproc = ga_nnodes()
      me = ga_nodeid()
c     we can also use GA_set_memory_limit BEFORE first ga_create call
c
      ma_heap = heap + fudge 
      call GA_set_memory_limit(util_mdtob(ma_heap))
c
      if(ga_nodeid().eq.0)then
         write(6,'(a)') ' GA initialized '
         call ffflush(6)
      endif
c
c***  Initialize the MA package
c     MA must be initialized before any global array is allocated
c
      status = ma_init(MT_DCPL, stack, ma_heap)
      if (.not. status) call ga_error('ma_init failed',-1) 
c
      if(me.eq.0)then
        write(6,'(a,i8,a)') 'using ', nproc, ' process(es)'
        call ffflush(6)
      endif
c
c   Test strided put, get, acc calls
c
      ndim = DIM
      dims(1) = nmax
      dims(2) = nmax
      maxval = 1
      do i = 1, ndim
        maxval = dims(i)*maxval
      end do
      maxval = maxval - 1
      do j = 1, nmax
        do i = 1, nmax
          a(i,j) = 0
        end do
      end do
c
      do i = 1, ndim
        width(i) = nwidth
        skip(i) = nskip + (i-1)
        chunk(i) = 1
        ld(i) = nmax
      end do

      type = MT_INT
      status = nga_create_ghosts(type, ndim, dims, width,
     +                           "test_array", chunk, g_a)
      if (status.and.me.eq.0) then
        write(6,'(a)') '*'
        write(6,'(a)') '* Global array creation was successful'
        write(6,'(a)') '*'
      elseif (.not.status) then
        write(6,'(a,i8)') 'Global array creation failure on ',me
      endif
      call ga_zero(g_a)
c
c     initialize g_a
c
      call ga_sync
      call nga_distribution(g_a, me, lo, hi)
c      write(6,'(i4,a,i8,a,i8,a,i8,a,i8)') me, ' (dist) lo(1): ',lo(1),
c     +  ' hi(1): ',hi(1),' lo(2): ',lo(2),' hi(2): ',hi(2)
c
c     fill local array a with unique values for each element
c
      do i = 1, hi(1) - lo(1) + 1
        do j = 1, hi(2) - lo(2) + 1
          a(i,j) = (j + lo(2) - 2)*dims(1) + (i + lo(1) - 2) + 1
        end do
      end do
      call nga_put(g_a, lo, hi, a, ld)
      call ga_sync
c
c     reset all values in local array a to zero
c
      do i = 1, hi(1) - lo(1) + 1
        do j = 1, hi(2) - lo(2) + 1
          a(i,j) = 0
        end do
      end do
      do i = 1, ndim
        lo(i) = 1
        hi(i) = dims(i)
      end do
      call random_range(lo,hi,lop,hip,ndim)
c      write(6,'(i4,a,i8,a,i8,a,i8,a,i8)') me, ' lop(1): ',lop(1),
c     +  ' hip(1): ',hip(1),' lop(2): ',lop(2),' hip(2): ',hip(2)
      call nga_strided_get(g_a, lop, hip, skip, a, ld)
c
c  check results
c
      do i = 1, ndim
        dim_strided(i) = hip(i) - lop(i)
        dim_strided(i) = dim_strided(i) - mod(dim_strided(i),skip(i))
        dim_strided(i) = dim_strided(i)/skip(i) + 1
      end do
      ichk = 0
      do j = 1, dim_strided(2)
        jj = (j-1)*skip(2) + lop(2)
        do i = 1, dim_strided(1)
          ii = (i-1)*skip(1) + lop(1)
          ival = (ii-1) + (jj-1)*dims(1) + 1
          if (a(i,j).ne.ival.and.me.eq.0) then
            write(6,'(a,i4,a,i4,a,i8,a,i8,a,i8)') 'a(',ii,',',jj,'): ',
     +                        a(i,j),' ival: ',ival,' dims(1): ',dims(1)
            ichk = ichk + 1
          endif
        end do
      end do
      if (ichk.eq.0.and.me.eq.0) then
        write(6,'(a)') '*'
        write(6,'(a)') '* nga_strided_get is 0K'
        write(6,'(a)') '*'
      endif
c
#if 1
      call ga_sync
      call ga_zero(g_a)
      if (me.eq.1) then
        call nga_strided_put(g_a, lop, hip, skip, a, ld)
        one = 1
        call nga_strided_acc(g_a, lop, hip, skip, a, ld, one)
      endif
      call ga_sync
      do i = 1, hi(1) - lo(1) + 1
        do j = 1, hi(2) - lo(2) + 1
          a(i,j) = 0
        end do
      end do
      call nga_get(g_a, lop, hip, a, ld)
      do i = 1, ndim
        dim_strided(i) = hip(i) - lop(i)
      end do
      ichk = 0
      do j = 1, dim_strided(2)
        jj = (j-1) + lop(2)
        do i = 1, dim_strided(1)
          ii = (i-1) + lop(1)
          ival = (ii-1) + (jj-1)*dims(1) + 1
          ival = 2*ival
          if (mod(i-1,skip(1)).eq.0.and.mod(j-1,skip(2)).eq.0) then
            if (a(i,j).ne.ival.and.me.eq.0) then
              write(6,'(a,i4,a,i4,a,i8,a,i8,a,i8)') 'a2(',ii,',',jj,
     +             '): ',a(i,j),' ival: ',ival,' dims(1): ',dims(1)
              ichk = ichk + 1
            endif
          else
            if (a(i,j).ne.0.0d00) ichk = ichk + 1
          endif
        end do
      end do
      if (ichk.eq.0.and.me.eq.0) then
        write(6,'(a)') '*'
        write(6,'(a)') '* nga_strided_put and nga_strided_acc is 0K'
        write(6,'(a)') '*'
      endif
#endif
c
c   print out values of a
c
c     call ga_print(g_a)
c
c***  Tidy up the GA package
c
      call ga_terminate()
c
c***  Tidy up after message-passing library
c
      call MP_FINALIZE()

c
      stop
      end
c
      subroutine aprint(a,nrow,ncol,ld,has_data)
#include "global.fh"
      integer ld
      integer a(ld,*)
      integer i, j, k, nproc
      logical has_data(0:1999)
       
      nproc = ga_nnodes()
      do k = 1, nproc
        call ga_sync
        if (k-1.eq.ga_nodeid().and.has_data(k-1)) then
          write(6,*) '*'
          write(6,*) '*   Data on processor ',k-1
          write(6,*) '*'
          do i = 1, min(nrow,12)
            write (6,102) (a(i,j), j = 1, min(ncol,12))
  102       format(14i5)
          end do
        endif
        call ffflush(6)
      enddo
c
      return
      end
c
      subroutine atest(a,nrow,ncol,ld,b,ld2,has_data)
#include "global.fh"
      integer ld
      integer a(ld,*), b(ld2,*)
      integer i, j, k, nproc
      logical has_data(0:1999), check_data

      nproc = ga_nnodes()
      check_data = .true.
      do k = 1, nproc
        call ga_sync
        if (k-1.eq.ga_nodeid().and.has_data(k-1)) then
          do i = 1, nrow
            do j = 1, ncol
              if (a(i,j).ne.b(i,j)) check_data = .false.
            end do
          end do
          if (check_data) then
            write(6,*) '*'
            write(6,*) '*   Data from nga_access_ghosts and'
            write(6,*) '*   nga_periodic_get is the same on'
            write(6,*) '*   processor ',k-1
            write(6,*) '*'
          else
            write(6,*) '*'
            write(6,*) '*   Data from nga_access_ghosts and'
            write(6,*) '*   nga_periodic_get is NOT the same on'
            write(6,*) '*   processor ',k-1
            write(6,*) '*'
          endif
        endif
        call ffflush(6)
      enddo
c
      return
      end
c
      subroutine factor(p,ndim,dims)
      implicit none
      integer i,j,p,ndim,dims(7),imin,mdim
      integer ip,ifac,pmax,prime(1000)
      integer fac(1000)
c
      i = 1
      ip = p
      do i = 1, ndim
        dims(i) = 1
      end do
c
c    factor p completely
c    first, find all prime numbers less than or equal to p
c
      pmax = 0
      do i = 2, p
        do j = 1, pmax
          if (mod(i,prime(j)).eq.0) go to 100
        end do
        pmax = pmax + 1
        prime(pmax) = i
  100   continue
      end do
c
c    find all prime factors of p
c
      ifac = 0
      do i = 1, pmax
  200   if (mod(ip,prime(i)).eq.0) then
          ifac = ifac + 1
          fac(ifac) = prime(i)
          ip = ip/prime(i)
          go to 200
        endif
      end do
c
c    determine dimensions of processor grid
c
      do i = ifac, 1, -1
c
c    find dimension with minimum value
c
        imin = dims(1)
        mdim = 1
        do j = 2, ndim
          if (dims(j).lt.imin) then
            imin = dims(j)
            mdim = j
          endif
        end do
        dims(mdim) = dims(mdim)*fac(i)
      end do
c
      return
      end
c
c    create a random patch
c
      subroutine random_range(lo,hi,lop,hip,ndim)
      implicit none
#include "testutil.fh"
      integer lo(1),hi(1),lop(1),hip(1),ndim
      integer i, range, swap, val,iran
      iran(range) = int(drand(0)*dble(range)) + 1
      do i = 1, ndim
         range = hi(i)-lo(i)+1
         val = iran(range)
         lop(i) = lo(i) + val
         val = iran(range)
         hip(i) = hi(i) - val
         if(hip(i) .lt. lop(i))then
            swap =hip(i)
            hip(i)=lop(i)
            lop(i)=swap
         endif
         hip(i)=MIN(hip(i),hi(i))
         lop(i)=MAX(lop(i),lo(i))
      enddo
      end
