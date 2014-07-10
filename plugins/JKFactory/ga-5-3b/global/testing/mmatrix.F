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

#define USE_FAST_MERGE 0

      program main
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
      integer heap, stack, fudge, ma_heap, me
      integer nmax, DIM, nwidth, MAXPROC, nloop
      parameter (nmax = 1025, DIM = 2, MAXPROC = 2000)
      parameter (nloop = 50)
      integer ndim, nprocs, tprocs, pdims(7), type, dcnt, maxval
      integer g_a, g_b, g_c, g_aa, g_bb, g_cc, iproc
      integer inode, tnodes, p_m, ix, iy, xinc, yinc
      integer i, j, k, l, ii, jj, dims(7), mdim, ndims(7)
      integer lo(7), hi(7), ld(7)
      integer lo2(7), hi2(7)
      integer lo3(7), hi3(7)
      integer lo4(7), hi4(7)
      integer lo5(7), hi5(7)
      integer chunk(7), imax, jmax, iloop, afail, bfail, cfail
      double precision a(nmax, nmax), b(nmax,nmax), c(nmax,nmax)
      double precision ta(nmax, nmax), tb(nmax,nmax),tc(nmax,nmax)
      double precision tat(nmax, nmax)
      double precision start,start1,t1,t2,t3,t4,t5,t6,t7
      double precision alpha, beta
#if USE_FAST_MERGE
      integer mblock,iblock
#endif
      logical status, ltest
      parameter (heap=400*400*4, fudge=100, stack=400*400*4)
c
c***  Intitialize a message passing library
c
#include "mp3.fh"
      mdim = nmax
c
c***  Initialize GA
c
c     There are 2 choices: ga_initialize or ga_initialize_ltd.
c     In the first case, there is no explicit limit on memory usage.
c     In the second, user can set limit (per processor) in bytes.
c
      call ga_initialize()
      me = ga_nodeid()
c     we can also use GA_set_memory_limit BEFORE first ga_create call
c
      ma_heap = heap + fudge 
c      call GA_set_memory_limit(util_mdtob(ma_heap))
c
      if(ga_nodeid().eq.0)then
         print *,' GA initialized '
         call ffflush(6)
      endif
c
c***  Initialize the MA package
c     MA must be initialized before any global array is allocated
c
      status = ma_init(MT_DBL, stack, ma_heap)
      if (.not. status) call ga_error('ma_init failed',-1) 
c
      if(me.eq.0)then
        print *, 'using ', ga_nnodes(), ' process(es)'
        call ffflush(6)
      endif
#if USE_FAST_MERGE
      if(me.eq.0)then
        print *,'*'
        print *,'* Using fast merge algorithm'
        print *,'*'
      endif
#endif
c
c   Test matrix multiply for mirrored arrays
c
      ndim = 2
c
c   processor distribution on node
c
      inode = ga_cluster_nodeid()
      nprocs = ga_cluster_nprocs(inode)
      tnodes = ga_cluster_nnodes()
      tprocs = ga_nnodes()
      iproc = mod(me,nprocs)
      call factor(tnodes,ndim,ndims)
      call factor(nprocs,ndim,pdims)
c
      if (me.eq.0) then
        write(6,*) '*'
        write(6,*) '*  Number of processors: ',tprocs
        write(6,*) '*  Number of processors/node : ',nprocs
        write(6,*) '*  Matrix dimension : ',mdim
        write(6,*) '*'
      endif
c
c   create global arrays
c
      do i = 1, ndim
        chunk(i) = 0
        dims(i) = mdim
        lo(i) = 1
        hi(i) = mdim
        ld(i) = nmax
      end do
      type = MT_DBL
      p_m = ga_pgroup_get_mirror()
      status = nga_create(type, ndim, dims,"aa", chunk, g_aa)
      if (status.and.me.eq.0) then
        write(6,*) '*'
        write(6,*) '* Global array AA creation was successful'
        write(6,*) '*'
      elseif (.not.status) then
        write(6,*) 'Global array creation failure on ',me
      endif
      status = nga_create(type, ndim, dims,"bb", chunk, g_bb)
      if (status.and.me.eq.0) then
        write(6,*) '*'
        write(6,*) '* Global array BB creation was successful'
        write(6,*) '*'
      elseif (.not.status) then
        write(6,*) 'Global array creation failure on ',me
      endif
      status = nga_create(type, ndim, dims,"cc", chunk, g_cc)
      if (status.and.me.eq.0) then
        write(6,*) '*'
        write(6,*) '* Global array CC creation was successful'
        write(6,*) '*'
      elseif (.not.status) then
        write(6,*) 'Global array creation failure on ',me
      endif
#define MIRROR 1
#if MIRROR
      status = nga_create_config(type, ndim, dims,"a", chunk, p_m, g_a)
      if (status.and.me.eq.0) then
        write(6,*) '*'
        write(6,*) '* Global array A creation was successful'
        write(6,*) '*'
      elseif (.not.status) then
        write(6,*) 'Global array creation failure on ',me
      endif
      status = nga_create_config(type, ndim, dims,"b", chunk, p_m, g_b)
      if (status.and.me.eq.0) then
        write(6,*) '*'
        write(6,*) '* Global array B creation was successful'
        write(6,*) '*'
      elseif (.not.status) then
        write(6,*) 'Global array creation failure on ',me
      endif
      status = nga_create_config(type, ndim, dims,"c", chunk, p_m, g_c)
      if (status.and.me.eq.0) then
        write(6,*) '*'
        write(6,*) '* Global array C creation was successful'
        write(6,*) '*'
      elseif (.not.status) then
        write(6,*) 'Global array creation failure on ',me
      endif
#endif
c
c      initialize a and b matrices
c
      k = 0
      l = 7
      do j = 1, mdim
        do i = 1, mdim
          k =  k+1
          k = mod(k,29)
          a(i,j) = dble(k)
          l = l + 1
          l = mod(l,37)
          b(i,j) = dble(l)
        end do
      end do
c
c      initialize global arrays
c
#if MIRROR
      if (iproc.eq.0) then
        call nga_put(g_a,lo,hi,a,ld)
        call nga_put(g_b,lo,hi,b,ld)
      endif
      call ga_zero(g_c)
#endif
      if (me.eq.0) then
        call nga_put(g_aa,lo,hi,a,ld)
        call nga_put(g_bb,lo,hi,b,ld)
      endif
      call ga_zero(g_cc)
      if (me.eq.0) then
        write(6,*) '*'
        write(6,*) '* Initialized global arrays with data'
        write(6,*) '*'
      endif
      call ga_sync()
c
c      work out which chunk of matrix c this node is responsible
c      for
c
      t1 = 0.0d00
      t2 = 0.0d00
      t3 = 0.0d00
      t4 = 0.0d00
      t5 = 0.0d00
      t6 = 0.0d00
      t7 = 0.0d00
      do iloop = 1, nloop
        call ga_zero(g_cc)
#if MIRROR
        call ga_zero(g_c)
        start = util_timer()
        start1 = util_timer()
#if !USE_FAST_MERGE
        ix = mod(inode,ndims(1))
        iy = (inode - ix)/ndims(1)
        xinc = dims(1)/ndims(1)
        yinc = dims(2)/ndims(2)
        lo2(1) = ix*xinc+1
        lo2(2) = iy*yinc+1
        if (ix.ne.ndims(1)-1) then
          hi2(1) = (ix+1)*xinc
        else
          hi2(1) = dims(1)
        endif
        if (iy.ne.ndims(2)-1) then
          hi2(2) = (iy+1)*yinc
        else
          hi2(2) = dims(2)
        endif
#else
        mblock = ga_num_mirrored_seg(g_c)
        do iblock = 1, mblock
        call ga_get_mirrored_block(g_c,iblock,lo2,hi2)
#endif
c
c      work out which chunk of matrix c this processor is responsible
c      for
c
        ix = mod(iproc,pdims(1))
        iy = (iproc - ix)/(hi2(1)-lo2(1)+1)
        xinc = (hi2(1)-lo2(1)+1)/pdims(1)
        yinc = (hi2(2)-lo2(2)+1)/pdims(2)
        lo3(1) = ix*xinc+lo2(1)
        lo3(2) = iy*yinc+lo2(2)
        if (ix.ne.pdims(1)-1) then
          hi3(1) = (ix+1)*xinc+lo2(1)-1
        else
          hi3(1) = hi2(1)
        endif
        if (iy.ne.pdims(2)-1) then
          hi3(2) = (iy+1)*yinc+lo2(2)-1
        else
          hi3(2) = hi2(2)
        endif
        t3 = t3 + util_timer() - start1
c
        start1 = util_timer()
#define DGEMM2 1
#if (DGEMM2 && !USE_FAST_MERGE)
        lo4(1) = lo2(1)
        hi4(1) = hi2(1)
        lo4(2) = 1
        hi4(2) = dims(2)
        lo5(1) = 1
        hi5(1) = dims(1)
        lo5(2) = lo2(2)
        hi5(2) = hi2(2)
c
c    perform matrix patch multiplies
c
        alpha = 1.0
        beta = 0.0
        start1 = util_timer()
        call ga_matmul_patch('n','n',alpha,beta,
     +                       g_a,lo4(1),hi4(1),lo4(2),hi4(2),
     +                       g_b,lo5(1),hi5(1),lo5(2),hi5(2),
     +                       g_c,lo2(1),hi2(1),lo2(2),hi2(2))

        t5 = t5 + util_timer() - start1
#else
        lo4(1) = lo3(1)
        hi4(1) = hi3(1)
        lo4(2) = 1
        hi4(2) = dims(2)
        call nga_get(g_a,lo4,hi4,ta,ld)
        lo4(1) = 1
        hi4(1) = dims(1)
        lo4(2) = lo3(2)
        hi4(2) = hi3(2)
        call nga_get(g_b,lo4,hi4,tb,ld)
        t4 = t4 + util_timer() - start1
        start1 = util_timer()
        imax = hi3(1) - lo3(1) + 1
        do i = 1, imax
          do j = 1, mdim
            tat(j,i) = ta(i,j)
          end do
        end do
        imax = hi3(1) - lo3(1) + 1
        jmax = hi3(2) - lo3(2) + 1
        do j = 1, jmax
          do i = 1, imax
            tc(i,j) = 0.0d00
            do k = 1, mdim
              tc(i,j) = tc(i,j) + tat(k,i)*tb(k,j)
            end do
          end do
        end do
        t5 = t5 + util_timer() - start1
        start1 = util_timer()
        call nga_put(g_c,lo3,hi3,tc,ld)
        t6 = t6 + util_timer() - start1
#endif
#if USE_FAST_MERGE
        end do
#endif
        call ga_mask_sync(.true.,.false.)
        start1 = util_timer()
#if USE_FAST_MERGE
        call ga_fast_merge_mirrored(g_c) 
#else
        call ga_merge_mirrored(g_c) 
#endif
        t7 = t7 + util_timer() - start1
        t1 = t1 + util_timer() - start
#else
        do i = 1, mdim
          do j = 1, mdim
            tat(j,i) = a(i,j)
          end do
        end do
        do i = 1, mdim
          do j = 1, mdim
            c(i,j) = 0.0d00
            do k = 1, mdim
              c(i,j) = c(i,j) + tat(k,i)*b(k,j)
            end do
          end do
        end do
#endif
c
c   do a conventional ga matrix multiply
c
        alpha = 1.0
        beta = 0.0
        call ga_sync()
        start = util_timer()
#define DGEMM 1
#if DGEMM
        call ga_dgemm('n','n',mdim,mdim,mdim,alpha,g_aa,g_bb,beta,g_cc)
#else
        do i = 1, mdim
          do j = 1, mdim
            tat(j,i) = a(i,j)
          end do
        end do
        do i = 1, mdim
          do j = 1, mdim
            tc(i,j) = 0.0d00
            do k = 1, mdim
              tc(i,j) = tc(i,j) + tat(k,i)*b(k,j)
            end do
          end do
        end do
#endif
        t2 = t2 + util_timer() - start
        call ga_sync()
c
c   check to see if matrix operations agree
c
#if MIRROR
        call nga_get(g_c,lo,hi,c,ld)
#endif
#if DGEMM
        call nga_get(g_cc,lo,hi,tc,ld)
#endif
        call ga_sync()
        ltest = .true.
        afail = 0
        bfail = 0
        cfail = 0
        do i = 1, mdim
          do j = 1, mdim
            if (c(i,j).ne.tc(i,j)) then
              cfail = cfail+1
              ltest = .false.
              if (c(i,j).eq.0.0d00) then
                afail = afail + 1
              endif
              if (tc(i,j).eq.0.0d00) then
                bfail = bfail + 1
              endif
            end if
          end do
        end do
        if (ltest) then
          if (me.eq.0) write(6,*) 'Matrices agree for iteration ',iloop
        else
          if (me.eq.0) then
            write(6,*) 'Matrices disagree for iteration',iloop
            write(6,*) '  Total number of disagreements',cfail
            write(6,*) '  Total number of zeros in mirrored array',afail
            write(6,*) '  Total number of zeros in dgemm array',bfail
          endif
        endif
      end do
c
      call ga_dgop(MT_DBL,t1,1,'+')
      call ga_dgop(MT_DBL,t2,1,'+')
      call ga_dgop(MT_DBL,t3,1,'+')
      call ga_dgop(MT_DBL,t4,1,'+')
      call ga_dgop(MT_DBL,t5,1,'+')
      call ga_dgop(MT_DBL,t6,1,'+')
      call ga_dgop(MT_DBL,t7,1,'+')
      t1 = t1/dble(nloop*tprocs)
      t2 = t2/dble(nloop*tprocs)
      t3 = t3/dble(nloop*tprocs)
      t4 = t4/dble(nloop*tprocs)
      t5 = t5/dble(nloop*tprocs)
      t6 = t6/dble(nloop*tprocs)
      t7 = t7/dble(nloop*tprocs)
      if (me.eq.0) then
        write(6,*) 'Elapsed time for Mirrored matrix multiply ',t1
        write(6,*) '  Elapsed time for initial setup ',t3
        write(6,*) '  Elapsed time for copy to buffers ',t4
        write(6,*) '  Elapsed time for partial multiply ',t5
        write(6,*) '  Elapsed time for copy from buffer',t6
        write(6,*) '  Elapsed time for merge',t7
        write(6,*) 'Elapsed time for Conventional matrix multiply ',t2
      endif
c
c***  Tidy up the GA package
c
      if(ga_nodeid().eq.0)print *,' All tests successful'
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
