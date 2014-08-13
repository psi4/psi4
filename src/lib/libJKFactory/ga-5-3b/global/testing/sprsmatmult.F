#if HAVE_CONFIG_H
#   include "config.fh"
#endif
#define USE_PETSC 0
#define TEST_MULT 0
#if TEST_MULT
c#   define LMAX 972
#else
c#   define LMAX 94173031
c#   define LMAX 9173031
#   define LMAX 973031
c#   define LMAX 97031
c#   define LMAX 9031
#endif
#define MAXVEC 50000000

#define bb_a(ib) bb_v(bb_i + (ib))
#define cc_a(ib) cc_v(cc_i + (ib))

      program main
      implicit none
#include "mafdecls.fh"
#include "global.fh"

#if USE_PETSC
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscsys.h"
#endif

      integer heap, stack, fudge, ma_heap
      integer nmax, ndim, nprocs, me
      integer g_a_data, g_a_i, g_a_j, g_a_sbit, g_tmp, isize
      integer g_b, g_c, g_tmp2
      integer i, j, k, one, icnt
      integer lo, hi, chunk, kp1, ld
      GA_ACCESS_INDEX_TYPE i8, idx, id_tmp
      integer g_chk(MAXVEC)
      double precision d_chk(MAXVEC)
      double precision t_beg, t_crt, t_mv, t_gat, t_elem, t_pack
      double precision t_shft, t_scnad
#if TEST_MULT
      double precision ddata(MAXVEC),bb(LMAX),cc(LMAX),gc(LMAX)
      integer ii(LMAX+1),jj(MAXVEC),sbit(MAXVEC)
      double precision rdot
#endif
      logical status
      integer idim, jdim, rmin, rmax, idum
      integer jmin, jmax, lsize
      parameter (heap=2000000, fudge=100, stack=2000000)
      integer i_chk(MAXVEC)
      double precision ran1
c
c   PETSc declarations
c
#if USE_PETSC
      PetscErrorCode p_ierr
      PetscReal bb_v(1),cc_v(1),zero,p_ij
      PetscOffset bb_i, cc_i
      Vec p_b, p_c
      PetscInt p_idim, p_nloc
      PetscReal, pointer :: val(:)
      PetscInt, pointer :: loc(:), iloc(:)
      PetscInt, pointer :: jloc(:)
      PetscInt p_i, p_j
      Mat p_m
      double precision t_petsc, t_assmbl, prdot
      double precision, pointer :: bbuf(:),cbuf(:)
#endif
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
      call ga_initialize()
      t_crt = 0.0d00
      t_mv = 0.0d00
      t_gat = 0.0d00
      t_elem = 0.0d00
      t_pack = 0.0d00
      t_shft = 0.0d00
      t_scnad = 0.0d00
#if USE_PETSC
      t_petsc = 0.0d00
      t_assmbl = 0.0d00
      call PetscInitialize(PETSC_NULL_CHARACTER,p_ierr)
#endif
      me = ga_nodeid()
      nprocs = ga_nnodes()
      one = 1
c     we can also use GA_set_memory_limit BEFORE first ga_create call
c
      ma_heap = heap + fudge 
c      call GA_set_memory_limit(util_mdtob(ma_heap))
c
      if(ga_nodeid().eq.0)then
        print *
        print *,'GA initialized '
        print *
#if TEST_MULT
        print *
        print *,'Comparing distributed multiply with serial multiply'
        print *
#endif
        call flush(6)
      endif
c
c***  Initialize the MA package
c     MA must be initialized before any global array is allocated
c
      status = ma_init(MT_DBL, stack, ma_heap)
      if (.not. status) call ga_error('ma_init failed',-1) 
c
c    create a sparse LMAX x LMAX matrix and two vectors of length
c    LMAX. The matrix is stored in compressed row format.
c    One of the vectors is filled with random data and the other
c    is filled with zeros.
c
      idim = LMAX
      jdim = LMAX
      if (me.eq.0) then
        print *
        print *, 'Nominal dimension of matrix: ',idim
        print *
      endif
      rmin = 5
      rmax = 10
      rmin = 3
      rmax = 7
      t_beg = ga_wtime()
      call create_rand_sprs_mat(idim,jdim,rmin,rmax,
     +                          g_a_data,g_a_i,g_a_j,g_a_sbit,isize)
      t_crt = ga_wtime() - t_beg
      if (me.eq.0) then
        print *
        print *, 'Number of non-zero elements in compressed matrix: ',
     +           isize
        print *
      endif
#if TEST_MULT
      if (me.eq.0) then
        call nga_get(g_a_data,one,isize,ddata,one)
        call nga_get(g_a_i,one,idim+1,ii,one)
        call nga_get(g_a_j,one,isize,jj,one)
      endif
      call ga_sync
#endif
      g_tmp = ga_create_handle()
      call ga_set_data(g_tmp,one,isize,MT_DBL)
      status = ga_allocate(g_tmp)
c
      g_tmp2 = ga_create_handle()
      call ga_set_data(g_tmp2,one,isize,MT_DBL)
      status = ga_allocate(g_tmp2)
c
      g_b = ga_create_handle()
      call ga_set_data(g_b,one,idim,MT_DBL)
      status = ga_allocate(g_b)
c
c   fill g_b with random values
c
      call nga_distribution(g_b,me,lo,hi)
      call nga_access(g_b,lo,hi,idx,ld)
      ld = hi-lo
      do i8 = 0, ld
        dbl_mb(i8+idx) = ran1(idum)
      end do
      call nga_release(g_b,lo,hi)
#if TEST_MULT
      call ga_sync
      if (me.eq.0) then
        call nga_get(g_b,one,idim,bb,one)
      endif
      call ga_sync
#endif
c
      g_c = ga_create_handle()
      call ga_set_data(g_c,one,idim,MT_DBL)
      status = ga_allocate(g_c)
      call ga_zero(g_c)
c
#if USE_PETSC
c
c   Create corresponding PETSC vectors
c
      p_idim = idim
      zero = 0.0d00
      p_nloc = ld+1
      call VecCreateMPI(PETSC_COMM_WORLD,p_nloc,PETSC_DECIDE,p_b,p_ierr)
      call VecDuplicate(p_b,p_c,p_ierr)
c
c   Copy values from global array to PETSC vector
c
      allocate(val(ld+1))
      allocate(loc(ld+1))
      call nga_access(g_b,lo,hi,idx,ld)
      ld = hi-lo
      do i8 = 0, ld
        i = i8+1
        val(i) = dbl_mb(i8+idx)
        loc(i) = lo+i8-1
        if (loc(i).lt.0.or.loc(i).ge.idim) then
          write(6,'(i4,a,i8,a,i12)') me,' loc(',i,'): ',loc(i)
        endif
      end do 
      call nga_release(g_b,lo,hi)
      call VecSetValues(p_b,p_nloc,loc,val,INSERT_VALUES,p_ierr)
      t_beg = ga_wtime()
      call VecAssemblyBegin(p_b,p_ierr)
      call VecAssemblyEnd(p_b,p_ierr)
      t_assmbl = ga_wtime() - t_beg
      call VecSet(p_c,zero,p_ierr)
c
cdbg      call VecGetArray(p_b,bb_v,bb_i,p_ierr)
cdbg      do i = 1, hi-lo+1
cdbg        write(6,'(i4,a,i4,a,2f12.6)') me,' b(',i,'): ',val(i),
cdbg     +                                bb_a(i)
cdbg      end do
      deallocate(val)
      deallocate(loc)
c
c   Create corresponding PETSC sparse matrix
c
      call MatCreate(PETSC_COMM_WORLD,p_m,p_ierr)
      call MatSetSizes(p_m,PETSC_DECIDE,PETSC_DECIDE,
     +                 p_idim,p_idim,p_ierr)
c      call MatSetType(p_m,MATMPIAIJ,p_ierr)
      call MatSetType(p_m,MATAIJ,p_ierr)
c      call MatSetFromOptions(p_m,p_ierr)
      ld = hi-lo+1
      allocate(iloc(ld+1))
      call nga_get(g_a_i,lo,hi,iloc,ld)
c      write(6,'(i4,a,i12)') me,' ld: ',ld
c      do i = 1, ld
c         if (iloc(i).lt.1) then
c           write(6,'(i4,a,i6,a,i12)') me,' iloc(',i,'): ',iloc(i)
c         endif
c      end do
      jmin = iloc(1)
      if (hi.lt.idim) then
        call nga_get(g_a_i,hi+1,hi+1,jmax,one)
        jmax = jmax - 1
      else
        jmax = isize
      endif
      iloc(ld+1) = jmax + 1
      lsize = jmax - jmin + 1
c      write(6,'(i4,a,i12)') me,' lsize: ',lsize
      allocate(jloc(lsize))
      allocate(val(lsize))
      call nga_get(g_a_j,jmin,jmax,jloc,one)
c      do i = 1, lsize
c         if (jloc(i).lt.1) then
c           write(6,'(i4,a,i6,a,i12)') me,' jloc(',i,'): ',jloc(i)
c         endif
c      end do
      call nga_get(g_a_data,jmin,jmax,val,one)
#if 1
      do i = 1, ld
        p_i = i+lo-2
        if (p_i.lt.0.or.p_i.ge.idim) then
          write(6,'(i4,a,i8)') me,' p_i out of bounds: ',p_i
        endif
        do j = iloc(i), iloc(i+1) - 1
          p_j = jloc(j-jmin+1) - 1
          if (p_j.lt.0.or.p_j.ge.idim) then
            write(6,'(i4,a,i8)') me,' p_j out of bounds: ',p_j
          endif
          p_ij = val(j-jmin+1)
cbjp          write(6,'(i4,a,2i8,f12.6)') me,' m(i,j): ',p_i, p_j,p_ij
          call MatSetValues(p_m,one,p_i,one,p_j,p_ij,INSERT_VALUES,
     +                      p_ierr)
        end do
      end do
      t_beg = ga_wtime()
      call MatAssemblyBegin(p_m,MAT_FINAL_ASSEMBLY,p_ierr)
      call MatAssemblyEnd(p_m,MAT_FINAL_ASSEMBLY,p_ierr)
      t_assmbl = ga_wtime() - t_beg
c      call MatView(p_m,PETSC_VIEWER_DRAW_WORLD,p_ierr)
#endif
#endif
c
c   Find out what section of j-index vector I own and use
c   this information to get a pointer to to the corresponding
c   parts of the j-index vector and g_tmp array.
c
      t_beg = ga_wtime()
      call nga_distribution(g_a_j,me,lo,hi)
      call nga_access(g_a_j,lo,hi,idx,ld)
      call nga_access(g_tmp,lo,hi,id_tmp,ld)
c
c   scatter elements of g_b into a temporary array that is the
c   same size as g_a_data.
c
      ld = hi-lo+1
      call nga_gather(g_b,dbl_mb(id_tmp),int_mb(idx),ld)
      call ga_sync
      call nga_release(g_a_j,lo,hi)
      call nga_release(g_tmp,lo,hi)
      t_gat = ga_wtime() - t_beg
c
      t_beg = ga_wtime()
      call ga_elem_multiply(g_tmp,g_a_data,g_tmp2)
      t_elem = ga_wtime() - t_beg
c      if (me.eq.0) then
c        do i = 0, 1
c          call nga_distribution(g_a_sbit,i,lo,hi)
c          call nga_get(g_a_sbit,lo,hi,i_chk,one)
c          write(6,'(25i5)') (i_chk(j),j=1,hi-lo+1)
c        end do
c      endif
      t_beg = ga_wtime()
      call ga_scan_add(g_tmp2,g_tmp,g_a_sbit,one,isize,0)
      t_scnad = ga_wtime() - t_beg
c      if (me.eq.0) then
c        do i = 0, 1
c          call nga_distribution(g_tmp,i,lo,hi)
c          call nga_get(g_tmp,lo,hi,d_chk,one)
c          write(6,'(10f8.2)') (d_chk(j),j=1,hi-lo+1)
c        end do
c      endif
#if TEST_MULT
      if (me.eq.0) then
        call nga_get(g_a_sbit,one,isize,sbit,one)
        icnt = 0
        if (isize.gt.MAXVEC) call ga_error('MAXVEC exceeded in sbit',
     +                                     isize) 
        do i = 1, isize
          if (sbit(i).eq.1) then
            icnt = icnt + 1
          endif
        end do
      endif
#endif
      t_beg = ga_wtime()
      call shift_minus(g_a_sbit,isize)
      t_shft = ga_wtime() - t_beg
c
      t_beg = ga_wtime()
      call ga_pack(g_tmp,g_c,g_a_sbit,one,isize,icnt)
      t_pack = ga_wtime() - t_beg
c
#if USE_PETSC
      t_beg = ga_wtime()
      call MatMult(p_m,p_b,p_c,p_ierr)
      t_petsc = ga_wtime() - t_beg
      call nga_distribution(g_c,me,lo,hi)
      allocate(cbuf(hi-lo+1))
      call nga_get(g_c,lo,hi,cbuf,one)
      call VecGetArray(p_c,cc_v,cc_i,p_ierr)
      prdot = 0.0d00
      do i = 1, hi-lo+1
cbjp        write(6,'(i4,a,i4,a,2f12.6)') me,' c(',i,'): ',cbuf(i),
cbjp     +                                cc_a(i)
        prdot = prdot + (cc_a(i)-cbuf(i))**2
      end do
      call ga_dgop(1,prdot,1,'+')
      if (me.eq.0) then
        write(6,'(a,e16.8)')
     +          'Difference between GA and PETSC results: ',
     +           sqrt(prdot)
      endif
#endif
#if TEST_MULT
      if (me.eq.0) then
        call nga_get(g_c,one,idim,gc,one)
        do i = 1, idim
          cc(i) = 0.0d00
          do j = ii(i), ii(i+1)-1
cbjp            write(6,'(i4,a,2i8,f12.6)') me,' gm(i,j): ',
cbjp     +                          i-1,jj(j)-1,ddata(j)
            cc(i) = cc(i) + ddata(j)*bb(jj(j))
          end do
          gc(i) = gc(i) - cc(i)
        end do
        rdot = 0.0d00
        do i = 1, idim
          rdot = rdot + gc(i)**2
        end do
        write(6,*) 'Difference between parallel and serial results: ',
     +             rdot
      endif
#endif
c
c   Clean up arrays
c
      status = ga_destroy(g_tmp)
      status = ga_destroy(g_tmp2)
      status = ga_destroy(g_b)
      status = ga_destroy(g_c)
      status = ga_destroy(g_a_data)
      status = ga_destroy(g_a_i)
      status = ga_destroy(g_a_j)
      status = ga_destroy(g_a_sbit)
#if USE_PETSC
      call VecDestroy(p_b,p_ierr)
      call VecDestroy(p_c,p_ierr)
      call MatDestroy(p_m,p_ierr)
#endif
c
      t_mv = t_gat + t_scnad + t_elem + t_shft + t_pack
      call ga_dgop(1,t_crt,1,'+')
      call ga_dgop(2,t_gat,1,'+')
      call ga_dgop(3,t_scnad,1,'+')
      call ga_dgop(4,t_elem,1,'+')
      call ga_dgop(5,t_shft,1,'+')
      call ga_dgop(6,t_pack,1,'+')
      call ga_dgop(7,t_mv,1,'+')
#if USE_PETSC
      call ga_dgop(8,t_petsc,1,'+')
      call ga_dgop(9,t_assmbl,1,'+')
#endif
c
      if (me.eq.0) then
        write(6,100) t_crt/dble(nprocs)
        write(6,200) t_gat/dble(nprocs)
        write(6,300) t_scnad/dble(nprocs)
        write(6,400) t_elem/dble(nprocs)
        write(6,500) t_shft/dble(nprocs)
        write(6,600) t_pack/dble(nprocs)
        write(6,700) t_mv/dble(nprocs)
#if USE_PETSC
        write(6,800) t_petsc/dble(nprocs)
        write(6,900) t_assmbl/dble(nprocs)
#endif
  100   format('Time to create sparse matrix:                  ',f12.4)
  200   format('Time to expand dense vector using gather:      ',f12.4)
  300   format('Time in segmented scan-add:                    ',f12.4)
  400   format('Time for element-wise multiplication:          ',f12.4)
  500   format('Time to shift bit-mask vector:                 ',f12.4)
  600   format('Time to pack result vector:                    ',f12.4)
  700   format('Total time for sparse matrix-vector multiply:  ',f12.4)
#if USE_PETSC
  800   format('Total time for PETSC matrix-vector multiply:   ',f12.4)
  900   format('Total time for PETSC assemble:                 ',f12.4)
#endif
      endif
c
#if USE_PETSC
      call PetscFinalize(p_ierr)
#endif
      call ga_terminate()
c
c***  Tidy up after message-passing library
c
      call MP_FINALIZE()

c
      stop
      end
c
c   create a random sparse matrix in compressed row form
c
      subroutine create_rand_sprs_mat(idim,jdim,rmin,rmax,
     +                                g_data,g_i,g_j,g_mask,isize)
c
c   idim: i-dimension of array
c   jdim: j-dimension of array
c   rmin, rmax: minimum and maximum number of row elements in array
c   g_data: global array of values
c   g_j: global array containing j indices
c   g_i: global array containing starting location of each row in g_j
c   g_mask: global array of 1s and 0s containing a 1 at the start of
c           every new row in g_data
c   isize: total size of compressed matrix
c
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer idim, jdim, rmin, rmax, g_data, g_i, g_j, g_tmp
      integer g_sbit, g_mask, isize
      integer nprocs, me, imin, imax, icnt, jcnt, jmin, jmax, ntot
      double precision rval(MAXVEC), x, dr
      integer jval(MAXVEC), inc(MAXVEC), ones(MAXVEC)
      integer jdx(MAXVEC), i, j, k, itmp, itot, one, lo, hi, ld
      integer idum
      logical status
      integer i_chk(1000000)
      double precision d_chk(1000000)
      double precision ran1
      GA_ACCESS_INDEX_TYPE idx
c
      me = ga_nodeid()
      nprocs = ga_nnodes()
      x = ran1(12345+me)
      dr = dble(rmax-rmin)
      one = 1
c
c  figure out min and max indices of rows that I am responsible for
c
      imin = nint(dble(me)*(dble(idim)/dble(nprocs))) + 1
      imax = nint(dble(me+1)*(dble(idim)/dble(nprocs)))
      if (me.eq.0) then
        write(6,*) 
        write(6,*) 'Minimum number of entries per row: ',
     +                     rmin
        write(6,*) 
        write(6,*) 
        write(6,*) 'Maximum number of entries per row: ',
     +                     rmin+nint(dr)
        write(6,*) 
      endif
      icnt = 0
      do i = imin, imax
        jmax = rmin + nint(dr*ran1(idum)) 
c        do j = 1, jmax
c          jdx(j) = int(dble(jdim)*ran1(idum))+1
c          if (jdx(j).gt.jdim) jdx(j) = jdim
c        end do
        j = 0
        do while (j.lt.jmax)
          j = j + 1
          jdx(j) = int(dble(jdim)*ran1(idum))+1
          if (jdx(j).gt.jdim) jdx(j) = jdim
          k = 0
          do while(k.lt.j-1)
            k = k + 1
            if (jdx(j).eq.jdx(k)) then
              k = j
              j = j - 1
            endif
          end do
        end do
c
c  sort j indices
c
        do j = 1, jmax
          do k = j+1, jmax
            if (jdx(j).gt.jdx(k)) then
              itmp = jdx(j)
              jdx(j) = jdx(k)
              jdx(k) = itmp
            endif
          end do
        end do
        do j = 1, jmax
          if (jdx(j).lt.1.or.jdx(j).gt.jdim) then
            write(6,'(i4,a,2i8)') me,' invalid j: ',i,jdx(j)
          endif
        end do
c
c  create array elements
c
        inc(i-imin+1) = jmax
        do j = 1, jmax
          icnt = icnt + 1
          rval(icnt) = ran1(idum)
          jval(icnt) = jdx(j)
        end do
      end do
      itot = icnt
      isize = 0
      do i = 1, imax-imin+1
        isize = isize + inc(i)
      end do
      if (isize.gt.MAXVEC)
     +   call ga_error('ISIZE exceeds MAXVEC in local arrays ',isize)
c
c  create global arrays to hold sparse matrix
c
      call ga_igop(1,itot,1,'+')
      isize = itot
c      write(6,*) me,' (create_rand) isize: ',isize
c
      g_data = ga_create_handle();
      call ga_set_data(g_data,one,itot,MT_DBL)
      status = ga_allocate(g_data)
c
      g_j = ga_create_handle();
      call ga_set_data(g_j,one,itot,MT_INT)
      status = ga_allocate(g_j)
c
      g_i = ga_create_handle();
      call ga_set_data(g_i,one,idim+1,MT_INT)
      status = ga_allocate(g_i)
      status = ga_duplicate(g_i, g_tmp, 'temp array')
c
      g_sbit = ga_create_handle();
      call ga_set_data(g_sbit,one,idim+1,MT_INT)
      status = ga_allocate(g_sbit)
c
      call ga_zero(g_sbit)
      if (me.eq.0) then
        call nga_put(g_sbit,one,one,one,one)
      endif
      call ga_sync
c
c  create offset vector
c
      lo = imin+1
      hi = imax+1
      call ga_zero(g_tmp)
      if (me.eq.0) then
        call nga_put(g_tmp,one,one,one,one)
      endif
      call nga_put(g_tmp,lo,hi,inc,one)
      call ga_sync
c      if (me.eq.0) then
c        do i = 0, nprocs-1
c          call nga_distribution(g_tmp,i,lo,hi)
c          call nga_get(g_tmp,lo,hi,i_chk,one)
c          write(6,'(25i5)') (i_chk(j),j=1,hi-lo+1)
c        end do
c      endif
      itot = idim+1
      call ga_scan_add(g_tmp,g_i,g_sbit,one,itot,0)
      call ga_sync
c      if (me.eq.0) then
c        do i = 0, nprocs-1
c          call nga_distribution(g_i,i,lo,hi)
c          call nga_get(g_i,lo,hi,i_chk,one)
c          write(6,'(25i5)') (i_chk(j),j=1,hi-lo+1)
c        end do
c      endif

      status = ga_destroy(g_tmp)
      status = ga_destroy(g_sbit)
c
c  create mask vector
c
      g_mask = ga_create_handle();
      call ga_set_data(g_mask,one,isize,MT_INT)
      status = ga_allocate(g_mask)
      call ga_zero(g_mask)
c
      call nga_distribution(g_i,me,lo,hi)
      call nga_access(g_i,lo,hi,idx,ld)
      ntot = hi - lo + 1
      if (ntot.gt.MAXVEC)
     +   call ga_error('NTOT exceeds MAXVEC in bit vector ',ntot)
c
c  decrease count by one on last processor because this
c  element doesn't exist
c
      if (me.eq.nprocs-1) then
        ntot = ntot - 1
      endif
      do i = 1, ntot
        ones(i) = 1
      end do
      call nga_scatter(g_mask,ones,int_mb(idx),ntot)
      call ga_sync
      call nga_release(g_i,lo,hi)
c      if (me.eq.0) then
c        do i = 0, 1
c          call nga_distribution(g_mask,i,lo,hi)
c          call nga_get(g_mask,lo,hi,i_chk,one)
c          write(6,'(25i5)') (i_chk(j),j=1,hi-lo+1)
c        end do
c      endif
c
c store values of j in g_j and data values in g_data
c
      call nga_get(g_i,imin,imin,jmin,one)
      call nga_get(g_i,imax+1,imax+1,jmax,one)
      jmax = jmax - 1
      call nga_put(g_j,jmin,jmax,jval,one)
      call nga_put(g_data,jmin,jmax,rval,one)
      call ga_sync
      return
      end
c
c   shift all the elements in a linear global array one element to
c   the left. Wrap elements around end if necessary
c
      subroutine shift_minus(g_a, isize)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer g_a, isize
      integer lo, hi, ld, me, itmp
      integer isav
      GA_ACCESS_INDEX_TYPE idx, ii
      me = ga_nodeid()
      call ga_sync
      call nga_distribution(g_a,me,lo,hi)
      call nga_access(g_a,lo,hi,idx,ld)
      ld = hi-lo
      isav = int_mb(idx)
      do ii = 1, ld
        int_mb(idx + ii - 1) = int_mb(idx + ii)
      end do
      if (lo.eq.1) then
        itmp = isize
      else
        itmp = lo-1
      endif
      call nga_release(g_a,lo,hi)
      call ga_sync
      call nga_put(g_a,itmp,itmp,isav,ld)
      call ga_sync
      return
      end
c
c   Wrapper for random number generator
c
      double precision function ran1(idum)
      implicit none
#include "testutil.fh"
      integer idum,iff
      data iff /1/
      save iff
      if (iff.eq.1) then
        ran1 = drand(iabs(idum))
        iff = 0
      endif
      ran1 = drand(0)
      return
      end
