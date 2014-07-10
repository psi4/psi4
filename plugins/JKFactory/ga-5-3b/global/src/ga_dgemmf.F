#if HAVE_CONFIG_H
#   include "config.fh"
#endif
#define xx_dgemm dgemm
      subroutine ga_dgemm(transa, transb, m, n, k, alpha, g_a, 
     $     g_b, beta, g_c)
C$Id: ga_dgemm.F,v 1.29 2000/11/04 01:46:31 d3h325 Exp $
      implicit none
      Character*1        transa, transb
      Integer            m, n, k
      Double precision   alpha, beta
      Integer            g_a, g_b, g_c
#include "mafdecls.fh"
#include "global.fh"
c     
c     GA_DGEMM  performs one of the matrix-matrix operations:
c           C := alpha*op( A )*op( B ) + beta*C,
c     where  op( X ) is one of
c           op( X ) = X   or   op( X ) = X`,
c     
c     alpha and beta are scalars, and A, B and C are matrices, with op( A )
c     an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
c     
c     On entry, TRANSA specifies the form of op( A ) to be used in
c     the matrix multiplication as follows:
c           transa = 'N' or 'n',  op( A ) = A.
c           transa = 'T' or 't',  op( A ) = A`.
c     
c     M      - On entry,  M  specifies  the number  of rows  of the  matrix
c              op( A )  and of the  matrix  C.  M  must  be at least  zero.
c     N      - On entry,  N  specifies the number  of columns of the matrix
c              op( B ) and the number of columns of the matrix C. N must be
c              at least zero.
c     K      - On entry,  K  specifies  the number of columns of the matrix
c              op( A ) and the number of rows of the matrix op( B ). K must
c              be at least  zero.
c
      integer ilo, ihi, jlo, jhi, klo, khi, ichunk, jchunk, kchunk
      integer idim, jdim, kdim, adim, bdim, cdim, ijk, me, nproc
      integer l_a, k_a, l_b, k_b
      logical status
C
      Logical Get_New_B ! Allow reuse of B patch when possible
C
      Double Precision Chunk_cube
      Integer Min_Tasks, Max_Chunk, Mem_Avail
      integer l_mxn,k_mxn,i0,i1,j0,j1,ldc,adrc
      integer an1, an2, bn1, bn2, cn1, cn2
      integer ilor, ihir,jlor,jhir,klor,khir,itipo,ijmax
      double precision t0,t1,gflop
      external MPI_Wtime
      double precision MPI_Wtime
      Parameter ( Min_Tasks = 10) ! Minimum acceptable tasks per node
c
C     Set defaults -- platform dependent
#ifdef GATIME
      t0=MPI_Wtime()
#endif
      ichunk = 512
      jchunk = 512
      kchunk = 512
C
      me = ga_nodeid()
      nproc = ga_nnodes()
c      if(me.eq.0)
c     W     write(6,*) ' transa, transb ', transa, transb
C
C     Make an estimate of how large patches can be and still insure
C     enough tasks per processor that loads will be reasonably balanced.
C     
C     Patches per dimension are M/chunk, N/chunk, K/chunk so total tasks
C     is roughly (K*M*N)/(chunk**3).  Assume all chunk sizes are the 
C     same and solve for the one that provides the minimum acceptable 
C     number of tasks.
C
C        Find out how much memory we can grab.  It will be used in
C        three chunks, and the result includes only the first one.
C
         Mem_Avail = MA_Inquire_Avail( MT_DBL )
     $      - 2 * MA_SizeOf_Overhead( MT_DBL )
	 Mem_Avail = 0.9 * Mem_Avail ! Do not use every last drop!
c         Call GA_IGOp(42, Mem_Avail, 1, 'min')
C

c
      if (beta .eq. 0.0d0) then
         call ga_zero(g_c)
      else
         call ga_scale(g_c, beta)
      endif
c
      call ga_distribution(g_c,
     .     ga_nodeid(), i0, i1, j0, j1)
      call ga_inquire(g_a,
     .     itipo, an1, an2)
      call ga_inquire(g_b,
     .     itipo, bn1, bn2)
      call ga_inquire(g_c,
     .     itipo, cn1, cn2)
      if (i0.gt.0 .and. i0.le.i1) then
         ilo=i0
         ihi=i1
         idim = ihi - ilo + 1
         jlo=j0
         jhi=j1
         jdim = jhi - jlo + 1
#if 0
         write(6,'(I4,A,4I6)') ga_nodeid(),' IJ ',i0,i1,j0,j1
      if(ga_nodeid().eq.0) call ffflush(6)
      if(ga_nodeid().eq.0) call ffflush(0)
#endif
         ijmax=max(idim,jdim)
         KChunk =  Int((DBLE(Mem_Avail/(2*ijmax)))) 
         kchunk=min(kchunk,ijmax)
      status = .true.
      status = ma_push_get(MT_DBL, idim*kchunk, 'ga_dgemm:a', l_a,k_a)
     $     .and. status
      status = ma_push_get(MT_DBL, kchunk*jdim, 'ga_dgemm:b', l_b,k_b)
     $     .and. status
      if (.not. status) call ga_error('ga_dgemm: insufficent memory?',
     A idim*kchunk+kchunk*jdim)
         call ga_access(g_c, i0, i1, j0, j1, adrc, ldc)
      ijk = 0
         do klo = 1, k, kchunk
            khi = min(k, klo+kchunk-1)
            kdim = khi - klo + 1
C
C           Each pass through the outer two loops means we need a
C           different patch of B.
C
            Get_New_B = .TRUE.
C
                  cdim = idim
                  if (transa.eq.'n' .or. transa.eq.'N') then
                     ilor=min(an1,ilo)
                     ihir=min(an1,ihi)
                     klor=min(an2,klo)
                     khir=min(an2,khi)
                     kdim=khir-klor+1
                     idim=ihir-ilor+1
                     adim = idim
                     cdim = idim
                     call ga_get(g_a, ilor, ihir, klor, khir,
     $                  dbl_mb(k_a), adim)
                  else
                     klor=min(an1,klo)
                     khir=min(an1,khi)
                     ilor=min(an2,ilo)
                     ihir=min(an2,ihi)
                     kdim=khir-klor+1
                     idim=ihir-ilor+1
                     adim = kdim
                     cdim=idim
                     call ga_get(g_a, klor, khir, ilor, ihir, 
     $                  dbl_mb(k_a), adim)
                  endif
C
C                 Avoid rereading B if it is the same patch as last time.
C
                  If ( Get_New_B ) then
                     if (transb.eq.'n' .or. transb.eq.'N') then
                        klor=min(bn1,klo)
                        khir=min(bn1,khi)
                        jlor=min(bn2,jlo)
                        jhir=min(bn2,jhi)
                        kdim=khir-klor+1
                        idim=ihir-ilor+1
                        bdim = kdim
                        call ga_get(g_b, klor, khir, jlor, jhir, 
     $                     dbl_mb(k_b), bdim)
                     else
                        jlor=min(bn1,jlo)
                        jhir=min(bn1,jhi)
                        klor=min(bn2,klo)
                        khir=min(bn2,khi)
                        kdim=khir-klor+1
                        jdim=jhir-jlor+1
                        bdim = jdim
                        call ga_get(g_b, jlor, jhir, klor, khir, 
     $                     dbl_mb(k_b), bdim)
                     endif
                     Get_New_B = .FALSE. ! Until J or K change again
                  EndIf
C
                  call xx_dgemm(transa, transb, idim, jdim, kdim,
     $                 alpha, dbl_mb(k_a), adim, dbl_mb(k_b), bdim,
     $                 1.0d0, dbl_mb(adrc), cdim)
         enddo
      status = ma_chop_stack(l_a)
      if (.not. status)call ga_error('ga_dgemm: pop of stack failed', 1)
      endif
      call ga_release_update(g_c, i0, i1, j0, j1)
      call ga_sync()
#ifdef GATIME
      if(ga_nodeid().eq.0) then
        call ffflush(6)
          t1=MPI_Wtime()-t0
          gflop=2d0*n*m*k/(t1*1d9)
      write(6,'(I4,A,3F14.6)') 
     G     ga_nodeid(),' dgemm done in ',t1,
     ,     gflop,gflop/ga_nnodes()
      endif
#endif
c
      end

      subroutine lga_acc(out, dim1,dim2,
     I     ilo, ihi, jlo, jhi, buf, 
     $     ld, alpha)
      implicit none
      integer dim1,dim2
      integer ilo, ihi, jlo, jhi,ld
      integer i,j
      double precision alpha
      double precision out(1:dim1,1:dim2)
      double precision buf(1:ld, 1:*)

      do j=jlo,jhi
         do i=ilo,ihi
            out(i, j) = out(i, j) + 
     +           alpha*buf(i-ilo+1, j-jlo+1)
         enddo
      enddo
      return
      end
      subroutine lga_accbrd(g_c,m,n,out)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer m,n
      double precision out(1:m,1:n)
      integer g_c
c
      integer ilo,ihi,jlo,jhi,numi,numj,k_in,l_in,i,j
      logical status
c
      call ga_distribution(g_c,
     .     ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. ilo.le.ihi) then
         numi =  ihi-ilo+1
         numj =  jhi-jlo+1
         if (numi.gt.0 .and. numj.gt.0) then
            if (.not.MA_Push_Get(MT_Dbl,numi*numj,'MxN',l_in,k_in))
     &           call ga_error('dft_scf: cannot allocate eval',0)
            call ga_get(g_c,ilo,ihi,jlo,jhi,
     .           dbl_mb(k_in),numi)
            do j=jlo,jhi
               do i=ilo,ihi
                  call daxpy(numi,1d0,dbl_mb(k_in+(j-jlo)*numi),1,
     1                 out,1)
               enddo
            enddo
            call ga_put(g_c,ilo,ihi,jlo,jhi,
     .           out(ilo,jlo),m)
            status = ma_pop_stack(l_in)
            if (.not. status)call ga_error('gad: pop of stack failed', 
     3           numi*numj)
         endif
      endif
      return
      end
