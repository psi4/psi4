#if HAVE_CONFIG_H
#   include "config.fh"
#endif
#ifdef SCALAPACK_I8
#define INTGR4 integer*8
#else
#define INTGR4 integer*4
#endif
cedo#define HAVE_ELPA 1
C     $Id: scalapack.F,v 1.11 2006/10/13 16:15:17 d3g293 Exp $
****************************************************************************
* GA to/from ScaLAPACK (square block scattered decomposition) interface    * 
*                                                                          *
* common variables:                                                        *
*     nnodes        - number of processors                                 *
*     iam           - my processor number                                  *
*     nprow/ npcol  - number of processor rows/cols in virtual proc. grid  *
*     myrow/mycol   - cordinates of my processor in virtual proc. grid     *
*                                                                          *
c***************************************************************************
c* 04/12/96  GVT  Changed the code to adapt to a new version of ScaLAPACK
c*           Giuseppe Vitillaro peppe@unipg.it
c***************************************************************************
      block data SLblock
      implicit none
c
#include "scalapack.fh"
      data init /.false./
      data init2 /.false./
      data init3 /.false./
      end


      subroutine SLinit
      implicit none
      INTGR4  icontxt
c
#include "scalapack.fh"
#include "global.fh"
      INTGR4 zero4
      parameter(zero4=0)
c
      if(init)return
c
c**** call ga_sync before to enter in BLACS and after
c     
      call blacs_pinfo(iam, nnodes)
c
      call FindGrid(nnodes, nprow, npcol)
c
      call blacs_get (zero4, zero4, icontxt)
      call blacs_gridinit(icontxt, 'R', nprow, npcol)
      iSLctxt = icontxt
c
      call blacs_gridinfo(iSLctxt, nprow, npcol, myrow, mycol)
      
c
      init=.true.
c
      call ga_sync()
c    
      end
c
      subroutine slinit2(n)
      implicit none
#include "scalapack.fh"
#include "global.fh"
      INTGR4 n
c
      INTGR4 zero4
      parameter(zero4=0)
      INTGR4 slgetmxproc
      external slgetmxproc
c
      if(init2)return
c
c**** call ga_sync before to enter in BLACS and after
      call ga_sync()
c     
      call blacs_pinfo(iam, nnodes)
c
c     determine optimal nprocs for eigensolvers based on matrix size n
c
      maxproc=slgetmxproc(n,nnodes)
      call FindGrid(maxproc, nprow2, npcol2)
c
      call blacs_get( zero4, zero4, islctxt2 )
      call blacs_gridinit(islctxt2, 'R', nprow2, npcol2)


c
      if(iam.lt.maxproc) then
         call blacs_gridinfo(iSLctxt2, nprow2, npcol2, myrow2, mycol2)
      else
         nprow2=0
         npcol2=0
         myrow2=0
         mycol2=0
      endif
      init2=.true.
c
      call ga_sync()
      
c    
      end
c
      subroutine SLinit3(g_a)
      implicit none
      INTGR4  icontxt
c
#include "scalapack.fh"
#include "global.fh"
      integer g_a, nblocks(2), ndim, dims(2), type
      INTGR4 zero4
      parameter(zero4=0)
c
      if(init3)return
c
c   find dimensions of processor grid for global array with handle g_a
c     
      call nga_inquire(g_a, type, ndim, dims)
      if (ndim.ne.2) then
        call ga_error('SCALAPAC array dimension is not 2',ndim)
      endif
      if (ga_uses_proc_grid(g_a)) then
        call ga_get_proc_grid(g_a, nblocks)
        nprow = nblocks(1)
        npcol = nblocks(2)
      else
        call ga_error('Array does not use processor grid SCALAPACK',0)
      endif

c
      call sl_init(icontxt, nprow, npcol)
      iSLctxt = icontxt
c
      call blacs_gridinfo(iSLctxt, nprow, npcol, myrow, mycol)
c
      init3=.true.
c
      call ga_sync()
c    
      end
c
      subroutine SLinit4(g_a)
      implicit none
#include "scalapack.fh"
#include "global.fh"
c
      INTGR4 icontxt, n
      INTGR4 slgetmxproc
      external slgetmxproc
      integer g_a, nblocks(2), ndim, dims(2), type
c
      if(init4)return
c
c**** call ga_sync before to enter in BLACS and after
      call nga_inquire(g_a, type, ndim, dims)
      if (ndim.ne.2) then
        call ga_error('SCALAPAC array dimension is not 2',ndim)
      endif
      if (ga_uses_proc_grid(g_a)) then
        call ga_get_proc_grid(g_a, nblocks)
        nprow = nblocks(1)
        npcol = nblocks(2)
      else
        call ga_error('Array does not use processor grid SCALAPACK',0)
      endif
c
      call blacs_pinfo(iam, nnodes)
      maxproc=slgetmxproc(n,nnodes)
      call sl_init(icontxt, nprow, npcol)
      islctxt2 = icontxt
c
      call blacs_gridinfo(iSLctxt2, nprow2, npcol2, myrow2, mycol2)
      
      init4 = .true.
c
      call ga_sync()
      
c    
      end

      INTGR4 function slgetmxproc(n,nnodes)
      implicit none
      INTGR4 nnodes
      INTGR4 n
      INTGR4 i
      double precision fact
      INTGR4 nmax,nprocs,twoi
      double precision nprocs0
      double precision otto
      parameter(nmax=11,fact=((7108d0*7108d0)/1024d0),otto=8d0)
cnew      parameter(nmax=11,fact=((7108d0*7108d0)/512d0),otto=8d0)
c     lower bound of 8 procs
      nprocs0=max((n*n)/fact,otto)
c
c     try to get powers of two
c
      do i = nmax, 0, -1 
         if(nint(nprocs0/(2d0**i)).eq.1) goto 1
      enddo
      i=4
1     twoi=2**i
      slgetmxproc=min(nnodes,twoi)
      return
      end
      subroutine slexit
      implicit none
c
#include "scalapack.fh"
#include "global.fh"
c
      if(.not.init)return
      call ga_sync()
      call blacs_gridexit(islctxt)
      init=.false.
      return
      end
c
      subroutine slexit2
      implicit none
      INTGR4 nprocs
c
#include "scalapack.fh"
#include "global.fh"
c
      if(.not.init2)return
      call ga_sync()
      if(iam.lt.maxproc)call blacs_gridexit(islctxt2)
      init2=.false.
      return
c    
      end
c
      subroutine slexit3
      implicit none
      INTGR4 nprocs
c
#include "scalapack.fh"
#include "global.fh"
c
      if(.not.init3)return
      call ga_sync()
      call blacs_gridexit(islctxt)
      init3=.false.
      return
c    
      end
c    
      subroutine slexit4
      implicit none
      INTGR4 nprocs
c
#include "scalapack.fh"
#include "global.fh"
c
      if(.not.init4)return
      call ga_sync()
      call blacs_gridexit(islctxt2)
      init4=.false.
      return
c    
      end

      subroutine FindGrid(nnodes, nprow, npcol)
c
c***  determine nprow, npcol from nnodes
c***  solution is searched in the neighborhood of the square grid 
c
      implicit none
      INTGR4 nnodes, nprow, npcol,i 
c
c     try to get the 1:4 ratio
c
      npcol = 2*int(sqrt(dble(nnodes)))
      do i = npcol, 1, -1 
         if(mod(nnodes,i).eq.0) goto 1
      enddo
1     continue
      npcol = i
      nprow = nnodes/npcol
      if(nprow.gt.npcol) then
         i=npcol
         npcol=nprow
         nprow=i
      endif
      end


      subroutine ga_to_SL(g_a, dim1, dim2, nbr, nbc, s_a, lda4, np, nq)
c
c***  transforms a GA to SL format
c***  reference: Dongarra et al, 'A look at scalable dense lin. alg. libs'
c
c   g_a         - handle for Global Array that is being converted SL format
c   dim1, dim2  - dimensions of Global Array
c   nbr         - number of block rows?
c   nbc         - number of block columns?
c   s_a         - local array holding SL formatted data
c   lda4        - leading dimension of s_a
c   np          - number of rows held by processor
c   nq          - number of columns held by processor
      implicit none
#include "scalapack.fh"
      INTGR4  nbr, nbc, lda4, np, nq, dim1, dim2
      integer g_a
      double  precision  s_a(lda4,*)
      integer row, col, tcol, trow, rbase, cbase
      INTGR4 pcol, prow,one4
      parameter(one4=1)
      integer lda,rowl,coll
      integer marg1,marg2

c***  Synchronize at the beginning
c    
      lda=lda4
      rbase = 1
      cbase = 1
c
      Trow = nbr * nprow
      Tcol = nbc * npcol
c     
      do col  = 1, dim2, nbc
         pcol = mod(col,Tcol)/nbc        ! processor column that holds "col"
         if(mycol.eq.pcol) then
            do row  = 1, dim1, nbr
               prow = mod(row,Trow)/nbr  ! processor row that holds "row"
               if(myrow.eq.prow) then
                  if(rbase.gt.np .or. cbase .gt. nq) then
                     call ga_error(' ga_to_SL: rows/cols error ',1)
                  endif
                  marg1=row+nbr-one4
                  marg2=dim1
                  rowl=min(marg1,marg2)
                  marg1=col+nbc-one4
                  marg2=dim2
                  coll=min(marg1,marg2)
                  call ga_get(g_a,row,rowl,col,coll,
     &                 s_a(rbase,cbase), lda) 
                  rbase = rbase + nbr
               endif
            enddo
            rbase = 1
            cbase = cbase + nbc
         endif
      enddo

c**** ... and at the end
      end



      subroutine ga_from_SL(g_a,dim1,dim2,nbr, nbc, s_a, lda4, np, nq)
c
c***  transforms a matrix from SL to GA format
c***  reference: Dongarra et al, 'A look at scalable dense lin. alg. libs'
c
      implicit none
#include "scalapack.fh"
      integer g_a,lda
      INTGR4  nbr, nbc, lda4, np, nq, dim1, dim2
      double  precision s_a(lda4,*)
      integer  tcol, trow
      integer rbase, cbase
      INTGR4 pcol, prow,one4
      parameter(one4=1)
      integer row,col,rowl,coll
      integer marg1,marg2
c
      lda=lda4
c
c**** Syncronize at the beginning
c   
      rbase = 1
      cbase = 1
c
      Trow = nbr * nprow
      Tcol = nbc * npcol
c
      do col  = 1, dim2, nbc
         pcol = mod(col,Tcol)/nbc        ! processor column that holds "col"
         if(mycol.eq.pcol) then
            do row  = 1, dim1, nbr
               prow = mod(row,Trow)/nbr  ! processor row that holds "row"
               if(myrow.eq.prow) then
                  if(rbase.gt.np .or. cbase .gt. nq)
     &               call ga_error(' ga_from_SL: rows/cols error ',1)
                  marg1=row+nbr-one4
                  marg2=dim1
                  rowl=min(marg1,marg2)
                  marg1=col+nbc-one4
                  marg2=dim2
                  coll=min(marg1,marg2)
                  call ga_put(g_a,row,rowl,col,
     &                    coll, s_a(rbase,cbase), lda) 
                  rbase = rbase + nbr
               endif
            enddo
            rbase = 1
            cbase = cbase + nbc
         endif
      enddo

c**** ... and at the end
      end

      subroutine ga_to_SL2(g_a, dim1, dim2, nbr, nbc, s_a, lda4, np, nq)
c
c***  transforms a GA to SL format
c***  reference: Dongarra et al, 'A look at scalable dense lin. alg. libs'
c
c   g_a         - handle for Global Array that is being converted SL format
c   dim1, dim2  - dimensions of Global Array
c   nbr         - number of block rows?
c   nbc         - number of block columns?
c   s_a         - local array holding SL formatted data
c   lda4        - leading dimension of s_a
c   np          - number of rows held by processor
c   nq          - number of columns held by processor
      implicit none
#include "scalapack.fh"
#include "global.fh"
      INTGR4  nbr, nbc, lda4, np, nq, dim1, dim2
      integer g_a
      double  precision  s_a(lda4,*)
      integer row, col, tcol, trow, rbase, cbase
      INTGR4 pcol, prow
      integer lda
      integer r0i,r1i,r0im1,r1im1,rowl,coll
      integer marg1,marg2
      logical putpending

c***  Synchronize at the beginning
c    
      lda=lda4
      rbase = 1
      cbase = 1
c
      Trow = nbr * nprow2
      Tcol = nbc * npcol2
      do col  = 1, dim2, nbc
         pcol = mod(col,Tcol)/nbc        ! processor column that holds "col"
         if(mycol2.eq.pcol) then
            r0im1=-9999
            r1im1=0
            putpending=.false.
            marg1=col+nbc-1
            marg2=dim2
            coll=min(marg1,marg2)

            do row  = 1, dim1, nbr
               prow = mod(row,Trow)/nbr  ! processor row that holds "row"
               if(myrow2.eq.prow) then
                  if(.not.putpending) then
                     r0im1=row
                     marg1=row+nbr-1
                     marg2=dim1
                     r1im1=min(marg1,marg2)
                  endif
                  r0i=row
                  marg1=row+nbr-1
                  marg2=dim1
                  r1i=min(marg1,marg2)
                  
                  if(r0i.eq.(r1im1+1).and.r1i.ne.dim1) then
                     r1im1=r1i
                     putpending=.true.
                  else
                  rowl=r1i
                  call ga_get(g_a,r0im1,rowl,col,
     &                    coll, s_a(rbase,cbase), lda) 
                     putpending=.false.
                     rbase = rbase + rowl - r0im1 +1
                     r0im1=-1
                     r1im1=-1
                  endif
               endif
            enddo
            rbase = 1
            cbase = cbase + nbc
         endif
      enddo

c**** ... and at the end
      end



      subroutine ga_from_SL2(g_a,dim1,dim2,nbr, nbc, s_a, lda4, np, nq)
c
c***  transforms a matrix from SL to GA format
c***  reference: Dongarra et al, 'A look at scalable dense lin. alg. libs'
c
      implicit none
#include "scalapack.fh"
#include "global.fh"
      integer g_a,lda
      INTGR4  nbr, nbc, lda4, np, nq, dim1, dim2
      double  precision s_a(lda4,*)
      integer  tcol, trow
      integer rbase, cbase
      INTGR4 pcol, prow
      integer row,col,rowl,coll
      integer r0i,r1i,r0im1,r1im1
      integer marg1,marg2
      logical putpending
c
      lda=lda4
c
c**** Syncronize at the beginning
c   
      rbase = 1
      cbase = 1
c
      Trow = nbr * nprow2
      Tcol = nbc * npcol2
c
      do col  = 1, dim2, nbc
         pcol = mod(col,Tcol)/nbc        ! processor column that holds "col"
         if(mycol2.eq.pcol) then
            r0im1=-9999
            r1im1=0
            putpending=.false.
            marg1=col+nbc-1
            marg2=dim2
            coll=min(marg1,marg2)
            do row  = 1, dim1, nbr
               prow = mod(row,Trow)/nbr  ! processor row that holds "row"
               if(myrow2.eq.prow) then
                  if(.not.putpending) then
                     r0im1=row
                     marg1=row+nbr-1
                     marg2=dim1
                     r1im1=min(marg1,marg2)
                  endif
                  r0i=row
                  marg1=row+nbr-1
                  marg2=dim1
                  r1i=min(marg1,marg2)
                  
                  if(r0i.eq.(r1im1+1).and.r1i.ne.dim1) then
!                  if(.false.) then
                     r1im1=r1i
                     putpending=.true.
                  else
                  rowl=r1i
                  call ga_put(g_a,r0im1,rowl,col,
     &                    coll, s_a(rbase,cbase), lda) 
                     putpending=.false.
                     rbase = rbase + rowl - r0im1 +1
                     r0im1=-1
                     r1im1=-1
                  endif
               endif
            enddo
            rbase = 1
            cbase = cbase + nbc
         endif
      enddo

c**** ... and at the end
      end
      subroutine ga_to_sl2cpl(
     G     g_a, dim1, dim2, nbr, nbc, s_a, lda4, np, nq)
c
c***  transforms a GA to SL format
c***  reference: Dongarra et al, 'A look at scalable dense lin. alg. libs'
c
c   g_a(1)      - handle for Global Array Real that is being converted SL format
c   g_a(2)      - handle for Global Array Img. that is being converted SL format
c   dim1, dim2  - dimensions of Global Array
c   nbr         - number of block rows?
c   nbc         - number of block columns?
c   s_a         - local array holding SL formatted data
c   lda4        - leading dimension of s_a
c   np          - number of rows held by processor
c   nq          - number of columns held by processor
      implicit none
#include "scalapack.fh"
#include "mafdecls.fh"
      INTGR4  nbr, nbc, lda4, np, nq, dim1, dim2
      integer g_a(2)
      complex*16  s_a(lda4,nq)
      integer row, col, tcol, trow, rbase, cbase
      INTGR4 pcol, prow
      integer lda
      integer rowl,coll
      integer ha1,adra1,ha2,adra2
      integer nbr8,nbc8
      integer marg1,marg2

c***  Synchronize at the beginning
c    
      lda=lda4
      nbr8=nbr
      nbc8=nbc
      rbase = 1
      cbase = 1
c
      Trow = nbr * nprow2
      Tcol = nbc * npcol2
c     
      if(.not.ma_push_get(MT_DBL,nbr8*nbc8,'a1',ha1,adra1))
     &       call ga_error('ga2sl2cpl: mem alloc failed A1 ', -1)
      if(.not.ma_push_get(MT_DBL,nbr8*nbc8,'a2',ha2,adra2))
     &       call ga_error('ga2sl2cpl: mem alloc failed A2 ', -1)
      do col  = 1, dim2, nbc
         pcol = mod(col,Tcol)/nbc        ! processor column that holds "col"
         if(mycol2.eq.pcol) then
            do row  = 1, dim1, nbr
               prow = mod(row,Trow)/nbr  ! processor row that holds "row"
               if(myrow2.eq.prow) then
                  if(rbase.gt.np .or. cbase .gt. nq) then
                     call ga_error(' gatoSLcpl: rows/cols error ',1)
                  endif
                  marg1=row+nbr -1
                  marg2=dim1
                  rowl=min(marg1,marg2)
                  marg1=col+nbc -1
                  marg2=dim2
                  coll=min(marg1,marg2)
                  call ga_get(g_a(1),row,rowl,col,
     &                    coll, dbl_mb(adra1), nbr8) 
                  call ga_get(g_a(2),row,rowl,col,
     &                    coll, dbl_mb(adra2), nbr8) 

                  call ga_cpydbl2cpl(s_a,lda,rbase,cbase,
     I                 dbl_mb(adra1),dbl_mb(adra2),nbr8,
     L                 rowl-row+1,coll-col+1)

                  rbase = rbase + nbr
               endif
            enddo
            rbase = 1
            cbase = cbase + nbc
         endif
      enddo
        if(.not.ma_chop_stack(ha1))
     G     call ga_error(' gatoSLcpl: chop_stack failed ',1)

      end
      subroutine ga_cpydbl2cpl(s_a,lda_s,rbase,cbase,
     I     a1,a2,lda_a,
     L     rowl,coll)
      implicit none
      integer lda_s,lda_a
      complex*16 s_a(lda_s,*)
      double precision a1(lda_a,*)
      double precision a2(lda_a,*)
      integer rbase,cbase,rowl,coll
      integer ccol,rrow
      do ccol=0,coll-1
         do rrow=0,rowl-1
            s_a(rbase+rrow , cbase+ccol )= dcmplx(
     D           a1(rrow+1,ccol+1),a2(rrow+1,ccol+1))
         enddo
      enddo
      return
      end
c
      subroutine ga_from_sl2cpl(
     G     g_a,dim1,dim2,nbr, nbc, s_a, lda4, np, nq)
c
c***  transforms a matrix from SL to GA format
c***  reference: Dongarra et al, 'A look at scalable dense lin. alg. libs'
c
      implicit none
#include "scalapack.fh"
#include "mafdecls.fh"
      integer g_a(2),lda
      INTGR4  nbr, nbc, lda4, np, nq, dim1, dim2
      complex*16 s_a(lda4,*)
      integer  tcol, trow
      integer rbase, cbase
      INTGR4 pcol, prow
      integer row,col,rowl,coll
      integer ha1,adra1,ha2,adra2
      integer nbr8,nbc8
      integer marg1,marg2
c
      lda=lda4
      nbr8=nbr
      nbc8=nbc
c
c**** Syncronize at the beginning
c   
      rbase = 1
      cbase = 1
c
      Trow = nbr * nprow2
      Tcol = nbc * npcol2
      if(.not.ma_push_get(MT_DBL,nbr8*nbc8,'a1',ha1,adra1))
     &       call ga_error('sl2gacpl: mem alloc failed A1 ', -1)
      if(.not.ma_push_get(MT_DBL,nbr8*nbc8,'a2',ha2,adra2))
     &       call ga_error('sl2gacpl: mem alloc failed A2 ', -1)
c
      do col  = 1, dim2, nbc
         pcol = mod(col,Tcol)/nbc        ! processor column that holds "col"
         if(mycol2.eq.pcol) then
            do row  = 1, dim1, nbr
               prow = mod(row,Trow)/nbr  ! processor row that holds "row"
               if(myrow2.eq.prow) then
                  if(rbase.gt.np .or. cbase .gt. nq)
     &               call ga_error(' ga_from_SL: rows/cols error ',1)
                  marg1=row+nbr-1
                  marg2=dim1
                  rowl=min(marg1,marg2)
                  marg1=col+nbc-1
                  marg2=dim2
                  coll=min(marg1,marg2)
c                  call ga_cpycpl2dbl(s_a,lda4,rbase,cbase,
c     I                 dbl_mb(adra1),dbl_mb(adra2),nbr8,
c     L                 rowl,coll)
                  call ga_cpycpl2dbl(s_a,lda,rbase,cbase,
     I                 dbl_mb(adra1),dbl_mb(adra2),nbr8,
     L                 rowl-row+1,coll-col+1)
                  call ga_put(g_a(1),row,rowl,col,
     &                    coll, dbl_mb(adra1), nbr8) 

                  call ga_put(g_a(2),row,rowl,col,
     &                    coll, dbl_mb(adra2), nbr8) 

                  rbase = rbase + nbr
               endif
            enddo
            rbase = 1
            cbase = cbase + nbc
         endif
      enddo
      if(.not.ma_chop_stack(ha1))
     G     call ga_error(' SL2gacpl: chop_stack failed ',1)

c**** ... and at the end
      end
      subroutine ga_cpycpl2dbl(s_a,lda_s,rbase,cbase,
     I     a1,a2,lda_a,
     L     lrow,lcol)
      implicit none
#include "mafdecls.fh"
      integer lda_s,lda_a
      complex*16 s_a(lda_s,*)
      double precision a1(lda_a,*)
      double precision a2(lda_a,*)
      integer rbase,cbase,lrow,lcol
      integer ccol,rrow
      do ccol=0,lcol-1
         do rrow=0,lrow-1
            a1(rrow+1,ccol+1)=dreal(s_a(rbase+rrow , cbase+ccol ))
            a2(rrow+1,ccol+1)=dimag(s_a(rbase+rrow , cbase+ccol ))
         enddo
      enddo
      return
      end


c**********************************************************************
c* Source : gai_lu_solve
c* Scope  : LU factorization
c*
c* 02/16/94  Jarek Nieplocha 
c* 04/12/96  GVT  Adapted the code to the new version of ScaLAPACK
c*           Giuseppe Vitillaro peppe@unipg.it
c**********************************************************************      
c
c length of scalapack array descriptor
# define DSCLEN 10
c
      subroutine gai_lu_solve(trans,g_a, g_b)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "scalapack.fh"
c
      character*1 trans         ! transpose or not  
      integer g_a               ! coefficient matrix A 
      integer g_b               ! rhs matrix, B, overwritten on exit by
                                ! the solution vector, X
c
c     solve the set of linear equations 
c
c           AX = B
c
c     with possibly multiple rhs stored as columns of matrix B
c
c     the matrix A is not destroyed
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      logical status
c
      integer ha, adra          !A
      integer hb, adrb          !B
      integer hi, adri          !ipiv
c
c
      logical oactive           ! true iff this process participates
      integer dimA18, dimA28, typeA
      integer dimB18, dimB28, typeB
      INTGR4 dimA1, dimA2
      INTGR4 dimB1, dimB2
      logical uses_sl_A, uses_sl_B
      integer alen, blen
c
      INTGR4 mpA, nqA          ! number of rows/cols of A held by the processor
      INTGR4 mpB, nqB          ! number of rows/cols of B held by the processor
      INTGR4 info
      integer me
      INTGR4 n,lda, ldb
      integer elemA, elemB
      INTGR4 numroc
      integer lipiv
      integer block_dims_A(2),block_dims_B(2),blocks(2)
      integer gridA(2), gridB(2)

      INTGR4 one4,zero4
      parameter(zero4=0,one4=1)
      INTGR4 descA(DSCLEN), descB(DSCLEN) ! descriptors for distr. matrixes A,B
      INTGR4 nb                           ! block size for A and rhs
      data nb /64/
      integer i
      logical use_direct
c
c***  check environment
c
      me     = ga_nodeid()
c
c***  check GA info for input arrays
c
      call ga_check_handle(g_a, 'gai_lu_solve: a')
      call ga_check_handle(g_b, 'gai_lu_solve: b')
      call ga_inquire(g_a, typeA, dimA18, dimA28)
      call ga_inquire(g_b, typeB, dimB18, dimB28)
      dima1=dima18
      dima2=dima28
      dimb1=dimb18
      dimb2=dimb28
c
      if (dimA1 .ne. dimA2) then
         call ga_error('gai_lu_solve: g_a must be square matrix ', 1)
      else if(dimA1 .ne. dimB1) then
         call ga_error('gai_lu_solve: dims of A and B do not match ',1)
      endif
c
c   Determine whether or not to use GA data directly
c
      use_direct = .true.
      uses_sl_A = ga_uses_proc_grid(g_a)
      uses_sl_B = .false.
      if (uses_sl_A) then
        uses_sl_B = ga_uses_proc_grid(g_b)
      endif
      if ((.not.uses_sl_A).or.(.not.uses_sl_B)) then
        use_direct = .false.
      endif
      if (uses_sl_A) then
        call ga_get_block_info(g_a,blocks,block_dims_A)
        if (block_dims_A(1).ne.block_dims_A(2)) then
          use_direct = .false.
        endif
      endif
      if (uses_sl_B) then
        call ga_get_block_info(g_b,blocks,block_dims_B)
        if (block_dims_B(1).ne.block_dims_B(2)) then
          use_direct = .false.
        endif
      endif
      if (use_direct.and.(block_dims_A(1).ne.block_dims_B(1))) then
        use_direct = .false.
      endif
      if (use_direct) then
        call ga_get_proc_grid(g_a,gridA)
        call ga_get_proc_grid(g_b,gridB)
        if (gridA(1).ne.gridB(1).or.gridA(2).ne.gridB(2)) then
          use_direct = .false.
        endif
      endif
      if (use_direct) then
        nb = block_dims_A(2)
      endif
c 
      n = dimA1
c
c***  initialize SL interface
c 
      if (.not.use_direct) then
        call SLinit()
      else
        call SLinit3(g_a)
      endif
c
c***  find SBS format parameters
c
      mpA = numroc(n, nb,myrow, zero4, nprow)
      nqA = numroc(n, nb,mycol, zero4, npcol)
      mpB = numroc(dimB1, nb,myrow, zero4, nprow)
      nqB = numroc(dimB2, nb,mycol, zero4, npcol) 
      lda = max(one4,mpA)
      ldb = max(one4,mpB)

      oactive = .true.
      if (oactive) then
c
c***     allocate A,B and ipiv arrays
c     
         elemA= mpA*nqA
         status = .true.
         if (use_direct) then
           call nga_access_block_segment(g_a,me,adra,alen)
         else
           if(elemA.ne.0) status = ma_push_get(MT_DBL,elemA,'a',ha,adra)
           if(.not.status) 
     &       call ga_error('gai_lu_solve: mem alloc failed A ', -1)
         endif
c
         elemB= mpB*nqB
         if (use_direct) then
           call nga_access_block_segment(g_b,me,adrb,blen)
         else
           if(elemB.ne.0)status = ma_push_get(MT_DBL,elemB,'b',hb,adrb)
           if(.not.status) 
     &       call ga_error('gai_lu_solve: mem alloc failed B ', -1)
         endif
c
         lipiv=mpa+nb
         status = ma_push_get(MT_INT,lipiv,'ipiv',hi,adri)
         if(.not.status) 
     &       call ga_error('gai_lu_solve: mem alloc failed ipiv ', -1)
c
c
c***     copy g_a to A and g_b to B using the SBS SL format 
c
         if (.not.use_direct) then
           call ga_to_SL(g_a,n,n,nb,nb,dbl_mb(adrA),lda,mpA,nqA)
           call ga_to_SL(g_b,n,dimB2,nb,nb,dbl_mb(adrB),ldb,mpB,nqB)
         endif
c     
c***     fill SCALAPACK matrices descriptors
c
         call descinit(descA, dimA1, dimA2, nb, nb, zero4, zero4, 
     &                 iSLctxt, lda, info)
         if(info.ne.0) call ga_error('gai_lu_solve: descinit A failed ',
     &                               -info)
         call descinit(descB, dimB1, dimB2, nb, nb, zero4, zero4, 
     &                 iSLctxt, ldb, info)
         if(info.ne.0) call ga_error('gai_lu_solve: descinit B failed ',
     &                               -info)
c
c***     LU factorization
c
         call ga_sync()
         call pdgetrf(n, n, dbl_mb(adrA),one4,one4, descA, 
     &                int_mb(adrI), info)
         call ga_sync()
c
c
         if(info.eq.0)then
c
c***        solve: forward/backward substitution
c
            call pdgetrs(trans, n, dimB2, dbl_mb(adrA),one4,one4, 
     A           descA, int_mb(adrI), dbl_mb(adrB),
     1           one4, one4, descB, info)
            if(info.eq.0)then
              if (.not.use_direct) then
c
c***           copy solution matrix back to g_b
c
                call ga_from_SL(g_b, dimB1, dimB2, nb, nb, dbl_mb(adrb),
     &                          mpB, ldb, nqB)
              endif
            else
               call ga_error('gai_lu_solve: pdgetrs failed ', -info)
            endif
c
         else
            call ga_error('gai_lu_solve: pdgetrf failed ', -info)
         endif
c
c***     deallocate work/SL arrays
c
         status = ma_pop_stack(hi)
         if (.not.use_direct) then
           if(elemB.ne.0) status = status .and. ma_pop_stack(hb)
           if(elemA.ne.0) status = status .and. ma_pop_stack(ha)
         endif
         if(.not.status)
     &       call ga_error('gai_lu_solve:ma_pop_stack failed',0)
      endif
      if (use_direct) then
        call slexit3
      else
        call slexit
      endif
c
      call ga_sync()
      end

c     alt entry point for C that avoids char args
c
      subroutine gai_lu_solve_alt(tr, g_a, g_b)
      implicit none
        integer g_a, g_b, tr
        if(tr.eq.0) then
           call gai_lu_solve('n',g_a, g_b) 
        else
           call gai_lu_solve('y',g_a, g_b) 
        endif
        end


c**********************************************************************
c* Source : ga_llt
c* Scope  : interface between GA and LLT SCALAPACK routines
c*
c* 04/12/96  GVT  First Implementation
c*           Giuseppe Vitillaro peppe@unipg.it
c* 08/28/96  Jarek Nieplocha, some cleanup and optimizations 
c**********************************************************************

#define NB_ 64
#define DSCLEN 10
#define ga_dnormF(g_a) sqrt(ga_ddot(g_a, g_a)
      
c**********************************************************************
c* subroutine: ga_zeroUL
c*
c*             Set to zero the L/U triangle part of a NxN
c*             double precision global array A
c***********************************************************************
      subroutine ga_zeroUL(uplo, g_A)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
c****
      character*1 uplo         ! (input) 'U' or 'L'
      integer     g_A          ! (input/output)
      logical status
      integer dimA1, dimA2, typeA
      integer me, nproc
      integer n
      integer i, j, hBUF, adrBUF

c**** Check Environment
      nproc = ga_nnodes()
      me    = ga_nodeid()
      
c**** Check GA info for input array A
      call ga_check_handle(g_A, 'ga_zeroUL: A')
      call ga_inquire(g_A, typeA, dimA1, dimA2)
c****
      if (dimA1.ne.dimA2) then
         call ga_error('ga_zeroUL: g_A must be a square matrix ', 1)
      endif
      if (typeA.ne.MT_DBL) then
         call ga_error('ga_zeroUL: g_A must be double precision ', 1)
      endif
c****      
      n = dimA1
      
c**** Allocate BUF
      status = ma_push_get(MT_DBL, n, 'BUF', hBUF, adrBUF)
      if (.not.status)
     &     call ga_error(' ga_zeroUL: mem alloc failed BUF ', -1)
      
      call ga_sync()

      do i = me+1, n, nproc
         call ga_get(g_A, 1, n, i, i, dbl_mb(adrBUF), n)
         if (uplo.eq.'L') then
c****       case L: make zero the upper triangle            
            call dcopy(i-1,0.0d0,0, dbl_mb(adrBUF),1)
         elseif (uplo.eq.'U') then
c****       case U: make zero the lower triangle            
            call dcopy(n-i,0.0d0,0, dbl_mb(adrBUF+i),1)
         else
            call ga_error('ga_symUL: uplo must be L or U ', 1)
         endif
         call ga_put(g_A, 1, n, i, i, dbl_mb(adrBUF), n)
      end do    !i
c
      status = ma_pop_stack(hBUF)
      call ga_sync()
      end


c**********************************************************************
c* subroutine: ga_symUL
c*
c*             Make a symmetric square matrix from
c*             double precision global array A in L/U triangle format
c***********************************************************************
      subroutine ga_symUL(uplo, g_A)

      implicit none
      
#include "mafdecls.fh"
#include "global.fh"
c****      
      character*1 uplo         ! (input) 'U' or 'L'
      integer     g_A          ! (input/output)
      logical status
      integer dimA1, dimA2, typeA
      integer me, nproc
      integer n
      integer i, j, hBUF, adrBUF, idx      

c**** Check Environment
      nproc = ga_nnodes()
      me    = ga_nodeid()
      
c**** Check GA info for input array A
      call ga_check_handle(g_A, 'ga_symUL: A')
      call ga_inquire(g_A, typeA, dimA1, dimA2)
c****
      if (dimA1.ne.dimA2) then
         call ga_error('ga_symUL: g_A must be a square matrix ', 1)
      endif
      if (typeA.ne.MT_DBL) then
         call ga_error('ga_symUL: g_A must be double precision ', 1)
      endif
c****      
      n = dimA1
c**** Allocate BUF
      status = .true.
      status = ma_push_get(MT_DBL, n, 'BUF', hBUF, adrBUF)
      if (.not.status)
     &     call ga_error(' ga_symUL: mem alloc failed BUF ', -1)
      
      call ga_sync()

      do i = me+1, n, nproc
         call ga_get(g_A, 1, n, i, i, dbl_mb(adrBUF), n)
         if (uplo.eq.'L') then
c****       case L : make sure upper triangle is zero
c****                and under the diag. is 2 times                 
            do j = 1, n
               idx = adrBUF + j - 1
               if (j.lt.i) then
                  dbl_mb(idx) = 0.d0
               elseif (j.gt.i) then
                  dbl_mb(idx) = 2.d0*dbl_mb(idx)
               endif
            end do ! j
         elseif (uplo.eq.'U') then
c****       case U : make sure lower triangle is zero
c****                and over the diag is 2 times               
            do j = 1, n
               idx = adrBUF + j - 1
               if (j.gt.i) then
                  dbl_mb(idx) = 0.d0
               elseif (j.lt.i) then
                  dbl_mb(idx) = 2.d0*dbl_mb(idx)
               endif
            end do ! j
         else
            call ga_error('ga_symUL: uplo must be L or U ', 1)
         endif
         call ga_put(g_A, 1, n, i, i, dbl_mb(adrBUF), n)
      end do    ! i
      
      status = ma_pop_stack(hBUF)
      call ga_symmetrize(g_A)
      end

      
c***********************************************************************
c* function :  ga_cholesky
c*
c*             Compute the Cholesky factorization of an NxN
c*             double precision symmetric positive definite matrix.
c*
c*             On succesful exit A will contain the L/U factor
c*             on the lower/upper triangular part of the matrix
c*                   
c*             It calls the PDPOTRF ScaLAPACK routine.
c*
c*             It returns
c*                   = 0 : successful exit
c*                   > 0 : the leading minor of this order
c*                         is not positive definite and the
c*                         factorization could not be completed
c*
c***********************************************************************
      integer function ga_cholesky(uplo, g_a)
      implicit none
      character*1 uplo         ! (input) 'U' or 'L'
      integer     g_A          ! (input/output)
      integer ga_llt_f
      ga_cholesky = ga_llt_f(uplo, g_A, -1)
      end



c***********************************************************************
c* function :  ga_llt_f
c*
c*             Internal function to compute the Cholesky factorization of 
c*             an NxN double precision symmetric positive definite GA.
c*
c*             Note: applications should use ga_cholesky
c*
c*             On succesful exit A will contain the L/U factor
c*             on the lower/upper triangular part of the matrix
c*
c*             This if (hsA.eq.-1), otherwise internal state
c*             is saved for future references. This is
c*             reserved for internal use. Users should
c*             not try to set hsA to something
c*             different from -1 if they do not
c*             know what they are doing! On exit hsA will
c*             contain the local A output array in SLS format.      
c*
c*             It calls the PDPOTRF ScaLAPACK routine.
c*      
c*             It returns
c*                   = 0 : successful exit
c*                   > 0 : the leading minor of this order
c*                         is not positive definite and the
c*                         factorization could not be completed      
c***********************************************************************
      integer function ga_llt_f(uplo, g_A, hsA)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "scalapack.fh"
c****
      character*1 uplo         ! (input) 'U' or 'L'
      integer     g_A          ! (input/output)
      integer     hsA          ! (input/output)
c****
      logical status
      integer hA, adrA
      integer dimA18, dimA28, typeA
      INTGR4 dimA1, dimA2
      INTGR4 mpA, nqA
      INTGR4 info
      integer me
      INTGR4 n, ldA
      integer elemA
      INTGR4 numroc
      INTGR4 nb
      INTGR4 descA(DSCLEN)
      integer intsize,info8
      INTGR4 one4,zero4
      parameter(zero4=0,one4=1)
      logical oactive
      integer alen,block_dims(2),blocks(2)
      logical use_direct
      data nb /NB_/

c**** Check Environment
      me = ga_nodeid()

c**** Check GA info for input array A
      call ga_check_handle(g_A, 'ga_llt_f: A')
      call ga_inquire(g_A, typeA, dimA18, dimA28)
      dima1=dima18
      dima2=dima28
c****
      if (dimA1.ne.dimA2) then
         call ga_error('ga_llt_f: g_A must be a square matrix ', 1)
      endif
c
c   Determine whether or not to use GA data directly
c
      use_direct = ga_uses_proc_grid(g_A)
      if (use_direct) then
        call ga_get_block_info(g_A,blocks,block_dims)
        if (block_dims(1).ne.block_dims(2)) then
          use_direct = .false.
        endif
      endif
      if (use_direct) then
        nb = block_dims(1)
      endif
c****
      n = dimA1
c**** Initialize SL Interface
      if (use_direct) then
        call SLinit3(g_A)
      else
        call SLinit()
      endif
c      oactive=iam.lt.maxproc
      oactive=.true.
      call ga_sync

c**** Find SBS format parameters
      if(oactive) then
         mpA = numroc(n, nb, myrow, zero4, nprow)
         nqA = numroc(n, nb, mycol, zero4, npcol)
         ldA = max(one4,mpA)
         
         
c**** Allocate A
         elemA = mpA * nqA
         status = .true.
         if (use_direct) then
           call nga_access_block_segment(g_A,me,adrA,alen)
         else
           if (elemA.ne.0)
     &          status = ma_push_get(MT_DBL, elemA, 'A', hA, adrA)
           if (.not.status)
     &          call ga_error('ga_llt_f: mem alloc failed A ', -1)
         endif
         
c**** Copy ga to A using SBS ScaLAPACK format      
         if (.not.use_direct) then
           call ga_to_SL(g_a, n, n, nb, nb, dbl_mb(adrA), ldA, mpA, nqA)
         endif
         endif
         call ga_sync()
         if(oactive) then
         
c**** Fill ScaLAPACK matrix descriptor
         call descinit(descA, dimA1, dimA2, nb, nb, zero4, zero4,
     &        iSLctxt, ldA, info)
         if(info.ne.0) call ga_error('ga_llt_f: descinit A failed',
     I        -info)
         
c**** CALL ScaLAPACK PDPOTRF LLT factorization routine *******
         call PDPOTRF(uplo, n, dbl_mb(adrA), one4, one4, descA, info)
c*************************************************************
         
         if (info.eq.0) then
c**** Copy solution matrix back to A if hsA==-1
c**** and zero the L/U triangle part according to uplo         
            if (hsA.eq.-1) then
              if (.not.use_direct) then
                call ga_from_SL(g_A, dimA1, dimA2, nb, nb, dbl_mb(adrA),
     &             mpA, ldA, nqA)
              endif
            endif
            
c
c**** If the SL A array was allocated
            if (elemA.ne.0) then
c**** and hsA==-1 or info>0 (i.e. fact. cannot be completed)
c**** then deallocate the SL A MA array
               if ((hsA.eq.-1).or.(info.gt.0)) then
                  if (.not.use_direct) then
                    status = ma_pop_stack(hA)
                  endif
c**** otherwise just save the hA MA handle            
               else
                  hsA = hA
               endif
            endif
         endif
      endif
c
      call ga_sync()
#if 0
      info8=info
      if(maxproc.lt.nnodes) then
         intsize=ma_sizeof(MT_INT,1,MT_BYTE)
         call ga_brdcst(1688,info8,intsize,0)
      endif
      info=info8
#endif
      if (info.eq.0) then
         if (hsA.eq.-1) then
            call ga_zeroUL(uplo, g_A)
         endif
      elseif (info.lt.0) then
         call ga_error('ga_llt_f: PDPOTRF failed ', -info)
      endif
      ga_llt_f = info
c
      if (use_direct) then
        call slexit3
      else
        call slexit
      endif
      end


c***********************************************************************
c* subroutine: ga_llt_s
c*
c*             Solves a system of linear equations
c*             
c*             A * X = B
c*             
c*             where A is the lower triangle L or the upper
c*             triangular Cholesky factor U of a NxN double precision
c*             symmetric positive definite global array (LL' or U'U)
c*             obtained from ga_llt_f routine.
c*
c*             If (hsA.eq.-1) then you "must" give in input
c*             the lower or upper triangular Cholesky factor.
c*             For internal use, if (hsA.ne.-1) then hsA is
c*             the MA handle of the Cholesky L/U factor already
c*             in SBS SL format.
c*
c*             On successful exit B will contain the solution X.      
c*      
c*             It calls the PDPOTRS ScaLAPACK routine.
c***********************************************************************
      subroutine ga_llt_s(uplo, g_A, g_B, hsA)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "scalapack.fh"
c****
      character*1 uplo         ! (input) 'U' or 'L'
      integer     g_A          ! (input)
      integer     g_B          ! (input/ouput)
      integer     hsA          ! (input)
c****
      logical status
      integer hA, adrA
      integer hB, adrB
      integer dimA18, dimA28, typeA
      integer dimB18, dimB28, typeB
      INTGR4 dimA1, dimA2
      INTGR4 dimB1, dimB2
      INTGR4 mpA, nqA
      INTGR4 mpB, nqB
      INTGR4 info
      integer me
      INTGR4 n, ldA, ldB
      integer elemA, elemB
      INTGR4 numroc
      INTGR4 nb
      INTGR4 descA(DSCLEN), descB(DSCLEN)
      INTGR4 one4,zero4
      parameter(zero4=0,one4=1)
c
      integer block_dims_A(2),block_dims_B(2),blocks(2)
      integer gridA(2), gridB(2)
      logical uses_sl_A, uses_sl_B
      integer alen, blen
      logical use_direct
c
      data nb /NB_/
c
c**** Check the Environment
      me = ga_nodeid()
      
c**** Check GA info for input array F, B
      call ga_check_handle(g_A, 'ga_llt_s: A')
      call ga_check_handle(g_B, 'ga_llt_s: B')
      call ga_inquire(g_A, typeA, dimA18, dimA28)
      call ga_inquire(g_B, typeB, dimB18, dimB28)
      dima1=dima18
      dima2=dima28
      dimb1=dimb18
      dimb2=dimb28
c****
      if (dimA1.ne.dimA2) then
         call ga_error('ga_llt_s: g_A must be a square matrix ', 1)
      else if (dimA1.ne.dimB1) then
         call ga_error('ga_llt_s: dims of A and B do not match ', 1)
      endif
c
c   Determine whether or not to use GA data directly
c
      use_direct = .true.
      uses_sl_A = ga_uses_proc_grid(g_a)
      uses_sl_B = .false.
      if (uses_sl_A) then
        uses_sl_B = ga_uses_proc_grid(g_b)
      endif
      if ((.not.uses_sl_A).or.(.not.uses_sl_B)) then
        use_direct = .false.
      endif
      if (uses_sl_A) then
        call ga_get_block_info(g_a,blocks,block_dims_A)
        if (block_dims_A(1).ne.block_dims_A(2)) then
          use_direct = .false.
        endif
      endif
      if (uses_sl_B) then
        call ga_get_block_info(g_b,blocks,block_dims_B)
        if (block_dims_B(1).ne.block_dims_B(2)) then
          use_direct = .false.
        endif
      endif
      if (use_direct.and.(block_dims_A(1).ne.block_dims_B(1))) then
        use_direct = .false.
      endif
      if (use_direct) then
        call ga_get_proc_grid(g_a,gridA)
        call ga_get_proc_grid(g_b,gridB)
        if (gridA(1).ne.gridB(1).or.gridA(2).ne.gridB(2)) then
          use_direct = .false.
        endif
      endif
      if (use_direct) then
        nb = block_dims_A(2)
      endif
c
      n = dimA1
      
c**** Initialize SL Interface
      if (use_direct) then
        call SLinit3(g_A)
      else
        call SLinit()
      endif
#if 1

c**** Find SBS format parameters
      mpA = numroc(n, nb, myrow, zero4, nprow)
      nqA = numroc(n, nb, mycol, zero4, npcol)
      mpB = numroc(n, nb, myrow, zero4, nprow)
      nqB = numroc(n, nb, mycol, zero4, npcol)
      ldA = max(one4,mpA)
      ldB = max(one4,mpB)

      call ga_sync()

c**** Allocate A or retrieve it from the ga_llt_f calls
      elemA = mpA * nqA
      status = .true.
      if (elemA.ne.0) then
        if (use_direct) then
          call nga_access_block_segment(g_A,me,adrA,alen)
        else
          if (hsA.eq.-1) then
            status = ma_push_get(MT_DBL, elemA, 'A', hA, adrA)
            if (.not.status)
     &        call ga_error('ga_llt_s: mem alloc failed A ', -1)
c****       copy g_A to A using SBS SL format            
            call ga_to_SL(g_A, n, n, nb, nb, dbl_mb(adrA),
     &                    ldA, mpA, nqA)
          else
            hA = hsA
            status = ma_get_index(hA, adrA)
            if (.not.status)
     &           call ga_error('ga_llt_s: get index failed A ', -1)
c****       if hsA was given the SBS SL conversion is avoided            
          endif
        endif
      endif
      
c**** Allocate B and copy g_B to B in SBS SL format
      elemB = mpB * nqB
      if (use_direct) then
        call nga_access_block_segment(g_B,me,adrB,blen)
      else
        if (elemB.ne.0) status = ma_push_get(MT_DBL,elemB,'B',hB,adrB)
        if (.not.status)
     &     call ga_error('ga_llt_s: mem alloc failed B ', -1)
        call ga_to_SL(g_B, n, dimB2, nb, nb, dbl_mb(adrB), ldB, mpB,nqB)
      endif
c
c**** Fill ScaLAPACK matrix descriptors for A and B
      call descinit(descA, dimA1, dimA2, nb, nb, zero4, zero4,
     &     iSLctxt, ldA, info)
      if(info.ne.0) call ga_error(' ga_llt_s: descinit A failed ',
     &     -info)
      
      call descinit(descB, dimB1, dimB2, nb, nb, zero4,zero4,
     D     iSLctxt, ldB, info)
      if(info.ne.0) call ga_error('ga_llt_s: descinit B failed', -info)
     
c**** CALL ScaLAPACK PDPOTRS solver routine ***********************
      call ga_sync()
      call PDPOTRS(uplo, n, dimB2,
     &     dbl_mb(adrA), one4, one4, descA,
     &     dbl_mb(adrB), one4, one4, descB,
     &     info)
c******************************************************************

      if (info.eq.0) then
c****    copy solution matrix back to g_B
         if (.not.use_direct) then
           call ga_from_SL(g_B, dimB1, dimB2, nb, nb, dbl_mb(adrB),
     &                   mpB, ldB, nqB)
         endif
      else
         call ga_error(' ga_llt_s: PDPOTRS failed: ', -info)
      endif

c**** deallocate work/SL arrays
c**** note that should not be others MA allocations
c**** between ga_llt_f and ga_llt_s beside A and B      
      if (.not.use_direct) then
        if (elemB.ne.0) status = ma_pop_stack(hB)
        if (elemA.ne.0) status = ma_pop_stack(hA)
      endif
c
#endif
      if (use_direct) then
        call slexit3
      else
        call slexit
      endif
c
      call ga_sync()
      end


c***********************************************************************
c* function :  ga_llt_i
c*
c*             It computes the inverse of a global array
c*             that is the lower triangle L or the upper
c*             triangular Cholesky factor U of a NxN double precision
c*             symmetric positive definite global array (LL' or U'U)
c*             obtained from the ga_llt_f routine.
c*
c*             If (hsA.eq.-1) then you "must" give in input
c*             the lower or upper triangular Cholesky factor.
c*             For internal use, if (hsA.ne.-1) then hsA is
c*             the MA handle of the Cholesky L/U factor already
c*             in SBS SL format.
c*
c*             On successful exit A will contain the inverse.
c*      
c*             It calls the PDPOTRI ScaLAPACK routine.
c*      
c*             It returns
c*                   = 0 : successful exit
c*                   > 0 : it returns the index i of the (i,i)
c*                         element of the factor L/U that is zero and,
c*                         so, the inverse could not be computed        
c***********************************************************************
      integer function ga_llt_i(uplo, g_A, hsA)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "scalapack.fh"
c****
      character*1 uplo         ! (input) 'U' or 'L'
      integer     g_A          ! (input/output)
      integer     hsA          ! (input/output)
c****
      logical status
      integer hA, adrA
      integer dimA18, dimA28, typeA
      INTGR4 dimA1, dimA2
      INTGR4 mpA, nqA
      INTGR4 info
      integer info8
      integer nproc, me
      INTGR4 n, ldA
      integer elemA
      INTGR4 numroc
      INTGR4 nb
      INTGR4 descA(DSCLEN)
      integer intsize
      logical oactive
      INTGR4 one4,zero4
      parameter(zero4=0,one4=1)
      integer alen, block_dims(2),blocks(2)
      logical use_direct
c****      
      data nb /NB_/
c
c**** Check Environment
      nproc = ga_nnodes()
      me = ga_nodeid()

c**** Check GA info for input array A
      call ga_check_handle(g_A, 'ga_llt_i: A')
      call ga_inquire(g_A, typeA, dimA18, dimA28)
      dima1=dima18
      dima2=dima28

c**** Check that is actually a square matrix
      if (dimA1.ne.dimA2) then
         call ga_error('ga_llt_i: g_A must be a square matrix ', 1)
      endif
c
c   Determine whether or not to use GA data directly
c
      use_direct = ga_uses_proc_grid(g_A)
      if (use_direct) then
        call ga_get_block_info(g_A,blocks,block_dims)
        if (block_dims(1).ne.block_dims(2)) then
          use_direct = .false.
        endif
      endif
      if (use_direct) then
        nb = block_dims(1)
      endif
c
      n = dimA1

c**** Initialize SL Interface
      if (use_direct) then
        call SLinit3(g_A)
      else
        call SLinit()
      endif
c      oactive=iam.lt.maxproc
      oactive=.true.
      call ga_sync()
      if(oactive) then
c**** Find SBS format parameters
        mpA = numroc(n, nb, myrow, zero4, nprow)
        nqA = numroc(n, nb, mycol, zero4, npcol)
        ldA = max(one4, mpA)


c**** Allocate A or retrieve it from ga_llt_f call
        elemA = mpA * nqA
           if (use_direct) then
             call nga_access_block_segment(g_A,me,adrA,alen)
           else
             if (elemA.ne.0) then
                status = .true.
                if (hsA.eq.-1) then
                   status = ma_push_get(MT_DBL, elemA, 'A', hA, adrA)
                   if (.not.status)
     &               call ga_error(' ga_llt_i: mem alloc failed A ', -1)
c****       copy g_A to A using SBS SL format
                   call ga_to_SL(g_A, n, n, nb, nb, dbl_mb(adrA),
     &                           ldA, mpA, nqA)
                else
                   hA = hsA
                   status = ma_get_index(hA, adrA)
                   if (.not.status)
     &               call ga_error(' ga_llt_i: get index failed A ', -1)
c****         if hsA was given the SBS SL conversion is avoided
                endif
             endif
           endif
        endif
        call ga_sync()
        if(oactive) then


c**** Fill ScaLAPACK matrix descriptor for A
        call descinit(descA, dimA1, dimA2, nb, nb, zero4, zero4,
     &       iSLctxt, ldA, info)
        if (info.ne.0) call ga_error(' ga_llt_i: descinit A failed ',
     &       -info)

c**** CALL ScaLAPACK PDPOTRI matrix inverter **********************
        call PDPOTRI(uplo, n, dbl_mb(adrA), one4, one4, descA, info)
c******************************************************************

        if (.not.use_direct) then
          if (info.eq.0) then
c****    Copy the inverse matrix back to A
c****    and symmetrize it         
             call ga_from_SL(g_A, dimA1, dimA2, nb, nb, dbl_mb(adrA),
     &            mpA, ldA, nqA)
         
          endif
          if(elemA.ne.0) status = ma_pop_stack(hA)
        endif
      endif
      call ga_sync()
#if 0
      info8=info
      if(maxproc.lt.nnodes) then
         intsize=ma_sizeof(MT_INT,1,MT_BYTE)
         call ga_brdcst(1688,info8,intsize,0)
      endif
      info=info8
#endif
      if (info.eq.0) then
         call ga_symUL(uplo, g_A)
      elseif (info.lt.0) then
         call ga_error(' ga_llt_i: PDPOTRI failed ', -info)
      endif

      ga_llt_i = info

c**** deallocate work/SL arrays
c**** note that should not be others MA allocations
c**** between ga_llt_f and ga_llt_i
c
      end

      
c***********************************************************************
c* function :  gai_llt_solve
c*
c*             Solves a system of linear equations
c*             
c*             A * X = B
c*
c*             using the Cholesky factorization of an NxN
c*             double precision symmetric positive definite
c*             global array A.
c*      
c*             On successful exit B will contain the solution X.      
c*      
c*             It calls the ga_llt_f and ga_llt_s and so it
c*             actually refer to PDPORTF/PDPOTRS ScaLAPACK
c*             routines.
c*      
c*             It returns
c*                   = 0 : successful exit
c*                   > 0 : the leading minor of this order
c*                         is not positive definite and the
c*                         factorization could not be completed      
c***********************************************************************      
      integer function gai_llt_solve(g_A, g_B)
      implicit none
      character*1 uplo         ! (internal parameter) 'U' or 'L'
      integer     g_A          ! (input)
      integer     g_B          ! (input/output)
      integer     hsA
      integer     irc
      integer     ga_llt_f
c
c**** call the Cholesky factorization routine
      hsA  = 0
      uplo = 'L'
      irc  = ga_llt_f(uplo, g_A, hsA)

c**** check if the factorization is OK
      if (irc.eq.0) then
c****    if the factorization is OK then solve the system
c****    call the ga_llt_s internal interface
         call ga_llt_s(uplo, g_A, g_B, hsA)
         gai_llt_solve = 0
      else
c****    if the factorization is not OK just return the error
         gai_llt_solve = irc
      endif
      end


c***********************************************************************
c* function :  gai_spd_invert
c*
c*             It computes the inverse of a  double precision
c*             using the Cholesky factorization of a NxN double
c*             precision symmetric positive definite global array A.
c*
c*             On successful exit A will contain the inverse.
c*      
c*             It calls the ga_llt_f and ga_llt_i and so it
c*             actually refer to PDPORTF/PDPOTRI ScaLAPACK
c*             routines.      
c*      
c*             It returns
c*                   = 0 : successful exit
c*                   > 0 : the leading minor of this order
c*                         is not positive definite and the
c*                         factorization could not be completed            
c*                   < 0 : it returns the index i of the (i,i)
c*                         element of the factor L/U that is zero and,
c*                         so, the inverse could not be computed        
c***********************************************************************
      integer function gai_spd_invert(g_A)
      implicit none
      character*1 uplo         ! (internal parameter) 'U' or 'L'
      integer     g_A          ! (input)
      integer     hsA
      integer     ircF, ircI
      integer     ga_llt_f
      integer     ga_llt_i
c
c**** call the Cholesky factorization routine
      hsA  = 0
      uplo = 'L' 
      ircF = ga_llt_f(uplo, g_A, hsA)
c
c**** check if the factorization is OK
      if (ircF.eq.0) then
c****    if the factorization is OK then invert the matrix
c****    call the ga_llt_i internal interface
         ircI = ga_llt_i(uplo, g_A, hsA)
         if (ircI.eq.0) then
            gai_spd_invert = 0
         else
            gai_spd_invert = -ircI
         endif
      else
c****    if the factorization is not OK just return the error
         gai_spd_invert = ircF
      endif
      end

c***********************************************************************
c* function: gai_solve
c*
c*             Solves a system of linear equations
c*             
c*             A * X = B
c*
c*             It first will call the Cholesky factorization
c*             routine and, if sucessfully, will solve the system
c*             with the Cholesky solver. If Cholesky will be not
c*             be able to factorize A, then it will call the LU
c*             factorization routine and will solve the system
c*             with forward/backward substitution.      
c*      
c*             On B will contain the solution X.      
c*      
c*             It calls the ga_llt_f and ga_llt_s and gai_lu_solve
c*             and so it actually refer to PDGETRF/PDGETRS and
c*             PDPORTF/PDPOTRS ScaLAPACK routines.
c*
c*             It returns
c*                   = 0 : Cholesky factoriztion was succesful
c*                   > 0 : the leading minor of this order
c*                         is not positive definite and Cholesky 
c*                         factorization could not be completed
c*
c*    Note: High overhead due to multiple conversions to/from scalapack
c*          format can be avoided if gai_llt_solve and gai_lu_solve
c*          were not used as building blocks - possible optimization
c*          for the future (JN) 
c***********************************************************************
      integer function gai_solve(g_A, g_B)
      implicit none
      integer     g_A          ! (input)
      integer     g_B          ! (input/output)
      integer     irc
      integer     gai_llt_solve
c
c**** Call Cholesky solver as the first try      
      gai_solve = gai_llt_solve(g_A, g_B)
      if (gai_solve.eq.0) then
c****    if Cholesky was succcesully just return 
      else
c****    else if Cholesky failed
         call gai_lu_solve('N', g_A, g_B)
      endif
      end
c******************************************************************
c     
c     ga_pdsyev
c
c     interface into scalapack's Hermitian eigensolver
c     
c     ga interface to scalapack
c     
c     g_a .g_b(*,i) = eval(i).g_b(*,i)
c     
c     assume equal size blocks ...
c     
c******************************************************************
c     
      subroutine ga_pdsyev(g_a, g_b, eval, nb8)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "scalapack.fh"
c     
      integer g_a               ! matrix A
      integer g_b               ! matrix B
      integer nb8                ! block size
c
c all eigenvalues are real
c
      double precision eval(*)
c     
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      logical status
c
      character*1 jobz,  uplo
c
      integer ha, adra          !A
      integer hb, adrb          !B
c
c
      logical oactive           ! true iff this process participates
      integer dimA18, dimA28, typeA
      integer dimB18, dimB28, typeB
      INTGR4 dimA1, dimA2
      INTGR4 dimB1, dimB2
c
      INTGR4 mpA, nqA          ! rows/cols of A held by the processor
      INTGR4 mpB, nqB          ! rows/cols of B held by the processor
c
      integer me
      INTGR4 lda, ldb
      integer elemA,elemB
      INTGR4 numroc
      
      INTGR4 nb                ! block size
      INTGR4 descA(9), descB(9) !descriptor for scalapack
      
c     
      integer lcwork,hcwork, adrcwork
      INTGR4 lcwork4
c
      INTGR4 nn,mq0, np0
      INTGR4 n
      INTGR4 info
      integer info8
      integer dblsize,ldc
      INTGR4 zero4,one4,two4
      parameter(zero4=0,one4=1,two4=2)

      double precision pdlamch,dumm
      INTGR4 iceil
      logical uses_sl_A, uses_sl_B
      integer alen, blen
      integer block_dims_A(2),block_dims_B(2),blocks(2)
      integer gridA(2), gridB(2)
      logical use_direct

      external pdlamch
      integer i,j

c     
c     processor dependent; machine dependent
c     
      if(nb8.eq.0) then
         nb=16
      else
         nb=nb8
      endif
c     
c***  check environment
c     
      me     = ga_nodeid()
c     
c***  check GA info for input arrays
c     
      call ga_check_handle(g_a, 'ga_pdsyev: a')
      call ga_check_handle(g_b, 'ga_pdsyev: b')
c     
      call ga_inquire(g_a, typeA, dimA18, dimA28)
      call ga_inquire(g_b, typeB, dimB18, dimB28)
      dima1=dima18
      dima2=dima28
      dimb1=dimb18
      dimb2=dimb28
      n=dima1
      if(nb.lt.1) nb=1
      if(dimA1.ne.dima2) call ga_error(
     %     'ga_pdsyev: matrix A not square ',0)
      if(dimb1.ne.dimb2) call ga_error(
     %     'ga_pdsyev: matrix B not square ',0)
      if(dimb1.ne.n) call ga_error(
     %     'ga_pdsyev: size matrix A and B differ ',0)
      use_direct = .true.
      uses_sl_A = ga_uses_proc_grid(g_a)
      uses_sl_B = .false.
      if (uses_sl_A) then
        uses_sl_B = ga_uses_proc_grid(g_b)
      endif
      if ((.not.uses_sl_A).or.(.not.uses_sl_B)) then
        use_direct = .false.
      endif
      if (uses_sl_A) then
        call ga_get_block_info(g_a,blocks,block_dims_A)
        if (block_dims_A(1).ne.block_dims_A(2)) then
          use_direct = .false.
        endif
      endif
      if (uses_sl_B) then
        call ga_get_block_info(g_b,blocks,block_dims_B)
        if (block_dims_B(1).ne.block_dims_B(2)) then
          use_direct = .false.
        endif
      endif
      if (use_direct.and.(block_dims_A(1).ne.block_dims_B(1))) then
        use_direct = .false.
      endif
      if (use_direct) then
        call ga_get_proc_grid(g_a,gridA)
        call ga_get_proc_grid(g_b,gridB)
        if (gridA(1).ne.gridB(1).or.gridA(2).ne.gridB(2)) then
          use_direct = .false.
        endif
      endif
      if (use_direct) then
        nb = block_dims_A(2)
      endif
c     
c     
c***  initialize SL interface
c     
      if (use_direct) then
        call SLinit4(g_a)
      else
        call SLinit2(n)
      endif
      oactive=iam.lt.maxproc
      if (oactive) then
c     
c***  find SBS format parameters
c     
c     
         mpA = numroc(dimA1, nb, myrow2, zero4, nprow2)
         nqA = numroc(dimA2, nb, mycol2, zero4, npcol2)
c     
         mpB = numroc(dimB1, nb, myrow2, zero4, nprow2)
         nqB = numroc(dimB2, nb, mycol2, zero4, npcol2) 
c     
c     
         lda = max(one4,mpA)
         ldb = max(one4,mpB)
c     
c     
c     let scalapack check for errors
c     
c     should check to see if this is a compute node
c     check to see how this works in the new data server model
c     
c     
         if (use_direct) then
            call nga_access_block_segment(g_a,me,adra,alen)
            call nga_access_block_segment(g_b,me,adrb,blen)
         else
            elemA= mpA*nqA
            status = .true.
            if (elemA.ne.0)status =
     $           ma_push_get(MT_DBL,elemA,'a',ha,adra)
            if (.not.status) 
     &           call ga_error('ga_pdsyev: mem alloc failed A ', -1)
c     
            elemB= mpB*nqB
c     
            if (elemB.ne.0)status =
     $           ma_push_get(MT_DBL,elemB,'b',hb,adrb)
            if (.not.status) 
     &           call ga_error('ga_pdsyev: mem alloc failed B ', -1)
c     
c     
c***  copy g_a to A using the block cyclic scalapack format 
c     
            call ga_to_SL2(g_a, dimA1, dimA2, nb, nb,
     $           dbl_mb(adrA), lda, mpA, nqA)
         endif
      endif
      call ga_sync()
      if(oactive) then
c     
c***  fill SCALAPACK matrix descriptors
c
         call descinit(descA, dimA1, dimA2, nb, nb, zero4, zero4,
     $        islctxt2, lda, info)
         info8=info
         if (info8.ne.0)
     $        call ga_error(' ga_pdsyev: descinit A failed ',-info8)
c     
         call descinit(descB, dimB1, dimB2, nb, nb, zero4, zero4,
     .        islctxt2, ldb, info)
         info8=info
         if (info8.ne.0)
     $        call ga_error(' ga_pdsyev: descinit B failed ',-info8)
c     
         jobz = 'V'
         uplo = 'L'
c     
         nn = max(n, nb, two4)
         np0 = numroc(nn, nb, zero4, zero4, nprow2)
         mq0 = numroc(nn, nb, zero4, zero4, npcol2)
c     
c     get lcwork
c
#if 1
         lcwork4=-1
         call pdsyev(jobz,  uplo,
     1        n, dbl_mb(adrA), one4, one4, descA,
     2        eval,dbl_mb(adrB), one4, one4, 
     D        descB, dumm, lcwork4,
     3        info)
         lcwork=dumm
#else
         lcwork = 5*n +MAX(3*NB,NB*(NP0+1))
c     lwmin
         ldc=numroc(n,nb,myrow2,zero4,nnodes)
         ldc=max(1,ldc)
         lcwork=max(lcwork,5*n+n*ldc+1+max(2*n-2,8192))
#endif
c     
         if (lcwork.ne.0)
     $        status = ma_push_get(MT_DBL, lcwork ,
     $        'cwork',hcwork,adrcwork)
         if (.not.status) 
     &      call ga_error('ga_pdsyev: mem alloc failed cwork  ', -1)
c     
c     
         lcwork4=lcwork
         call pdsyev(jobz,  uplo,
     1        n, dbl_mb(adrA), one4, one4, descA,
     2        eval,dbl_mb(adrB), one4, one4, 
     D        descB, dbl_mb(adrcwork), lcwork4,
     3        info)
c     

         if ( info .ne. 0 ) then
            if ( info .gt. 0 ) then
               call ga_error(' ga_pdsyev: argument is illegal ', info)
            else
               call ga_error(' ga_pdsyev: eigenvecs not converged ',
     $              info)
            endif
         endif
c     
c***  copy solution matrix back to g_c
c     
         if (.not.use_direct) then
           call ga_from_SL2(g_b, dimA1, dimB2,
     $        nb, nb, dbl_mb(adrB),
     &        ldb, mpb, nqB)
         endif
c     
c***  deallocate work/SL arrays
c     
         if ( lcwork .ne. 0 ) status = ma_pop_stack(hcwork)
         if (.not.use_direct) then
           if ( elemB .ne. 0 ) status = ma_pop_stack(hb)
           if ( elemA .ne. 0 ) status = ma_pop_stack(ha)
        endif
      endif
c     
      call ga_sync()
      if(maxproc.lt.nnodes.or.dima1.le.nb) then
c     broadcast evals
         dblsize=ma_sizeof(MT_DBL,1,MT_BYTE)*dima1
         call ga_brdcst(1688,eval,dblsize,0)
      endif
      if (use_direct) then
        call SLexit4
      else
        call SLexit2
      endif
      return
      end
c******************************************************************
c     
c     ga_pdsyevx
c
c******************************************************************
c     
      subroutine ga_pdsyevx(g_a, g_b, eval, nb8)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "scalapack.fh"
c     
      integer g_a               ! matrix A
      integer g_b               ! matrix B
      integer nb8                ! block size
c
c all eigenvalues are real
c
      double precision eval(*)
c     
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      logical status
c
      character*1 jobz, range, uplo
c
      integer ha, adra          !A
      integer hb, adrb          !B
c
c
      logical oactive           ! true iff this process participates
      integer dimA18, dimA28, typeA
      integer dimB18, dimB28, typeB
      INTGR4 dimA1, dimA2
      INTGR4 dimB1, dimB2
c
      INTGR4 mpA, nqA          ! rows/cols of A held by the processor
      INTGR4 mpB, nqB          ! rows/cols of B held by the processor
c
      integer me
      INTGR4 lda, ldb
      integer elemA,elemB
      INTGR4 numroc
      
      INTGR4 nb                ! block size
      INTGR4 descA(9), descB(9) !descriptor for scalapack
      
c     
      integer ngaps, hgap, adrgaps
      integer iclu, hclustr, adrclustr
      integer if, hfail, adrfail
      integer liwork,hiwork, adriwork
      integer lcwork,hcwork, adrcwork
      INTGR4 lcwork4
      INTGR4 liwork4
c
      INTGR4 nn,mq0, np0
      double precision vl, vu, abstol, orfac
      INTGR4 il, iu
      INTGR4 m, nz
      INTGR4 n
      INTGR4 info
      integer info8
      integer dblsize
      INTGR4 zero4,one4,two4,four4
      integer two4n
      parameter(zero4=0,one4=1,two4=2,four4=4)

      double precision pdlamch,dumm
      external pdlamch
      INTGR4 iceil

c     
c     processor dependent; machine dependent
c     
      if(nb8.eq.0) then
         nb=16
      else
         nb=nb8
      endif
c     
c***  check environment
c     
      me     = ga_nodeid()
c     
c***  check GA info for input arrays
c     
      call ga_check_handle(g_a, 'ga_pdsyevx: a')
      call ga_check_handle(g_b, 'ga_pdsyevx: b')
c     
      call ga_inquire(g_a, typeA, dimA18, dimA28)
      call ga_inquire(g_b, typeB, dimB18, dimB28)
      dima1=dima18
      dima2=dima28
      dimb1=dimb18
      dimb2=dimb28
      n=dima1
      if(nb.lt.1) nb=1
      if(dimA1.ne.dima2) call ga_error(
     %     'ga_pdsyevx: matrix A not square ',0)
      if(dimb1.ne.dimb2) call ga_error(
     %     'ga_pdsyevx: matrix B not square ',0)
      if(dimb1.ne.n) call ga_error(
     %     'ga_pdsyevx: size matrix A and B differ ',0)
      

c     
c     
c***  initialize SL interface
c     
      call SLinit2(n)
      oactive=iam.lt.maxproc
      call ga_sync
      if (oactive) then
c     
c***  find SBS format parameters
c     
c     
         mpA = numroc(dimA1, nb, myrow2, zero4, nprow2)
         nqA = numroc(dimA2, nb, mycol2, zero4, npcol2)
c     
         mpB = numroc(dimB1, nb, myrow2, zero4, nprow2)
         nqB = numroc(dimB2, nb, mycol2, zero4, npcol2) 
c     
c     
         lda = max(one4,mpA)
         ldb = max(one4,mpB)
c     
c     
c     let scalapack check for errors
c     
c     should check to see if this is a compute node
c     check to see how this works in the new data server model
c     
c     
         elemA= mpA*nqA
         status = .true.
         if(elemA.ne.0)status =
     $        ma_push_get(MT_DBL,elemA,'a',ha,adra)
         if(.not.status) 
     &        call ga_error('ga_pdsyevx: mem alloc failed A ', -1)
c     
c***  copy g_a to A using the block cyclic scalapack format 
c     
         
         call ga_to_SL2(g_a, dimA1, dimA2, nb, nb,
     $        dbl_mb(adrA), lda, mpA, nqA)
c     
         elemB= mpB*nqB
c     
         if(elemB.ne.0)status =
     $        ma_push_get(MT_DBL,elemB,'b',hb,adrb)
         if(.not.status) 
     &        call ga_error('ga_pdsyevx: mem alloc failed B ', -1)
c     
         ngaps = nprow2*npcol2
         if(ngaps.ne.0)status =
     $        ma_push_get(MT_DBL, ngaps ,'gap',hgap,adrgaps)
         if(.not.status) 
     &        call ga_error('ga_pdsyevx: mem alloc failed gaps  ', -1)
c     
         iclu = 2*nprow2*npcol2
         two4n=two4*n
         iclu = max(two4n,iclu)
         if(iclu.ne.0)status =
     $        ma_push_get(MT_INT, iclu ,'iclustr',hclustr,adrclustr)
         if(.not.status) 
     &        call ga_error('ga_pdsyevx: mem alloc failed iclustr',-1)
c     
         if = n
         if(if.ne.0) then
            status =
     $           ma_push_get(MT_INT, if ,'ifail',hfail,adrfail)
         endif
         if(.not.status)
     &        call ga_error('ga_pdsyevx: mem alloc failed ifail  ', -1)
      endif
      call ga_sync()
      if(oactive) then
         
c     
c***  fill SCALAPACK matrix descriptors
c     
         call descinit(descA, dimA1, dimA2, nb, nb, zero4, zero4,
     $        islctxt2, lda, info)
         info8=info
         if(info8.ne.0) call ga_error(' ga_pdsyevx: descinit A failed ',
     .        -info8)
c     
         call descinit(descB, dimB1, dimB2, nb, nb, zero4, zero4,
     .        islctxt2, ldb, info)
         info8=info
         if(info8.ne.0) call ga_error(' ga_pdsyevx: descinit B failed ',
     .        -info8)
c     
         jobz = 'V'
         range = 'A'
         range = 'I'
         uplo = 'L'
         vl = 0.d0
         vu = 0.d0
         il = 0
         iu = 0
         il = 1
!         iu = neigen
         iu = n
         nz = 0
c     
c     ability to deal with orthonormality ; let's just
c     have the regular scalapack stuff for the moment
c     
         liwork = 6*max(n, nprow2*npcol2+one4, four4) 
         liwork=liwork
         if(liwork.ne.0)status =
     $        ma_push_get(MT_INT, liwork ,'iwork',hiwork,adriwork)
         if(.not.status) 
     &        call ga_error('ga_pdsyevx: mem alloc failed iwork  ', -1)
c     
         nn = max(n, nb, two4)
         np0 = numroc(nn, nb, zero4, zero4, nprow2)
         mq0 = numroc(nn, nb, zero4, zero4, npcol2)
c     
c     
         orfac = 1.d-3
c     
c     
#if 0
         lcwork4=-1

         call pdsyevx(jobz, range, uplo,
     1        n, dbl_mb(adrA), one4, one4, descA,vl, 
     2        vu, il, iu, abstol, m, nz, eval, orfac, dbl_mb(adrB), 
     4        one4, one4, descB, dumm, lcwork4,
     3        int_mb(adriwork), liwork4, int_mb(adrfail),
     4        int_mb(adrclustr), dbl_mb(adrgaps), info)
         lcwork=dumm
#else
         lcwork = 5*n +MAX(5*NN,(NP0*MQ0 + 2*nb*nb))+
     +          ICEIL( N, NPROW2*NPCOL2)*NN+1
#endif

c     
         if(lcwork.ne.0)
     $        status = ma_push_get(MT_DBL, lcwork ,
     $        'cwork',hcwork,adrcwork)
         if(.not.status) 
     &        call ga_error('ga_pdsyevx: mem alloc failed cwork  ', -1)
c     
c     
         abstol=pdlamch(islctxt2, 'U')
c     
c     
         liwork4=liwork
         lcwork4=lcwork
         call pdsyevx(jobz, range, uplo,
     1        n, dbl_mb(adrA), one4, one4, descA,vl, 
     2        vu, il, iu, abstol, m, nz, eval, orfac, dbl_mb(adrB), 
     4        one4, one4, descB, dbl_mb(adrcwork), lcwork4,
     3        int_mb(adriwork), liwork4, int_mb(adrfail),
     4        int_mb(adrclustr), dbl_mb(adrgaps), info)
c     

         if (nz .ne. n ) then
            if ( info .ne. 0 ) then
               if ( info .gt. 0 ) then
         call ga_error(' ga_pdsyevx: argument is illegal ', info)
               else
         call ga_error(' ga_pdsyevx: eigenvectors failed to converge ',
     $                 info)
               endif
            endif
         endif
c     
c     
c     
c***  copy solution matrix back to g_c
c     
         call ga_from_SL2(g_b, dimA1, dimB2,
     $        nb, nb, dbl_mb(adrB),
     &        ldb, mpb, nqB)
c     
c     
c     
c***  deallocate work/SL arrays
c     
         if ( lcwork .ne. 0 ) status = ma_pop_stack(hcwork)
         if ( liwork .ne. 0 ) status = ma_pop_stack(hiwork)
         if ( if .ne. 0 ) status = ma_pop_stack(hfail)
         if ( iclu .ne. 0 ) status = ma_pop_stack(hclustr)
         if ( ngaps.ne.0 ) status = ma_pop_stack(hgap)
         if ( elemB .ne. 0 ) status = ma_pop_stack(hb)
         if ( elemA .ne. 0 ) status = ma_pop_stack(ha)
      endif
c     
      call ga_sync()
      if(maxproc.lt.nnodes.or.dima1.le.nb) then
c     broadcast evals
         dblsize=ma_sizeof(MT_DBL,1,MT_BYTE)*dima1
         call ga_brdcst(1688,eval,dblsize,0)
      endif
      return
      end
c******************************************************************
c     
c     ga_pdsygv
c
c     interface into scalapack's Hermitian eigensolver
c     
c     ga interface to scalapack
c     
c     g_a .g_c(*,i) = eval(i).g_s(*,i) g_b(*,i)
c     
c     assume equal size blocks ...
c     
c******************************************************************
      subroutine ga_pdsygv(g_a, g_s, g_b, eval)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "scalapack.fh"
c     
      integer g_a               ! matrix A
      integer g_b               ! matrix B (eigenvectors)
      integer g_s               ! metric
c
c all eigenvalues are real
c
      double precision eval(*)
c     
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      logical status
c
      character*1 jobz, range, uplo
c
      integer ha, adra          !A (input)
      integer hb, adrb          !B (output)
      integer hs, adrs          !S (input)
c
c
      logical oactive           ! true if this process participates
      integer dimA18, dimA28, typeA
      integer dimB18, dimB28, typeB
      integer dims18, dims28, types
      INTGR4 dimA1, dimA2
      INTGR4 dimB1, dimB2
      INTGR4 dims1, dims2
c
      INTGR4 mpA, nqA          ! rows/cols of A held by the processor
      INTGR4 mpB, nqB          ! rows/cols of B held by the processor
      INTGR4 mps, nqs          ! rows/cols of B held by the processor
c
      integer me
      INTGR4 lda, ldb, lds
      integer elemA, elemB, elems
      INTGR4 numroc
      
      INTGR4 nb                ! block size
      INTGR4 descA(9), descB(9) !descriptor for scalapack
      INTGR4 descs(9)
      
c     
      integer ngaps, hgap, adrgaps
      integer iclu, hclustr, adrclustr
      integer if, hfail, adrfail
      integer liwork, hiwork, adriwork
      integer lcwork, hcwork, adrcwork
      INTGR4 lcwork4
      INTGR4 liwork4
c
      INTGR4 nn, mq0, np0
      double precision vl, vu, abstol, orfac
      INTGR4 il, iu
      INTGR4 m, nz
      INTGR4 n
      INTGR4 info
      integer info8,dblsize
      INTGR4 zero4,one4,two4,four4
      parameter(zero4=0,one4=1,two4=2,four4=4)

      INTGR4 iceil
      double precision pdlamch
      integer i
      logical use_direct
      logical uses_sl_A, uses_sl_B, uses_sl_S
      integer alen, blen, slen
      integer block_dims_A(2),block_dims_B(2),block_dims_S(2),blocks(2)
      integer gridA(2), gridB(2), gridS(2)

      external pdlamch
      external iceil
      integer j
!      intrinsic max
c     
c     processor dependent; machine dependent
c     
      data nb/16/
c     
c***  check environment
c     
      me     = ga_nodeid()
c     
c***  check GA info for input arrays
c     
      call ga_check_handle(g_a, 'ga_pdsygv: a')
      call ga_check_handle(g_b, 'ga_pdsygv: b')
      call ga_check_handle(g_s, 'ga_pdsygv: s')
c     
      call ga_inquire(g_a, typeA, dimA18, dimA28)
      call ga_inquire(g_b, typeB, dimB18, dimB28)
      call ga_inquire(g_s, types, dims18, dims28)
      dima1=dima18
      dima2=dima28
      dimb1=dimb18
      dimb2=dimb28
      dims1=dims18
      dims2=dims28
      if(dimA1.ne.dima2) call ga_error(
     %     'ga_pdsygv: matrix A not square ',0)
      n=dima1
      if(nb.lt.1) nb=1
      if(dimb1.ne.dimb2) call ga_error(
     %     'ga_pdsygv: matrix B not square ',0)
      if(dimb1.ne.n) call ga_error(
     %     'ga_pdsygv: size matrix B and A differ ',0)
      if(dims1.ne.dims2) call ga_error(
     %     'ga_pdsygv: matrix S not square ',0)
      if(dims1.ne.n) call ga_error(
     %     'ga_pdsygv: size matrix A and S differ ',0)
c     
      use_direct = .true.
      uses_sl_A = ga_uses_proc_grid(g_a)
      uses_sl_B = .false.
      uses_sl_S = .false.
      if (uses_sl_A) then
        uses_sl_B = ga_uses_proc_grid(g_b)
        uses_sl_S = ga_uses_proc_grid(g_s)
      endif
      if ((.not.uses_sl_A).or.(.not.uses_sl_B)) then
        use_direct = .false.
      endif
      if (uses_sl_A) then
        call ga_get_block_info(g_a,blocks,block_dims_A)
        if (block_dims_A(1).ne.block_dims_A(2)) then
          use_direct = .false.
        endif
      endif
      if (uses_sl_B) then
        call ga_get_block_info(g_b,blocks,block_dims_B)
        if (block_dims_B(1).ne.block_dims_B(2)) then
          use_direct = .false.
        endif
      endif
      if (uses_sl_S) then
        call ga_get_block_info(g_s,blocks,block_dims_S)
        if (block_dims_S(1).ne.block_dims_S(2)) then
          use_direct = .false.
        endif
      endif
      if (use_direct.and.((block_dims_A(1).ne.block_dims_B(1)).or.
     %                   (block_dims_B(1).ne.block_dims_S(1)))) then
        use_direct = .false.
      endif
      if (use_direct) then
        call ga_get_proc_grid(g_a,gridA)
        call ga_get_proc_grid(g_b,gridB)
        call ga_get_proc_grid(g_s,gridS)
        if ((gridA(1).ne.gridB(1).or.gridA(2).ne.gridB(2)).or.
     %      (gridA(1).ne.gridS(1).or.gridA(2).ne.gridS(2))) then
          use_direct = .false.
        endif
      endif
      if (use_direct) then
        nb = block_dims_A(2)
      endif
c     
c***  initialize SL interface
c     
!      call SLinit()
      if (use_direct) then
        call SLinit4(g_a)
      else
        call SLinit2(n)
      endif
      oactive=iam.lt.maxproc
      call ga_sync
      if (oactive) then
c     
c***  find SBS format parameters
c     
c     
         mpA = numroc(dimA1, nb, myrow2, zero4, nprow2)
         nqA = numroc(dimA2, nb, mycol2, zero4, npcol2)
c     
         mpB = numroc(dimB1, nb, myrow2, zero4, nprow2)
         nqB = numroc(dimB2, nb, mycol2, zero4, npcol2) 
c     
         mps = numroc(dims1, nb, myrow2, zero4, nprow2)
         nqs = numroc(dims2, nb, mycol2, zero4, npcol2) 
c     
c     
         lda = max(one4,mpA)
         ldb = max(one4,mpB)
         lds = max(one4,mps)
c     
c     
c     let scalapack check for errors
c     
         if (use_direct) then
           call nga_access_block_segment(g_a, me, adra, alen)
           call nga_access_block_segment(g_b, me, adrb, blen)
           call nga_access_block_segment(g_s, me, adrs, slen)
         else
           elemA= mpA*nqA
           status = .true.
           if(elemA.ne.0)            status =
     $        ma_push_get(MT_DBL,elemA,'a',ha,adra)
           if(.not.status) 
     &        call ga_error('ga_pdsygv: mem alloc failed A ', elema)
c     
           elemB= mpB*nqB
           if(elemB.ne.0)status =
     $        ma_push_get(MT_DBL,elemB,'b',hb,adrb)
           if(.not.status) 
     &        call ga_error('ga_pdsygv: mem alloc failed B ', -1)
           elems= mps*nqs
           if(elems.ne.0)status =
     $        ma_push_get(MT_DBL,elems,'s',hs,adrs)
           if(.not.status) 
     &        call ga_error('ga_pdsygv: mem alloc failed S ', -1)
         endif
c     
         ngaps = nprow2*npcol2
         if(ngaps.ne.0)status =
     $        ma_push_get(MT_DBL, ngaps ,'gap',hgap,adrgaps)
         if(.not.status) 
     &        call ga_error('ga_pdsygv: mem alloc failed gaps  ', -1)
c     
         iclu = 2*nprow2*npcol2
         iclu = 2*n
         if(iclu.ne.0)status =
     $        ma_push_get(MT_INT, iclu ,'iclustr',hclustr,adrclustr)
         if(.not.status) 
     &        call ga_error('ga_pdsygv: mem alloc failed iclustr  ', -1)
c     
         if = n
         if(if.ne.0) then
            status =
     $           ma_push_get(MT_INT, if ,'ifail',hfail,adrfail)
         endif
         if(.not.status)
     &        call ga_error('ga_pdsygv: mem alloc failed ifail  ', -1)
         
c     
c***  copy g_a to A and g_b to B using the block cyclic scalapack format 
c     
         
         if (.not.use_direct) then
           call ga_to_SL2(g_a, dimA1, dimA2, nb, nb,
     $          dbl_mb(adrA), lda, mpA, nqA)
c     
           call ga_to_SL2(g_b, dimB1, dimB2, nb,nb,
     $          dbl_mb(adrB), ldb, mpB, nqB)
           call ga_to_SL2(g_s, dims1, dims2, nb,nb,
     $          dbl_mb(adrs), ldb, mps, nqs)
         endif
c     
c     
c***  fill SCALAPACK matrix descriptors
c     
         call descinit(descA, dimA1, dimA2, nb, nb, zero4, zero4,
     $        islctxt2, lda, info)
         info8=info
         if(info8.ne.0) call ga_error(' ga_pdsygv: descinit A failed ',
     .        -info8)
c     
         call descinit(descB, dimB1, dimB2, nb, nb, zero4, zero4,
     .        islctxt2, ldb, info)
         info8=info
         if(info8.ne.0) call ga_error(' ga_pdsygv: descinit B failed ',
     .        -info8)
         call descinit(descs, dims1, dims2, nb, nb, zero4, zero4,
     .        islctxt2, lds, info)
         info8=info
         if(info8.ne.0) call ga_error(' ga_pdsygv: descinit S failed ',
     .        -info8)
c     
c     call scalapack
c     
c     
         jobz = 'V'
         range = 'A'
         uplo = 'L'
         vl = 0.d0
         vu = 0.d0
         il = 0
         iu = 0
         nz = 0
c     
c     ability to deal with orthonormality ; let's just
c     have the regular scalapack stuff for the moment
c     
         liwork = 6*max(n, nprow2*npcol2+one4, four4) 
         if(liwork.ne.0)status =
     $        ma_push_get(MT_INT, liwork ,'iwork',hiwork,adriwork)
         if(.not.status) 
     &        call ga_error('ga_pdsygv: mem alloc failed iwork  ', -1)
c     
         nn = max(n, nb, two4)
         np0 = numroc(nn, nb, zero4, zero4, nprow2)
         mq0 = numroc(nn, nb, zero4, zero4, npcol2)
c     
c     
         orfac = 1.d-3
c     
c     
         lcwork = 5*n +MAX(5*NN,(NP0*MQ0 + 2*nb*nb))+
     +          ICEIL( N, NPROW2*NPCOL2)*NN+1 

c     
         if(lcwork.ne.0)
     $        status = ma_push_get(MT_DBL, lcwork ,
     $        'cwork',hcwork,adrcwork)
         if(.not.status) 
     &        call ga_error('ga_pdsygv: mem alloc failed cwork  ', -1)
c     
         abstol=pdlamch(islctxt2, 'U')
c     
c     
c     
         liwork4=liwork
         lcwork4=lcwork
         call pdsygvx(one4,jobz, range, uplo,
     1        n, dbl_mb(adrA), one4, one4, 
     2        descA, dbl_mb(adrS), one4, one4,descs, vl,vu, il, iu, 
     3        abstol, m, nz, eval, orfac, dbl_mb(adrB), 
     4        one4, one4, descB, 
     4        dbl_mb(adrcwork), lcwork4,int_mb(adriwork), liwork4, 
     4        int_mb(adrfail), int_mb(adrclustr), 
     5        dbl_mb(adrgaps), info)
c     
c     
c
         
c     
         if (nz .ne. n ) then
            if ( info .ne. 0 ) then
               if ( info .gt. 0 ) then
         call ga_error(' ga_pdsygv: argument is illegal ', info)
               else
         call ga_error(' ga_pdsygv: eigenvectors failed to converge ',
     $                 info)
               endif
            endif
         endif
c     
c     
c     
c***  copy solution matrix back to g_c
c     
c     
         if (.not.use_direct) then
           call ga_from_SL2(g_b, dimA1, dimB2,
     $          nb, nb, dbl_mb(adrB),
     &          ldb, mpb, nqB)
         endif
c     
     
c     
c***  deallocate work/SL arrays
c     
         if ( lcwork .ne. 0 ) status = ma_pop_stack(hcwork)
         if ( liwork .ne. 0 ) status = ma_pop_stack(hiwork)
         if ( if .ne. 0 ) status = ma_pop_stack(hfail)
         if ( iclu .ne. 0 ) status = ma_pop_stack(hclustr)
         if ( ngaps.ne.0 ) status = ma_pop_stack(hgap)
         if (.not.use_direct) then
           if ( elemS .ne. 0 ) status = ma_pop_stack(hs)
           if ( elemB .ne. 0 ) status = ma_pop_stack(hb)
           if ( elemA .ne. 0 ) status = ma_pop_stack(ha)
         endif
      endif
c     
      call ga_sync()
      if(maxproc.lt.nnodes.or.dima1.le.nb) then
c     broadcast evals
         dblsize=ma_sizeof(MT_DBL,1,MT_BYTE)*dima1
         call ga_brdcst(1688,eval,dblsize,0)
      endif
      if (use_direct) then
        call SLexit4
      else
        call SLexit2
      endif
      return
      end
c
c     interface to pdsyevd (divide&conquer)
c
      subroutine ga_pdsyevd(g_a, g_b, eval, nb8)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "scalapack.fh"
c     
      integer g_a               ! matrix A
      integer g_b               ! matrix B
      integer nb8
c
c all eigenvalues are real
c
      double precision eval(*)
c     
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      logical status
c
      character*1 jobz, uplo
c
      integer ha, adra          !A
      integer hb, adrb          !B
c
c
      logical oactive           ! true iff this process participates
      integer dimA18, dimA28, typeA
      integer dimB18, dimB28, typeB
      INTGR4 dimA1, dimA2
      INTGR4 dimB1, dimB2
c
      INTGR4 mpA, nqA          ! rows/cols of A held by the processor
      INTGR4 mpB, nqB          ! rows/cols of B held by the processor
c
      integer me
      INTGR4 lda, ldb
      integer elemA,elemB
      INTGR4 numroc
      
      INTGR4 nb                ! block size
      INTGR4 descA(9), descB(9) !descriptor for scalapack
      
c     
      integer liwork,hiwork, adriwork
      integer lcwork,hcwork, adrcwork
      INTGR4 liwork4,lcwork4
c
      INTGR4 nn,nq, np,trilwmin,iarow,iacol
      INTGR4 n
      INTGR4 info,one4,zero4,two4
      parameter(zero4=0,one4=1,two4=2)
      integer info8,dblsize

      double precision pdlamch,dum
      INTGR4 iceil,indxg2p
      logical uses_sl_A, uses_sl_B
      integer alen, blen
      integer block_dims_A(2),block_dims_B(2),blocks(2)
      integer gridA(2), gridB(2)
      logical use_direct

      external pdlamch,iceil,indxg2p
c     
c     processor dependent; machine dependent
c     
      if(nb8.eq.0) then
         nb=16
      else
         nb=nb8
      endif
c     
c***  check environment
c     
      me     = ga_nodeid()
c     
c***  check GA info for input arrays
c     
      call ga_check_handle(g_a, 'ga_pdsyevd: a')
      call ga_check_handle(g_b, 'ga_pdsyevd: b')
c     
      call ga_inquire(g_a, typeA, dimA18, dimA28)
      call ga_inquire(g_b, typeB, dimB18, dimB28)
      dima1=dima18
      dima2=dima28
      dimb1=dimb18
      dimb2=dimb28
      n=dima1
      if(nb.lt.1) nb=1

      if(dimA1.ne.dima2) call ga_error(
     %     'ga_pdsyevd: matrix A not square ',0)
      if(dimb1.ne.dimb2) call ga_error(
     %     'ga_pdsyevd: matrix B not square ',0)
      if(dimb1.ne.n) call ga_error(
     %     'ga_pdsyevd: size matrix B and B differ ',0)
      use_direct = .true.
      uses_sl_A = ga_uses_proc_grid(g_a)
      uses_sl_B = .false.
      if (uses_sl_A) then
        uses_sl_B = ga_uses_proc_grid(g_b)
      endif
      if ((.not.uses_sl_A).or.(.not.uses_sl_B)) then
        use_direct = .false.
      endif
      if (uses_sl_A) then
        call ga_get_block_info(g_a,blocks,block_dims_A)
        if (block_dims_A(1).ne.block_dims_A(2)) then
          use_direct = .false.
        endif
      endif
      if (uses_sl_B) then
        call ga_get_block_info(g_b,blocks,block_dims_B)
        if (block_dims_B(1).ne.block_dims_B(2)) then
          use_direct = .false.
        endif
      endif
      if (use_direct.and.(block_dims_A(1).ne.block_dims_B(1))) then
        use_direct = .false.
      endif
      if (use_direct) then
        call ga_get_proc_grid(g_a,gridA)
        call ga_get_proc_grid(g_b,gridB)
        if (gridA(1).ne.gridB(1).or.gridA(2).ne.gridB(2)) then
          use_direct = .false.
        endif
      endif
      if (use_direct) then
        nb = block_dims_A(2)
      endif
c     
c     
c***  initialize SL interface
c     
!      call SLinit()
      if (use_direct) then
        call SLinit4(g_a)
      else
        call SLinit2(n)
      endif
      oactive=iam.lt.maxproc
      call ga_sync
      if (oactive) then
c     
c***  find SBS format parameters
c     
c     
         mpA = numroc(dimA1, nb, myrow2, zero4, nprow2)
         nqA = numroc(dimA2, nb, mycol2, zero4, npcol2)
c     
         mpB = numroc(dimB1, nb, myrow2, zero4, nprow2)
         nqB = numroc(dimB2, nb, mycol2, zero4, npcol2) 
c     
c     
         lda = max(one4,mpA)
         ldb = max(one4,mpB)
c     
c     
c     let scalapack check for errors
c     
c     should check to see if this is a compute node
c     check to see how this works in the new data server model
c     
c     
c     
         if (use_direct) then
            call nga_access_block_segment(g_a,me,adra,alen)
            call nga_access_block_segment(g_b,me,adrb,blen)
         else
            elemA= mpA*nqA
            status = .true.
            if(elemA.ne.0)status =
     $           ma_push_get(MT_DBL,elemA,'a',ha,adra)
            if(.not.status) 
     &           call ga_error('ga_pdsyevd: mem alloc failed A ', -1)
c     
c***  copy g_a to A and g_b to B using the block cyclic scalapack format 
c     
         
            call ga_to_SL2(g_a, dimA1, dimA2, nb, nb,
     $           dbl_mb(adrA), lda, mpA, nqA)
c     
            elemB= mpB*nqB
c     
            if(elemB.ne.0)status =
     $           ma_push_get(MT_DBL,elemB,'b',hb,adrb)
            if(.not.status) 
     &           call ga_error('ga_pdsyevd: mem alloc failed B ', -1)
         endif
      endif
      call ga_sync()
      if(oactive) then
c     
c     
c     
c***  fill SCALAPACK matrix descriptors
c     
      call descinit(descA, dimA1, dimA2, nb, nb, zero4, zero4,
     $     islctxt2, lda, info)
      info8=info
      if(info8.ne.0) call ga_error(' ga_pdsyevd: descinit A failed ',
     .     -info8)
c     
      call descinit(descB, dimB1, dimB2, nb, nb, zero4, zero4,
     .     islctxt2, ldb, info)
      info8=info
      if(info8.ne.0) call ga_error(' ga_pdsyevd: descinit B failed ',
     .     -info8)
c     
      jobz = 'V'
      uplo = 'L'
c     
c     ability to deal with orthonormality ; let's just
c     have the regular scalapack stuff for the moment
c     
      liwork = 7*n + 8*npcol2 +2
      if(liwork.ne.0)status =
     $     ma_push_get(MT_INT, liwork ,'iwork',hiwork,adriwork)
      if(.not.status) 
     &     call ga_error('ga_pdsyevd: mem alloc failed iwork  ', -1)
c     
      nn = max(n, nb, two4)
      IAROW = INDXG2P(one4, nb, MYROW2, DESCA( 7 ), NPROW2 )
      IACOL = INDXG2P(one4, nb, MYCOL2, DESCA( 8 ), NPCOL2 )
            
      np = numroc(n, nb, myrow2, iarow, nprow2)
      nq = numroc(n, nb, mycol2, iacol, npcol2)
c     
c     
c     
      liwork4=liwork
#if 0
      lcwork4=-1
          call pdsyevd(jobz, uplo,
     1         n, dbl_mb(adrA), one4, one4, descA,
     1         eval, dbl_mb(adrB), one4, one4, 
     2         descB, dbl_mb(adrcwork), lcwork4,
     2         dum, liwork4, info)
          lcwork=dum
#else

          lcwork = MAX( 1+6*N+2*NP*NQ, 
     =         3*N + MAX( NB*( NP+1 ), 3*NB ))+ 2*N
          lcwork=max(lcwork,16384)
#endif

c     
      if(lcwork.ne.0)
     $     status = ma_push_get(MT_DBL, lcwork ,
     $     'cwork',hcwork,adrcwork)
      if(.not.status) 
     &     call ga_error('ga_pdsyevd: mem alloc failed cwork  ', -1)
c     
c     
c     
      lcwork4=lcwork
      call pdsyevd(jobz, uplo,
     1     n, dbl_mb(adrA), one4, one4, descA,
     1     eval, dbl_mb(adrB), one4, one4, 
     2     descB, dbl_mb(adrcwork), lcwork4,
     2     int_mb(adriwork), liwork4, info)
c     
      if ( info .ne. 0 ) then
         if ( info .gt. 0 ) then
            call ga_error(' ga_pdsyevd: argument is illegal ', info)
         else
            call ga_error(
     $           ' ga_pdsyevd: eigenvectors failed to converge ',
     $           info)
         endif
      endif
c     
c     
c     
c***  copy solution matrix back to g_b
c     
      if (.not.use_direct) then
        call ga_from_SL2(g_b, dimA1, dimB2,
     $       nb, nb, dbl_mb(adrB),
     &       ldb, mpb, nqB)
      endif
c     
c     
c     
c***  deallocate work/SL arrays
c     
      if ( lcwork .ne. 0 ) status = ma_pop_stack(hcwork)
      if ( liwork .ne. 0 ) status = ma_pop_stack(hiwork)
      if (.not.use_direct) then
        if ( elemB .ne. 0 ) status = ma_pop_stack(hb)
        if ( elemA .ne. 0 ) status = ma_pop_stack(ha)
      endif
c     
      endif
      call ga_sync()
      if(maxproc.lt.nnodes.or.dima1.le.nb) then
c     broadcast evals
         dblsize=ma_sizeof(MT_DBL,1,MT_BYTE)*dima1
         call ga_brdcst(1688,eval,dblsize,0)
      endif
      return
      end
#if HAVE_PDSYEVR
      subroutine ga_pdsyevr(g_a, g_b, eval, nb8, mout)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "scalapack.fh"
c     
      integer g_a               ! matrix A
      integer g_b               ! matrix B
      integer nb8
      integer mout ! no. of evals/evecs to compute
c
c all eigenvalues are real
c
      double precision eval(*)
c     
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      logical status
c
      character*1 jobz, range, uplo
c
      integer ha, adra          !A
      integer hb, adrb          !B
c
c
      logical oactive           ! true if this process participates
      integer dimA18, dimA28, typeA
      integer dimB18, dimB28, typeB
      INTGR4 dimA1, dimA2,mout4
      INTGR4 dimB1, dimB2
c
      INTGR4 mpA, nqA          ! rows/cols of A held by the processor
      INTGR4 mpB, nqB          ! rows/cols of B held by the processor
c
      integer me
      INTGR4 lda, ldb
      integer elemA,elemB
      INTGR4 numroc
      
      INTGR4 nb                ! block size
      INTGR4 descA(9), descB(9) !descriptor for scalapack
      
c     

      integer n8, hfail, adrfail
      integer liwork,hiwork, adriwork
      integer lcwork,hcwork, adrcwork
      INTGR4 liwork4
      INTGR4 lcwork4
c
      INTGR4 nn,mq0, np0
      double precision vl, vu
      INTGR4 il, iu
      INTGR4 m, nz
      INTGR4 n
      INTGR4 info,one4,zero4,two4,four4
      parameter(one4=1,zero4=0,two4=4,four4=4)
      integer info8,dblsize

      double precision pdlamch,dumm
      INTGR4 dumm2
      external pdlamch
      INTGR4 iceil
c     
c     processor dependent; machine dependent
c     
      if(nb8.eq.0) then
         nb=16
      else
         nb=nb8
      endif
c     
c***  check environment
c     
      me     = ga_nodeid()
c     
c***  check GA info for input arrays
c     
      call ga_check_handle(g_a, 'ga_pdsyevr: a')
      call ga_check_handle(g_b, 'ga_pdsyevr: b')
c     
      call ga_inquire(g_a, typeA, dimA18, dimA28)
      call ga_inquire(g_b, typeB, dimB18, dimB28)
      dima1=dima18
      dima2=dima28
      dimb1=dimb18
      dimb2=dimb28
      if(dimA1.ne.dima2) call ga_error(
     %     'ga_pdsyevr: matrix A not square ',0)
      n=dima1
      if(dimb1.ne.dimb2) call ga_error(
     %     'ga_pdsyevr: matrix B not square ',0)
      if(dimb1.ne.n) call ga_error(
     %     'ga_pdsyevr: size matrix B and B differ ',0)
      
      if(mout.le.0.or.mout.gt.n) mout=n

c     
c     
c***  initialize SL interface
c     
      call SLinit2(n)
      oactive=iam.lt.maxproc
      call ga_zero(g_b)
      if (oactive) then
c     
c***  find SBS format parameters
c     
c     
      mpA = numroc(dimA1, nb, myrow2, zero4, nprow2)
      nqA = numroc(dimA2, nb, mycol2, zero4, npcol2)
c     
      mpB = numroc(dimB1, nb, myrow2, zero4, nprow2)
      nqB = numroc(dimB2, nb, mycol2, zero4, npcol2) 
c     
c     
      lda = max(one4,mpA)
      ldb = max(one4,mpB)
c     
c     
c     let scalapack check for errors
c     
c     should check to see if this is a compute node
c     check to see how this works in the new data server model
c     
c     
         elemA= mpA*nqA
         status = .true.
         if(elemA.ne.0)status =
     $        ma_push_get(MT_DBL,elemA,'a',ha,adra)
         if(.not.status) 
     &        call ga_error('ga_pdsyevr: mem alloc failed A ', -1)
c     
c***  copy g_a to A and g_b to B using the block cyclic scalapack format 
c     
         
         call ga_to_SL2(g_a, dimA1, dimA2, nb, nb,
     $        dbl_mb(adrA), lda, mpA, nqA)
         endif
         call ga_sync()
         if(oactive) then
c     
         elemB= mpB*nqB
c     
         if(elemB.ne.0)status =
     $        ma_push_get(MT_DBL,elemB,'b',hb,adrb)
         if(.not.status) 
     &        call ga_error('ga_pdsyevr: mem alloc failed B ', -1)
c     
         n8 = n
         if(n8.ne.0) then
            status =
     $           ma_push_get(MT_INT, n8 ,'ifail',hfail,adrfail)
         endif
         if(.not.status)
     &        call ga_error('ga_pdsyevr: mem alloc failed ifail  ', n)
         
c     
c***  fill SCALAPACK matrix descriptors
c     
         call descinit(descA, dimA1, dimA2, nb, nb, zero4, zero4,
     $        islctxt2, lda, info)
         info8=info
         if(info.ne.0) call ga_error(' ga_pdsyevr: descinit A failed ',
     .        -info8)
c     
         call descinit(descB, dimB1, dimB2, nb, nb, zero4, zero4,
     .        islctxt2, ldb, info)
         info8=info
         if(info.ne.0) call ga_error(' ga_pdsyevr: descinit B failed ',
     .        -info8)
c     
         jobz = 'V'
         uplo = 'L'
         vl = 0.d0
         vu = 0.d0
c     
         nn = max(n, nb, two4)
         np0 = numroc(nn, nb, zero4, zero4, nprow2)
         mq0 = numroc(nn, nb, zero4, zero4, npcol2)
c     
         il=1
         iu=mout
         range='i'
         lcwork4=-1
         liwork4=-1
         call pdsyevr(jobz, range, uplo,
     1        n, dbl_mb(adrA), one4, one4, descA,vl, 
     2        vu, il, iu, m, nz, eval,  dbl_mb(adrB), one4, 
     3        one4, descB, dumm, lcwork4,
     3        dumm2, liwork4, info)

         lcwork=dumm
         liwork=dumm2
c        lcwork = 5*n +MAX(18*NN,(NP0*MQ0 + 2*nb*nb))+
c     +          (2 + ICEIL( N, NPROW2*NPCOL2))*NN+1
c        liwork = 12*max(n, nprow2*npcol2+one4, four4) +2 *n

c     
         if(liwork.ne.0)status =
     $        ma_push_get(MT_INT, liwork ,'iwork',hiwork,adriwork)
         if(.not.status) 
     &        call ga_error('ga_pdsyevr: mem alloc failed iwork  ', -1)
c     
         if(lcwork.ne.0)
     $        status = ma_push_get(MT_DBL, lcwork ,
     $        'cwork',hcwork,adrcwork)
         if(.not.status) 
     &        call ga_error('ga_pdsyevr: mem alloc failed cwork  ', -1)
c     
         liwork4=liwork
         lcwork4=lcwork

         call pdsyevr(jobz, range, uplo,
     1        n, dbl_mb(adrA), one4, one4, descA,vl, 
     2        vu, il, iu, m, nz, eval,  dbl_mb(adrB), one4, 
     3        one4, descB, dbl_mb(adrcwork), lcwork4,
     3        int_mb(adriwork), liwork4, info)
c     
         mout=m
         if(iam.eq.0.and.mout.ne.n) 
     W        write(0,*) ' computed eval ',mout,' out of ',
     A        dima1


         if (nz .ne. n ) then
            if ( info .ne. 0 ) then
               if ( info .gt. 0 ) then
         call ga_error(' ga_pdsyevr: argument is illegal ', info)
               else
         call ga_error(' ga_pdsyevr: eigenvectors failed to converge ',
     $                 info)
               endif
            endif
         endif
c     
c     
c     
c***  copy solution matrix back to g_b
c     
         mout4=mout
         call ga_from_SL2(g_b, dimA1, mout4,
     $        nb, nb, dbl_mb(adrB),
     &        ldb, mpb, nqB)
c     
c     
c     
c***  deallocate work/SL arrays
c     
         if ( lcwork .ne. 0 ) status = ma_pop_stack(hcwork)
         if ( liwork .ne. 0 ) status = ma_pop_stack(hiwork)
         if ( n8 .ne. 0 ) status = ma_pop_stack(hfail)
         if ( elemB .ne. 0 ) status = ma_pop_stack(hb)
         if ( elemA .ne. 0 ) status = ma_pop_stack(ha)
c     
      endif
      call ga_sync()
      if(maxproc.lt.nnodes.or.dima1.le.nb) then
c     broadcast evals
         dblsize=ma_sizeof(MT_DBL,1,MT_BYTE)*dima1
         call ga_brdcst(1688,eval,dblsize,0)
      endif
      
      return
      end
#endif
#ifdef HAVE_ELPA
      subroutine ga_evp_real(g_a, g_b, eval, nb8, mout)
      use ELPA1
      implicit none
#include "mpif.h"
#include "mafdecls.fh"
#include "global.fh"
#include "scalapack.fh"
c     
      integer g_a               ! matrix A
      integer g_b               ! matrix B
      integer nb8
      integer mout ! no. of evals/evecs to compute
c
c all eigenvalues are real
c
      double precision eval(*)
c     
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      logical status
c
      character*1 jobz, range, uplo
c
      integer ha, adra          !A
      integer hb, adrb          !B
c
c
      logical oactive           ! true iff this process participates
      integer dimA18, dimA28, typeA
      integer dimB18, dimB28, typeB
      INTGR4 dimA1, dimA2,mout4
      INTGR4 dimB1, dimB2
c
      INTGR4 mpA, nqA          ! rows/cols of A held by the processor
      INTGR4 mpB, nqB          ! rows/cols of B held by the processor
c
      integer me
      INTGR4 lda, ldb
      integer elemA,elemB
      INTGR4 numroc
      
      INTGR4 nb                ! block size
      INTGR4 descA(9), descB(9) !descriptor for scalapack
      INTGR4 elpa_comm_rows, elpa_comm_cols
c     
      integer n8
c
      INTGR4 m
      INTGR4 n
      INTGR4 info,one4,zero4,two4,four4
      INTGR4 mpierr,elpa_comm,elpa_group
      integer mxproc,i
      parameter(mxproc=100000)
      INTGR4 group_world, group_members(mxproc)

      parameter(one4=1,zero4=0,two4=4,four4=4)
      integer info8,dblsize
c     
c     processor dependent; machine dependent
c     
      if(maxproc.gt.mxproc) call ga_error(
     A     'ga-evp_real: increase mxproc to ',maxproc)
      if(nb8.eq.0) then
         nb=16
      else
         nb=nb8
      endif
c     
c***  check environment
c     
      me     = ga_nodeid()
c     
c***  check GA info for input arrays
c     
      call ga_check_handle(g_a, 'ga_evp_real: a')
      call ga_check_handle(g_b, 'ga_evp_real: b')
c     
      call ga_inquire(g_a, typeA, dimA18, dimA28)
      call ga_inquire(g_b, typeB, dimB18, dimB28)
      dima1=dima18
      dima2=dima28
      dimb1=dimb18
      dimb2=dimb28
      if(dimA1.ne.dima2) call ga_error(
     %     'ga_evp_real: matrix A not square ',0)
      n=dima1
      if(dimb1.ne.dimb2) call ga_error(
     %     'ga_evp_real: matrix B not square ',0)
      if(dimb1.ne.n) call ga_error(
     %     'ga_evp_real: size matrix B and B differ ',0)
      
      if(mout.le.0.or.mout.gt.n) mout=n
      mout4=mout
c     
c     
c***  initialize SL interface
c     
      call SLinit2(n)
      oactive=iam.lt.maxproc
      call ga_zero(g_b)
celpa
      call MPI_Comm_group(MPI_COMM_WORLD, group_world, mpierr)
      do i=1,maxproc
         group_members(i)=i-1
      enddo
      call MPI_Group_incl(group_world, maxproc, group_members, 
     C     elpa_group, mpierr)
      if(mpierr.ne.0) call ga_error(
     A     'ga-evp_real: mpigroupincl failed ',mpierr)
c          write(6,'(I3,A,3I12)') ga_nodeid(), 'comms ',elpa_group,
c     W         group_world,MPI_COMM_WORLD
      call mpi_comm_create(MPI_COMM_WORLD, elpa_group, elpa_comm,
     W mpierr)
      call ga_sync()
      if (oactive) then
c     
c***  find SBS format parameters
c     
c     
      mpA = numroc(dimA1, nb, myrow2, zero4, nprow2)
      nqA = numroc(dimA2, nb, mycol2, zero4, npcol2)
c     
      mpB = numroc(dimB1, nb, myrow2, zero4, nprow2)
      nqB = numroc(dimB2, nb, mycol2, zero4, npcol2) 
c     
c     
      lda = max(one4,mpA)
      ldb = max(one4,mpB)
c     
c     
c     let scalapack check for errors
c     
c     should check to see if this is a compute node
c     check to see how this works in the new data server model
c     
c     
         elemA= mpA*nqA
         status = .true.
         if(elemA.ne.0)status =
     $        ma_push_get(MT_DBL,elemA,'a',ha,adra)
         if(.not.status) 
     &        call ga_error('ga_evp_real: mem alloc failed A ', -1)
c     
         elemB= mpB*nqB
c         elemb=max(mpb,nqb)**2
c     
         if(elemB.ne.0)status =
     $        ma_push_get(MT_DBL,elemB,'b',hb,adrb)
         if(.not.status) 
     &        call ga_error('ga_evp_real: mem alloc failed B ', -1)
c     
c     
         n8 = n
         
c     
c***  copy g_a to A and g_b to B using the block cyclic scalapack format 
c     
         call ga_to_SL2(g_a, dimA1, dimA2, nb, nb,
     $        dbl_mb(adrA), lda, mpA, nqA)
c     
c***  fill SCALAPACK matrix descriptors
c     
         call descinit(descA, dimA1, dimA2, nb, nb, zero4, zero4,
     $        islctxt2, lda, info)
         info8=info
         if(info.ne.0) call ga_error(' ga_evp_real: descinit A failed ',
     .        -info8)
c     
         call descinit(descB, dimB1, dimB2, nb, nb, zero4, zero4,
     .        islctxt2, ldb, info)
         info8=info
         if(info.ne.0) call ga_error(' ga_evp_real: descinit B failed ',
     .        -info8)
c     
c     

          call get_elpa_row_col_comms(elpa_comm,
     V        myrow2, mycol2, 
     V        elpa_comm_rows, elpa_comm_cols)
 

         call solve_evp_real(n, mout4, 
     A        dbl_mb(adra), lda, 
     E        eval, dbl_mb(adrb), ldb, nb, 
     A        elpa_comm_rows, elpa_comm_cols)

         if(iam.eq.0.and.mout4.ne.n) 
     W        write(6,*) ' computed eval ',mout,' out of ',
     A        dima1


c     
c     
c     
c***  copy solution matrix back to g_b
c     
         call ga_from_SL2(g_b, dimA1, mout4,
     $        nb, nb, dbl_mb(adrB),
     &        ldb, mpb, nqB)

c     
c     
c     
c***  deallocate work/SL arrays
c     
         if ( elemB .ne. 0 ) status = ma_pop_stack(hb)

         if ( elemA .ne. 0 ) status = ma_pop_stack(ha)
           call MPI_Group_free( elpa_group, mpierr)
         if(mpierr.ne.0) call ga_error(
     A     'ga-evp_real: mpigroupfree failed ',mpierr)
           call MPI_comm_free( elpa_comm, mpierr)
         if(mpierr.ne.0) call ga_error(
     A     'ga-evp_real: mpicommfree failed ',mpierr)
c     
      endif
      if(maxproc.lt.nnodes.or.dima1.le.nb) then
c     broadcast evals
         dblsize=ma_sizeof(MT_DBL,1,MT_BYTE)*dima1
         call ga_brdcst(1688,eval,dblsize,0)
      endif
      
      return
      end
#endif
      subroutine ga_pzheevd(g_a, g_b, eval, nb8)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "scalapack.fh"
c     
      integer g_a(2)               ! matrix A
      integer g_b(2)               ! matrix B
      integer nb8
c
c  real eigenvalues 
c
      double precision eval(*)
c     
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      logical status
c
      character*1 jobz, uplo
c
      integer ha, adra          !A
      integer hb, adrb          !B
c
c
      logical oactive           ! true iff this process participates
      integer dimA18, dimA28, typeA
      integer dimB18, dimB28, typeB
      INTGR4 dimA1, dimA2,mout4
      INTGR4 dimB1, dimB2
c
      INTGR4 mpA, nqA          ! rows/cols of A held by the processor
      INTGR4 mpB, nqB          ! rows/cols of B held by the processor
c
      integer me
      INTGR4 lda, ldb
      integer elemA,elemB
      INTGR4 numroc
      
      INTGR4 nb                ! block size
      INTGR4 descA(9), descB(9) !descriptor for scalapack
      
c     
      integer n8, hfail, adrfail
      integer liwork,hiwork, adriwork
      integer lwork,hwork, adrwork
      integer lrwork,hrwork, adrrwork
      INTGR4 liwork4
      INTGR4 lwork4
      INTGR4 lrwork4
c
      INTGR4 m, nz
      INTGR4 n
      INTGR4 info,one4,zero4,two4,four4
      parameter(one4=1,zero4=0,two4=4,four4=4)
      integer info8,dblsize

      double precision pdlamch
      complex*16 dumm
      double precision dumm1
      INTGR4 dumm2
      external pdlamch
      INTGR4 iceil
c     
c     processor dependent; machine dependent
c     
      if(nb8.eq.0) then
         nb=16
      else
         nb=nb8
      endif
c     
c***  check environment
c     
      me     = ga_nodeid()
c     
c***  check GA info for input arrays
c     
      call ga_check_handle(g_a(1), 'ga_pzheevd: a')
      call ga_check_handle(g_b(1), 'ga_pzheevd: b')
c     
      call ga_inquire(g_a(1), typeA, dimA18, dimA28)
      call ga_inquire(g_b(1), typeB, dimB18, dimB28)
      dima1=dima18
      dima2=dima28
      dimb1=dimb18
      dimb2=dimb28
      if(dimA1.ne.dima2) call ga_error(
     %     'ga_pzheevd: matrix A not square ',0)
      n=dima1
c      nb=min(nb,n/2)
      if(dimb1.ne.dimb2) call ga_error(
     %     'ga_pzheevd: matrix B not square ',0)
      if(dimb1.ne.n) call ga_error(
     %     'ga_pzheevd: size matrix B and B differ ',0)
      
c     
c***  initialize SL interface
c     
      call SLinit2(n)
      oactive=iam.lt.maxproc
      call ga_sync
      if (oactive) then
c     
c***  find SBS format parameters
c     
c     
      mpA = numroc(dimA1, nb, myrow2, zero4, nprow2)
      nqA = numroc(dimA2, nb, mycol2, zero4, npcol2)
c     
      mpB = numroc(dimB1, nb, myrow2, zero4, nprow2)
      nqB = numroc(dimB2, nb, mycol2, zero4, npcol2) 
c     
c     
      lda = max(one4,mpA)
      ldb = max(one4,mpB)
c     
c     
c     let scalapack check for errors
c     
c     should check to see if this is a compute node
c     check to see how this works in the new data server model
c     
c     
         elemA= mpA*nqA
         status = .true.
         if(elemA.ne.0)status =
     $        ma_push_get(MT_dcpL,elemA,'a',ha,adra)
         if(.not.status) 
     &        call ga_error('ga_pzheevd: mem alloc failed A ', -1)
         
c     
c***  copy g_a to A and g_b to B using the block cyclic scalapack format 
c     
         
         call ga_to_SL2cpl(g_a, dimA1, dimA2, nb, nb,
     $        dcpl_mb(adrA), lda, mpA, nqA)
         endif
         call ga_sync()
         if(oactive) then
c     
         elemB= mpB*nqB
c     
         if(elemB.ne.0)status =
     $        ma_push_get(MT_dcpL,elemB,'b',hb,adrb)
         if(.not.status) 
     &        call ga_error('ga_pzheevd: mem alloc failed B ', -1)
c     
         n8 = n
         if(n8.ne.0) then
            status =
     $           ma_push_get(MT_INT, n8 ,'ifail',hfail,adrfail)
         endif
         if(.not.status)
     &        call ga_error('ga_pzheevd: mem alloc failed ifail  ', n)
c     
c***  fill SCALAPACK matrix descriptors
c     
         call descinit(descA, dimA1, dimA2, nb, nb, zero4, zero4,
     $        islctxt2, lda, info)
         info8=info
         if(info.ne.0) call ga_error(' ga_pzheevd: descinit A failed ',
     .        -info8)
c     
         call descinit(descB, dimB1, dimB2, nb, nb, zero4, zero4,
     .        islctxt2, ldb, info)
         info8=info
         if(info.ne.0) call ga_error(' ga_pzheevd: descinit B failed ',
     .        -info8)
c     
         jobz = 'V'
         uplo = 'L'
c     
c     
         lwork4=-1
         lrwork4=-1
         liwork4=-1
         call pzheevd(jobz,  uplo, n, dcpl_mb(adrA), one4, one4, descA,
     2        eval,  dcpl_mb(adrB), one4, one4, 
     3        descB, dumm, lwork4, dumm1,lrwork4, dumm2, 
     4        liwork4, info)

         lwork=dumm
c         lwork=max(lwork,n)
c         lwork=2*max(lwork,mpa*(nb+nqa))
c         write(0,*) ga_nodeid(), ' lw orig new ',dumm,lwork
         lrwork=dumm1
         liwork=dumm2

c     
         if(liwork.ne.0)status =
     $        ma_push_get(MT_INT, liwork ,'iwork',hiwork,adriwork)
         if(.not.status) 
     &        call ga_error('ga_pzheevd: mem alloc failed iwork  ', -1)
c     
         if(lwork.ne.0)
     $        status = ma_push_get(MT_DCPL, lwork ,
     $        'work',hwork,adrwork)
         if(.not.status) 
     &        call ga_error('ga_pzheevd: mem alloc failed work  ', -1)
         if(lrwork.ne.0)
     $        status = ma_push_get(MT_DBL, lrwork ,
     $        'rwork',hrwork,adrrwork)
         if(.not.status) 
     &        call ga_error('ga_pzheevd: mem alloc failed rwork  ', -1)
c     
         liwork4=liwork
         lwork4=lwork
         lrwork4=lrwork

         call pzheevd(jobz,  uplo, n, dcpl_mb(adrA),one4,one4, descA,
     2        eval,  dcpl_mb(adrB), one4, one4, 
     3        descB, dcpl_mb(adrwork), lwork4, dbl_mb(adrrwork),
     4        lrwork4, int_mb(adriwork), liwork4, info)


         if (nz .ne. n ) then
            if ( info .ne. 0 ) then
               if ( info .gt. 0 ) then
         call ga_error(' ga_pzheevd: argument is illegal ', info)
               else
         call ga_error(' ga_pzheevd: eigenvectors failed to converge ',
     $                 info)
               endif
            endif
         endif
c     
c     
c     
c***  copy solution matrix back to g_b
c     

         call ga_from_SL2cpl(g_b, dimA1, dima1,
     $        nb, nb, dcpl_mb(adrB),
     &        ldb, mpb, nqB)
c     
c     
c     
c***  deallocate work/SL arrays
c     
         if ( lrwork .ne. 0 ) status = ma_pop_stack(hrwork)
         if ( lwork .ne. 0 ) status = ma_pop_stack(hwork)
         if ( liwork .ne. 0 ) status = ma_pop_stack(hiwork)
         if ( n8 .ne. 0 ) status = ma_pop_stack(hfail)
         if ( elemB .ne. 0 ) status = ma_pop_stack(hb)
         if ( elemA .ne. 0 ) status = ma_pop_stack(ha)
c     
      endif
      call ga_sync()
      if(maxproc.lt.nnodes.or.dima1.le.nb) then
c     broadcast evals
         dblsize=ma_sizeof(MT_DBL,1,MT_BYTE)*dima1
         call ga_brdcst(1688,eval,dblsize,0)
      endif
      
      return
      end
#if HAVE_PDSYEVR
      subroutine ga_pzheevr(g_a, g_b, eval, nb8, mout)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "scalapack.fh"
c     
      integer g_a(2)               ! matrix A
      integer g_b(2)               ! matrix B
      integer nb8
      integer mout ! no. of evals/evecs to compute
c
c  real eigenvalues 
c
      double precision eval(*)
c     
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      logical status
c
      character*1 jobz, range, uplo
c
      integer ha, adra          !A
      integer hb, adrb          !B
c
c
      logical oactive           ! true iff this process participates
      integer dimA18, dimA28, typeA
      integer dimB18, dimB28, typeB
      INTGR4 dimA1, dimA2,mout4
      INTGR4 dimB1, dimB2
c
      INTGR4 mpA, nqA          ! rows/cols of A held by the processor
      INTGR4 mpB, nqB          ! rows/cols of B held by the processor
c
      integer me
      INTGR4 lda, ldb
      integer elemA,elemB
      INTGR4 numroc
      
      INTGR4 nb                ! block size
      INTGR4 descA(9), descB(9) !descriptor for scalapack
      
c     
      integer n8, hfail, adrfail
      integer liwork,hiwork, adriwork
      integer lwork,hwork, adrwork
      integer lrwork,hrwork, adrrwork
      INTGR4 liwork4
      INTGR4 lwork4
      INTGR4 lrwork4
c
      double precision vl, vu
      INTGR4 il, iu
      INTGR4 m, nz
      INTGR4 n
      INTGR4 info,one4,zero4,two4,four4
      parameter(one4=1,zero4=0,two4=4,four4=4)
      integer info8,dblsize

      double precision pdlamch
      complex*16 dumm
      double precision dumm1
      INTGR4 dumm2
      external pdlamch
      INTGR4 iceil
c     
c     processor dependent; machine dependent
c     
      if(nb8.eq.0) then
         nb=16
      else
         nb=nb8
      endif
c     
c***  check environment
c     
      me     = ga_nodeid()
c     
c***  check GA info for input arrays
c     
      call ga_check_handle(g_a(1), 'ga_pzheevr: a')
      call ga_check_handle(g_b(1), 'ga_pzheevr: b')
c     
      call ga_inquire(g_a(1), typeA, dimA18, dimA28)
      call ga_inquire(g_b(1), typeB, dimB18, dimB28)
      dima1=dima18
      dima2=dima28
      dimb1=dimb18
      dimb2=dimb28
      if(dimA1.ne.dima2) call ga_error(
     %     'ga_pzheevr: matrix A not square ',0)
      n=dima1
c      nb=min(nb,n/2)
      if(dimb1.ne.dimb2) call ga_error(
     %     'ga_pzheevr: matrix B not square ',0)
      if(dimb1.ne.n) call ga_error(
     %     'ga_pzheevr: size matrix B and B differ ',0)
      
      if(mout.le.0.or.mout.gt.n) mout=n

c     
c     
c***  initialize SL interface
c     
      call SLinit2(n)
      oactive=iam.lt.maxproc
      call ga_zero(g_b(1))
      call ga_zero(g_b(2))
      if (oactive) then
c     
c***  find SBS format parameters
c     
c     
      mpA = numroc(dimA1, nb, myrow2, zero4, nprow2)
      nqA = numroc(dimA2, nb, mycol2, zero4, npcol2)
c     
      mpB = numroc(dimB1, nb, myrow2, zero4, nprow2)
      nqB = numroc(dimB2, nb, mycol2, zero4, npcol2) 
c     
c     
      lda = max(one4,mpA)
      ldb = max(one4,mpB)
c     
c     
c     let scalapack check for errors
c     
c     should check to see if this is a compute node
c     check to see how this works in the new data server model
c     
c     
         elemA= mpA*nqA
         status = .true.
         if(elemA.ne.0)status =
     $        ma_push_get(MT_dcpL,elemA,'a',ha,adra)
         if(.not.status) 
     &        call ga_error('ga_pzheevr: mem alloc failed A ', -1)
c     
c***  copy g_a to A and g_b to B using the block cyclic scalapack format 
c     
         
         call ga_to_SL2cpl(g_a, dimA1, dimA2, nb, nb,
     $        dcpl_mb(adrA), lda, mpA, nqA)
         endif
         call ga_sync()
         if(oactive) then
c     
         elemB= mpB*nqB
c     
         if(elemB.ne.0)status =
     $        ma_push_get(MT_dcpL,elemB,'b',hb,adrb)
         if(.not.status) 
     &        call ga_error('ga_pzheevr: mem alloc failed B ', -1)
c     
         n8 = n
         if(n8.ne.0) then
            status =
     $           ma_push_get(MT_INT, n8 ,'ifail',hfail,adrfail)
         endif
         if(.not.status)
     &        call ga_error('ga_pzheevr: mem alloc failed ifail  ', n)
         
c     
c***  fill SCALAPACK matrix descriptors
c     
         call descinit(descA, dimA1, dimA2, nb, nb, zero4, zero4,
     $        islctxt2, lda, info)
         info8=info
         if(info.ne.0) call ga_error(' ga_pzheevr: descinit A failed ',
     .        -info8)
c     
         call descinit(descB, dimB1, dimB2, nb, nb, zero4, zero4,
     .        islctxt2, ldb, info)
         info8=info
         if(info.ne.0) call ga_error(' ga_pzheevr: descinit B failed ',
     .        -info8)
c     
         jobz = 'V'
         uplo = 'L'
         vl = 0.d0
         vu = 0.d0
c     
         il=1
         iu=mout
         range='i'
         lwork4=-1
         lrwork4=-1
         liwork4=-1
         call pzheevr(jobz, range, uplo,
     1        n, dcpl_mb(adrA), one4, one4, 
     D        descA,vl,vu, il, iu, m, nz, eval,  dcpl_mb(adrB), one4, 
     3        one4, descB, 
     4        dumm, lwork4, dumm1,lrwork4,
     3        dumm2, liwork4, info)

         lwork=dumm
         lrwork=dumm1
         liwork=dumm2
c     
         if(liwork.ne.0)status =
     $        ma_push_get(MT_INT, liwork ,'iwork',hiwork,adriwork)
         if(.not.status) 
     &        call ga_error('ga_pzheevr: mem alloc failed iwork  ', -1)
c     
         if(lwork.ne.0)
     $        status = ma_push_get(MT_DCPL, lwork ,
     $        'work',hwork,adrwork)
         if(.not.status) 
     &        call ga_error('ga_pzheevr: mem alloc failed work  ', -1)
         if(lrwork.ne.0)
     $        status = ma_push_get(MT_DBL, lrwork ,
     $        'rwork',hrwork,adrrwork)
         if(.not.status) 
     &        call ga_error('ga_pzheevr: mem alloc failed rwork  ', -1)
c     
         liwork4=liwork
         lwork4=lwork
         lrwork4=lrwork

         call pzheevr(jobz, range, uplo,
     1        n, dcpl_mb(adrA), one4, one4, 
     D        descA,vl,vu, il, iu, m, nz, eval,  dcpl_mb(adrB), one4, 
     3        one4, descB, 
     4        dcpl_mb(adrwork), lwork4, dbl_mb(adrrwork),lrwork4,
     3        int_mb(adriwork), liwork4, info)
         mout=m
         if(iam.eq.0.and.mout.ne.n) 
     W        write(6,*) ' computed eval ',mout,' out of ',
     A        dima1

         if (nz .ne. n ) then
            if ( info .ne. 0 ) then
               if ( info .gt. 0 ) then
         call ga_error(' ga_pzheevr: argument is illegal ', info)
               else
         call ga_error(' ga_pzheevr: eigenvectors failed to converge ',
     $                 info)
               endif
            endif
         endif
c     
c     
c     
c***  copy solution matrix back to g_b
c     
         mout4=mout

         call ga_from_SL2cpl(g_b, dimA1, mout4,
     $        nb, nb, dcpl_mb(adrB),
     &        ldb, mpb, nqB)
c     
c     
c     
c***  deallocate work/SL arrays
c     
         if ( lrwork .ne. 0 ) status = ma_pop_stack(hrwork)
         if ( lwork .ne. 0 ) status = ma_pop_stack(hwork)
         if ( liwork .ne. 0 ) status = ma_pop_stack(hiwork)
         if ( n8 .ne. 0 ) status = ma_pop_stack(hfail)
         if ( elemB .ne. 0 ) status = ma_pop_stack(hb)
         if ( elemA .ne. 0 ) status = ma_pop_stack(ha)
c     
      endif
      call ga_sync()
      if(maxproc.lt.nnodes.or.dima1.le.nb) then
c     broadcast evals
         dblsize=ma_sizeof(MT_DBL,1,MT_BYTE)*dima1
         call ga_brdcst(1688,eval,dblsize,0)
      endif
      
      return
      end
#endif
      SUBROUTINE PXERBLA( ICTXT, SRNAME, INFO )
*
*  -- ScaLAPACK auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER*4            ICTXT, INFO
*     ..
*     .. Array Arguments ..
      CHARACTER*(*)      SRNAME
*     ..
*
*  Purpose
*  =======
*
*  PXERBLA is an error handler for the ScaLAPACK routines.  It is called
*  by a ScaLAPACK routine if an input parameter has an invalid value.  A
*  message is printed. Installers may consider modifying this routine in
*  order to call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  ICTXT   (local input) INTEGER*4
*          On entry,  ICTXT  specifies the BLACS context handle, indica-
*          ting the global  context of the operation. The context itself
*          is global, but the value of ICTXT is local.
*
*  SRNAME  (global input) CHARACTER*(*)
*          On entry, SRNAME specifies the name of the routine which cal-
*          ling PXERBLA.
*
*  INFO    (global input) INTEGER*4
*          On entry, INFO  specifies the position of the invalid parame-
*          ter in the parameter list of the calling routine.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER*4            MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO
*     ..
*     .. Executable Statements ..
*
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      WRITE( *, FMT = 9999 ) MYROW, MYCOL, SRNAME, INFO
*
 9999 FORMAT( '{', I5, ',', I5, '}:  On entry to ', A,
     $        ' parameter number ', I4, ' had an illegal value' )
*
      call ga_error(' pxerbla error: info equal ',
     I        info)
      RETURN
*
*     End of PXERBLA
*
      END

