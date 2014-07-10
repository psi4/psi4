#if HAVE_CONFIG_H
#   include "config.fh"
#endif
C
C     Test the minval, minloc, maxval, maxloc, and enum functions in GA.
C
      program main
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer heap, stack, fudge, ma_heap, me, nproc
      logical status
      parameter (heap=100*100*4, fudge=100, stack=100*100)
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
      nproc = ga_nnodes()
      me = ga_nodeid()
c     we can also use GA_set_memory_limit BEFORE first ga_create call
c
      if(ga_nodeid().eq.0)then
         print *,' GA initialized '
         call ffflush(6)
      endif
c
c***  Initialize the MA package
c     MA must be initialized before any global array is allocated
c
      status = ma_init(MT_DCPL, stack, heap)
      if (.not. status) call ga_error('ma_init failed',-1) 
c
      if(me.eq.0)then
        print *, 'using ', nproc, ' process(es)'
        call ffflush(6)
      endif
c
      call test_nga_bin()

c     if(me.eq.0) call ga_print_stats()
c
c***  Tidy up the GA package
c
      call ga_terminate()
c
c***  Tidy up after message-passing library
c
      call MP_FINALIZE()
c
      end

      logical function ga_create_bin(nelem, type, name, g_off, g_bin)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
      integer nelem, type, g_off, g_bin
      character*1 name
c
      integer lo,hi, map(100),i,me,nproc
c
      nproc = ga_nnodes()
      me = ga_nodeid()
c
      do i = 1,nproc
         map(i) = 0
      enddo
c
c     call ga_print_distribution(g_off)
      call nga_distribution(g_off,me,lo,hi)
      call ga_sync()
c
      call nga_get(g_off,lo,lo,map(me+1),1)
      if(me.gt.0. and. map(me+1).eq.0)then
         print *,'my lo/hi in offset array',lo,hi,map(me+1)
         call ga_error('got 0',0)
      endif
      map(me+1) = map(me+1)+1
      call ga_igop(1,map,nproc,'+')
c
      ga_create_bin=nga_create_irreg(type,1,nelem,name,map,nproc,g_bin)
      end
      
     
     
      subroutine test_nga_bin()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
C
      integer g_ielm1, g_result1                      ! handles to global arrays
      integer g_ibin, g_ibincnt, g_ibinoff            ! handles to global arrays
      integer g_ibinga, g_idiff                       ! handles to global arrays
      integer g_ielm2, g_sbit2, g_src2, g_sink2, g_result2 ! handles to global arrays
      integer g_irowmat, g_icolmat                    ! handles to global arrays
      integer g_velm4, g_velm5                        ! handles to global arrays
C
      integer ilo, ihi, num
      integer ilo4, ihi4, num4
      integer i, it, count
      integer imax_bin, imin_bin
      integer iboff1, iboff2
      integer me, nproc                       ! my processor & number of procs 
      integer ndim,dims(1),chunck(1)
C
      integer nelements_max, nen_max
      parameter (nelements_max=1000, nen_max=4)
      integer ipack(nelements_max), 
     *        jpack(nelements_max), 
     *        kpack(nelements_max), 
     *        ielm1(nelements_max), 
     *        ielm2(nelements_max), 
     *        iarray(nelements_max), 
     *        ielm4(nen_max*nelements_max)
      double precision velm4(nen_max*nelements_max), 
     *                 velm5(nen_max*nelements_max)
C
      integer nbins_max
      parameter (nbins_max=1000)
      integer ibincnt(nbins_max), ibinoff(nbins_max)
c
      integer ilocmax, ilocmin
      integer ilocmax_ga, ilocmin_ga
      integer imax_ga, imin_ga
C
      integer nelements, nen
      integer isum, icount, jcount, itoff
      double precision result
      integer i1, i2, ibin
      logical ga_create_bin
C
C
c
c***  check parallel environment
      me = ga_nodeid()
      nproc = ga_nnodes()
c
      nelements=80
      nen=4
c
c***  create a global 1-D array
      ndim=1
      dims(1)=nen*nelements
      chunck(1)=min(1,dims(1)/nproc)
      if (.not. nga_create(MT_INT, ndim, dims, 'ielm1', chunck, 
     $                     g_ielm1))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'result1', chunck, 
     $                     g_result1))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'sbit2', chunck, 
     $                     g_sbit2))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'src2', chunck, 
     $                     g_src2))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'sink2', chunck, 
     $                     g_sink2))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'result2', chunck, 
     $                     g_result2))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'irowmat', chunck, 
     $                     g_irowmat))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'icolmat', chunck, 
     $                     g_icolmat))
     $     call ga_error(' ga_create failed ',0)
C
      if (.not. nga_create(MT_DBL, ndim, dims, 'velm4', chunck, 
     $                     g_velm4))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_DBL, ndim, dims, 'velm5', chunck, 
     $                     g_velm5))
     $     call ga_error(' ga_create failed ',0)
C
c
c     Enumerate the sbit arrary to get a sequential vector.
      call ga_fill(g_ielm1,0)
      call ga_patch_enum(g_ielm1,1,nen*nelements,1,1)
C     *** ielm1: 1 2 3 4 5 ... nelements ... 0 0 0 0 0
c
c***  compute local ilo, ihi, num for each processor
      call nga_distribution(g_ielm1,me,ilo,ihi,chunck)
      num=ihi-ilo+1
      if(ihi.le.0) num=0
      if(ihi.gt.nen*nelements) then
         ihi=nelements
         num=max(0,ihi-ilo+1)
      endif
C*****print *,"Random: ",ilo,ihi,ndim,dims(1)
      if(num.gt.nelements_max) then
         print *,"Too many elements "
         goto 9999
      endif
C*****print *, 'me=',me, num,ilo,ihi
c
      itoff=nen*nelements
      do i=ilo,ihi
c        print *, 'rand=',10.0*drand(0),
c    *             nint(10.0*drand(0)), drand(0),drand(0)
         ielm1(i-ilo+1)=max(1,min(nbins_max,1+nint(20.0*drand(0))))
      enddo  
      call ssortii(ielm1,ielm1,num,1)
      if(num.gt.0) call nga_put(g_ielm1,ilo,ihi,ielm1,1)
c
      if(num.gt.0) call nga_get(g_ielm1,ilo,ihi,ielm1,1)
C*****print *,"ielm1: ",me,(ielm1(i-ilo+1),i=ilo,ihi)
      call ga_sync()
C*****call ga_print(g_ielm1)
c
      call nga_select_elem(g_ielm1,'max',imax_bin,ilocmax_ga)
      call nga_select_elem(g_ielm1,'min',imin_bin,ilocmin_ga)
C
      ndim=1
      dims(1)=imax_bin
      chunck(1)=min(1,dims(1)/nproc)
      if (.not. nga_create(MT_INT, ndim, dims, 'ibincnt', chunck, 
     $                     g_ibincnt))
     $     call ga_error(' ga_create failed ',0)
      if (.not. nga_create(MT_INT, ndim, dims, 'ibinoff', chunck, 
     $                     g_ibinoff))
     $     call ga_error(' ga_create failed ',0)
      ndim=1
      dims(1)=nen*nelements
      chunck(1)=min(1,dims(1)/nproc)
      if (.not. nga_create(MT_INT, ndim, dims, 'ibin', chunck, 
     $                     g_ibin))
     $     call ga_error(' ga_create failed ',0)
c
      ndim=1
      dims(1)=nen*nelements
      chunck(1)=min(1,dims(1)/nproc)
      if (.not. nga_create(MT_INT, ndim, dims, 'idiff', chunck, 
     $                     g_idiff))
     $     call ga_error(' ga_create failed ',0)
c
      iboff1=1
      iboff2=0
      call ga_fill(g_src2,0)
      call ga_patch_enum(g_src2,1,itoff,1,1)
      if(num.gt.0) call nga_get(g_src2,ilo,ihi,ielm2,1)
      do ibin=imin_bin,imax_bin
         icount=0
         do i=ilo,ihi
            if(ielm1(i-ilo+1).eq.ibin) then
               icount=icount+1
               ipack(i-ilo+1)=1
            else
               ipack(i-ilo+1)=0
            endif
         enddo
         if(num.gt.0) then
            call nga_put(g_sbit2,ilo,ihi,ipack,1)
         endif 
         if(icount.gt.0)  then
            jcount=0
            do i=ilo,ihi
               if(ipack(i-ilo+1).eq.1) then
                  jcount=jcount+1
                  jpack(jcount)=i
               endif
            enddo
            icount=0
            call ga_fill(g_sink2,0)
            call ga_pack(g_src2, g_sink2,g_sbit2,1,itoff,icount)
            iboff2=iboff1+icount-1
            if(me.eq.0) print *,"Pack bin: ",ibin," count ",icount,
     *                          iboff1,iboff2
            call nga_copy_patch('n',g_sink2,1,icount,
     *                              g_ibin,iboff1,iboff2)
            iboff1=iboff1+icount
         endif
         ibincnt(ibin)=icount
         call ga_sync()
      enddo
      if(me.eq.0) then
         call nga_put(g_ibincnt,1,imax_bin,ibincnt,1)
         icount=0
         do i=1,imax_bin
            ibinoff(i)=icount
            icount=icount+ibincnt(i)
            print *, 'bin off',i, ibinoff(i)
         enddo
         call nga_put(g_ibinoff,1,imax_bin,ibinoff,1)
      endif
C
C*****call ga_print(g_ibin)
C
#if 1
      if(.not.ga_create_bin(nen*nelements,MT_INT,'ibin',g_ibinoff,
     $     g_ibinga))
     $     call ga_error(' ga_create failed ',0)
#else
      ndim=1
      dims(1)=nen*nelements
      chunck(1)=min(1,dims(1)/nproc)
      if (.not. nga_create(MT_INT, ndim, dims, 'ibinga', chunck,
     $                     g_ibinga))
     $     call ga_error(' ga_create failed ',0)
#endif
c
c     call ga_print_distribution(g_ibinga)
c     call ga_print_distribution(g_ibincnt)
c     call ga_print(g_ibinoff)
      call GA_bin_index(g_ibinga,g_ibincnt,g_ibinoff,ielm2,ielm1,num,1)
*     call ga_fill(g_ibinga,0)
*     call ga_copy(g_ibin,g_ibinga)
C
c     call ga_print(g_ibinga)
c     call ga_print(g_ibincnt)
      if(num.gt.0) call nga_get(g_ibin  ,ilo,ihi,ielm1,1)
      if(num.gt.0) call nga_get(g_ibinga,ilo,ihi,ielm2,1)
      isum=0
      do i=1,num
         if(ielm1(i).ne.ielm2(i)) then
            isum=isum+1
         endif
      enddo
      call ga_sync()
C
      call ga_igop(MT_INT,isum,1,'+') 
      if(me.eq.0) then
         if(isum.ne.0) then
            print *,"Binning: failed: ",isum
         else
            print *,"Binning: successful"
         endif
      endif
C
      call ga_sync()
C
      goto 9999
 9999 continue
C
C     ..................................................................
C
      return
      end
*dk,ssortii
      subroutine ssortii(x,y,n,kflag)
C
C #####################################################################
C
C     PURPOSE -
C
C        None
C
C     INPUT ARGUMENTS -
C
C        None
C
C     OUTPUT ARGUMENTS -
C
C        None
C
C     CHANGE HISTORY -
C
C        $Log: not supported by cvs2svn $
CPVCS    
CPVCS       Rev 1.2   08/03/95 13:53:18   dcg
CPVCS    replace print * with writloga calls
CPVCS
CPVCS       Rev 1.1   04/20/95 11:28:46   ejl
CPVCS    Fixed message when number of items is not greater then zero.
CPVCS
CPVCS
CPVCS       Rev 1.0   02/14/95 14:39:22   dcg
CPVCS    Original version
C
C ######################################################################
C
      implicit real*8 (a-h,o-z)
C***BEGIN PROLOGUE  SSORT
C***DATE WRITTEN   761101   (YYMMDD)
C***REVISION DATE  861211   (YYMMDD)
C***CATEGORY NO.  N6A2B1
C***KEYWORDS  LIBRARY=SLATEC,
C             TYPE=SINGLE PRECISION(SSORT-S DSORT-D ISORT-I),QUICKSORT,
C             SINGLETON QUICKSORT,SORT,SORTING
C***AUTHOR  JONES, R. E., (SNLA)
C           WISNIEWSKI, J. A., (SNLA)
C***PURPOSE  SSORT sorts array X and optionally makes the same
C            interchanges in array Y.  The array X may be sorted in
C            increasing order or decreasing order.  A slightly modified
C            QUICKSORT algorithm is used.
C***DESCRIPTION
C
C     Written by Rondall E. Jones
C     Modified by John A. Wisniewski to use the Singleton quicksort
C     algorithm.  Date 18 November 1976.
C
C     Abstract
C         SSORT sorts array X and optionally makes the same
C         interchanges in array Y.  The array X may be sorted in
C         increasing order or decreasing order.  A slightly modified
C         quicksort algorithm is used.
C
C     Reference
C         Singleton, R. C., Algorithm 347, An Efficient Algorithm for
C         Sorting with Minimal Storage, CACM,12(3),1969,185-7.
C
C     Description of Parameters
C         X - array of values to be sorted   (usually abscissas)
C         Y - array to be (optionally) carried along
C         N - number of values in array X to be sorted
C         KFLAG - control parameter
C             =2  means sort X in increasing order and carry Y along.
C             =1  means sort X in increasing order (ignoring Y)
C             =-1 means sort X in decreasing order (ignoring Y)
C             =-2 means sort X in decreasing order and carry Y along.
C***REFERENCES  SINGLETON,R.C., ALGORITHM 347, AN EFFICIENT ALGORITHM
C                 FOR SORTING WITH MINIMAL STORAGE, CACM,12(3),1969,
C                 185-7.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  SSORT
      integer X(N),Y(N),IL(21),IU(21)
      real*4 r
      character*80 logmess
      integer TY, TTY
C***FIRST EXECUTABLE STATEMENT  SSORT
      NN = N
      IF (NN.LT.0) THEN
         write(*,'(a)')
     x           'SSORT- THE NUMBER OF VALUES TO BE SORTED IS NEGATIVE'
         RETURN
      ELSEIF (NN.EQ.0) THEN
         write(*,'(a)')
     x            'SSORT- THE NUMBER OF VALUES TO BE SORTED IS ZERO'
         RETURN
      ENDIF
C
      KK = IABS(KFLAG)
      IF ((KK.EQ.1).OR.(KK.EQ.2)) GO TO 15
      write(*,'(a)')
     x  'SSORT- THE SORT CONTROL PARAMETER, K, WAS NOT 2, 1, -1, OR -2.'
      RETURN
C
C ALTER ARRAY X TO GET DECREASING ORDER IF NEEDED
C
   15 IF (KFLAG.GE.1) GO TO 30
      DO 20 I=1,NN
   20 X(I) = -X(I)
   30 GO TO (100,200),KK
C
C SORT X ONLY
C
  100 CONTINUE
      M=1
      I=1
      J=NN
      R=.375
  110 IF (I .EQ. J) GO TO 155
  115 IF (R .GT. .5898437) GO TO 120
      R=R+3.90625E-2
      GO TO 125
  120 R=R-.21875
  125 K=I
C                                  SELECT A CENTRAL ELEMENT OF THE
C                                  ARRAY AND SAVE IT IN LOCATION T
C*****IJ = I + IFIX (FLOAT (J-I) * sngl(R))
      IJ = I + IFIX (FLOAT (J-I) * R)
      T=X(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN T, INTERCHANGE WITH T
      IF (X(I) .LE. T) GO TO 130
      X(IJ)=X(I)
      X(I)=T
      T=X(IJ)
  130 L=J
C                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
C                                  T, INTERCHANGE WITH T
      IF (X(J) .GE. T) GO TO 140
      X(IJ)=X(J)
      X(J)=T
      T=X(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN T, INTERCHANGE WITH T
      IF (X(I) .LE. T) GO TO 140
      X(IJ)=X(I)
      X(I)=T
      T=X(IJ)
      GO TO 140
  135 TT=X(L)
      X(L)=X(K)
      X(K)=TT
C                                  FIND AN ELEMENT IN THE SECOND HALF OF
C                                  THE ARRAY WHICH IS SMALLER THAN T
  140 L=L-1
      IF (X(L) .GT. T) GO TO 140
C                                  FIND AN ELEMENT IN THE FIRST HALF OF
C                                  THE ARRAY WHICH IS GREATER THAN T
  145 K=K+1
      IF (X(K) .LT. T) GO TO 145
C                                  INTERCHANGE THESE ELEMENTS
      IF (K .LE. L) GO TO 135
C                                  SAVE UPPER AND LOWER SUBSCRIPTS OF
C                                  THE ARRAY YET TO BE SORTED
      IF (L-I .LE. J-K) GO TO 150
      IL(M)=I
      IU(M)=L
      I=K
      M=M+1
      GO TO 160
  150 IL(M)=K
      IU(M)=J
      J=L
      M=M+1
      GO TO 160
C                                  BEGIN AGAIN ON ANOTHER PORTION OF
C                                  THE UNSORTED ARRAY
  155 M=M-1
      IF (M .EQ. 0) GO TO 300
      I=IL(M)
      J=IU(M)
  160 IF (J-I .GE. 1) GO TO 125
      IF (I .EQ. 1) GO TO 110
      I=I-1
  165 I=I+1
      IF (I .EQ. J) GO TO 155
      T=X(I+1)
      IF (X(I) .LE. T) GO TO 165
      K=I
  170 X(K+1)=X(K)
      K=K-1
      IF (T .LT. X(K)) GO TO 170
      X(K+1)=T
      GO TO 165
C
C SORT X AND CARRY Y ALONG
C
  200 CONTINUE
      M=1
      I=1
      J=NN
      R=.375
  210 IF (I .EQ. J) GO TO 255
  215 IF (R .GT. .5898437) GO TO 220
      R=R+3.90625E-2
      GO TO 225
  220 R=R-.21875
  225 K=I
C                                  SELECT A CENTRAL ELEMENT OF THE
C                                  ARRAY AND SAVE IT IN LOCATION T
C*****IJ = I + IFIX (FLOAT (J-I) *sngl(R))
      IJ = I + IFIX (FLOAT (J-I) * R)
      T=X(IJ)
      TY= Y(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN T, INTERCHANGE WITH T
      IF (X(I) .LE. T) GO TO 230
      X(IJ)=X(I)
      X(I)=T
      T=X(IJ)
       Y(IJ)= Y(I)
       Y(I)=TY
      TY= Y(IJ)
  230 L=J
C                                  IF LAST ELEMENT OF ARRAY IS LESS THAN
C                                  T, INTERCHANGE WITH T
      IF (X(J) .GE. T) GO TO 240
      X(IJ)=X(J)
      X(J)=T
      T=X(IJ)
       Y(IJ)= Y(J)
       Y(J)=TY
      TY= Y(IJ)
C                                  IF FIRST ELEMENT OF ARRAY IS GREATER
C                                  THAN T, INTERCHANGE WITH T
      IF (X(I) .LE. T) GO TO 240
      X(IJ)=X(I)
      X(I)=T
      T=X(IJ)
       Y(IJ)= Y(I)
       Y(I)=TY
      TY= Y(IJ)
      GO TO 240
  235 TT=X(L)
      X(L)=X(K)
      X(K)=TT
      TTY= Y(L)
       Y(L)= Y(K)
       Y(K)=TTY
C                                  FIND AN ELEMENT IN THE SECOND HALF OF
C                                  THE ARRAY WHICH IS SMALLER THAN T
  240 L=L-1
      IF (X(L) .GT. T) GO TO 240
C                                  FIND AN ELEMENT IN THE FIRST HALF OF
C                                  THE ARRAY WHICH IS GREATER THAN T
  245 K=K+1
      IF (X(K) .LT. T) GO TO 245
C                                  INTERCHANGE THESE ELEMENTS
      IF (K .LE. L) GO TO 235
C                                  SAVE UPPER AND LOWER SUBSCRIPTS OF
C                                  THE ARRAY YET TO BE SORTED
      IF (L-I .LE. J-K) GO TO 250
      IL(M)=I
      IU(M)=L
      I=K
      M=M+1
      GO TO 260
  250 IL(M)=K
      IU(M)=J
      J=L
      M=M+1
      GO TO 260
C                                  BEGIN AGAIN ON ANOTHER PORTION OF
C                                  THE UNSORTED ARRAY
  255 M=M-1
      IF (M .EQ. 0) GO TO 300
      I=IL(M)
      J=IU(M)
  260 IF (J-I .GE. 1) GO TO 225
      IF (I .EQ. 1) GO TO 210
      I=I-1
  265 I=I+1
      IF (I .EQ. J) GO TO 255
      T=X(I+1)
      TY= Y(I+1)
      IF (X(I) .LE. T) GO TO 265
      K=I
  270 X(K+1)=X(K)
       Y(K+1)= Y(K)
      K=K-1
      IF (T .LT. X(K)) GO TO 270
      X(K+1)=T
       Y(K+1)=TY
      GO TO 265
C
C CLEAN UP
C
  300 IF (KFLAG.GE.1) RETURN
      DO 310 I=1,NN
  310 X(I) = -X(I)
      RETURN
      END
