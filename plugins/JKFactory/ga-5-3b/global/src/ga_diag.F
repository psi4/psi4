#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine gai_diag_std(g_a, g_v, eval)
C$Id: ga_diag.F,v 1.26 2004-10-19 22:39:20 edo Exp $
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "diag.fh"
      integer g_a               ! matrix to diagonalize
      integer g_v               ! global matrix to return evecs
      double precision eval(*)  ! local array to return evals
c
c
c     solve the standard symmetric eigenvalue problem returning
c     all eigenvectors and values in ascending order
c
c     the input matrix is not destroyed
c
c     proc - # of processors
c     pans - panel size in columns
c     n    - problem size
c
c     
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      integer iproc, pan, istart, iend, k, n, npan, proc, myproc
      integer mypanel, myelem, mycol, elemz, ipan, ik, info
      logical status
c
      integer hma, adrma        !mapa
      integer hmz, adrmz        !mapz
      integer ga_diag_type
c
      integer ha, adra          !a
      integer hz, adrz          !z
c
      integer isize, rsize, ptr_size
      integer his, adris
      integer hsc, adrsc
      integer hip, adrip
c
      external  pdspev
c
      logical oactive           ! true iff this process participates
      integer dimA1, dimA2, typeA
      integer dimV1, dimV2, typeV
c
      integer hls, adrls        ! list of TCG nodes
c
      integer Elem
      Elem(istart,iend,n)=((iend-istart+1)*(2*n-istart-iend+2))/2
c
c
c***  check environment
c
      myproc = ga_nodeid()
      proc   = ga_nnodes()
c
c***  check GA info
c
      call ga_check_handle(g_a, 'gai_diag: a')
      call ga_check_handle(g_v, 'gai_diag: v')
      call ga_inquire(g_a,  typeA, dimA1, dimA2)
      call ga_inquire(g_v,  typeV, dimV1, dimV2)
c
c
      if (dimA1 .ne. dimA2) then
         call ga_error('gai_diag_std: can only diag square matrix ',1)
      else if (dimV1 .ne. dimV2) then
         call ga_error('gai_diag_std: can only diag square matrix ',3)
      else if(dimA1.ne.dimV1 )then
         call ga_error('gai_diag_std: A & V dimension mismatch ', 4)
      endif
c
      n = dimA1
      if (n .lt. 64) then
             call ga_diag_std_seq(g_a,  g_v, eval)
             return
      endif
c
      call ga_sync()
c
      pan = 1                   !panel size
      npan = n/pan
c
      proc = min(n/30, ga_nnodes())
c     
      oactive = myproc .lt. proc
      if (oactive) then
*     print *,myproc,' participate in eig: ',proc
*     call ffflush(6)
c     
c     allocate map arrays
c     
         status = .true.
         status = status .and. ma_push_get(mt_int,n,'mapa',hma,adrma)
         status = status .and. ma_push_get(mt_int,n,'mapz',hmz,adrmz)
         if (.not. status)
     $        call ga_error('gai_diag_std: insufficient memory 1', n)
c
c        allocate and setup list of TCG nodes
         if( .not. ma_push_get(mt_int, proc,'list',hls,adrls))
     $        call ga_error('gai_diag_std: insufficient memory 1b', n)
         call ga_list_nodeid(int_mb(adrls), proc)
*         print *,myproc,'list ', (int_mb(adrls+k),k=0,proc-1)
*         call ffflush(6)
c         
c        determine distribution of a,b, and z
c     
         mypanel = 0            ! number of panels at given processor
         mycol   = 0            ! number of columns at given processor
         myelem  = 0            ! number of elements at given processor
c     
c        allocate first half of the panels - forwards
c
         do ipan = 1, npan/2 
            iproc = mod(ipan -1, proc)
            istart = (ipan - 1)*pan + 1
            iend   = istart + pan - 1
            do k = istart, iend
               int_mb(-1+adrma+k) = int_mb(adrls+iproc)
               int_mb(-1+adrmz+k) = int_mb(adrls+iproc)
            enddo
            if(iproc .eq. myproc)then
               myelem = myelem + elem(istart,iend,n)
               mypanel = mypanel+1	
            endif
         enddo	   
c
c        allocate second half of the panels - backwards
c
         ik = 1                 !ik is used to forward number processors
         do ipan =   npan-1, npan/2 +1 , -1
            ik = ik+1
            iproc = mod(ik -1, proc)
            istart = (ipan - 1)*pan + 1
            iend   = istart + pan - 1
            do k = istart, iend
               int_mb(-1+adrma+k) = int_mb(adrls+iproc)
               int_mb(-1+adrmz+k) = int_mb(adrls+iproc)
            enddo
            if(iproc .eq. myproc)then
               myelem = myelem + elem(istart,iend,n)
               mypanel = mypanel+1	
            endif
         enddo	   
c
c***     actually, there is one more panel left for the mismatch 
c     
         iproc  = 0             !processor 0 gets the mismatch panel
         istart = (npan - 1)*pan + 1
         iend   = n
         do k = istart, iend
            int_mb(-1+adrma+k) = int_mb(adrls+iproc)
            int_mb(-1+adrmz+k) = int_mb(adrls+iproc)
         enddo
         mycol = mypanel * pan
         if(iproc .eq. myproc)then
            myelem = myelem + elem(istart,iend,n)
            mycol  = mycol + iend - istart + 1	
         endif
         elemz = mycol * n
c
c***     allocate a,b, and z arrays
c     
         status = .true.
         status = status .and.
     $        ma_push_get(mt_dbl,myelem,'gai_diag_std:a',ha,adra)
         status = status .and.
     $        ma_push_get(mt_dbl,elemz,'gai_diag_std:z',hz,adrz)
         if (.not. status)
     $        call ga_error('gai_diag_std: insufficient mem 2',
     $        myelem+elemz)
c
c***     load  matrix A 
c
         istart = 0
         do k = 1, n
            if( int_mb(-1+adrma+k) .eq. int_mb(adrls + myproc))then
               myelem = n-k+1
               call ga_get(g_a,k,n,k,k,dbl_mb(adra+istart),1)
               istart = istart + myelem
            endif
         enddo
c     
c***     allocate some memory for std eigensolver 
c
*
* temporary patch following george reccomendation to use 6n
*
         if (.not. ma_push_get(mt_int,6*n,'iscratch',his,adris))
     $        call ga_error('gai_diag_std: insufficient mem 3', 6*n)
         call fmemreq(1,n,int_mb(adrma),int_mb(adrma),int_mb(adrmz),
     &        isize, rsize, ptr_size,int_mb(adris))
         status = ma_pop_stack(his)
         if(.not. ma_push_get(mt_int,isize,'iscratch',his,adris))
     &        call ga_error('gai_diag_std allocator',1)
         if(.not. ma_push_get(mt_dbl,rsize,'scratch',hsc,adrsc))
     &        call ga_error('gai_diag_std allocator',2)
         if(.not. ma_push_get(mt_dbl,ptr_size,'iptr',hip,adrip))
     &        call ga_error('gai_diag_std allocator',3)
c     
      endif
c
c
      call ga_sync()            ! To reduce message buffering problems
c
c***  calling PEIGS parallel eigensolver
c
      if (oactive) then
*          print *,myproc,' participate in eig: ',proc
*          call ffflush(6)
*        if(myproc.eq.0)then
*          print *,myproc,'map a ', (int_mb(adrma +k),k=0,n-1)
*          print *,myproc,'a ', (dbl_mb(adra +k),k=0,5)
*          print *,myproc,'map z ', (int_mb(adrmz +k),k=0,n-1)
*          print *, ' '
*          call ffflush(6)
*        endif
*        call ga_sync()
*        if(myproc.eq.1)then
*          print *,myproc,'map a ', (int_mb(adrma +k),k=0,n-1)
*          print *,myproc,'a ', (dbl_mb(adra +k),k=0,3)
*          print *,myproc,'map z ', (int_mb(adrmz +k),k=0,n-1)
*          call ffflush(6)
*        endif
         call pdspev(  n, dbl_mb(adra), int_mb(adrma),
     &        dbl_mb(adrz), int_mb(adrmz),  
     &        eval,  
     &        int_mb(adris), isize,  
     &        dbl_mb(adrip), ptr_size,  
     &        dbl_mb(adrsc), rsize, info)  
*          print *,myproc,'out: participate in eig: ',proc
*          call ffflush(6)
         if (info .ne. 0)
     $        call ga_error('gai_diag_std: ... info non-zero ', myproc)
      endif
c     
      call ga_sync()            ! To reduce message buffering problems      
c
c***  store the eigenvector matrix in global array g_v
c
      if (oactive) then
         istart = 0
         do k = 1, n
            if( int_mb(-1+adrmz+k) .eq. int_mb(adrls + myproc))then
               call ga_put(g_v,1,n,k,k,dbl_mb(adrz+istart),1)
               istart = istart + n
            endif
         enddo
c     
c***     deallocate work arrays 
c     
         status = .true.
         status = status .and. ma_pop_stack(hip)
         status = status .and. ma_pop_stack(hsc)
         status = status .and. ma_pop_stack(his)
         status = status .and. ma_pop_stack(hz)
         status = status .and. ma_pop_stack(ha)
         status = status .and. ma_pop_stack(hls)
         status = status .and. ma_pop_stack(hmz)
         status = status .and. ma_pop_stack(hma)
         if (.not. status) call ga_error('gai_diag_std: ma_pop',0)
      endif
c
c
c     Make sure that even processes that did not participate in the
c     diagonalization have the eigen vectors
c
      ga_diag_type = 32760 - 30 
      call ga_brdcst(ga_diag_type, eval, 
     $               ma_sizeof(MT_DBL,n,MT_BYTE), 0)
c     if(myproc.eq.0)then
c       print *, 'evals after brdcast'
c       print *, (eval(k),k=1,n)
c       call ffflush(6)
c     endif

c
      call ga_sync()
c
      end


      subroutine gai_diag(g_a, g_s, g_v, eval)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "diag.fh"
      integer g_a               ! matrix to diagonalize
      integer g_s               ! metric
      integer g_v               ! global matrix to return evecs
      double precision eval(*)  ! local array to return evals
c
c
c     solve the generalized eigen-value problem returning
c     all eigen-vectors and values in ascending order
c
c     the input matrices are not destroyed
c
c     proc - # of processors
c     pans - panel size in columns
c     n    - problem size
c
c     j. nieplocha 08.22.93
c     
c     upgraded for shared memory: 01.14.94
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      integer iproc, pan, istart, iend, k, n, npan, proc, myproc
      integer mypanel, myelem, mycol, elemz, ipan, ik, info
      logical status
c
      integer hma, adrma        !mapa
      integer hmb, adrmb        !mapb
      integer hmz, adrmz        !mapz
      integer ga_diag_type
c
      integer ha, adra          !a
      integer hb, adrb          !b
      integer hz, adrz          !z
c
      integer isize, rsize, ptr_size
      integer his, adris
      integer hsc, adrsc
      integer hip, adrip
c
      double precision drand
      external pdspgv, drand
      intrinsic mod
c
      logical oactive           ! true iff this process participates
      integer dimA1, dimA2, typeA
      integer dimS1, dimS2, typeS
      integer dimV1, dimV2, typeV
c
      integer hls, adrls        ! list of TCG nodes
c
c
      integer Elem
      Elem(istart,iend,n)=((iend-istart+1)*(2*n-istart-iend+2))/2
c
c
c***  check environment
c
      myproc = ga_nodeid()
      proc   = ga_nnodes()
c
c***  check GA info
c
      call ga_check_handle(g_a, 'gai_diag: a')
      call ga_check_handle(g_s, 'gai_diag: s')
      call ga_check_handle(g_v, 'gai_diag: v')
      call ga_inquire(g_a, typeA, dimA1, dimA2)
      call ga_inquire(g_s, typeS, dimS1, dimS2)
      call ga_inquire(g_v, typeV, dimV1, dimV2)
c
c
      if (dimA1 .ne. dimA2) then
         call ga_error('gai_diag: can only diag square matrix ', 1)
      else if (dimS1 .ne. dimS2) then
         call ga_error('gai_diag: can only diag square matrix ', 2)
      else if (dimV1 .ne. dimV2) then
         call ga_error('gai_diag: can only diag square matrix ', 3)
      else if(dimA1.ne.dimS1 .or. dimS1 .ne. dimV1)then
         call ga_error('gai_diag: A, S, V dimension mismatch ', 4)
      endif
c
      n = dimA1
      if (n .lt. 64) then
         call ga_diag_seq(g_a, g_s, g_v, eval)
         return
      endif
c
      call ga_sync()
c
      pan = 1                   !panel size
      npan = n/pan
c
	proc = min(ga_nnodes(), n/30)
c
c         print *,myproc,' participate in eig: ',proc
c         call ffflush(6)
      oactive = myproc .lt. proc
      if (oactive) then
c     
c        allocate map arrays
c
         status = .true.
         status = status .and. ma_push_get(mt_int,n,'mapa',hma,adrma)
         status = status .and. ma_push_get(mt_int,n,'mapb',hmb,adrmb)
         status = status .and. ma_push_get(mt_int,n,'mapz',hmz,adrmz)
         if (.not. status)
     $        call ga_error('gai_diag: insufficient memory 1', n)
c
c        allocate and setup list of TCG nodes
         if( .not. ma_push_get(mt_int, proc,'list',hls,adrls))
     $        call ga_error('gai_diag: insufficient memory 1b', n)
         call ga_list_nodeid(int_mb(adrls), proc)
*         print *,myproc,'list ', (int_mb(adrls+k),k=0,proc-1)
*         call ffflush(6)
c         
c         
c        determine distribution of a,b, and z
c     
         mypanel = 0            ! number of panels at given processor
         mycol   = 0            ! number of columns at given processor
         myelem  = 0            ! number of elements at given processor
c     
c        allocate first half of the panels - forwards
c
         do ipan = 1, npan/2 
            iproc = mod(ipan -1, proc)
            istart = (ipan - 1)*pan + 1
            iend   = istart + pan - 1
            do k = istart, iend
               int_mb(-1+adrma+k) = int_mb(adrls+iproc)
               int_mb(-1+adrmb+k) = int_mb(adrls+iproc)
               int_mb(-1+adrmz+k) = int_mb(adrls+iproc)
            enddo
            if(iproc .eq. myproc)then
               myelem = myelem + elem(istart,iend,n)
               mypanel = mypanel+1	
            endif
         enddo	   
c
c        allocate second half of the panels - backwards
c
         ik = 1                 !ik is used to forward number processors
         do ipan =   npan-1, npan/2 +1 , -1
            ik = ik+1
            iproc = mod(ik -1, proc)
            istart = (ipan - 1)*pan + 1
            iend   = istart + pan - 1
            do k = istart, iend
               int_mb(-1+adrma+k) = int_mb(adrls+iproc)
               int_mb(-1+adrmb+k) = int_mb(adrls+iproc)
               int_mb(-1+adrmz+k) = int_mb(adrls+iproc)
            enddo
            if(iproc .eq. myproc)then
               myelem = myelem + elem(istart,iend,n)
               mypanel = mypanel+1	
            endif
         enddo	   
c
c        actually, there is one more panel left for the mismatch 
c     
         iproc  = 0             !processor 0 gets the mismatch panel
         istart = (npan - 1)*pan + 1
         iend   = n
         do k = istart, iend
            int_mb(-1+adrma+k) = int_mb(adrls+iproc)
            int_mb(-1+adrmb+k) = int_mb(adrls+iproc)
            int_mb(-1+adrmz+k) = int_mb(adrls+iproc)
         enddo
         mycol = mypanel * pan
         if(iproc .eq. myproc)then
            myelem = myelem + elem(istart,iend,n)
            mycol  = mycol + iend - istart + 1	
         endif
         elemz = mycol * n
c
c        allocate a,b, and z arrays
c     
         status = .true. 
         status = status .and.
     $        ma_push_get(mt_dbl,myelem,'gai_diag:a',ha,adra)
         status = status .and.
     $        ma_push_get(mt_dbl,myelem,'gai_diag:b',hb,adrb)
         status = status .and.
     $        ma_push_get(mt_dbl,elemz,'gai_diag:z',hz,adrz)
         if (.not. status) call ga_error('gai_diag: memory?', myelem)
c
c        load A and S matrices 
c
         istart = 0
         do k = 1, n
            if( int_mb(-1+adrma+k) .eq. int_mb(adrls + myproc))then
               myelem = n-k+1
               call ga_get(g_a,k,n,k,k,dbl_mb(adra+istart),1)
               call ga_get(g_s,k,n,k,k,dbl_mb(adrb+istart),1)
               istart = istart + myelem
            endif
         enddo
c     
c        allocate some memory for eigensolver 
c
*
* temporary patch following george reccomendation to use 6n
*
         if (.not. ma_push_get(mt_int,6*n,'iscratch',his,adris))
     $        call ga_error('gai_diag: memory 2?', 6*n)
         call fmemreq(0,n,int_mb(adrma),int_mb(adrmb),int_mb(adrmz),
     &        isize, rsize, ptr_size,int_mb(adris))
         if (.not. ma_pop_stack(his))
     $        call ga_error('gai_diag: ma pop his', 0)
         if(.not. ma_push_get(mt_int,isize,'iscratch',his,adris))
     &        call ga_error('eig allocator',1)
         if(.not. ma_push_get(mt_dbl,rsize,'scratch',hsc,adrsc))
     &        call ga_error('eig allocator',2)
         if(.not. ma_push_get(mt_dbl,ptr_size,'iptr',hip,adrip))
     &        call ga_error('eig allocator',3)
c     
      endif
c
c
      call ga_sync()            ! To reduce message buffering problems
c
c***  calling PEIGS parallel eigensolver
c
      if (oactive) then
*          print *,myproc,' participate in eig: ',proc
*          call ffflush(6)
         call pdspgv( 1, n, dbl_mb(adra), int_mb(adrma),
     &        dbl_mb(adrb), int_mb(adrmb),  
     &        dbl_mb(adrz), int_mb(adrmz),  
     &        eval,  
     &        int_mb(adris), isize,  
     &        dbl_mb(adrip), ptr_size,  
     &        dbl_mb(adrsc), rsize, info)  
*          print *,myproc,'out: participate in eig: ',proc
*          call ffflush(6)
         if (info .ne. 0)
     $        call ga_error('gai_diag: info non-zero ', info)
      endif
c     
      call ga_sync()            ! To reduce message buffering problems      
c
c***  store the eigenvector matrix in global array g_v
c
      if (oactive) then
         istart = 0
         do k = 1, n
            if( int_mb(-1+adrmz+k) .eq. int_mb(adrls + myproc))then
               call ga_put(g_v,1,n,k,k,dbl_mb(adrz+istart),1)
               istart = istart + n
            endif
         enddo
c     
c***     deallocate work arrays 
c     
         status = .true.
         status = status .and. ma_pop_stack(hip) 
         status = status .and. ma_pop_stack(hsc) 
         status = status .and. ma_pop_stack(his) 
         status = status .and. ma_pop_stack(hz) 
         status = status .and. ma_pop_stack(hb) 
         status = status .and. ma_pop_stack(ha) 
         status = status .and. ma_pop_stack(hls)
         status = status .and. ma_pop_stack(hmz) 
         status = status .and. ma_pop_stack(hmb) 
         status = status .and. ma_pop_stack(hma)
         if (.not. status) call ga_error('gai_diag: ma_pop?', 0)
      endif
c
c
c     Make sure that even processes that did not participate in the
c     diagonalization have the eigen vectors
c
c     if(myproc.eq.proc-1)then
c       print *, 'evals P-1'
c       print *, (eval(k),k=1,n)
c       call ffflush(6)
c     endif

      ga_diag_type = 32760 - 30 
      call ga_brdcst(ga_diag_type, eval, 
     $               ma_sizeof(MT_DBL,n,MT_BYTE), 0)
c
      call ga_sync()
c
      end
c
      subroutine gai_diag_reuse(reuse,g_a, g_s, g_v, eval)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "diag.fh"
      integer reuse             ! allows to reuse factorized g_s   
                                !     0  first time,
                                !    >0  following calls 
                                !    <0  only deletes factorized g_s 
      integer g_a               ! matrix to diagonalize
      integer g_s               ! metric
      integer g_v               ! global matrix to return evecs
      double precision eval(*)  ! local array to return evals
c
c
c     solve the generalized eigen-value problem returning
c     all eigen-vectors and values in ascending order
c
c     the input matrices are not destroyed
c
c     proc - # of processors
c     pans - panel size in columns
c     n    - problem size
c
c     j. nieplocha 02.08.94
c     
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      integer iproc, pan, istart, iend, k, n, npan, proc, myproc
      integer mypanel, myelem, mycol, elemz, ipan, ik, info
      logical status, first_time
      integer g_ps, eig_flag, thr_seq
      save g_ps, first_time
      
c
      integer hma, adrma        !mapa
      integer hmb, adrmb        !mapb
      integer hmz, adrmz        !mapz
      integer ga_diag_type
c
      integer ha, adra          !a
      integer hb, adrb          !b
      integer hz, adrz          !z
c
      integer isize, rsize, ptr_size
      integer his, adris
      integer hsc, adrsc
      integer hip, adrip
c
      double precision drand
      external pdspgv, drand
      intrinsic mod
c
      logical oactive           ! true iff this process participates
      integer dimA1, dimA2, typeA
      integer dimS1, dimS2, typeS
      integer dimV1, dimV2, typeV
c
      integer hls, adrls        ! list of TCG nodes
c
c
      integer Elem
      Elem(istart,iend,n)=((iend-istart+1)*(2*n-istart-iend+2))/2
      data first_time /.true./
      data thr_seq /64/         !min problem size to be run in parallel 
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c***  check the environment
c
      myproc = ga_nodeid()
      proc   = ga_nnodes()
c
c***  check GA handles, types and dimensions 
c
      call ga_check_handle(g_a, 'gai_diag: a')
      call ga_check_handle(g_s, 'gai_diag: s')
      call ga_check_handle(g_v, 'gai_diag: v')
c
c...  unless it is the first call, check if g_s has not changed 
c  
      if(first_time)then
         g_ps = g_s
      else
         if(g_ps .ne. g_s) 
     $    call ga_error('gai_diag:attempt to reuse invalid data ',g_ps) 
      endif
c
      call ga_inquire(g_a, typeA, dimA1, dimA2)
      call ga_inquire(g_s, typeS, dimS1, dimS2)
      call ga_inquire(g_v, typeV, dimV1, dimV2)
c
c...  verify if the dimensions match
c
      if (dimA1 .ne. dimA2) then
         call ga_error(' gai_diag: can only diag square matrix ', 1)
      else if (dimS1 .ne. dimS2) then
         call ga_error(' gai_diag: can only diag square matrix ', 2)
      else if (dimV1 .ne. dimV2) then
         call ga_error(' gai_diag: can only diag square matrix ', 3)
      else if(dimA1.ne.dimS1 .or. dimS1 .ne. dimV1)then
         call ga_error(' gai_diag: A, S, V dimension mismatch ', 4)
      endif
c
c***  for smallest problems use sequential eigensolver
c
      n = dimA1
      if (n .lt. thr_seq) then
         call ga_diag_seq(g_a, g_s, g_v, eval)
         return
      endif
c
      call ga_sync()
c
c***  test the reuse value and take the proper action
c
      if((reuse .ne. 0) .and. first_time)
     $ call ga_error(' gai_diag: cannot reuse uninitialized data',reuse)
c
c...  detect repeated "first-time" calls -- they are allowed if 
c...         storage allocated for B is deallocated first
c
      if((reuse .eq. 0) .and. .not. first_time)
     $ call ga_error(' gai_diag_reuse: repeated call, memory problem',0)
c
c...  case for < 0 is just a cleanup of storage for B
c
      if(reuse .lt. 0) then
       if(.not. first_time) then
         if (.not. ma_free_heap(hb))
     $         call ga_error('ga_reuse: ma_free hb', 0)
         if (.not. ma_free_heap(hmb))
     $         call ga_error('ga_reuse: ma_free hmb', 0)
         first_time = .true.
       endif
       return
      endif
c
c...  set the flag for George eigensolver (orthogonal to reuse)
c
      if(reuse.eq.0) then
         eig_flag = 1
      else
         eig_flag = 0
      endif
c
c***  for smaller problems use fewer processors (performance reasons)
c
      pan = 1                   !panel size
      npan = n/pan
c
	proc = min( ga_nnodes(), n/30)
c
      oactive = myproc .lt. proc
      if (oactive) then
c     
c***     allocate storage for map arrays: A & Z from stack and B from heap
c
         status = ma_push_get(mt_int,n,'mapa',hma,adrma)
         if(.not.status)call ga_error('gai_diag: mapa not allocated ',0)

         if(first_time)status = ma_alloc_get(mt_int,n,'mapb',hmb,adrmb)
         if(.not.status)call ga_error('gai_diag: mapb not allocated ',0)

         status = ma_push_get(mt_int,n,'mapz',hmz,adrmz)
         if(.not.status)call ga_error('gai_diag: mapc not allocated ',0)
c
c        allocate and setup list of TCG nodes
         if( .not. ma_push_get(mt_int, proc,'list',hls,adrls))
     $        call ga_error('gai_diag_std: insufficient memory 1b', n)
         call ga_list_nodeid(int_mb(adrls), proc)
c         
c         
c***     determine distribution of a,b, and z
c     
         mypanel = 0            ! number of panels at given processor
         mycol   = 0            ! number of columns at given processor
         myelem  = 0            ! number of elements at given processor
c     
c...     allocate first half of the panels - forwards
c
         do ipan = 1, npan/2 
            iproc = mod(ipan -1, proc)
            istart = (ipan - 1)*pan + 1
            iend   = istart + pan - 1
            do k = istart, iend
               int_mb(-1+adrma+k) = int_mb(adrls+iproc)
               if(first_time)int_mb(-1+adrmb+k) = int_mb(adrls+iproc)
               int_mb(-1+adrmz+k) = int_mb(adrls+iproc)
            enddo
            if(iproc .eq. myproc)then
               myelem = myelem + elem(istart,iend,n)
               mypanel = mypanel+1	
            endif
         enddo	   
c
c...     allocate second half of the panels - backwards
c
         ik = 1                 !ik is used to forward number processors
         do ipan =   npan-1, npan/2 +1 , -1
            ik = ik+1
            iproc = mod(ik -1, proc)
            istart = (ipan - 1)*pan + 1
            iend   = istart + pan - 1
            do k = istart, iend
               int_mb(-1+adrma+k) = int_mb(adrls+iproc)
               if(first_time)int_mb(-1+adrmb+k) = int_mb(adrls+iproc)
               int_mb(-1+adrmz+k) = int_mb(adrls+iproc)
            enddo
            if(iproc .eq. myproc)then
               myelem = myelem + elem(istart,iend,n)
               mypanel = mypanel+1	
            endif
         enddo	   
c
c...     actually, there is one more panel left for the mismatch 
c     
         iproc  = 0             !processor 0 gets the mismatch panel
         istart = (npan - 1)*pan + 1
         iend   = n
         do k = istart, iend
            int_mb(-1+adrma+k) = int_mb(adrls+iproc)
            if(first_time)int_mb(-1+adrmb+k) = int_mb(adrls+iproc)
            int_mb(-1+adrmz+k) = int_mb(adrls+iproc)
         enddo
         mycol = mypanel * pan
         if(iproc .eq. myproc)then
            myelem = myelem + elem(istart,iend,n)
            mycol  = mycol + iend - istart + 1	
         endif
         elemz = mycol * n
c
c***     allocate arrays: A & Z from stack and B from heap
c     
         status = ma_push_get(mt_dbl,myelem,'a',ha,adra)
         if(.not. status)call ga_error('gai_diag: a not allocated ',0)
         if(first_time)status = ma_alloc_get(mt_dbl,myelem,'b',hb,adrb)
         if(.not. status)call ga_error('gai_diag: b not allocated ',0)
         status = ma_push_get(mt_dbl,elemz,'z',hz,adrz)
         if(.not. status)call ga_error('gai_diag: c not allocated ',0)
c
c***     load A and S matrices 
c
         istart = 0
         do k = 1, n
            if( int_mb(-1+adrma+k) .eq. int_mb(adrls + myproc))then
               myelem = n-k+1
               call ga_get(g_a,k,n,k,k,dbl_mb(adra+istart),1)
               if(first_time)
     %            call ga_get(g_s,k,n,k,k,dbl_mb(adrb+istart),1)
               istart = istart + myelem
            endif
         enddo
c     
c***     allocate some scratch memory for eigensolver 
c
*
* temporary patch following george reccomendation to use 6n
*
         status = ma_push_get(mt_int,6*n,'iscratch',his,adris)
         if (.not. status) call ga_error('ga_reuse: 6*n', 6*n)
         call fmemreq(0,n,int_mb(adrma),int_mb(adrmb),int_mb(adrmz),
     &        isize, rsize, ptr_size,int_mb(adris))
         status = ma_pop_stack(his)
         if (.not. status) call ga_error('ga_reuse: ma pop his', 0)
         if(.not. ma_push_get(mt_int,isize,'iscratch',his,adris))
     &        call ga_error('eig allocator',1)
         if(.not. ma_push_get(mt_dbl,rsize,'scratch',hsc,adrsc))
     &        call ga_error('eig allocator',2)
         if(.not. ma_push_get(mt_dbl,ptr_size,'iptr',hip,adrip))
     &        call ga_error('eig allocator',3)
c     
      endif
c
c
      call ga_sync()            ! To reduce message buffering problems
c
c***  calling PEIGS parallel eigensolver
c
      if (oactive) then
         call pdspgv( eig_flag, n, dbl_mb(adra), int_mb(adrma),
     &        dbl_mb(adrb), int_mb(adrmb),  
     &        dbl_mb(adrz), int_mb(adrmz),  
     &        eval,  
     &        int_mb(adris), isize,  
     &        dbl_mb(adrip), ptr_size,  
     &        dbl_mb(adrsc), rsize, info)  
         if (info .ne. 0)
     $        call ga_error('gai_diag: info non-zero ', info)
      endif
c     
      call ga_sync()            ! To reduce message buffering problems      
c
c***  store the eigenvector matrix in global array g_v
c
      if (oactive) then
         istart = 0
         do k = 1, n
            if( int_mb(-1+adrmz+k) .eq. int_mb(adrls + myproc))then
               call ga_put(g_v,1,n,k,k,dbl_mb(adrz+istart),1)
               istart = istart + n
            endif
         enddo
c     
c***     deallocate storage for temporary arrays 
c     
         status = .true.
         status = status .and. ma_pop_stack(hip) 
         status = status .and. ma_pop_stack(hsc) 
         status = status .and. ma_pop_stack(his) 
         status = status .and. ma_pop_stack(hz) 
         status = status .and. ma_pop_stack(ha) 
         status = status .and. ma_pop_stack(hls)
         status = status .and. ma_pop_stack(hmz) 
         status = status .and. ma_pop_stack(hma) 
         if (.not. status) call ga_error('ga_reuse: ma_pop?', 0)
      endif
c
c
c***  Make sure that even processes that did not participate in the
c***  diagonalization have the eigen vectors
c
      ga_diag_type = 32760 - 30 
      call ga_brdcst(ga_diag_type, eval, 
     $               ma_sizeof(MT_DBL,n,MT_BYTE), 0)
c
      if(reuse.eq.0) first_time = .false.
c
      call ga_sync()
      end
