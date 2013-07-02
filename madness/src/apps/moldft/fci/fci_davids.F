*
* $Id: fci_davids.F,v 1.2 1999/07/28 00:23:39 d3e129 Exp $
*
      logical function fci_davids(n, nroot, maxiter, maxsub, tol,
     $     work, lenwrk, evecfile, diagfile, eval, ecore, h, g)
      implicit none
c
      integer n                 ! [input] Matrix dimension
      integer nroot             ! [input] No. of eigen vectors sought
      integer maxiter           ! [input] Maximum no. of iterations
      integer maxsub            ! [input] Maximum dimension of iterative subspace
      double precision tol      ! [input] Required norm of residual
      integer lenwrk            ! [input] Length of scratch space
      double precision work(lenwrk) ! [scratch] Work space
      integer evecfile          ! [input] Handle for eigenvector file
      integer diagfile          ! [input] Handle for file with diagonals
      double precision eval(nroot) ! [output] Eigenvalues
      double precision ecore, h(*), g(*) ! [input] Integrals for fci_hv
c
c     Solve eigenvalue problem using multi-root block Davidson with
c     eigenvectors on disk in direct access files.
c
c     On input fileevec should contain nroot orthonormal guess vectors.
c     On output it will contain the eigenvectors.
c
      integer workfile          ! Handle for workspace file
      double precision test_tol, denom
      logical converged
      integer maxdim
      parameter (maxdim = 100)   ! Maximum size of iterative subspace
      integer buflen            ! Buffer size for IO on scratch file
      double precision aa(maxdim,maxdim), ss(maxdim,maxdim), 
     $     yy(maxdim,maxdim), ee(maxdim), fv1(maxdim), fv2(maxdim),
     $     rnorm(maxdim)
      double precision xnorm, ab
      integer i, nsub, info, cur_root, iroot, ivec
      integer icbuf, isbuf, irbuf, itop, ic, is, iw, lenw
      integer ierr, iter
      double precision ddot, fci_time
      external ddot, fci_time
      logical odebug
      data odebug /.false./
c
      if (maxsub .gt. maxdim) call fci_err
     $     ('fci_davids: maxsub > maxdim', maxdim)
      if (maxsub .lt. nroot+2) call fci_err
     $     ('fci_davids: subspace too small', maxsub)
c
c     Scratchfile is organized as v1,H*v1,v2,H*v2,...
c
c     These buffers are used to process the scratch file and compute
c     the reduced space matrices in one pass through the file
c     
      buflen = min(4096,n)
      icbuf = 1
      isbuf = icbuf + maxsub*buflen
      irbuf = isbuf + maxsub*buflen
      itop  = irbuf + buflen
      if (itop .gt. lenwrk) call fci_err
     $     ('fci_davids: need more work', itop)
c
      workfile = 3
      call wopen(workfile, 'workfile', 0, 1, ierr)
      if (ierr .ne. 0) call fci_err
     $     ('fci_davidsnew: failed to open workfile', 0)
c
c     Form initial matrix-vector products
c     
      ic = 1
      is = ic + n
      iw = is + n
      lenw = lenwrk - iw
      if (lenw .le. 0) call fci_err('fci_davids: need more work', is)
      do i = 1, nroot
         call getwa(evecfile, work(ic), (i-1)*n + 1, n, ierr)
         if (ierr.ne.0) call fci_err('fci_davids: IO error (a)',ierr)
         xnorm = 1.0d0/sqrt(ddot(n,work(ic),1,work(ic),1)) ! Normalize it
         call dscal(n, xnorm, work(ic), 1)
c
         call fci_hv(ecore, h, g, work(iw), lenw, work(ic), work(is))
c
         call putwa(workfile, work(ic), 2*(i-1)*n + 1, n, ierr)
         if (ierr.ne.0) call fci_err('fci_davids: IO error (b)',ierr)
         call putwa(workfile, work(is), 2*(i-1)*n+n+1, n, ierr)
         if (ierr.ne.0) call fci_err('fci_davids: IO error (c)',ierr)
      enddo
c
      write(6,1) 
 1    format(
     $     7x, ' iter nsub  cur root   eigenvalue      residual',
     $     '    time'/
     $     7x, ' ---- ---- ---- ---- --------------- -----------',
     $     ' -------')
c     
      test_tol = 0.5d0
      cur_root  = 1             ! Current vector being updated
      nsub = nroot              ! No. of vectors within subspace
      do iter = 1, maxiter
c     
c     Now have nsub vectors and products in workfile.
c     
 555     call fci_davids_redmat(workfile, n, nsub, buflen, 
     $        work(icbuf), work(isbuf), aa, ss, maxdim)
c     
         if (odebug) then
            write(6,*) ' fci_davids: Reduced space H matrix'
            call fci_output(aa, 1, nsub, 1, nsub, maxdim, nsub, 1)
            write(6,*) ' fci_davids: Reduced space S matrix'
            call fci_output(ss, 1, nsub, 1, nsub, maxdim, nsub, 1)
         endif
c     
c     Solve the reduced space equations (Ay = Sye)
c     
         call rsg(maxdim, nsub, aa, ss, ee, 1, yy, fv1, fv2, info)
         if (info .ne. 0) call fci_err('fci_davids: rsg failed', info)
         if (odebug) then
            write(6,*) ' fci_davids: Reduced space eigenvectors'
            call fci_output(yy, 1, nsub, 1, nsub, maxdim, nsub, 1)
            write(6,*) ' fci_davids: Reduced space eigenvalues'
            call fci_output(ee, 1, nsub, 1, 1, nsub, 1, 1)
         endif
         call dcopy(nroot, ee, 1, eval, 1)
c
c     Resum if we must ... not optimal right now
c
         if (nsub .eq. maxsub) then 
            write(6,*) '   Resumming '
            call fci_davids_resum(workfile, 
     $           n, nsub, nroot, buflen, work(icbuf), work(isbuf), 
     $           work(irbuf), yy, maxdim)
            nsub = nroot
            goto 555            ! Regenerate reduced space matrices
         endif
c
c     Form the residuals in evecfile and their norms in rnorm
c
         call fci_davids_residual(evecfile, workfile, 
     $        n, nsub, nroot, buflen, work(icbuf), work(isbuf), 
     $        work(irbuf), ee, yy, maxdim, rnorm)
c     
c     Print
c     
         do iroot = 1, nroot
            write(6,3) iter, nsub, cur_root, iroot, eval(iroot),
     $           rnorm(iroot), fci_time()
 3          format(6x,4i5,2x,f15.8,f12.8,f8.1)
         enddo
*         call util_flush(6)
c     
c     Determine the next vector to update
c     
 20      converged = .true.
         cur_root = 0           ! Comment this out to cycle
         do iroot = 1, nroot
            cur_root = mod(cur_root,nroot) + 1
            if (rnorm(iroot) .gt. test_tol) then
               converged = .false.
               goto 10
            endif
         enddo
c         
 10      if (converged) then
            if (test_tol .le. tol) then
               goto 1000
            else
               test_tol = max(test_tol * 0.01d0, tol)
               goto 20
            endif
         endif
c     
         if (iter.eq.maxiter) then
            write(6,*) ' fci_davids: exceeded max iters ... ',
     $           'aborting solve'
            goto 1000
         endif
c     
c     Read diagonals into work(is) and form update for current root in work(ic)
c     
         call getwa(diagfile, work(is), 1, n, ierr)
         if (ierr.ne.0) call fci_err('fci_davids: IO error (diag)',ierr)
         call getwa(evecfile, work(ic), (cur_root-1)*n+1, n, ierr)
         if (ierr.ne.0) call fci_err('fci_davids: IO error (d)',ierr)
         do i = 1, n
            if (ierr.ne.0) call fci_err('fci_davids: IO error (e)',ierr)
            denom = eval(cur_root) - work(is+i-1)
            if (abs(denom).lt.0.1d0) denom = sign(0.1d0,denom)
            work(ic+i-1) = work(ic+i-1) / denom
         enddo
c
c     Approximately orthogonalize it to the previous vectors (assumed
c     normalized), force symmetries and then normalize it
c
         do ivec = 1, nsub
            call getwa(workfile, work(is), 2*(ivec-1)*n+1, n, ierr)
            ab = ddot(n, work(ic), 1, work(is), 1)
            call daxpy(n, -ab, work(ic), 1, work(is), 1)
         enddo
c
         xnorm = 1.0d0/sqrt(ddot(n,work(ic),1,work(ic),1)) ! Normalize it
         call dscal(n, xnorm, work(ic), 1)
c
         call fci_spin_force(work(ic))
c
c     Form new product and shove vectors them into the work file
c     
         call fci_hv(ecore, h, g, work(iw), lenw, work(ic), work(is))
*         write(6,*) ' NEW VECTOR '
*         call fci_output(work(ic),1,n,1,1,n,1,1)
*         write(6,*) ' NEW H*V '
*         call fci_output(work(is),1,n,1,1,n,1,1)
c
         nsub = nsub + 1
c
         call putwa(workfile, work(ic), 2*(nsub-1)*n+1, n, ierr)
         if (ierr.ne.0) call fci_err('fci_davids: IO error (f)',ierr)
         call putwa(workfile, work(is), 2*(nsub-1)*n+n+1, n, ierr)
         if (ierr.ne.0) call fci_err('fci_davids: IO error (g)',ierr)
c     
      enddo
c     
c     Should never actually fall out of the bottom
c     
      call fci_err('fci_davids: fell out of loop bottom?', 0)
c     
 1000 continue
c     
c     Resum current solution into the output file
c     
      call fci_davids_resum(workfile, 
     $     n, nsub, nroot, buflen, work(icbuf), work(isbuf), 
     $     work(irbuf), yy, maxdim)
c
      do ivec = 1, nroot
         call getwa(workfile, work, 2*(ivec-1)*n+1, n, ierr)
         if (ierr.ne.0) call fci_err
     $        ('fci_davids: failed getting final vector', ierr)
         call putwa(evecfile, work, (ivec-1)*n+1, n, ierr)
         if (ierr.ne.0) call fci_err
     $        ('fci_davids: failed putting final vector', ierr)
      enddo
c         
      call dcopy(nroot, ee, 1, eval, 1)
c
      call wclose(workfile, ierr)
c     
      end
      subroutine fci_davids_redmat(workfile, n, nvec, buflen, 
     $     cbuf, sbuf, aa, ss, maxdim)
      implicit none
c
      integer workfile          ! [input] Handle to workfile
      integer n                 ! [input] Vector length
      integer nvec              ! [input] No. of vectors in the file
      integer buflen            ! [input] IO buffer length
      double precision cbuf(buflen,nvec) ! [scratch] C buffer
      double precision sbuf(buflen,nvec) ! [scratch] S=HC buffer
      integer maxdim            ! [input] Dimension of matrices
      double precision aa(maxdim,maxdim) ! [output] Reduced H matrix
      double precision ss(maxdim,maxdim) ! [output] Reduced S matrix
c
c     Compute the Hamiltonian and Overlap matrices over the Davidson
c     subspace vectors and their products stored in the work file.
c
c     Scratchfile is organized as v1,H*v1,v2,H*v2,...
c
      integer ilo, ihi, ndo, ivec, jvec, ioff, ierr
      double precision ddot
      external ddot
c
      call dfill(0.0d0,aa,1,maxdim*maxdim)
      call dfill(0.0d0,ss,1,maxdim*maxdim)
c
c     Loop thru blocks of the vectors
c
*      write(6,*) ' NVEC,N ', nvec, n
c
      do ilo = 1, n, buflen
         ihi = min(ilo+buflen-1,n)
         ndo = ihi - ilo + 1
c
c     Read blocks in
c
         do ivec = 1, nvec
            ioff = ilo + 2*(ivec-1)*n
            call getwa(workfile,cbuf(1,ivec), ioff, ndo, ierr)
            if (ierr.ne.0) call fci_err('fci_davids: IO error (h)',ierr)
*            write(6,*) ' workfile v ', ivec, ilo
*            call fci_output(cbuf(1,ivec),1,ndo,1,1,ndo,1,1)
            ioff = ioff + n
            call getwa(workfile,sbuf(1,ivec), ioff, ndo, ierr)
*            write(6,*) ' workfile s ', ivec, ilo
*            call fci_output(sbuf(1,ivec),1,ndo,1,1,ndo,1,1)
         enddo
c
c     Accumlate the products
c
         call dgemm('t','n',nvec,nvec,ndo,1.0d0,cbuf,buflen,
     $        sbuf,buflen,1.0d0,aa,maxdim)
         call dgemm('t','n',nvec,nvec,ndo,1.0d0,cbuf,buflen,
     $        cbuf,buflen,1.0d0,ss,maxdim)
*         do ivec = 1, nvec
*            do jvec = 1, nvec
*               aa(ivec,jvec) = aa(ivec,jvec) + 
*     $              ddot(ndo, cbuf(1,ivec), 1, sbuf(1,jvec), 1)
*               ss(ivec,jvec) = ss(ivec,jvec) + 
*     $              ddot(ndo, cbuf(1,ivec), 1, cbuf(1,jvec), 1)
*            enddo
*         enddo
      enddo
c
      end
      subroutine fci_davids_residual(residfile, workfile, 
     $     n, nvec, nroot, buflen, cbuf, sbuf, rbuf,
     $     ee, yy, maxdim, rnorm)
      implicit none
c
      integer residfile         ! [input] Handle to file for residuals
      integer workfile          ! [input] Handle to work file
      integer n                 ! [input] Vector length
      integer nvec              ! [input] No. of vectors in the file
      integer nroot             ! [input] No. of roots 
      integer buflen            ! [input] IO buffer length
      double precision cbuf(buflen,nvec) ! [scratch] C buffer
      double precision sbuf(buflen,nvec) ! [scratch] S=HC buffer
      double precision rbuf(buflen) ! [scratch] To form residual
      integer maxdim            ! [input] Dimension of matrices
      double precision ee(maxdim) ! [input] Eigenvalues
      double precision yy(maxdim,nroot) ! [input] Reduced space eigenvectors
      double precision rnorm(nroot) ! [output] Residuals for each root
c
c     Compute the residual vectors for each root storing it in
c     residfile and also return the residual norms
c
c     r = (H - E)c = s - Ec = sum(i) (si - Eci)*yi
c
      integer ilo, ihi, ndo, iroot, ivec, jvec, ioff, i, ierr
c
      call dfill(0.0d0,rnorm,1,nroot)
c
c     Loop thru blocks of the vectors
c
      do ilo = 1, n, buflen
         ihi = min(ilo+buflen-1,n)
         ndo = ihi - ilo + 1
c
c     Read blocks in
c
         do ivec = 1, nvec
            ioff = ilo + 2*(ivec-1)*n
            call getwa(workfile,cbuf(1,ivec), ioff, ndo, ierr)
            if (ierr.ne.0) call fci_err('fci_davids: IO error (i)',ierr)
            ioff = ioff + n
            call getwa(workfile,sbuf(1,ivec), ioff, ndo, ierr)
         enddo
c
c     Form the residuals 
c
         do iroot = 1, nroot
            call dfill(0.0d0,rbuf,1,ndo)
            do jvec = 1, nvec
               do i = 1, ndo
                  rbuf(i) = rbuf(i) + 
     $                 (sbuf(i,jvec) - ee(iroot)*cbuf(i,jvec)) * 
     $                 yy(jvec,iroot)
               enddo
            enddo
            do i = 1, ndo
               rnorm(iroot) = rnorm(iroot) + rbuf(i)*rbuf(i)
            enddo
            ioff = ilo + (iroot-1)*n
            call putwa(residfile,rbuf, ioff, ndo, ierr)
            if (ierr.ne.0) call fci_err('fci_davids: IO error (j)',ierr)
         enddo
      enddo
c
      do iroot = 1, nroot
         rnorm(iroot) = sqrt(rnorm(iroot))
      enddo
c
*      write(6,*) ' RNORM '
*      call fci_output(rnorm,1,nroot,1,1,nroot,1,1)
c
      end
      subroutine fci_davids_resum(workfile, 
     $     n, nvec, nroot, buflen, cbuf, sbuf, rbuf,
     $     yy, maxdim)
      implicit none
c
      integer workfile          ! [input] Handle to work file
      integer n                 ! [input] Vector length
      integer nvec              ! [input] No. of vectors in the file
      integer nroot             ! [input] No. of roots 
      integer buflen            ! [input] IO buffer length
      double precision cbuf(buflen,nvec) ! [scratch] C buffer
      double precision sbuf(buflen,nvec) ! [scratch] S=HC buffer
      double precision rbuf(buflen) ! [scratch] To form residual
      integer maxdim            ! [input] Dimension of matrices
      double precision yy(maxdim,nroot) ! [input] Reduced space eigenvectors
c
c     Resum the nvec vectors and products on the workfile down 
c     to nroot vectors
c
      integer ilo, ihi, ndo, iroot, ivec, jvec, ioff, i, ierr
c
c     Loop thru blocks of the vectors
c
      do ilo = 1, n, buflen
         ihi = min(ilo+buflen-1,n)
         ndo = ihi - ilo + 1
c
c     Read blocks in
c
         do ivec = 1, nvec
            ioff = ilo + 2*(ivec-1)*n
            call getwa(workfile,cbuf(1,ivec), ioff, ndo, ierr)
            if (ierr.ne.0) call fci_err('fci_davids: IO error (k)',ierr)
            ioff = ioff + n
            call getwa(workfile,sbuf(1,ivec), ioff, ndo, ierr)
         enddo
c
         do iroot = 1, nroot
c
c     Resum CI vector
c
            call dfill(0.0d0,rbuf,1,ndo)
            do jvec = 1, nvec
               do i = 1, ndo
                  rbuf(i) = rbuf(i) + cbuf(i,jvec) * yy(jvec,iroot)
               enddo
            enddo
            ioff = ilo + 2*(iroot-1)*n
            call putwa(workfile,rbuf, ioff, ndo, ierr)
            if (ierr.ne.0) call fci_err('fci_davids: IO error (l)',ierr)
c
c     Resum product vector
c
            call dfill(0.0d0,rbuf,1,ndo)
            do jvec = 1, nvec
               do i = 1, ndo
                  rbuf(i) = rbuf(i) + sbuf(i,jvec) * yy(jvec,iroot)
               enddo
            enddo
            ioff = ilo + 2*(iroot-1)*n + n
            call putwa(workfile,rbuf, ioff, ndo, ierr)
            if (ierr.ne.0) call fci_err('fci_davids: IO error (m)',ierr)
         enddo
      enddo
c
      end
