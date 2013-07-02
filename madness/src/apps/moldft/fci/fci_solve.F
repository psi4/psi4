*
* $Id: fci_solve.F,v 1.2 1999/07/28 00:23:42 d3e129 Exp $
*
      subroutine fci_solve(nroot, thresh, maxvec, ecore, h, g, 
     $     work, lenwrk, roots, evecfilename, diagfilename)
      implicit none
#include "fciinfo.fh"
c     
c     Solve the CI to the requested accuracy.  Right now just uses
c     determinants but is readily extended to also accomodate CSF.
c     
c     The eigenvalues are returned in roots.  The eigenvectors are
c     returned in file evecfilename.
c
c     The minimum amount of memory required is independent of the
c     number of roots, and is slightly larger than 2*ndeterminants.
c     The more memory the better, but 2*ndeterminants+500000 will
c     probably always be enough.
c     
      integer nroot             ! [input] No. of roots to find
      double precision thresh   ! [input] Threshold for residual norm
      integer maxvec            ! [input] Davidson subspace dim PER ROOT
      double precision ecore    ! [input] Nuclear repulsion+frozen core
      double precision h(nacttr) ! [input] h(ij) = h_{ij}
      double precision g(nacttr,nacttr) ! [input] g(ij,kl) = (ij|kl)
      integer lenwrk            ! [input] Size of work space
      double precision work(lenwrk) ! [scratch] Work space
      double precision roots(nroot) ! [output] Returns the eigenvalues
      character*(*) evecfilename ! [input] Name of file for eigen vectors
      character*(*) diagfilename ! [input] Name of file for diagonals
c     
      integer maxv              ! Local copy of maxvec so we can reduce it
      logical oprint
      integer itermx            ! Maximum number of iterations
      integer iv, ihd, itop2,  ionel, kpt, mkpt, ihmod
      integer ierr, nonelt, nonemx
c
      logical fci_davids
      external fci_davids
c
      integer evecfile, diagfile
c
      data itermx / 200/
c     
      oprint = .true.
      call fci_tstart(9,oprint)
c
      call fci_intchk(h, g)
c
c     Open files
c
      evecfile = 1
      diagfile = 2
      call wopen(evecfile, evecfilename, 0, 0, ierr)
      if (ierr .ne. 0) call fci_err
     $     ('fci_davids: failed to open evecfile', ierr)
      call wopen(diagfile, diagfilename, 0, 0, ierr)
      if (ierr .ne. 0) call fci_err
     $     ('fci_davids: failed to open diagfile', ierr)
c
c     Allocate temporary work space for initial guess and diags
c
      call fci_cntone(nonelt, nonemx)
      ihd = 1
      ihmod = ihd + ndets             ! h + 1/2 sum(k) (ik|kj)
      ionel = ihmod + nacttr
      kpt   = ihmod + (nonemx*5)/nav + 1 ! For computing diagonals
      mkpt  = kpt + nonemx/nav + 1 ! For computing diagonals
      itop2 = mkpt + ndetmx/nav + 1
      if (itop2 .gt. lenwrk+1) call fci_err
     $     ('fci_davids: insufficient memory for diagonals',itop2)
c
c     Modify the input 1-e integrals 
c
      call fci_hmod(h,g,work(ihmod))
c     
c     make the diagonals of the hamiltonian
c     
      call fci_diags(work(ihd),ecore,work(ihmod),g,work(ionel),
     $     nonemx,work(kpt), work(mkpt))
c
      call putwa(diagfile, work(ihd), 1, ndets, ierr)
      if (ierr.ne.0) call fci_err
     $     ('fci_davids: failed to write diagonals',ierr)
c
c     Form initial guess from the lowest energy orbital occupations
c     of the correct symmetry
c
      iv = ihd + ndets
      itop2 = iv + ndets
      if (itop2 .gt. lenwrk+1) call fci_err
     $     ('fci_davids: insufficient memory for guess',itop2)
      call fci_guess(evecfile, work(iv), work(ihd), nroot)
c
c     No memory is in use now.
c
      maxv = maxvec
      maxv = min(maxv,itermx,ncsf+1) ! Logic for small problems
c
      if (.not. fci_davids(ndets, nroot, itermx, maxv, thresh,
     $     work, lenwrk, evecfile, diagfile, roots, ecore, h, g))
     $     call fci_err('fci_davids failed', 0)
c     
      call wclose(evecfile,ierr)
      call wclose(diagfile,ierr)
c
      call fci_tfinish(9,oprint)
c     
      end
      subroutine fci_guess(evecfile, vec, diags, nroot)
      implicit none
#include "fciinfo.fh"
      integer evecfile          ! [input] Handle to vector file
      integer nroot             ! [input]
      double precision vec(ndets) ! [scratch]
      double precision diags(ndets) ! [input] Diagonal energies
c
c     Form an initial guess for nroot CI vectors of the correct symmetry.
c
c     Loop thru orbital configurations and determine the nroot
c     lowest energy orbital occupations.  Set the lowest energy 
c     determinant of each orbital occupation to 1 and then
c     spin project to get nroot orthogonal vectors of the correct
c     symmetry.
c
      integer iroot, iwalk, idone(50), jroot, sym, i, ilowest
      integer iocc(nactmx), ns, ioff, idim, iwalklo, ierr
      double precision elowest
#include "bitops.fh"
c
      if (nroot .gt. 50) call fci_err('fci_guess: nroot too big',nroot)
c
      write(6,*)
      do iroot = 1, nroot
         ilowest = -1
         elowest = 9999999.0d0
         do 10 iwalk = 1, nwalko
            do jroot = 1, iroot-1
               if (idone(jroot).eq.iwalk) goto 10
            enddo
            call fci_owalk_info(iwalk, iocc, ns, ioff, idim, sym)
            if (sym .eq. statesym) then
               do i = 1, idim
                  if (diags(i+ioff).lt.elowest) then
                     ilowest = i+ioff
                     iwalklo = iwalk
                     elowest = diags(ilowest)
                  endif
               enddo
            endif
 10      continue
         if (ilowest .le. 0) call fci_err
     $        ('fci_guess: insufficient independent states',nroot)
         call dfill(0.0d0, vec, 1, ndets)
c
         write(6,11) iroot, ilowest, elowest
 11      format(' Guess for root ', i2, ' is determinant ', 
     $        i8, ' energy ', f12.6)
         vec(ilowest) = 1.0d0
         call fci_spin_force(vec)
         call putwa(evecfile, vec, (iroot-1)*ndets + 1, ndets, ierr)
         if (ierr.ne.0) call fci_err
     $        ('fci_guess: failed writing vector ', ierr)
         idone(iroot) = iwalklo
      enddo
      write(6,*)
c
      end
      subroutine fci_readevectors(ndets, nvec, v, filename)
      implicit none
      integer ndets             ! [input] Vector length
      integer nvec             ! [input] No. 
      double precision v(ndets,nvec)
      character*(*) filename
c
      integer ierr, file
c
      file = 1
      call wopen(file, filename, 0, 0, ierr)
      if (ierr .ne. 0) call fci_err
     $     ('fci_readvectors: failed to open file', ierr)
      call getwa(file, v, 1, nvec*ndets, ierr)
      if (ierr .ne. 0) call fci_err
     $     ('fci_readvectors: failed to read file', ierr)
      call wclose(file,ierr)
c
      end
