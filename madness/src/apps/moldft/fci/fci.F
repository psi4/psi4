*
* $Id: fci.F,v 1.2 1999/07/28 00:23:38 d3e129 Exp $
*
      program fci
      implicit none
      integer nactiv, nelec, multi
      integer nwalko, ndets, ncsf, nroot
      integer maxact, maxtri, maxroot, maxhami
      parameter (maxact=14, maxroot=3, maxhami=441)
      parameter (maxtri = maxact*(maxact+1)/2)
      double precision ecore, h(maxtri), g(maxtri*maxtri)
      double precision thresh, prthresh, weight(maxroot)
      integer lenwrk, iwork, iv, lwork, iroot, orbsyms(maxact),
     $     symstate, nwalkosym
      parameter (lenwrk = 16000000)
      double precision work(lenwrk)
c
      character*9 evecfilename, diagfilename
c
      double precision hami(maxhami*(maxhami+1)/2), evals(maxhami)
      integer pivot(maxhami), i
c
      double precision onepdm(maxact*maxact*maxroot*maxroot)
      double precision twopdm(maxact*maxact*maxact*maxact*
     $     maxroot*maxroot)
      double precision roots(maxroot)
      data weight/maxroot*1.0D0/
      data evecfilename, diagfilename /'evecs.dat', 'diags.dat'/
c
c     Simple test driver program for the FCI 
c
c     The first line of the input is
c     ..   nactive, nelectrons, multiplicity, symmetry, nroot, thresh, prtol
c
c     The second is the list of orbital symmetries
c
c     All subsequent lines are read as
c     ..   g, i, j, k, l 
c
c     if (k>0) then 
c     ..   g is (ij|kl)
c     else if (i>0) then
c     ..   g is hij (including frozen core contribution)
c     else
c     ..   g is the nuclear repulsion + frozen core energy
c     endif
c     
      read (5,*) nactiv,nelec,multi, symstate, nroot, thresh, prthresh
      read(5,'(20i2)') (orbsyms(i),i=1,nactiv)
      if (nroot .gt. maxroot) nroot = maxroot ! Dimension of roots and 1pdm
c
c initialise tables and common blocks for ci and onelcc
c
      call fci_setup(multi, nactiv, nelec, 
     $     orbsyms, symstate, nwalkosym, ndets, ncsf, nwalko)
c
      write(6,1) nactiv,nelec,multi,symstate,nwalkosym,ndets,ncsf,
     $     (orbsyms(i),i=1,nactiv)
1     format(/11x,'small full configuration interaction program'/
     &        11x,'--------------------------------------------'//
     &       ' no. of active orbitals .......... ',i8/
     &       ' no. of active electrons ......... ',i8/
     &       ' spin multiplicity ............... ',i8/
     &       ' state symmetry .................. ',i8/
     &       ' no. of orbital configurations ... ',i8/
     &       ' no. of determinants in full ci .. ',i8/
     &       ' no. of configurations in full ci  ',i8/
     &       ' orbital symmetries .............. ',20i2)
      write(6,*)
c
c     Read the integrals from the test ascii input file
c
      call fci_load(h, g, ecore)
c
c$$$      iv    = 1
c$$$      iwork = iv + nroot*ndets
c$$$      call dfill(0.0d0, work(iv), 1, ndets*nroot)
c$$$      lwork = lenwrk - nroot*ndets
c$$$c
c$$$c     With a zero input vector it will use a default guess
c$$$c
c$$$      call fci_solve(nroot, thresh, 32, ecore, h, g, 
c$$$     $     work(iwork), lwork, roots, work(iv))
c
      call fci_solve(nroot, thresh, 32, ecore, h, g, 
     $     work, lenwrk, roots, evecfilename, diagfilename)
c
c     Eigenvectors now on disk, eigenvalues in array roots()
c
c     Allocate memory to hold all of the vectors, read them in
c     and analyze them.
c
      iv    = 1
      iwork = iv + nroot*ndets
      lwork = lenwrk - nroot*ndets
c
      call fci_readevectors(ndets, nroot, work(iv), evecfilename)
c
      do iroot = 1, nroot
         write(6,*)
         write(6,*) ' Eigenvector number ', iroot
         call fci_print_vector(work(iv+(iroot-1)*ndets), 
     $        prthresh, nwalko, nactiv, symstate-1)
      enddo
c
      call fci_onepdm(work(iv), nroot, weight, onepdm)
      call fci_1pdm(work(iv), nroot, onepdm)
      call fci_2pdm(work(iv), nroot, twopdm, onepdm, 
     $     ecore, h, g)
*      call fci_3pdm(work(iv), nroot, work(iwork), twopdm)
c
c     For a small FCI make the full Hamiltonian and diagonalize it
c
      if (ndets .le. maxhami) then
         call fci_hami(hami, ecore, h, g, work, lenwrk)
         if (lenwrk.lt.maxhami*(maxhami+2)) 
     $        call fci_err('more memory for hami ', lenwrk)
         call yacobi(hami, work, ndets, evals, pivot, 1)
         write(6,*)
         write(6,*) ' Eigenvalues of the determinant Hamiltionian '
         write(6,1101) (i,evals(i),i=1,ndets)
 1101    format(1x,i5,2x,f20.8)
         write(6,*)
         write(6,*) ' Eigenvectors of the determinant Hamiltonian '
         do i = 1, min(10,ndets)
            write(6,*)
            write(6,*) ' Eigenvector number ', i, evals(i)
            write(6,*)
            call fci_print_vector(work(1+(i-1)*ndets), prthresh,
     $           nwalko, nactiv, symstate-1)
            write(6,*)
         enddo
      endif
c
      call fci_tprint()
c
      end
      subroutine fci_print_vector(c, thresh, nwalko, nactive, statesym)
      implicit none
#include "fciparam.fh"
c      
c     Print large components of the CI vector out
c
      double precision c(*)   ! [input] The determinant CI vector
      double precision thresh ! [input] Print coeffs > thresh
      integer nwalko          ! [input] No. of orbital occupations
      integer nactive         ! [input] No. of active orbitals
      integer statesym        ! [input] Symmetry of the CI state (0...7)
c
      integer iwalk, mwalk, ns, ioff, numdet, i, testo, testm, sym
      integer iocc(nactmx), mocc(nsmax)
      character*80 fmt
      character*2 spin(2)
      data spin/'+ ','- '/
c
      write(6,1)
 1    format(/
     $     '    Index     Owalk Mwalk      Coeff       Orbital and',
     $     ' Spin occupation'/
     $     ' --------    ------ -----    ---------     -----------',
     $     '----------------')
      do iwalk = 1, nwalko
         call fci_owalk_info(iwalk, iocc, ns, ioff, numdet,sym)
         if (sym .eq. statesym) then
            write(fmt,12) nactive, max(ns,1)
 12         format('(1x,i8,4x,i6,2x,i4,4x,f9.6,4x,',i2,
     $           'i2,4x,',i2,'a2)')
            do mwalk = 1, numdet
               if (abs(c(ioff+mwalk)).gt.thresh) then
                  call fci_mfind_occ(mwalk, mocc, ns)
                  write(6,fmt) ioff+mwalk,iwalk,mwalk,
     $                 c(ioff+mwalk),
     $                 (iocc(i)-1,i=1,nactive),
     $                 (spin(mocc(i)),i=1,ns)
c     
c     This just to check that the det_to_walk routine is working
c     
                  call fci_det_to_walk(ioff+mwalk, testo, testm)
                  if (testo.ne.iwalk .or. testm.ne.mwalk)
     $                 call fci_err('fci_det_to_walk failed',0)
               endif
            enddo
         endif
      enddo
c
      write(6,*)
c
      end
