*
* $Id: fci_setup.F,v 1.2 1999/07/28 00:23:42 d3e129 Exp $
*
      subroutine fci_setup(multiplicity, nactive, nelectrons,
     $     orbsyms, symstate,
     $     norbitalconfig, ndeterminants, nconfigurations,
     $     norbconfignosym)
      implicit none
#include "fciinfo.fh"
      integer multiplicity      ! [input] Spin multiplicity
      integer nactive           ! [input] No. of active orbitals
      integer nelectrons        ! [input] No. of active electrons
      integer orbsyms(nactive)  ! [input] Symmetry of orbitals (1...8)
      integer symstate          ! [input] Desired symmetry of state
      integer norbitalconfig    ! [output] No. of orb. conf. WITH symmetry
      integer ndeterminants     ! [output] No. of determinants in FCI
      integer nconfigurations   ! [output] No. of configurations in FCI
      integer norbconfignosym   ! [output] No. of orb. conf. WITHOUT symmetry
c     
c     initialise tables and common blocks for fullci of
c     nelectrons electrons in nactive orbitals and multiplicity=2s+1.
c
c     norbitalconfig, ndeteterminants and ncconfiguration return the 
c     no. of orbital configurations, the no. of determinants in the ci,
c     and the no. of csf in the ci respectively
c
c     symmetry information is accepted and enforced but is not currently
c     used to reduce computation or memory usage ... coming soon!
c     
      integer iocc(nactmx)
      integer i, iwalk, ns, junk, idim, isym
      logical oprint
c
      call fci_tinit()
c
      oprint = .false.
      call fci_tstart(7,oprint)
c
c     Put arguments into common and then initialize everthing else
c
      nactiv = nactive
      nelec  = nelectrons
      multi  = multiplicity
      statesym = symstate - 1
      do i = 1, nactiv
         orbsym(i) = orbsyms(i)-1
         if (orbsym(i).lt.0 .or. orbsym(i).gt.7) call fci_err
     $        ('fci_setup: orbital symmetry is bad ', orbsyms(i))
      enddo
c
      if (statesym.lt.0 .or. statesym.gt.7) call fci_err
     $     ('fci_setup: state symmetry is bad',symstate)
      if (nelec.lt.1 .or. nelec.gt.nelmax) call fci_err(
     & ' invalid no. of electrons in full ci space ',nelec)
      if (nactiv.lt.1 .or. nactiv.gt.nactmx) call fci_err(
     & ' invalid no. of orbitals in full ci space ',nactiv)
      if (multi.lt.1 .or. multi.gt.nelec+1) call fci_err(
     & ' invalid spin multiplicity in full ci space ',multi)
c
      nacttr = nactiv*(nactiv+1)/2
c     
      do i = 1,nactiv
         iky(i) = ((i-1)*i)/2
      enddo
c
c     make the spin graphs
c
      call fci_mgraph
c     
c     now the orbital graph for nelec and nelec-1 electrons
c
      call fci_ograph
c
      if (nwalko.lt.1 .or. nwalko.gt.nwkomx) call fci_err(
     &     ' invalid no. of n e- orbital configurations ',nwalko)
c
c     count the determinants and csf, make indexo vector
c
      if (oprint) write(6,*) ' orbital configurations '
      ncsf = 0
      ndets = 0
      nwalkosym = 0
c
      do 10 iwalk = 1,nwalko
         call fci_owalk_info(iwalk,iocc,ns,junk,idim,isym)
         if (isym .eq. statesym) then
            if (oprint) write(6,1) iwalk,mnsdet(ns),mnscsf(ns),
     &           (iocc(i)-1,i=1,nactiv)
 1          format(1x,i6,4x,2i4,4x,10i2)
            indexo(iwalk) = ndets
            ndets = ndets + mnsdet(ns)
            ncsf = ncsf + mnscsf(ns)
            nwalkosym = nwalkosym + 1
         else
            indexo(iwalk) = ndets ! Must do this for detwalk
         endif
 10   continue
      indexo(nwalko+1) = ndets
      if (oprint) then
         write(6,2) (indexo(i),i=1,nwalko+1)
 2       format(/' index from orbital configurations to ci vector'//
     &        100(1x,13i6/))
      endif
c     
c     Make the Eaj information
c
      call fci_make_all_eaj()
c     
      call fci_tfinish(7,.false.)
c
      norbitalconfig  = nwalkosym
      ndeterminants   = ndets
      nconfigurations = ncsf
      norbconfignosym = nwalko
c
      end
