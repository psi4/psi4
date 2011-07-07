*
* $Id: fci_couple1.F,v 1.2 1999/07/28 00:23:39 d3e129 Exp $
*
      subroutine fci_couple1(vec, nvec, p, iocc, idim, 
     $     numq, qlist, epq)
      implicit none
#include "fciinfo.fh"
c     
c     Compute <V|Epq|Iw> for given CI vectors V, fixed 
c     orbital p, fixed orbital occupancy I (iocc), and all orbitals q
c     and primitive spin functions w.
c     
c     numq returns the number of orbitals q giving a nonzero interaction
c     
c     qlist(1...numq) contains the values of q
c     
c     epq(1...idim,1...nvec,1...numq) contains the coupling coeffcients
c     
      integer nvec              ! [input] No. of vectors on the left hand side
      double precision vec(ndets,nvec) ! [input] The CI vectors
      integer p                 ! [input] Fixed label on Epq
      integer iocc(nactiv)      ! [input] Orbital occupation of I
      integer idim              ! [input] No. of spin functions for I
      integer numq              ! [output] Returns the number of q
      integer qlist(*)          ! [output] Returns the list of interacting q
      double precision epq(idim,nvec,*) ! [output] epq(w,i,q) = <Vi|Epq|Iw>
c     
c     <V|Epq|Iw> = sum(Jv) V(Jv) <Jv|Epq|Iw> 
c     .          = sum(Jv) V(Jv) sum(u) <Jv|Epa|Ku><Ku|Eaq|Iw>
c     .          = sum(Jv) V(Jv) sum(u) <Ku|Eap|Jv><Ku|Eaq|Iw>
c     
c     The strategy is then to
c     
c     1) Loop over q occupied in I 
c     2) Form <Ku|Eaq|Iw>
c     3) Determine K
c     4) Determine occupations in J 
c     5) Form <Ku|Eap|Jv>
c     6) For each vector V <V|Epq|Iw> = sum(u) <Ku|Eap|Jv><Ku|Eaq|Iw>
c     .  noting that v = v(u).  If result is non-zero then add it in.
c     
      integer indbar(nactmx)    ! Position amoung singly occupied for I/K
      integer jperm(ndetmx), iperm(ndetmx) ! Permutations of spin functions
      double precision jval(ndetmx), ival(ndetmx) ! Phase factors
c     
      integer nsi, nsj, nsk     ! No. of unpaired spins in I, J, and K
      integer jdim, kdim        ! No. of determinants for J and K
      integer np, nq            ! Occuption of p in J and q in I
      integer jwalk             ! Lexical index of J
      integer jdet              ! Offset of J in CI vector
      integer symi              ! Symmetry of I
      integer symq              ! Required symmetry for q
      integer q, w, u, v, pbar, qbar, deltaq, idummy, i, qlo, qhi, ivec
      double precision sum, value
#include "bitops.fh"      
c     
      call fci_tstart(1,.false.)
      numq = 1
c     
c     Find I occupancy, compute no. of open shells and I indbar vector
c     
      nsi = 0
      symi = 0
      do i = 1,nactiv
         if (iocc(i).eq.2) then
            nsi = nsi + 1
            symi = ieor(symi,orbsym(i))
         endif
         indbar(i) = nsi
      enddo
      if (idim .ne. mnsdet(nsi)) call fci_err
     $     ('fci_couple1: mismatch for idim ', idim)
c
      symq = ieor(ieor(symi,orbsym(p)),statesym)
c     
c     Loop thru occupied q.  Note that if p is doubly occupied in
c     I then the only non-zero value is for q=p
c     
      if (iocc(p) .eq. 3) then
         qlo = p
         qhi = p
      else
         qlo = 1
         qhi = nactiv
      endif
c     
      do q = qlo, qhi
         nq = iocc(q)-1
         if (nq.gt.0 .and. orbsym(q).eq.symq) then
            iocc(q) = iocc(q) - 1 ! iocc now refers to K
c
            np = iocc(p)        ! Since iocc refers to K and we want J
            if (np .lt. 3) then
               iocc(p) = iocc(p) + 1 ! iocc now refers to J
c
c     Form <Ku|Eaq|Iw>
c     
               qbar = indbar(q)
               nsk  = nsi
               deltaq = -1
               if (nq .eq. 2) then
                  nsk = nsk + 2
                  qbar = qbar + 1
                  deltaq = 1
               endif
               kdim = mnsdet(nsk)
               call fci_apply_eaj(qbar, nq, idim, kdim, 
     $              eajdata(eajptr(nsi)), eajphase(1,1,nsi), .true.,
     $              idummy, idummy, iperm, ival)
c     
c     Determine info about J orbital occupation
c     
               call fci_ofind_weight(jwalk,iocc)
               pbar = indbar(p)
               if (np .eq. 1) pbar = pbar + 1 ! Since indbar refers to K
               if (p .ge. q) pbar = pbar + deltaq
               nsj = nsk
               if (np .eq. 2) nsj = nsj - 2 
               jdim = mnsdet(nsj)
c     
c     Form <Ku|Eap|Jv>
c     
               call fci_apply_eaj(pbar, np, jdim, kdim, 
     $              eajdata(eajptr(nsj)), eajphase(1,1,nsj), .true.,
     $              idummy, idummy, jperm, jval)
c     
c     Now use <Jv|Epq|Iw> = <Ku|Eap|Jv>*<Ku|Eaq|Iw>
c     
c     noting w(u)=iperm(u), v(u)=jperm(u), <Ku|Eap|Jv>=jval(u) and
c     .      <Ku|Eaq|Iw>=ival(u)
c     
               do ivec = 1, nvec
                  do w = 1, idim
                     epq(w,ivec,numq) = 0.0d0
                  enddo
               enddo
c     
               jdet = indexo(jwalk)
               sum = 0.0d0
               do ivec = 1, nvec
                  do u = 1, kdim
                     v = jperm(u)
                     w = iperm(u)
                     value = vec(jdet+v,ivec)*jval(u)*ival(u)
                     epq(w,ivec,numq) = epq(w,ivec,numq) + value
                     sum = sum + value*value
                  enddo
               enddo
               if (abs(sum) .gt. 0.0d0) then
                  qlist(numq) = q
                  numq = numq + 1
               endif
c
               iocc(p) = iocc(p) - 1 ! iocc now refers to K
            endif
c     
            iocc(q) = iocc(q) + 1 ! iocc now refers to I
         endif
      enddo
c     
      numq = numq - 1
c     
      call fci_tfinish(1,.false.)
c
      end
      
