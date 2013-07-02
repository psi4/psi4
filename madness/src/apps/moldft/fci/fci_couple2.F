*
* $Id: fci_couple2.F,v 1.2 1999/07/28 00:23:39 d3e129 Exp $
*
      subroutine fci_couple2(vec, nvec, p, r, iocc, idim, 
     $     numqs, qslist, epqrs)
      implicit none
#include "fciinfo.fh"
c     
c     Compute <V|Epq,rs|Iw> for given CI vectors V, fixed orbitals
c     p and r, fixed orbital occupancy I (iocc), and all orbitals q and s,
c     and primitive spin functions w.
c
c     The constraint (pq) >= (rs) is applied, so we loop over all q but
c     restrict the s loop.  If couple2 is called with p>=r then all unique
c     coeffs are generated.
c
c     numqs returns the number of orbitals q, s giving a nonzero interaction
c     
c     qslist(1, 1...numqs) contains the values of q
c     qslist(2, 1...numqs) contains the values of s
c     
c     epqrs(1...idim,1...nvec,1...numqs) contains the coupling coeffcients
c     
      integer nvec              ! [input] No. of vectors on the left hand side
      double precision vec(ndets,nvec) ! [input] The CI vectors
      integer p, r              ! [input] Fixed labels on Epq,rs
      integer iocc(nactiv)      ! [input] Orbital occupation of I
      integer idim              ! [input] No. of spin functions for I
      integer numqs             ! [output] Returns the number of q, s
      integer qslist(2,*)       ! [output] Returns the list of interacting q,s
      double precision epqrs(idim,nvec,*) ! [output] epq(w,i,q) = <Vi|Epqrs|Iw>
c     
c     <V|Epq,rs|Iw> = sum(Jv) V(Jv) <Jv|EpaErbEbsEaq|Iw> 
c     .  = sum(Jx) V(Jv) <Jv|EpaErb|Mx><Mx|EbsEaq|Iw>
c     .  = sum(Jx) V(Jv) (<Jv|Epa|Ny><Ny|Erb|Mx>)(<Mx|Ebs|Ku><Ku|Eaq|Iw>)
c     .  = sum(Jx) V(Jv) (<Ny|Eap|Jv><Mx|Ebr|Ny>)(<Mx|Ebs|Ku><Ku|Eaq|Iw>)
c
c     Where we have v(x), y(x), u(x), w(x)
c
c     where a and b are fictitious orbitals with b > a
c     
c     The strategy is then to
c     
c     1) Loop over q occupied in I 
c     2) Form <Ku|Eaq|Iw>
c     3) Determine occupations in K
c     4) Loop over s occupied in K
c     5) Determine occupations in M
c     6) Form <Mx|Ebs|Ku><Ku|Eaq|Iw>
c     7) Determine occupations in N
c     8) Determine occupations in J 
c     9) Form <Ny|Eap|Jv>
c     10) Form <Ny|Eap|Jv><Mx|Ebr|Ny>
c     11) For each vector V complete the final summation over x
c     .   noting that we have v(x) and w(x).  If the result is non-zero
c     .   then add it in.
c
c     Defer as much computation as possible until we have satisfied
c     all possible occupation and sparsity tests.
c     
      integer indbar(nactmx)    ! Position amoung singly occupied I/K/M/N/J
      integer kperm(ndetmx), jperm(ndetmx), 
     $     nperm(ndetmx), iperm(ndetmx) ! Permutations of spin functions
      double precision kval(ndetmx), nval(ndetmx),
     $     jval(ndetmx), ival(ndetmx) ! Accumulated phase factors
c     
      integer nsi, nsj, nsk, nsm, nsn ! No. of unpaired spins in I/J/K/M/N
      integer jdim, kdim, mdim, ndim ! No. of determinants for J/K/M/N
      integer np                ! Occuption of p in J
      integer nq                ! Occupation of q in I
      integer ns                ! Occupation of s in K
      integer nr                ! Occupation of r in N
      integer jwalk             ! Lexical index of J
      integer jdet              ! Offset of J in CI vector
      integer sympr             ! Sym product  state*I*p*r
      integer syms              ! Required symmetry of s
      integer q, s, w, x, v, pbar, qbar, rbar, sbar, shi,
     $     deltaq, deltar, deltas, idummy, i, ivec
      double precision sum, dummy
      logical omade_ival        ! Used to defer computation
      double precision thresh
      parameter (thresh = 1d-12) ! Should be SQUARE of desired accuracy
#include "bitops.fh"      
c     
      call fci_tstart(2,.false.)
      numqs = 1
c     
c     Find I occupancy, compute no. of open shells and I indbar vector
c     
      nsi  = 0
      sympr = ieor(statesym,ieor(orbsym(p),orbsym(r)))
      do i = 1,nactiv
         if (iocc(i).eq.2) then
            nsi = nsi + 1
            sympr = ieor(sympr,orbsym(i))
         endif
         indbar(i) = nsi
      enddo
      if (idim .ne. mnsdet(nsi)) call fci_err
     $     ('fci_couple2: mismatch for idim ', idim)
c     
c     Loop thru occupied q.
c     
      do q = 1, nactiv
         nq = iocc(q)-1
         if (nq .gt. 0) then
            iocc(q) = iocc(q) - 1 ! iocc now refers to K
c
c     Get info about q and K
c     
            qbar = indbar(q)
            nsk  = nsi
            deltaq = -1
            if (nq .eq. 2) then
               deltaq = 1
               nsk = nsk + 2
               qbar = qbar + 1
            endif
            kdim = mnsdet(nsk)
c
c     Defer forming <Ku|Eaq|Iw> until we know it will be used
c
            omade_ival = .false.
c     
c     Loop thru s occupied in K with condition (pq) <= (rs)
c
            syms = ieor(orbsym(q),sympr)
            shi = nactiv
            if (p .eq. r) shi = q
            do s = 1, shi
               ns = iocc(s)-1   ! Occupation of s in K
               if (ns.gt.0 .and. orbsym(s).eq.syms) then
                  iocc(s) = iocc(s) - 1 ! iocc now refers to M
c
                  nr = iocc(r)  ! Since iocc refers to M and we want N
                  if (nr .lt. 3) then
                     iocc(r) = iocc(r) + 1 ! iocc now refers to N
c     
                     np = iocc(p) ! Since iocc refers to N and we want J
                     if (np .lt. 3) then
                        iocc(p) = iocc(p) + 1 ! iocc now refers to J
c
                        nsm  = nsk
                        if (ns .eq. 2) nsm = nsm + 2
c
                        nsn = nsm
                        if (nr .eq. 2) nsn = nsn - 2
c
                        nsj = nsn
                        if (np .eq. 2) nsj = nsj - 2 
c     
c     Determine the index of J in the CI vector and apply sparsity test
c
                        jdim = mnsdet(nsj)
                        call fci_ofind_weight(jwalk,iocc)
                        jdet = indexo(jwalk)
c
                        sum = 0.0d0
                        do ivec = 1, nvec
                           do v = 1, jdim
                              sum = sum + 
     $                             vec(jdet+v,ivec)*vec(jdet+v,ivec)
                           enddo
                        enddo
                        if (sum .gt. thresh) then
c
c     More detailed info about s and M
c
                           deltas = -1
                           sbar = indbar(s)
                           if (ns .eq. 2) then
                              sbar = sbar + 1
                              deltas = 1
                           endif
                           mdim = mnsdet(nsm)
                           if (s .ge. q) sbar = sbar + deltaq
c
c     More detailed info about r and N
c
                           deltar = 1
                           rbar = indbar(r)
                           if (nr .eq. 2) then
                              deltar = -1
                              rbar = rbar + 1
                           endif
                           ndim = mnsdet(nsn)
                           rbar = rbar + deltar
                           if (r .ge. q) rbar = rbar + deltaq
                           if (r .ge. s) rbar = rbar + deltas
c
c     More detailed info about p and J
c
                           pbar = indbar(p)
                           if (np .eq. 1) pbar = pbar + 1 ! Indbar refers to N
                           if (p .ge. q) pbar = pbar + deltaq
                           if (p .ge. s) pbar = pbar + deltas
                           if (p .ge. r) pbar = pbar + deltar
c
c     If not already done form <Ku|Eaq|Iw>
c
                           if (.not. omade_ival) then
                              call fci_apply_eaj(qbar, nq, idim, kdim, 
     $                             eajdata(eajptr(nsi)), 
     $                             eajphase(1,1,nsi), .true.,
     $                             idummy, dummy, iperm, ival)
                              omade_ival = .true.
                           endif
c     
c     Form <Mx|EbsEaq|Iw> = <Mx|Ebs|Ku><Ku|Eaq|Iw>
c     
                           call fci_apply_eaj(sbar, ns, kdim, mdim, 
     $                          eajdata(eajptr(nsk)), eajphase(1,1,nsk),
     $                          .false., iperm, ival, kperm, kval)
c     
c     Form <Ny|Eap|Jv>
c     
                           call fci_apply_eaj(pbar, np, jdim, ndim, 
     $                          eajdata(eajptr(nsj)), eajphase(1,1,nsj),
     $                          .true., idummy, dummy, jperm, jval)
c     
c     Form <Jv|EpaErb|Mx> = <Ny|Eap|Jv><Mx|Ebr|Ny>
c     
                           call fci_apply_eaj(rbar, nr, ndim, mdim, 
     $                          eajdata(eajptr(nsn)), eajphase(1,1,nsn),
     $                          .false., jperm, jval, nperm, nval)
c     
c     Now use <V|Epq,rs|Iw> = sum(x) V(Jv) <Mx|EapEbr|Ny><Mx|EbsEaq|Iw>
c     .
c     .   noting that we have v(x) and w(x)
                           
c     
                           do ivec = 1, nvec
                              do w = 1, idim
                                 epqrs(w,ivec,numqs) = 0.0d0
                              enddo
                           enddo
c     
                           do ivec = 1, nvec
                              do x = 1, mdim
                                 v = nperm(x)
                                 w = kperm(x)
                                 epqrs(w,ivec,numqs) = 
     $                                epqrs(w,ivec,numqs) + 
     $                                vec(jdet+v,ivec)*nval(x)*kval(x)
                              enddo
                           enddo
                           qslist(1,numqs) = q
                           qslist(2,numqs) = s
                           numqs = numqs + 1
                        endif   ! sum.gt.thresh
c
                        iocc(p) = iocc(p) - 1 ! iocc now refers to N
                     endif      ! np.lt.3
c
                     iocc(r) = iocc(r) - 1 ! iocc now refers to M
                  endif         ! nr.lt.3
c
                  iocc(s) = iocc(s) + 1 ! iocc now refers to K
               endif            ! ns.gt.0
            enddo               ! do s
c     
            iocc(q) = iocc(q) + 1 ! iocc now refers to I
         endif                  ! nq.gt.0
      enddo                     ! do q
c     
      numqs = numqs - 1
c     
      call fci_tfinish(2,.false.)
c
      end
      
