*
* $Id: fci_couple3.F,v 1.2 1999/07/28 00:23:39 d3e129 Exp $
*
      subroutine fci_couple3(vec, nvec, p, q, r, t, iocc, idim, 
     $     numsu, sulist, epqrstu)
      implicit none
#include "fciinfo.fh"
c     
c     Compute <V|Epq,rs,tu|Iw> for given CI vectors V, fixed orbitals
c     p,q,r,t and fixed orbital occupancy I (iocc), and all orbitals s,u
c     and primitive spin functions w.
c
c     Constrain the loops over s and u so that (pq)>=(rs)>=(tu).  If the
c     calling program externally constrains p>=r>=t (all q) then all
c     unique coupling coefficients will be generated.
c     
c     numsu returns the number of orbitals s,u giving a nonzero interaction
c     
c     sulist(1, 1...numsu) contains the values of s
c     sulist(2, 1...numsu) contains the values of u
c     
c     epqrstu(1...idim,1...nvec,1...numsu) contains the coupling coeffcients
c     
      integer nvec              ! [input] No. of vectors on the left hand side
      double precision vec(ndets,nvec) ! [input] The CI vectors
      integer p, q, r, t        ! [input] Fixed labels on Epq,rs,tu
      integer iocc(nactiv)      ! [input] Orbital occupation of I
      integer idim              ! [input] No. of spin functions for I
      integer numsu             ! [output] Returns the number of s,u
      integer sulist(2,*)       ! [output] Returns the list of interacting s,u
      double precision epqrstu(idim,nvec,*) ! [output] <Vi|Epqrstu|Iw>
c     
c     <V|Epq,rs,tu|Iw> = sum(Jv) V(Jv) <Jv|EpaErbEtcEcuEbsEaq|Iw> 
c     .  = sum(Jz) V(Jv) <Jv|EpaErbEtc|Oz><Oz|EcuEbsEaq|Iw>
c     .  = sum(Jz) V(Jv) (<Jv|Epa|Ny><Ny|Erb|Hg><Hg|Etc|Oz>)*
c     .                  (<Oz|Ecu|Mx><Mx|Ebs|Kf><Kf|Eaq|Iw>)
c     .  = sum(Jz) V(Jv) (<Ny|Eap|Jv><Hg|Ebr|Ny><Oz|Ect|Hg>)*
c     .                  (<Oz|Ecu|Mx><Mx|Ebs|Kf><Kf|Eaq|Iw>)
c     
c     where a,b,c are fictitious orbitals with c > b > a
c     and we have v(z), y(z), f(z), w(z), x(z), g(z)
c     
c     .      -q   -s   -u   +t   +r   +p
c     .    I -> K -> M -> O -> H -> N -> J
c     
c     Defer as much computation as possible until we have satisfied
c     all possible occupation and sparsity tests.
c     
      integer indbar(nactmx)    ! Posn amoung sing. occ. I/K/M/O/H/N/J
      integer kperm(ndetmx), mperm(ndetmx),hperm(ndetmx),
     $     nperm(ndetmx), iperm(ndetmx) ! Permutations of spin functions
      double precision kval(ndetmx), nval(ndetmx), mval(ndetmx),
     $     hval(ndetmx), ival(ndetmx) ! Accumulated phase factors
*     integer jperm(ndetmx)          ! Not used by optimized code
*     double precision jval(ndetmx)  ! Not used by optimized code
c     
      integer nsi, nsj, nsk, nsm, nsn, nsh, nso ! No. of sing. occ in I/...
      integer jdim, kdim, mdim, ndim, odim, hdim ! No. of dets for J/...
      integer np                ! Occuption of p in J
      integer nq                ! Occupation of q in I
      integer ns                ! Occupation of s in K
      integer nr                ! Occupation of r in N
      integer nt                ! Occupation of t in H
      integer nu                ! Occupation of u in M
      integer jwalk             ! Lexical index of J
      integer jdet              ! Offset of J in CI vector
      integer sympqrt           ! Sym product state*I*p*q*r*t
      integer symu              ! Required symmetry for u
      integer u, s, w, z, v, pbar, qbar, rbar, sbar, tbar, ubar,
     $     deltat, deltau, deltaq, deltar, deltas, idummy, i, ivec
      integer shi, uhi
      double precision sum, dummy
      logical omade_kval        ! Used to defer computation of kperm/kval
      logical omade_ival        ! Used to defer computaiton of iperm/ival
      double precision thresh
      parameter (thresh = 1d-12) ! Should be SQUARE of desired accuracy
#include "bitops.fh"      
c     
c     Check to see if q is really occupied, compute no. of open shells
c     and I indbar vector
c     
      nq = iocc(q)-1
      if (nq .le. 0) then
         numsu = 0
         return
      endif
c
      call fci_tstart(3,.false.)
c     
      nsi = 0
      sympqrt = ieor(ieor(orbsym(p),orbsym(r)),statesym)
      sympqrt = ieor(ieor(sympqrt,orbsym(q)),orbsym(t))
      do i = 1,nactiv
         if (iocc(i).eq.2) then 
            nsi = nsi + 1
            sympqrt = ieor(sympqrt,orbsym(i))
         endif
         indbar(i) = nsi
      enddo
      if (idim .ne. mnsdet(nsi)) call fci_err
     $     ('fci_couple3: mismatch for idim ', idim)
c     
      numsu = 1
c     
      iocc(q) = iocc(q) - 1     ! iocc now refers to K
c
c     More detailed info about q and K
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
c    Defer computation of <Ku|Eaq|Iw> until there is a non-zero term
c     
      omade_ival = .false.
c     
c     Loop thru s occupied in K
c     
      shi = nactiv
      if (p .eq. r) shi = q
      do s = 1, shi
         ns = iocc(s)-1         ! Occupation of s in K
         if (ns .gt. 0) then
            iocc(s) = iocc(s) - 1 ! iocc now refers to M
c     
            nsm  = nsk
            if (ns .eq. 2) nsm = nsm + 2
c     
c     Defer computation of <Mx|EbsEaq|Iw> = <Mx|Ebs|Ku><Ku|Eaq|Iw>
c     until we are sure we will need it
c     
            omade_kval = .false.
c     
c     Loop thru u occupied in M
c     
            uhi = nactiv
            if (t .eq. r) uhi = s
            symu = ieor(sympqrt,orbsym(s))
            do u = 1, uhi
               nu = iocc(u)-1   ! Occupation of u in M
               if (nu.gt.0 .and. orbsym(u).eq.symu) then
                  iocc(u) = iocc(u) - 1 ! iocc now refers to O
c     
                  nt = iocc(t)  ! Since iocc refers to O and we want H
                  if (nt .lt. 3) then
                     iocc(t) = iocc(t) + 1 ! iocc now refers to H
c     
                     nr = iocc(r) ! Since iocc refers to H and we want N
                     if (nr .lt. 3) then
                        iocc(r) = iocc(r) + 1 ! iocc now refers to N
c     
                        np = iocc(p) ! Since iocc refers to N and we want J
                        if (np .lt. 3) then
                           iocc(p) = iocc(p) + 1 ! iocc now refers to J
c     
                           nso = nsm
                           if (nu .eq. 2) nso = nso + 2
c     
                           nsh = nso
                           if (nt .eq. 2) nsh = nsh - 2
c     
                           nsn = nsh
                           if (nr .eq. 2) nsn = nsn - 2
c     
                           nsj = nsn
                           if (np .eq. 2) nsj = nsj - 2 
c     
c     Determine the index of J in the CI vector and apply sparsity test
c     
                           jdim = mnsdet(nsj)
*                           call fci_ofind_weight(jwalk,iocc)
*     manually inlined version of fci_ofind_weight
*
                           z = 1
                           jwalk = 1
                           do i = 1,nactiv
                              jwalk = jwalk + karcwt(iocc(i),z,i)
                              z = karc(iocc(i),z,i)
                           enddo
*
                           jdet = indexo(jwalk)
c     
                           sum = 0.0d0
                           do ivec = 1, nvec
                              do v = 1, jdim
                                 sum = sum + 
     $                                vec(jdet+v,ivec)*vec(jdet+v,ivec)
                              enddo
                           enddo
                           if (sum .gt. thresh) then
c
c     Make <Ku|Eaq|Iw> if have not already done so
c
                              if (.not. omade_ival) then
                                 call fci_apply_eaj(qbar, nq, idim,kdim, 
     $                                eajdata(eajptr(nsi)), 
     $                                eajphase(1,1,nsi), .true.,
     $                                idummy, dummy, iperm, ival)
                                 do i = q, nactiv ! Adjust indbar to K
                                    indbar(i) = indbar(i) + deltaq
                                 enddo
                                 omade_ival = .true.
                              endif
c
c     If not already done form <Mx|EbsEaq|Iw> = <Mx|Ebs|Ku><Ku|Eaq|Iw>
c
                              if (.not. omade_kval) then
                                 sbar = indbar(s)
                                 deltas = -1
                                 if (ns .eq. 2) then
                                    sbar = sbar + 1
                                    deltas = 1
                                 endif
                                 do i = s, nactiv ! Adjust indbar to M
                                    indbar(i) = indbar(i) + deltas
                                 enddo
                                 mdim = mnsdet(nsm)
                                 call fci_apply_eaj(sbar, ns,kdim,mdim, 
     $                                eajdata(eajptr(nsk)),
     $                                eajphase(1,1,nsk),
     $                                .false., iperm, ival, kperm, kval)
                                 omade_kval = .true.
                              endif
c     
c     More detailed info about u and O
c     
                              deltau = -1
                              ubar = indbar(u)
                              if (nu .eq. 2) then
                                 ubar = ubar + 1
                                 deltau = 1
                              endif
                              odim = mnsdet(nso)
c     
c     More detailed info about t and H
c     
                              deltat = 1
                              tbar = indbar(t)
                              if (nt .eq. 2) then
                                 deltat = -1
                                 tbar = tbar + 1
                              endif
                              hdim = mnsdet(nsh)
                              tbar = tbar + deltat
                              if (t .ge. u) tbar = tbar + deltau
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
                              if (r .ge. u) rbar = rbar + deltau
                              if (r .ge. t) rbar = rbar + deltat
c     
c     More detailed info about p and J
c     
                              pbar = indbar(p)
                              if (np .eq. 1) pbar = pbar + 1
                              if (p .ge. r) pbar = pbar + deltar
                              if (p .ge. u) pbar = pbar + deltau
                              if (p .ge. t) pbar = pbar + deltat
c     
c     Form <Oz|Ecu|Mx><Mx|Ebs|Kf><Kf|Eaq|Iw>
c     
                              call fci_apply_eaj(ubar, nu, mdim, odim, 
     $                             eajdata(eajptr(nsm)), 
     $                             eajphase(1,1,nsm),
     $                             .false., kperm, kval, mperm, mval)
c
                              call fci_apply_eaj3(
     $                             pbar, np, eajdata(eajptr(nsj)), 
     $                             eajphase(1,1,nsj),
     $                             rbar, nr, eajdata(eajptr(nsn)), 
     $                             eajphase(1,1,nsn),
     $                             tbar, nt, eajdata(eajptr(nsh)), 
     $                             eajphase(1,1,nsh),
     $                             jdim, ndim, hdim, odim,
     $                             nperm, nval, hperm, hval)
c$$$
c$$$c     
c$$$c     Form <Ny|Eap|Jv>
c$$$c     
c$$$                              call fci_apply_eaj(pbar, np, jdim, ndim, 
c$$$     $                             eajdata(eajptr(nsj)), 
c$$$     $                             eajphase(1,1,nsj),
c$$$     $                             .true., idummy, dummy, jperm, jval)
c$$$c     
c$$$c     Form <Jv|EpaErb|Hg> = <Ny|Eap|Jv><Hg|Ebr|Ny>
c$$$c     
c$$$                              call fci_apply_eaj(rbar, nr, ndim, hdim, 
c$$$     $                             eajdata(eajptr(nsn)), 
c$$$     $                             eajphase(1,1,nsn),
c$$$     $                             .false., jperm, jval, nperm, nval)
c$$$c     
c$$$c     Form <Ny|Eap|Jv><Hg|Ebr|Ny><Oz|Ect|Hg>
c$$$c     
c$$$                              call fci_apply_eaj(tbar, nt, hdim, odim, 
c$$$     $                             eajdata(eajptr(nsh)), 
c$$$     $                             eajphase(1,1,nsh),
c$$$     $                             .false., nperm, nval, hperm, hval)
c     
c     Now sum(Jz) V(Jv) <Jv|EpaErbEtc|Oz><Oz|EcuEbsEaq|Iw>
c     .
c     .   noting that we have v(z) and w(z)
c     
                              do ivec = 1, nvec
                                 do w = 1, idim
                                    epqrstu(w,ivec,numsu) = 0.0d0
                                 enddo
                              enddo
c     
                              do ivec = 1, nvec
                                 do z = 1, odim
                                    v = hperm(z)
                                    w = mperm(z)
                                    epqrstu(w,ivec,numsu) = 
     $                                   epqrstu(w,ivec,numsu) + 
     $                                   vec(jdet+v,ivec)*
     $                                   hval(z)*mval(z)
                                 enddo
                              enddo
                              sulist(1,numsu) = s
                              sulist(2,numsu) = u
                              numsu = numsu + 1
                           endif ! sum.gt.thresh
c     
                           iocc(p) = iocc(p) - 1 ! iocc now refers to N
                        endif   ! np.lt.3
c     
                        iocc(r) = iocc(r) - 1 ! iocc now refers to H
                     endif      ! nr.lt.3
c     
                     iocc(t) = iocc(t) - 1 ! iocc now refers to O
                  endif         ! nt.lt.3
c     
                  iocc(u) = iocc(u) + 1 ! iocc now refers to M
               endif            ! nu.gt.0
            enddo               ! do u
c     
            if (omade_kval) then
               do i = s, nactiv ! Adjust indbar to refer to K
                  indbar(i) = indbar(i) - deltas
               enddo
            endif
            iocc(s) = iocc(s) + 1 ! iocc now refers to K
         endif                  ! ns.gt.0
      enddo                     ! do s
c     
      iocc(q) = iocc(q) + 1     ! iocc now refers to I
c
      numsu = numsu - 1
c
      call fci_tfinish(3,.false.)
c     
      end
      
