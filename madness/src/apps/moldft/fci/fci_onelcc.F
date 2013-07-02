*
* $Id: fci_onelcc.F,v 1.2 1999/07/28 00:23:41 d3e129 Exp $
*
      subroutine fci_onelcc(iocc,idet,idim,isym,
     $     ionel,nonel,nonemx)
      implicit none
#include "fciinfo.fh"
c     
c     compute all non-zero one-electron coupling coefficents
c     g(i,j;I,J) = <I|Eij|J> for all orbitals ij and for all
c     determinants generated from orbital configuration iwalk
c     
c     ionel(1,*) = i
c     ionel(2,*) = j
c     ionel(3,*) = ideterminant
c     ionel(4,*) = jdeterminant
c     ionel(5,*) = value (only -1,+1,2 are valid)
c     
c     <Iv|Eij|Jw> = <Iv|Eia|Ku><Ku|Eaj|Jw>
c     .           = <Ku|Eai|Iv><Ku|Eaj|Jw>
c
c     Orbital configuration I may be of any symmetry (specified 
c     only by its occupation etc.) but J is constrained to be
c     of the same symmetry as the CI state.
c     
      integer iocc(nactiv),idet,idim,isym ! [input] Info about I
      integer nonemx            ! [input] last dimension of ionel
      integer ionel(5,nonemx)   ! [output] returns coupling coeff info
      integer nonel             ! [output] no. of coeff returned
c     
      integer indbar(nactmx)
      integer jperm(ndetmx), iperm(ndetmx)
      double precision jval(ndetmx), ival(ndetmx)
c     
      integer nsi, nsj, nsk, jdet, jwalk, ni, nj,
     $     jdim, kdim, i, j, u, v, ibar, jbar, deltai, value,
     $     idummy, ii, jsym
      logical oprint
#include "bitops.fh"
c     
      oprint = .false.
      nonel = 0
c     
c     Form I indbar vector
c     
      if (oprint) then
         write(6,*) ' iocc ', (iocc(i),i=1,nactiv)
         write(6,*) ' isym ', isym, statesym
      endif
c
      nsi = 0
      do i = 1,nactiv
         if (iocc(i).eq.2) nsi = nsi + 1
         indbar(i) = nsi
      enddo
      if (idim .ne. mnsdet(nsi)) call fci_err('fci_onelcc: nsi?',nsi)
c     
c     Loop thru i
c     
      do i = 1, nactiv
         ni = iocc(i)-1
         if (ni .gt. 0) then
            if (oprint) write(6,*) ' i = ', i, ' ni = ', ni
c     
c     Handle Eii special case for efficiency
c     
            if (isym .eq. statesym) then
               call fci_ofind_weight(jwalk,iocc) ! idet may not be absolute
               jdet = indexo(jwalk)
               do v = 1,idim
                  nonel = nonel + 1
                  ionel(1,nonel) = i
                  ionel(2,nonel) = i
                  ionel(3,nonel) = v + idet
                  ionel(4,nonel) = v + jdet
                  ionel(5,nonel) = ni
               enddo
            endif
c     
c     Form <Ku|Eai|Iv>
c     
            ibar = indbar(i)
            nsk  = nsi
            if (ni .eq. 2) then
               nsk = nsk + 2
               ibar = ibar + 1
            endif
            kdim = mnsdet(nsk)
            call fci_apply_eaj(ibar, ni, idim, kdim, 
     $           eajdata(eajptr(nsi)), eajphase(1,1,nsi), .true.,
     $           idummy, idummy, iperm, ival)
c
c
            if (oprint) write(6,17) i, (iperm(u),ival(u),u=1,kdim)
 17         format(i5,20(i3,f5.1,2x))
c     
c     Remove the electron in i from occupation vectors
c     
            if(oprint)write(6,88) ' I ', (iocc(u),indbar(u),u=1,nactiv)
 88         format(a,2x,20(i2,'(',i2,')',2x))
            iocc(i) = iocc(i) - 1
            deltai = -1
            if (ni .eq. 2) deltai = 1
            do ii = i,nactiv
               indbar(ii) = indbar(ii) + deltai
            enddo
            if(oprint)write(6,88) ' K ', (iocc(u),indbar(u),u=1,nactiv)
c
            jsym = ieor(ieor(isym,orbsym(i)),statesym) ! Sym. of orbital j
c     
c     indbar and iocc now refer to the occupation of K
c     
c     Now loop thru j!=i, and not doubly occupied in K, and form <Ku|Eaj|Jw>
c     
            do j = 1, nactiv
               nj = iocc(j)     ! Since iocc refers to K and we want J
               if ((nj .lt. 3) .and. (i .ne. j) .and. 
     $              (orbsym(j).eq.jsym)) then
                  if(oprint)write(6,*) ' j = ', j, ' nj = ', nj
c     
c     Determine index of J orbital occupation
c     
                  iocc(j) = iocc(j) + 1
                  if(oprint)write(6,11) i, j, (iocc(v),v=1,nactiv)
 11               format(2i5,5x,20i3)
                  call fci_ofind_weight(jwalk,iocc)
                  iocc(j) = iocc(j) - 1
c     
                  jbar = indbar(j)
                  if (nj .eq. 1) jbar = jbar + 1 ! Since indbar refers to K
                  nsj = nsk
                  if (nj .eq. 2) nsj = nsj - 2 
                  jdim = mnsdet(nsj)
                  if (oprint) write(6,*) 
     $                 ' jbar ', jbar, ' nsj ', nsj, ' nsk ', nsk
                  call fci_apply_eaj(jbar, nj, jdim, kdim, 
     $                 eajdata(eajptr(nsj)), eajphase(1,1,nsj), .true.,
     $                 idummy, idummy, jperm, jval)
                  if(oprint) write(6,17) j, (jperm(u),jval(u),u=1,kdim)
c     
c     Now use <Iv|Eij|Jw> = <Ku|Eai|Iv>*<Ku|Eaj|Jv>
c     
                  jdet = indexo(jwalk)
                  do u = 1, kdim
                     value = nint(ival(u)*jval(u))
                     if (value .ne. 0) then
                        nonel = nonel + 1
                        ionel(1,nonel) = i
                        ionel(2,nonel) = j
                        ionel(3,nonel) = iperm(u)+idet
                        ionel(4,nonel) = jperm(u)+jdet
                        ionel(5,nonel) = value
                     endif
                  enddo
               endif
            enddo
c     
c     Restore the electron in i to the occupation vectors
c     
            iocc(i) = iocc(i) + 1
            do ii = i,nactiv
               indbar(ii) = indbar(ii) - deltai
            enddo
c     
c     indbar and iocc now again refer to iocc
c     
         endif
      enddo
c
      if (oprint) call fci_prteop(ionel,nonel)
c
*      stop 1
c
      end
