*
* $Id: fci_mgraph.F,v 1.2 1999/07/28 00:23:41 d3e129 Exp $
*
      subroutine fci_mgraph
      implicit none
#include "fciinfo.fh"
c     
c     generate determinant spin graphs for up to nsmax electrons
c     with multi = 2sz+1. storage of graphs is inefficient but
c     easy to use ... they are small anyway!
c     
      logical oprint
      integer nsoc, is, irow, k, nspt, krow, i, mdead, iwsum
c      
      oprint = .false.
c     
      if(multi.lt.1 .or. multi.gt.(2*nsmax+1)) 
     $     call fci_err('multi error',multi)
      call ifill(0,marc,1,2*nsmax1*nsmax1*nsmax2)
      call ifill(0,marcwt,1,2*nsmax1*nsmax1*nsmax2)
      call ifill(0,mwt,1,nsmax1*nsmax1*nsmax2)
      call ifill(0,mnsdet,1,nsmax)
      call ifill(0,mnscsf,1,nsmax)
c     
c     loop over no. of unpaired electrons and generate arrays
c     
      do 10 nsoc = multi-1,nsmax,2
         nspt = nsoc/2+1
c     fill in arc and wt ignoring sz
         mwt(1,1,nspt)=1
         do 20 is = 1,nsoc
            do 30 irow = 1,is
               do 40 k = 1,2
                  krow = irow + (2-k)
                  marc(k,irow,is,nspt) = krow
                  mwt(krow,is+1,nspt)=mwt(krow,is+1,nspt) +
     &                 mwt(irow,is  ,nspt)
 40            continue
 30         continue
 20      continue
c     eliminate unwanted arcs/rows
         do 50 irow = 1,nsoc+1
            if(irow.ne.multi+(nsoc-multi+1)/2) mwt(irow,nsoc+1,nspt)=0
 50      continue
c     
c     set no. of dets and csf for no. of open shells. assume that sz = s
c     
         mnsdet(nsoc) = mwt(multi+(nsoc-multi+1)/2,nsoc+1,nspt)
         ncsf = multi
         do 52 i=(nsoc+1+multi)/2+1,nsoc+1
            ncsf = ncsf * i
 52      continue
         do 53 i=1,(nsoc-multi+1)/2
            ncsf = ncsf / i
 53      continue
         mnscsf(nsoc) = ncsf/(nsoc+1)
         do 60 is = nsoc,1,-1
            do 70 irow = 1,is
               mdead = 0
               do 80 k=1,2
                  krow = irow + (2-k)
                  if(mwt(krow,is+1,nspt).eq.0) then
                     marc(k,irow,is,nspt)=0
                     mdead = mdead+1
                  endif
 80            continue
               if(mdead.eq.2) mwt(irow,is,nspt)=0
 70         continue
 60      continue
c     generate arc weight array
         do 90 is = nsoc+1,2,-1
            do 100 irow = 1,is
               if(mwt(irow,is,nspt).gt.0) then
                  iwsum = 0
                  do 110 k = 2,1,-1
                     krow = irow-(2-k)
                     if(krow.gt.0) then
                        if(marc(k,krow,is-1,nspt).gt.0) then
                           marcwt(k,krow,is-1,nspt)=iwsum
                           iwsum = iwsum + mwt(krow,is-1,nspt)
                        endif
                     endif
 110              continue
               endif
 100        continue
 90      continue
         if (oprint) then
            write(6,*)
            write(6,*) ' mwt,marc,marcwt arrays for nsoc = ',nsoc
            write(6,*) ' no. of determinants = ',mnsdet(nsoc)
            write(6,*) ' no. of configurations = ',mnscsf(nsoc)
            write(6,1101) (irow,irow=1,nsoc+1)
 1101       format('           ',11(5x,i2,5x)/
     &           '           ',11(2x,'----------',2x)/)
            do 120 is=1,nsoc+1
               write(6,1201) is,(mwt(irow,is,nspt),irow=1,is)
               write(6,1202) ((marc(k,irow,is,nspt),k=1,2),irow=1,is)
               write(6,1202) ((marcwt(k,irow,is,nspt),k=1,2),irow=1,is)
 120        continue
 1201       format(' level ',i3/'           ',11(5x,i2,5x))
 1202       format('           ',11(2x,i3,2x,i3,2x))
         endif
c     
 10   continue
c
      end
