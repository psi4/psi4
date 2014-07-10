#if HAVE_CONFIG_H
#   include "config.fh"
#endif
       subroutine output (z,rowlow,rowhi,collow,colhi,rowdim,coldim,
     $     nctl)
c.......................................................................
c output prints a real*8 matrix in formatted form with numbered rows
c and columns.  the input is as follows;
c        matrix(*,*).........matrix to be output
c        rowlow..............row number at which output is to begin
c        rowhi...............row number at which output is to end
c        collow..............column number at which output is to begin
c        colhi...............column number at which output is to end
c        rowdim..............row dimension of matrix(*,*)
c        coldim..............column dimension of matrix(*,*)
c        nctl................carriage control flag; 1 for single space
c                                                   2 for double space
c                                                   3 for triple space
c the parameters that follow matrix are all of type integer*4.  the
c program is set up to handle 5 columns/page with a 1p5d24.15 format for
c the columns.  if a different number of columns is required, change
c formats 1000 and 2000, and initialize kcol with the new number of
c columns.
c author;  nelson h.f. beebe, quantum theory project, university of
c          florida, gainesville
c.......................................................................
C$Id: output.F,v 1.1 2005/03/08 23:58:03 d3g293 Exp $
      implicit double precision (a-h,o-z)
#include "mafdecls.fh"
#include "global.fh"
#include "mp3def.fh"
      integer rowlow,rowhi,collow,colhi,rowdim,coldim,begin,kcol
      dimension z(rowdim,coldim)
      dimension asa(3)
      data column/8hcolumn   /,asa/8h          ,8h00000000  ,
     1     8h--------  /,blank/8h          /
      data kcol/8/
      data zero/0.d00/
      if (ga_nodeid().ne.0) return
      do 11 i=rowlow,rowhi
         do 10 j=collow,colhi
            if (z(i,j).ne.zero) go to 15
 10      continue
 11   continue
      write (6,3000)
 3000 format (/' zero matrix'/)
      go to 3
 15   continue
      if (rowhi.lt.rowlow) go to 3
      if (colhi.lt.collow) go to 3
      last = min0(colhi,collow+kcol-1)
      do 2 begin = collow,colhi,kcol
         write (6,1000) (column,i,i = begin,last)
         do 1 k = rowlow,rowhi
            do 4 i=begin,last
               if (z(k,i).ne.zero) go to 5
 4          continue
            go to 1
c
5       write (6,2000) k,(z(k,i), i = begin,last)               

 1       continue
         last = min0(last+kcol,colhi)
 2    continue
 3    return
 1000 format (/11x,7(a3,i3,3x),(a3,i3)) 
 2000 format ('row',i4,1x,8f9.4) 
      end
