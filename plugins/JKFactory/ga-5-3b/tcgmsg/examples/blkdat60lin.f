      block data
C$Id: blkdat60lin.f,v 1.2 1995-02-02 23:23:52 d3g681 Exp $
      implicit double precision (a-h, o-z)
      include 'cscf.h'
c
c     initalize data in common ... clumsy but avoids code to read in data
c
c     line of four be atoms 4.0 a.u. apart with 60 orbitals
c
c     have 9s functions on each center and simulate p's by having
c     s function at +- 1 in each of x, y, z
c
      data ax /   0.0d0, 4.0d0, 8.0d0, 12.0d0/
      data ay /   0.0d0, 0.0d0, 0.0d0,  0.0d0/
      data az /   0.0d0, 0.0d0, 0.0d0,  0.0d0/
      data  q /4*4.0d0/, enrep/17.3333333333333d0/
c
      data x /9*0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,
     $        9*4.0d0, 5.6d0,  2.4d0, 4.0d0, 4.0d0, 4.0d0, 4.0d0,
     $        9*8.0d0, 9.6d0,  6.4d0, 8.0d0, 8.0d0, 8.0d0, 8.0d0,
     $       9*12.0d0,13.6d0, 10.4d0,12.0d0,12.0d0,12.0d0,12.0d0/
      data y /9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0/
      data z /9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0,
     $        9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0/
      data expnt /1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0,
     $            1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0/
      end
