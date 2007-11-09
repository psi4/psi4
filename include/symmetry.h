/*-------------------------------------------------------
  symmetry.h : Definitions of various symmetry constants

  Edward Valeev, Oct. 1999
 -------------------------------------------------------*/

#ifndef _psi_include_symmetry_h_
#define _psi_include_symmetry_h_

/*--------------------------------------------------------------
  Point groups at this moment are limited to D2h and sumbgroups
  Therefore the number of symmetry elements is at most 8. Each
  symmetry operation corresponds to a bit of a byte word. The
  correspondence is hardwired via defines "GFLAG" where G is
  the operation symbol.

  To describe nuclear stabilizers or subgroups in general I use
  a byte in which bits corresponding to the symmetry operations
  that constitute the group are set. The result is that each
  operation's contribution to the byte equals "GCODE".
 --------------------------------------------------------------*/
#define EFLAG 0
#define C2ZFLAG 1
#define C2YFLAG 2
#define C2XFLAG 3
#define IFLAG 4
#define SIGXYFLAG 5
#define SIGXZFLAG 6
#define SIGYZFLAG 7
#define C2XCODE 1<<C2XFLAG
#define C2YCODE 1<<C2YFLAG
#define C2ZCODE 1<<C2ZFLAG
#define ICODE 1<<IFLAG
#define SIGXYCODE 1<<SIGXYFLAG
#define SIGXZCODE 1<<SIGXZFLAG
#define SIGYZCODE 1<<SIGYZFLAG
#define ECODE 1<<EFLAG

/*-----------------
  Indices for axes
 -----------------*/
#define XAXIS 0
#define YAXIS 1
#define ZAXIS 2

#endif /* header guard */
