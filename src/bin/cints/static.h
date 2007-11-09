#ifndef _psi_src_bin_cints_static_h
#define _psi_src_bin_cints_static_h

/*! \file static.h
    \ingroup (CINTS)
    \brief Static variables.
*/
//namespace psi {
//  namespace CINTS {
  /*------------------------------------------------
    \internal Useful definitions for some arrays which should
    be static to make functions work well
    ------------------------------------------------*/
#ifdef STATIC_OO2NP1
  static double oo2np1[] = {1.0,  1.0/3.0,  1.0/5.0,  1.0/7.0,  1.0/9.0,
			    1.0/11.0, 1.0/13.0, 1.0/15.0, 1.0/17.0, 1.0/19.0,
			    1.0/21.0, 1.0/23.0, 1.0/25.0, 1.0/27.0, 1.0/29.0,
			    1.0/31.0, 1.0/33.0, 1.0/35.0, 1.0/37.0, 1.0/39.0,
			    1.0/41.0, 1.0/43.0, 1.0/45.0, 1.0/47.0, 1.0/49.0,
			    1.0/51.0, 1.0/53.0, 1.0/55.0, 1.0/57.0, 1.0/59.0,
			    1.0/61.0, 1.0/63.0, 1.0/65.0, 1.0/67.0, 1.0/69.0,
			    1.0/71.0, 1.0/73.0, 1.0/75.0, 1.0/77.0, 1.0/79.0};
#endif
  
#ifdef STATIC_OON
  static double oon[] = {0.0, 1.0, 1.0/2.0, 1.0/3.0, 1.0/4.0, 1.0/5.0};
#endif
  
//};};
#endif
