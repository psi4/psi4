/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/

#ifndef _psi3_DBOC_defines_h_
#define _psi3_DBOC_defines_h_

namespace psi{ namespace DBOC {

const int MAX_NUM_DISP=4;
typedef enum { MinusDelta = 0, PlusDelta = 1, Minus2Delta = 2, Plus2Delta = 3} DisplacementIndex;

}} /* namespace psi::DBOC */

//#define LONG_DOUBLE 1   // Use long doubles for FLOAT

#define USE_MOINFO 0  // set to 1 if using MOInfo

#endif
