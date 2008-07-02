/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include "input.h"
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

void init_globals()
{
  int i, j, errcod;

  disp_num = 0;
  
  /*Make ioff and stuff*/
  setup(MAXIOFF);

  /*No rotation has been performed yet -
    hence canonical and reference frames are equivalent */
  Rref = block_matrix(3,3);
  canon_eq_ref_frame();

  /*GTO data has not been initialized yes */
  GTOs.max_angmom = 0;

  return;
}

}} // namespace psi::input
