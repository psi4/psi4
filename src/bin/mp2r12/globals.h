/*! \file
    \ingroup MP2R12
    \brief Enter brief description of file here 
*/
#ifndef _psi_src_bin_mp2r12_globals_h_
#define _psi_src_bin_mp2r12_globals_h_

extern "C" {
  extern FILE *infile, *outfile;
}
namespace psi{
  namespace mp2r12{
  extern int *ioff;
  extern struct MOInfo moinfo;
  extern struct Params params;
  } /* Namespace mp2r12 */
} /* Namespace psi */


#endif /* Header guard */
