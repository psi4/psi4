#ifndef _psi_src_bin_cints_Tools_molecule_h
#define _psi_src_bin_cints_Tools_molecule_h

/*! \file molecule.h
    \ingroup CINTS
*/

namespace psi { namespace CINTS {
void init_molecule();
void cleanup_molecule();
void compute_enuc(void);

};}
#endif
