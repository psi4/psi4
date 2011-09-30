/*!
  \file
  \ingroup CHKPT
*/

#include <cstdlib>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

int* Chkpt::rd_atom_dummy(void)
{
        int num_allatoms;
        int *atom_dummy;
        char *keyword;
        keyword = build_keyword("Dummy atom flags");

        num_allatoms = rd_nallatom();
        atom_dummy = new int[num_allatoms];
        //atom_dummy = (int *) malloc(sizeof(int)*num_allatoms);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) atom_dummy,
                num_allatoms*sizeof(int));

        free(keyword);
        return atom_dummy;
}

void Chkpt::wt_atom_dummy(int* atom_dummy)
{
        int num_allatoms = rd_nallatom();
        char *keyword;
        keyword = build_keyword("Dummy atom flags");

        psio->write_entry(PSIF_CHKPT, keyword, (char *) atom_dummy,
                num_allatoms*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_atom_dummy()
**
** Reads the array of flags which indicate whether the atom in full_geom
** is dummy
**
** Parameters: none
**
** Returns: atom_dummy = array of integers nallatom long.
** \ingroup CHKPT
*/
        int* chkpt_rd_atom_dummy(void)
        {
                return _default_chkpt_lib_->rd_atom_dummy();
        }


/*!
** chkpt_wt_atom_dummy()
**
** Writes the array of flags which indicate whether the atom in full_geom
** is dummy
**
** \param atom_dummy = array of integers nallatom long.
**
** Returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_atom_dummy(int* atom_dummy)
        {
                _default_chkpt_lib_->wt_atom_dummy(atom_dummy);
        }
}

