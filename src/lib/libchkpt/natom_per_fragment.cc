/*!
  \file
  \ingroup CHKPT
*/

#include <stdio.h>
#include <stdlib.h>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

int *Chkpt::rd_natom_per_fragment(void)
{
        int nfragment;
        int *natom_per_fragment;
        char *keyword;
        keyword = build_keyword("Num. atoms per fragment");

        nfragment = rd_nfragment();
        natom_per_fragment = array<int>(nfragment);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) natom_per_fragment,
          nfragment*sizeof(int));

        free(keyword);
        return natom_per_fragment;
}

void Chkpt::wt_natom_per_fragment(int *natom_per_fragment)
{
        int nfragment;
        char *keyword;
        keyword = build_keyword("Num. atoms per fragment");

        nfragment = rd_nfragment();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) natom_per_fragment,
          nfragment*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_natom_per_fragment():  Reads in the number of frozen doubly occupied molecular
**   orbitals in each irrep.
**
**   takes no arguments.
**
**   returns:
**     int *natom_per_fragment  an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each
**                 element contains the number of frozen doubly occupied
**                 molecular orbitals for
**                 that irrep. Also, see chkpt_rd_sopi().
** \ingroup CHKPT
*/
        int *chkpt_rd_natom_per_fragment(void)
        {
                int *natom_per_fragment;
                natom_per_fragment = _default_chkpt_lib_->rd_natom_per_fragment();
                return natom_per_fragment;
        }


/*!
** chkpt_wt_natom_per_fragment():  Writes the number of frozen doubly occupied molecular
**   orbitals in each irrep
**
** \param natom_per_fragment = an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each
**                 element contains the number of frozen doubly occupied
**                 molecular orbitals for that irrep.  See also
**                 chkpt_rd_sopi().
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_natom_per_fragment(int *natom_per_fragment)
        {
                _default_chkpt_lib_->wt_natom_per_fragment(natom_per_fragment);
        }
}
