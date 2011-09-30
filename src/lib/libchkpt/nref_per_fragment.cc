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

int *Chkpt::rd_nref_per_fragment(void)
{
        int nfragment;
        int *nref_per_fragment;
        char *keyword;
        keyword = build_keyword("Ref. atoms per fragment");

        nfragment = rd_nfragment();
        nref_per_fragment = array<int>(nfragment);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) nref_per_fragment,
          nfragment*sizeof(int));

        free(keyword);
        return nref_per_fragment;
}

void Chkpt::wt_nref_per_fragment(int *nref_per_fragment)
{
        int nfragment;
        char *keyword;
        keyword = build_keyword("Ref. atoms per fragment");

        nfragment = rd_nfragment();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) nref_per_fragment,
          nfragment*sizeof(int));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_nref_per_fragment():  Reads in the number of frozen doubly occupied molecular
**   orbitals in each irrep.
**
**   takes no arguments.
**
**   returns:
**     int *nref_per_fragment  an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each
**                 element contains the number of frozen doubly occupied
**                 molecular orbitals for
**                 that irrep. Also, see chkpt_rd_sopi().
** \ingroup CHKPT
*/
        int *chkpt_rd_nref_per_fragment(void)
        {
                int *nref_per_fragment;
                nref_per_fragment = _default_chkpt_lib_->rd_nref_per_fragment();
                return nref_per_fragment;
        }


/*!
** chkpt_wt_nref_per_fragment():  Writes the number of frozen doubly occupied molecular
**   orbitals in each irrep
**
** \param nref_per_fragment = an array which has an element for each irrep of the
**                 point group of the molecule (n.b. not just the ones
**                 with a non-zero number of basis functions). each
**                 element contains the number of frozen doubly occupied
**                 molecular orbitals for that irrep.  See also
**                 chkpt_rd_sopi().
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_nref_per_fragment(int *nref_per_fragment)
        {
                _default_chkpt_lib_->wt_nref_per_fragment(nref_per_fragment);
        }
}
