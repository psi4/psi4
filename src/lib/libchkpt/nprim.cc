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

int Chkpt::rd_nprim(const char *key2)
{
    int nprim;
    char *keyword;
    keyword = build_keyword("Num. primitives", key2);

    psio->read_entry(PSIF_CHKPT, keyword, (char *) &nprim, sizeof(int));

    free(keyword);
    return nprim;
}

void Chkpt::wt_nprim(int nprim, const char *key2)
{
    char *keyword;
    keyword = build_keyword("Num. primitives", key2);

    psio->write_entry(PSIF_CHKPT, keyword, (char *) &nprim, sizeof(int));

    free(keyword);
}

extern "C" {
/*!
** int chkpt_rd_nprim()
** Reads in the total number of primitive Gaussian functions
** (only primitives of symmetry independent atoms are taken into account!).
**
** returns: nprim = total number of primitive Gaussian functions.
** \ingroup CHKPT
*/
int chkpt_rd_nprim(void)
{
    return _default_chkpt_lib_->rd_nprim();
}

/*!
** void chkpt_wt_nprim(int)
** Writes out the total number of primitive Gaussian functions
** (only primitives of symmetry independent atoms are taken into account!).
**
** \param nprim = total number of primitive Gaussian functions.
**
** returns: none
** \ingroup CHKPT
*/
void chkpt_wt_nprim(int nprim, const char *key2)
{
    _default_chkpt_lib_->wt_nprim(nprim, key2);
}
}
