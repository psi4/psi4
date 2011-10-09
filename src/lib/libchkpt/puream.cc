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

bool Chkpt::rd_puream(const char *key2)
{
        int puream;
        char *keyword;
        keyword = build_keyword("Pure Harmonics?", key2);

        psio->read_entry(PSIF_CHKPT, keyword, (char *) &puream, sizeof(int));

        free(keyword);
        return (puream == 1);
}

void Chkpt::wt_puream(bool puream, const char *key2)
{
        char *keyword;
        keyword = build_keyword("Pure Harmonics?", key2);

        int p = puream ? 1 : 0;
        psio->write_entry(PSIF_CHKPT, keyword, (char *) &p, sizeof(int));

        free(keyword);
}

extern "C" {
/*!
**  int chkpt_rd_puream()
**  Reads whether cartesian or spherical harmonics are used (Psi is currently
**  limited to only using one type of functions at a time)
**
**  returns: 1 (harmonics) or 0 (cartesian)
**  \ingroup CHKPT
*/
        int chkpt_rd_puream(void)
        {
                return _default_chkpt_lib_->rd_puream() ? 1 : 0;
        }

/*!
**  void chkpt_wt_puream(int)
**  Writes whether cartesian or spherical harmonics are used (Psi is currently
**  limited to only using one type of functions at a time)
**
**  \param 1 (harmonics) or 0 (cartesian)
**
**  returns: none
**  \ingroup CHKPT
*/
        void chkpt_wt_puream(int puream, const char *key2)
        {
                _default_chkpt_lib_->wt_puream(puream == 1, key2);
        }
}
