/*!
  \file
  \ingroup CHKPT
*/

#include <cstdio>
#include <cstdlib>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

double Chkpt::rd_e_t()
{
        double energy;
        char *keyword;
        keyword = build_keyword("(T) Energy");

        // Read the energy in
        psio->read_entry(PSIF_CHKPT, keyword, (char*)&energy, sizeof(double));

        // Return the value to the user
        return energy;
}

void Chkpt::wt_e_t(double e_t)
{
        char *keyword;
        keyword = build_keyword("(T) Energy");

        psio->write_entry(PSIF_CHKPT, keyword, (char*)&e_t, sizeof(double));

        free(keyword);
}


extern "C" {
/*!
** chkpt_rd_e_t(): Reads in the (T) contribution to total energy.
**
**   takes no arguments.
**
**   returns: double e_t  the (T) energy.
** \ingroup CHKPT
*/
        double chkpt_rd_e_t(void)
        {
                double e_t;
                e_t = _default_chkpt_lib_->rd_e_t();
                return e_t;
        }

/*!
** chkpt_wt_e_t(): Writes out the (T) contribution to total energy.
**
** \param e_t = the (T) energy.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_e_t(double e_t)
        {
                _default_chkpt_lib_->wt_e_t(e_t);
        }
}
