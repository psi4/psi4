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

double Chkpt::rd_emp2()
{
        double energy;
        char *keyword;
        keyword = build_keyword("MP2 Energy");

        // Read the energy in
        psio->read_entry(PSIF_CHKPT, keyword, (char*)&energy, sizeof(double));

        // Return the value to the user
        return energy;
}

void Chkpt::wt_emp2(double emp2)
{
        char *keyword;
        keyword = build_keyword("MP2 Energy");

        psio->write_entry(PSIF_CHKPT, keyword, (char*)&emp2, sizeof(double));

        free(keyword);
}


extern "C" {
/*!
** chkpt_rd_emp2(): Reads in the MP2 contribution to total energy.
**
**   takes no arguments.
**
**   returns: double emp2  the MP2 energy.
** \ingroup CHKPT
*/
        double chkpt_rd_emp2(void)
        {
                double emp2;
                emp2 = _default_chkpt_lib_->rd_emp2();
                return emp2;
        }

/*!
** chkpt_wt_emp2(): Writes out the MP2 contribution to total energy.
**
** \param emp2 = the MP2 energy.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_emp2(double emp2)
        {
                _default_chkpt_lib_->wt_emp2(emp2);
        }
}
