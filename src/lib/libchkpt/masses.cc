/*!
    \file
*/

#include <cstdio>
#include <cstdlib>
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libchkpt/chkpt.h>

namespace psi {

double *Chkpt::rd_masses(void)
{
	int natom;
	double *masses;
	char *keyword;
	keyword = build_keyword("Atomic masses");

	natom = rd_natom();
	masses = array<double>(natom);

	psio->read_entry(PSIF_CHKPT, keyword, (char *) masses, 
		natom*sizeof(double));

	free(keyword);
	return masses;
}

void Chkpt::wt_masses(double *masses)
{
	int natom;
	char *keyword;
	keyword = build_keyword("Atomic masses");

	natom = rd_natom();

	psio->write_entry(PSIF_CHKPT, keyword, (char *) masses, 
		natom*sizeof(double));

	free(keyword);
}

/*!
** chkpt_rd_masses()
** Reads the atomic masses from the checkpoint file.
**
** arguments: none
**
** returns: 
**   double *masses: An array of the masses
*/

/*!
** chkpt_wt_masses()
** Writes the atomic masses to the checkpoint file.
**
** \param double *masses: An array of the masses
**
** returns: nothing
*/

}

