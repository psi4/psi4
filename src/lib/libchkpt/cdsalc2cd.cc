/*!
  \file cdsalc2cd.cc
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include <psifiles.h>
#include <libpsio/psio.hpp>
extern "C" {
	#include <libchkpt/chkpt.h>
}
#include <libchkpt/chkpt.hpp>

using namespace psi;

double **Chkpt::rd_cdsalc2cd(void)
{
	const int num_cd = 3*rd_natom();
	double **cdsalc2cd = matrix<double>(num_cd,num_cd);
	psio_address ptr = PSIO_ZERO;
	char *keyword = build_keyword("cartdisp->SALC matrix");

	psio_read(PSIF_CHKPT, keyword, (char *) cdsalc2cd[0], num_cd*num_cd*sizeof(double), ptr, &ptr);

	free(keyword);
	return cdsalc2cd;
}

void Chkpt::wt_cdsalc2cd(const double **cdsalc2cd)
{
	const int num_cd = 3*rd_natom();
	psio_address ptr = PSIO_ZERO;
	char *keyword = build_keyword("cartdisp->SALC matrix");

	psio_write(PSIF_CHKPT, keyword, (char *) cdsalc2cd[0], num_cd*num_cd*sizeof(double), ptr, &ptr);

	free(keyword);
}

extern "C" {
/*!
** chkpt_rd_cdsalc2cd(): Read in (normalized) SALCs of cartesian displacements
**
** takes no arguments.
**
** returns: cdsalc2cd = A natom*3 by natom*3 blocked matrix of doubles. columnts correpond to symmetry-blocked SALCs
** 
** \ingroup (CHKPT)
*/
	double **chkpt_rd_cdsalc2cd(void)
	{
		double **cdsalc2cd;
		cdsalc2cd = _default_chkpt_lib_->rd_cdsalc2cd();
		return cdsalc2cd;
	}


/*!
** chkpt_wt_cdsalc2cd(): Writes out (normalized) SALCs of cartesian displacements
**
** \param cdsalc2cd = A natom*3 by natom*3 blocked matrix of doubles. columnts correpond to symmetry-blocked SALCs
**
** returns: none
**
** \ingroup (CHKPT)
*/
	void chkpt_wt_cdsalc2cd(const double **cdsalc2cd)
	{
		_default_chkpt_lib_->wt_cdsalc2cd(cdsalc2cd);
	}
}
