/*!
  \file irr_labs.cc
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

char **Chkpt::rd_irr_labs(void)
{
	int i,nirreps;
	char **irr_labs;
	psio_address ptr;
	char *keyword;
	keyword = build_keyword("Irrep labels");

	nirreps = rd_nirreps();

	ptr = PSIO_ZERO;
	irr_labs = (char **)malloc(sizeof(char *)*nirreps);
	for(i=0;i<nirreps;i++) {
		irr_labs[i] = (char *) malloc(4*sizeof(char));
		psio->read(PSIF_CHKPT, keyword, (char *) irr_labs[i],  4*sizeof(char), ptr, &ptr);
		irr_labs[i][3] = '\0';
	}

	free(keyword);
	return irr_labs;
}

void Chkpt::wt_irr_labs(char **irr_labs)
{
	int i,nirreps;
	psio_address ptr;
	char *keyword;
	keyword = build_keyword("Irrep labels");

	nirreps = rd_nirreps();

	ptr = PSIO_ZERO;
	for(i=0;i<nirreps;i++)
		psio->write(PSIF_CHKPT, keyword, (char *) irr_labs[i], 4*sizeof(char), ptr, &ptr);

	free(keyword);
}

extern "C" {
/*!
** chkpt_rd_irr_labs(): Read in the symmetry labels for all irreps in the 
**   point group in which the molecule is considered.
**
**   takes no arguments.
**
**   returns: irr_labs =  an array of labels (strings) which denote
**      the irreps for the point group	in which the molecule is considered,
**      _regardless_ of whether there exist any symmetry orbitals which 
**      transform as that irrep.  
** \ingroup (CHKPT)
*/
	char **chkpt_rd_irr_labs(void)
	{
		return _default_chkpt_lib_->rd_irr_labs();
	}

/*!
** chkpt_wt_irr_labs(): Write out the symmetry labels for all irreps in the 
**   point group in which the molecule is considered.
**
**  arguemnts: 
**  \param irr_labs = an array of labels (strings) which denote
**      the irreps for the point group	in which the molecule is considered,
**      _regardless_ of whether there exist any symmetry orbitals which 
**      transform as that irrep.  
** \ingroup (CHKPT)
*/
	void chkpt_wt_irr_labs(char **irr_labs)
	{
		_default_chkpt_lib_->wt_irr_labs(irr_labs);
	}
}

