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

char **Chkpt::rd_felement(void)
{
        char **label;
        int nallatom, i;
        psio_address ptr;
        char *keyword;
        keyword = build_keyword("Full atom labels");

        nallatom = rd_nallatom();

        label = (char **)malloc(nallatom*sizeof(const char*));
        for(i=0; i < nallatom; i++)
                label[i] = (char *) malloc(MAX_ELEMNAME*sizeof(char));

        ptr = PSIO_ZERO;
        for(i=0; i < nallatom; i++)
                psio->read(PSIF_CHKPT, keyword, (char *) label[i],
                                        MAX_ELEMNAME*sizeof(char), ptr, &ptr);

        free(keyword);
        return label;
}

void Chkpt::wt_felement(char ** const label)
{
        int nallatom, i;
        psio_address ptr;
        char *keyword;
        keyword = build_keyword("Full atom labels");

        nallatom = rd_nallatom();

        ptr = PSIO_ZERO;
        for(i=0; i < nallatom; i++)
                psio->write(PSIF_CHKPT, keyword, (char *) label[i],
                                        MAX_ELEMNAME*sizeof(char), ptr, &ptr);

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_felement():  Reads in element labels including dummy atoms
**
**   takes no arguments.
**
**   returns: char **label element label matrix
** \ingroup CHKPT
*/
        char **chkpt_rd_felement(void)
        {
                return _default_chkpt_lib_->rd_felement();
        }

/*!
** chkpt_wt_felement():  Writes out element labels including dummy atoms
**
** arguments:
**   \param label = element label matrix.
**
** returns: none
** \ingroup CHKPT
*/
        void chkpt_wt_felement(char ** const label)
        {
                _default_chkpt_lib_->wt_felement(label);
        }
}
