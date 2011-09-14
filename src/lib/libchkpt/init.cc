/*!
  \file
  \ingroup CHKPT
*/

#include <cstdio>
#include <cstdlib>
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>
#include <boost/shared_ptr.hpp>

using namespace psi;

/* Definition of global data */
boost::shared_ptr<Chkpt> psi::_default_chkpt_lib_;

Chkpt::Chkpt(psi::PSIO *psioObject, int status, std::string use_prefix) : psio(psioObject)
{
    char *prefix;

    psio->open(PSIF_CHKPT, status);

    if (use_prefix.empty()) {
        if(psio->tocscan(PSIF_CHKPT, "Default prefix") != NULL) {
            prefix = rd_prefix();
            set_prefix(prefix);
            free(prefix);
        }
        else {
            set_prefix("");
            commit_prefix();  /* assume no default prefix existed in PSIF_CHKPT */
        }
    }
    else {
        set_prefix(use_prefix.c_str());
    }
}

Chkpt::Chkpt(psi::PSIO& psioObject, int status, std::string use_prefix) : psio(&psioObject)
{
    char *prefix;

    psio->open(PSIF_CHKPT, status);

    if (use_prefix.empty()) {
        if(psio->tocscan(PSIF_CHKPT, "Default prefix") != NULL) {
            prefix = rd_prefix();
            set_prefix(prefix);
            free(prefix);
        }
        else {
            set_prefix("");
            commit_prefix();  /* assume no default prefix existed in PSIF_CHKPT */
        }
    }
    else {
        set_prefix(use_prefix.c_str());
    }
}

Chkpt::Chkpt(boost::shared_ptr<PSIO> psioObject, int status, std::string use_prefix) : psio(psioObject.get())
{
    char *prefix;

    psio->open(PSIF_CHKPT, status);

    if (use_prefix.empty()) {
        if(psio->tocscan(PSIF_CHKPT, "Default prefix") != NULL) {
            prefix = rd_prefix();
            set_prefix(prefix);
            free(prefix);
        }
        else {
            set_prefix("");
            commit_prefix();  /* assume no default prefix existed in PSIF_CHKPT */
        }
    }
    else {
        set_prefix(use_prefix.c_str());
    }
}

boost::shared_ptr<Chkpt> Chkpt::shared_object()
{
    return _default_chkpt_lib_;
}

void Chkpt::rehash() {
  psio->rehash(PSIF_CHKPT);
}

extern "C" {
/*!
**  chkpt_init()  Initializes the checkpoint file for other chkpt_
**    functions to perform their duties.
**
**  arguments:
**    int status: boolean indicating if the chkpt file should be
**                initialized (PSIO_OPEN_NEW) or the old chkpt
**                file should be used (PSIO_OPEN_OLD).
**
**  returns: zero.  Perhaps this will change some day.
**  \ingroup CHKPT
*/
    int chkpt_init(int status)
    {
        if (_default_chkpt_lib_.get() == 0) {
                        boost::shared_ptr<Chkpt> temp(new Chkpt(_default_psio_lib_, status));
            _default_chkpt_lib_ = temp;
            if (_default_chkpt_lib_ == 0) {
                fprintf(stderr, "LIBCHKPT::init() -- failed to allocate memory\n");
                exit(PSI_RETURN_FAILURE);
            }
        }
        return 0;
    }
}

