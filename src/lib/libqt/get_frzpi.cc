/*!
  \file
  \brief Get frozen core/virtuals per irrep
  \ingroup QT
*/

#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libchkpt/chkpt.h>
#include <psifiles.h>
#include "qt.h"

namespace psi {

#define PSIO_INIT if (!psio_state()) { \
    psio_init(); psio_ipv1_config(); \
    need_to_init_psio = 1; \
  }

#define CHKPT_INIT(n) if (!psio_open_check(PSIF_CHKPT)) { \
    chkpt_init(n); \
    need_to_init_chkpt = 1; \
  }

#define PSIO_DONE if (need_to_init_psio) \
    psio_done();

#define CHKPT_DONE if (need_to_init_chkpt) \
    chkpt_close();

/*!
** get_frzcpi(): Get frozen core per irrep array.
**
** Parameters: none
**
** Returns: pointer to int array
**
** \ingroup QT
*/
int* get_frzcpi()
{
  int errcod;
  int nirreps, nirr;
  int need_to_init_psio = 0;
  int need_to_init_chkpt = 0;
  int if_exists;
  int* frzcpi;

PSIO_INIT
CHKPT_INIT(PSIO_OPEN_OLD);
  nirreps = chkpt_rd_nirreps();
  
  frzcpi = init_int_array(nirreps);
  
  if_exists = ip_exist("FROZEN_DOCC",0);
  if (if_exists) {
    errcod = ip_int_array("FROZEN_DOCC",frzcpi,nirreps);
  }
  else {
    free(frzcpi);
    frzcpi = chkpt_rd_frzcpi();
  }

CHKPT_DONE
PSIO_DONE

  return frzcpi;
}


/*!
** get_frzvpi(): Get frozen virtuals per irrep array.
**
** Parameters: none
**
** Returns: pointer to int array
**
** \ingroup QT
*/
int* get_frzvpi()
{
  int errcod;
  int nirreps, nirr;
  int need_to_init_psio = 0;
  int need_to_init_chkpt = 0;
  int if_exists;
  int* frzvpi;

PSIO_INIT
CHKPT_INIT(PSIO_OPEN_OLD);
  nirreps = chkpt_rd_nirreps();
  
  frzvpi = init_int_array(nirreps);
  
  if_exists = ip_exist("FROZEN_UOCC",0);
  if (if_exists) {
    errcod = ip_int_array("FROZEN_UOCC",frzvpi,nirreps);
  }
  else {
    free(frzvpi);
    frzvpi = chkpt_rd_frzvpi();
  }

CHKPT_DONE
PSIO_DONE

  return frzvpi;
}

}

