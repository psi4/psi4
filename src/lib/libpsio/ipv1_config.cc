#include <cstdlib>
#include <cstring>
#include <psifiles.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <psi4-dec.h>

namespace psi {

  /*!
   This function initializes the libpsio library using the .psirc file + input file
   */
  int psiopp_ipv1_config(PSIO *psio_obj) {
    int i, unit, nkwds;
    char* ip_token;
    char *name;
    char* userhome;
    FILE* psirc;
    char* filename;
    int errcod;

    PSIO::_error_exit_code_ = PSI_RETURN_FAILURE;

    /* Open user's general .psirc file, if exists */
    userhome = getenv("HOME");
    filename = (char*) malloc((strlen(userhome)+8)*sizeof(char));
    sprintf(filename, "%s%s", userhome, "/.psirc");
    psirc = fopen(filename, "r");
    if (psirc != NULL) {
      ip_append(psirc, stdout);
      fclose(psirc);
    }
    free(filename);

    /*
     implement some default PSI3 behavior:
     1) checkpoint file should by default be in "./"
     2) all other files should go to "/tmp/"
     3) default name is psi_file_prefix
     4) 1 volume
     */
    for (i=1; i<=PSIO_MAXVOL; ++i) {
      char kwd[20];
      sprintf(kwd, "VOLUME%u", i);
      psio_obj->filecfg_kwd("DEFAULT", kwd, PSIF_CHKPT, "./");
      psio_obj->filecfg_kwd("DEFAULT", kwd, -1, "/tmp/");
    }
    psio_obj->filecfg_kwd("DEFAULT", "NAME", -1, psi_file_prefix);
    psio_obj->filecfg_kwd("DEFAULT", "NVOLUME", -1, "1");

    /*
     transfer the necessary keywords from IPV1 to PSIO **LAST** (so that the defaults are overridden as necessary)
     */
    /* need "NAME", "NVOLUME", and "VOLUMEX" -- total of 2+PSIO_MAXVOL keywords */
    nkwds = 2+PSIO_MAXVOL;
    char** kwds = (char**) malloc(nkwds*sizeof(char*));
    kwds[0] = strdup("NAME");
    kwds[1] = strdup("NVOLUME");
    for (i=1; i<=PSIO_MAXVOL; ++i) {
      char kwd[20];
      sprintf(kwd, "VOLUME%u", i);
      kwds[1+i] = strdup(kwd);
    }

    /* allocate ip_token
     conservative estimate for its length = strlen(gprgid())+80 */
    char *module_name = module.gprgid();
    ip_token = (char*) malloc( (strlen(module_name)+80)*sizeof(char));
    name = (char*) malloc( 80 * sizeof(char));

    for (i=0; i<nkwds; ++i) {
      const char* kwd = kwds[i];

      /* unit and program specific */
      for (unit=0; unit<PSIO_MAXUNIT; ++unit) {
        sprintf(ip_token, ":%s:FILES:FILE%u:%s", module_name, unit, kwd);
        errcod = ip_data(ip_token, "%s", name, 0);
        if (errcod == IPE_OK) {
          psio_obj->filecfg_kwd(module_name, kwd, unit, name);
        }
      }

      /* program specific */
      sprintf(ip_token, ":%s:FILES:DEFAULT:%s", module_name, kwd);
      errcod = ip_data(ip_token, "%s", name, 0);
      if (errcod == IPE_OK) {
        psio_obj->filecfg_kwd(module_name, kwd, -1, name);
      }

      /* unit specific in PSI section */
      for (unit=0; unit<PSIO_MAXUNIT; ++unit) {
        sprintf(ip_token, ":PSI:FILES:FILE%u:%s", unit, kwd);
        errcod = ip_data(ip_token, "%s", name, 0);
        if (errcod == IPE_OK) {
          psio_obj->filecfg_kwd("PSI", kwd, unit, name);
        }
      }

      /* in PSI section */
      sprintf(ip_token, ":PSI:FILES:DEFAULT:%s", kwd);
      errcod = ip_data(ip_token, "%s", name, 0);
      if (errcod == IPE_OK) {
        psio_obj->filecfg_kwd("PSI", kwd, -1, name);
      }

      /* unit specific in DEFAULT section */
      for (unit=0; unit<PSIO_MAXUNIT; ++unit) {
        sprintf(ip_token, ":DEFAULT:FILES:FILE%u:%s", unit, kwd);
        errcod = ip_data(ip_token, "%s", name, 0);
        if (errcod == IPE_OK) {
          psio_obj->filecfg_kwd("DEFAULT", kwd, unit, name);
        }
      }

      /* in DEFAULT section */
      sprintf(ip_token, ":DEFAULT:FILES:DEFAULT:%s", kwd);
      errcod = ip_data(ip_token, "%s", name, 0);
      if (errcod == IPE_OK) {
        psio_obj->filecfg_kwd("DEFAULT", kwd, -1, name);
      }
    }
    delete [] module_name;

    for (i=0; i<nkwds; ++i) {
      free(kwds[i]);
    }
    free(kwds);
    free(ip_token);
    free(name);

    return 0;
  }

  int psio_ipv1_config() {
      return psiopp_ipv1_config(_default_psio_lib_);
  }

}

