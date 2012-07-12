/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
/*
**  CCDENSITY: Program to calculate the coupled-cluster one- and
**             two-particle densities.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libiwl/iwl.h>
#include <liboptions/liboptions.h>
#include <psi4-dec.h>
#include <cmath>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#include "globals.h"

namespace psi { namespace ccdensity {

void init_io(void);
void title(void);
void get_moinfo(void);
void get_frozen(void);
void get_params(Options& options);
void exit_io(void);
void onepdm(struct RHO_Params);
void sortone(struct RHO_Params);
void twopdm(void);
void energy(struct RHO_Params);
void resort_tei(void);
void resort_gamma(void);
void lag(struct RHO_Params rho_params);
void build_X(void);
void build_A(void);
void build_Z(void);
void relax_I(void);
void relax_D(struct RHO_Params rho_params);
void sortI(void);
void fold(struct RHO_Params rho_params);
void deanti(struct RHO_Params rho_params);
void add_ref_RHF(struct iwlbuf *);
void add_ref_ROHF(struct iwlbuf *);
void add_ref_UHF(struct iwlbuf *, struct iwlbuf *, struct iwlbuf *);
void add_core_ROHF(struct iwlbuf *);
void add_core_UHF(struct iwlbuf *, struct iwlbuf *, struct iwlbuf *);
void dump_RHF(struct iwlbuf *, struct RHO_Params rho_params);
void dump_ROHF(struct iwlbuf *, struct RHO_Params rho_params);
void dump_UHF(struct iwlbuf *, struct iwlbuf *, struct iwlbuf *, struct RHO_Params rho_params);
void kinetic(void);
void dipole(void);
void probable(void);
int **cacheprep_rhf(int level, int *cachefiles);
int **cacheprep_uhf(int level, int *cachefiles);
void cachedone_rhf(int **cachelist);
void cachedone_uhf(int **cachelist);
void setup_LR(struct RHO_Params);
void G_build(void);
void x_oe_intermediates(struct RHO_Params);
void x_onepdm(struct RHO_Params);
void x_te_intermediates(void);
void x_Gijkl(void);
void x_Gabcd(void);
void x_Gibja(void);
void x_Gijka(void);
void x_Gijab(void);
void x_Gciab(void);
void V_build_x(void);
void x_xi1(void);
void x_xi_zero(void);
void x_xi2(void);
void x_xi_oe_intermediates(void);
void G_norm(void);
void zero_onepdm(struct RHO_Params rho_params);
void zero_twopdm(void);
void get_rho_params(Options& options);
void get_td_params(Options& options);
void td_setup(struct TD_Params S);
void tdensity(struct TD_Params S);
void td_print(void);
void oscillator_strength(struct TD_Params *S);
void rotational_strength(struct TD_Params *S);
void ael(struct RHO_Params *rho_params);
void cleanup(void);
void td_cleanup(void);
void x_oe_intermediates_rhf(struct RHO_Params rho_params);
void x_te_intermediates_rhf(void);
void x_xi_intermediates(void);
void V_build(void);
void densgrid_RHF(void);

PsiReturnType ccdensity(Options& options)
{
  int i;
  int **cachelist, *cachefiles;
  struct iwlbuf OutBuf;
  struct iwlbuf OutBuf_AA, OutBuf_BB, OutBuf_AB;
  dpdfile2 D;
  double tval;

  init_io();
  title();
  /*  get_frozen(); */
  get_params( options );
  get_moinfo();
  get_rho_params(options);

  if ((moinfo.nfzc || moinfo.nfzv) && params.relax_opdm) {
    fprintf(outfile, "\n\tGradients/orbital relaxation involving frozen orbitals not yet available.\n");
    throw PsiException("ccdensity: error", __FILE__, __LINE__);
  }

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/
    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL,
             2, moinfo.occpi, moinfo.occ_sym, moinfo.virtpi, moinfo.vir_sym);

  }
  else if(params.ref == 2) { /** UHF **/
    cachelist = cacheprep_uhf(params.cachelev, cachefiles);

    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles,
             cachelist, NULL, 4, moinfo.aoccpi, moinfo.aocc_sym, moinfo.avirtpi,
             moinfo.avir_sym, moinfo.boccpi, moinfo.bocc_sym, moinfo.bvirtpi, moinfo.bvir_sym);
  }

  for (i=0; i<params.nstates; ++i) {

    /* CC_GLG will contain L, or R0*L + Zeta, if relaxed and zeta is available */
    /* CC_GL will contain L */
    setup_LR(rho_params[i]);

    /* Calculate Xi, put Xi in EOM_XI, and quit */
    if ( params.calc_xi ) {
      /* these intermediates go into EOM_TMP and are used to compute Xi;
         they may be reused to compute the excited-state density matrix */
      if (params.ref == 0) {
        x_oe_intermediates_rhf(rho_params[i]);
        x_te_intermediates_rhf();
      }
      else {
        x_oe_intermediates(rho_params[i]);
        x_te_intermediates();
      }
      x_xi_intermediates(); /*Xi intermediates put in EOM_TMP_XI */
      x_xi_zero(); /* make blank Xi */
      x_xi1();
      x_xi2();
      dpd_close(0);
      if(params.ref == 2) cachedone_uhf(cachelist);
      else cachedone_rhf(cachelist);
      free(cachefiles);
      cleanup();
      psio_close(EOM_TMP_XI,0); /* delete EOM_TMP_XI */
      psio_open(EOM_TMP_XI,PSIO_OPEN_NEW);
      exit_io();
      return Success;
    }

    /* compute ground state parts of onepdm or put zeroes there */
    if ( ((rho_params[i].L_irr == rho_params[i].G_irr) || (params.use_zeta)) ) {
      zero_onepdm(rho_params[i]);
      onepdm(rho_params[i]);
    }
    else
      zero_onepdm(rho_params[i]);

    /* if the one-electron excited-state intermediates are not already on disk (from a Xi
       calculation, compute them.  They are nearly all necessary to compute the excited-state
       onepdm. Then complete excited-state onepdm.*/
    if (!rho_params[i].R_ground) {
      x_oe_intermediates(rho_params[i]); /* change to x_oe_intermediates_rhf() when rho gets spin-adapted */
      x_onepdm(rho_params[i]);
    }

    /* begin construction of twopdm */
    if (!params.onepdm) {

      /* Compute intermediates for construction of ground-state twopdm */
      if ( (params.L_irr == params.G_irr) || (params.use_zeta) ) {
        V_build(); /* uses CC_GLG, writes tau2*L2 to CC_MISC */
        G_build(); /* uses CC_GLG, writes t2*L2 to CC_GLG */
      }

      /* Compute ground-state twopdm or ground-state-like contributions to the excited twodpm */
      if ( (params.L_irr == params.G_irr) || (params.use_zeta) )
        twopdm();
      else
        zero_twopdm();

      /* Compute intermediates for construction of excited-state twopdm */
      if (!params.ground) {
        x_te_intermediates(); /* change to x_te_intermediates_rhf() when rho gets spin-adapted */
        V_build_x(); /* uses CC_GL, writes t2*L2 to EOM_TMP */

        /* add in non-R0 parts of onepdm and twopdm */
        x_Gijkl();
        x_Gabcd();
        x_Gibja();
        x_Gijka();
        x_Gciab();
        x_Gijab();
      }
    }

    sortone(rho_params[i]); /* puts full 1-pdm into moinfo.opdm */
    if (!params.onepdm) {
      if(!params.aobasis) energy(rho_params[i]);

      kinetic(); /* puts kinetic energy integrals into MO basis */

      lag(rho_params[i]); /* builds the orbital lagrangian pieces, I */

      /* dpd_init(1, moinfo.nirreps, params.memory, 2, frozen.occpi, frozen.occ_sym,
         frozen.virtpi, frozen.vir_sym); */

      /*  if(moinfo.nfzc || moinfo.nfzv) {
          resort_gamma();
          resort_tei();
          } */

      build_X(); /* builds orbital rotation gradient X */
      build_A(); /* construct MO Hessian A */
      build_Z(); /* solves the orbital Z-vector equations */

      relax_I(); /* adds orbital response contributions to Lagrangian */

      if (params.relax_opdm) {
        relax_D(rho_params[i]); /* adds orbital response contributions to onepdm */
      }
      sortone(rho_params[i]); /* builds large moinfo.opdm matrix */
      sortI(); /* builds large lagrangian matrix I */
      fold(rho_params[i]);
      deanti(rho_params[i]);
    }

    /*  dpd_close(0); dpd_close(1); */

    if(params.ref == 0) { /** RHF **/

      iwl_buf_init(&OutBuf, PSIF_MO_TPDM, params.tolerance, 0, 0);

      add_core_ROHF(&OutBuf);
      add_ref_RHF(&OutBuf);

      // ==> One-Electron Properties <== //
      fprintf(outfile, "  ==> Properties: Root %d <==\n\n", i);
      dipole();

      if(params.onepdm_grid_dump) densgrid_RHF();
      dump_RHF(&OutBuf, rho_params[i]);

      iwl_buf_flush(&OutBuf, 1);
      iwl_buf_close(&OutBuf, 1);
    }
    else if(params.ref == 1) { /** ROHF **/

      iwl_buf_init(&OutBuf, PSIF_MO_TPDM, params.tolerance, 0, 0);

      add_core_ROHF(&OutBuf);
      add_ref_ROHF(&OutBuf);
      dump_ROHF(&OutBuf, rho_params[i]);

      iwl_buf_flush(&OutBuf, 1);
      iwl_buf_close(&OutBuf, 1);
    }
    else if(params.ref == 2) { /** UHF **/

      iwl_buf_init(&OutBuf_AA, PSIF_MO_AA_TPDM, params.tolerance, 0, 0);
      iwl_buf_init(&OutBuf_BB, PSIF_MO_BB_TPDM, params.tolerance, 0, 0);
      iwl_buf_init(&OutBuf_AB, PSIF_MO_AB_TPDM, params.tolerance, 0, 0);

      /*    add_core_UHF(&OutBuf_AA, &OutBuf_BB, &OutBuf_AB); */
      add_ref_UHF(&OutBuf_AA, &OutBuf_BB, &OutBuf_AB);
      dump_UHF(&OutBuf_AA, &OutBuf_BB, &OutBuf_AB, rho_params[i]);

      iwl_buf_flush(&OutBuf_AA, 1);
      iwl_buf_flush(&OutBuf_BB, 1);
      iwl_buf_flush(&OutBuf_AB, 1);
      iwl_buf_close(&OutBuf_AA, 1);
      iwl_buf_close(&OutBuf_BB, 1);
      iwl_buf_close(&OutBuf_AB, 1);
    }

    free_block(moinfo.opdm);

    psio_close(CC_TMP,0);   psio_open(CC_TMP,PSIO_OPEN_NEW);
    psio_close(EOM_TMP0,0); psio_open(EOM_TMP0,PSIO_OPEN_NEW);
    psio_close(EOM_TMP1,0); psio_open(EOM_TMP1,PSIO_OPEN_NEW);
    psio_close(CC_GLG,0);   psio_open(CC_GLG,PSIO_OPEN_NEW);
    psio_close(CC_GL,0);    psio_open(CC_GL,PSIO_OPEN_NEW);
    psio_close(CC_GR,0);    psio_open(CC_GR,PSIO_OPEN_NEW);
    if (!params.calc_xi) {
      psio_close(EOM_TMP,0);
      psio_open(EOM_TMP,PSIO_OPEN_NEW);
    }
  }

/*
  if ( params.ael && (params.nstates > 1) )
    ael(rho_params);
*/

  if(params.transition) {

    get_td_params(options);
    for(i=0; i < params.nstates; i++) {
      td_setup(td_params[i]);
      tdensity(td_params[i]);
      oscillator_strength(&(td_params[i]));
      if(params.ref == 0) {
        rotational_strength(&(td_params[i]));
      }
      td_cleanup();
    }
    td_print();
  }

  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup();
  exit_io();
  return Success;
}


// must be fixed with options for excited state densities
void init_io(void)
{
  int i, num_unparsed;
  char *argv_unparsed[100];;

  params.onepdm = 0;
  params.prop_all = 0;
  params.calc_xi = 0;
  params.restart = 0;
  params.use_zeta = 0;
  params.transition = 0;

/*
  for (i=1, num_unparsed=0; i<argc; ++i) {
    if(!strcmp(argv[i], "--onepdm")) {
      params.onepdm = 1; // generate ONLY the onepdm (for one-electron properties)
    }
    else if (!strcmp(argv[i],"--use_zeta")) {
      params.use_zeta = 1;
      params.ground = 0;
      params.restart = 1;
    }
    else if (!strcmp(argv[i],"--calc_xi")) {
      params.calc_xi = 1;
      params.ground = 0;
      params.restart = 0;
    }
    else if (!strcmp(argv[i],"--prop_all")) {
      params.prop_all = 1;
    }
    else if (!strcmp(argv[i],"--transition")) {
      params.transition = 1;
      params.relax_opdm = 0;
      params.ground = 0;
    }
    else {
      argv_unparsed[num_unparsed++] = argv[i];
    }
  }
*/

  tstart();

  for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i,PSIO_OPEN_OLD);
  // erase old files
  psio_close(CC_GR,0);
  psio_close(CC_GL,0);
  psio_close(EOM_TMP0,0);
  psio_open(CC_GR,PSIO_OPEN_NEW);
  psio_open(CC_GL,PSIO_OPEN_NEW);
  psio_open(EOM_TMP0,PSIO_OPEN_NEW);
}

void title(void)
{
  fprintf(outfile, "\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*        CCDENSITY       *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\n");
}

void exit_io(void)
{
  int i;

  /* delete temporary EOM files */
  psio_close(EOM_TMP0,0);
  psio_close(EOM_TMP1,0);
  psio_close(CC_GLG,0);
  psio_open(EOM_TMP0,PSIO_OPEN_NEW);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);
  psio_open(CC_GLG,PSIO_OPEN_NEW);
  if (!params.calc_xi) {
    psio_close(EOM_TMP,0);
    psio_open(EOM_TMP,PSIO_OPEN_NEW);
  }
  if (params.use_zeta) { /* we're done with Xi amplitudes */
    psio_close(EOM_XI,0);
    psio_open(EOM_XI,PSIO_OPEN_NEW);
  }

  /* Close all dpd data files here */
  for(i=CC_MIN; i <= CC_MAX; i++) psio_close(i,1);

  tstop();
}

}} // namespace psi::ccdensity
