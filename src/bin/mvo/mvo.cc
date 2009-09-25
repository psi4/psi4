/*! \defgroup MVO mvo: Compute modified virtual orbitals, etc */

/*! 
** \file
** \ingroup MVO
** \brief Compute modified virtual orbitals, etc.
**
** This program reads the frozen core operator from TRANSQT 
** output ** and diagonalizes the virtual-virtual block to make 
** Modified Virtual Orbitals (MVO's) as described in 
** C. W. Bauschlicher, J. Chem. Phys. 72, 880 (1980)
**
** C. David Sherrill
** Georgia Institute of Technology
** March 9, 2001
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <libiwl/iwl.h>
#include "MOInfo.h"
#include "params.h"
#include "globals.h"


/* First definitions of globals */
extern "C" {
  FILE *infile, *outfile;
  char *psi_file_prefix;
}

namespace psi { namespace mvo {

int *ioff;
struct MOInfo moinfo;
struct Params params;

/* Max length of ioff array */
#define IOFF_MAX 32641


#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))


/* Function prototypes */
void init_io(int argc, char *argv[]);
void title();
void init_ioff();
void get_parameters();
void print_parameters();
void get_moinfo();
void cleanup();
void exit_io();
void get_reorder_array(void);
void get_fzc_operator(void);
void get_mvos(void);
void get_canonical(void);
void get_mp2nos(void);
/* void get_unos(void); */
}} // end namespace psi::mvo

using namespace psi::mvo;

int main(int argc, char *argv[])
{
  params.print_lvl = 1;
  init_io(argc, argv);
  title();
  get_parameters();
  get_moinfo();
  print_parameters();
  get_reorder_array();
  init_ioff();
  if (params.mp2nos) {
    get_mp2nos();
  }
  else if (params.unos) {
    /* get_unos(); */
  }
  else if (params.canonical) {
    get_canonical();
  }
  else {
    get_fzc_operator();
    get_mvos();
  }
  cleanup();
  exit_io();
  return(0);
}

namespace psi { namespace mvo {

void init_io(int argc, char *argv[])
{
  int i;
  int num_extra_args = 0;
  char **extra_args;

  extra_args = (char **) malloc(argc*sizeof(char *));

   for (i=1; i<argc; i++) {
     if (strcmp(argv[i], "--quiet") == 0) {
       params.print_lvl = 0;
     }
     else {
       extra_args[num_extra_args++] = argv[i];
     }
   }
 
  psi_start(&infile,&outfile,&psi_file_prefix,num_extra_args,extra_args,0); 
  if (params.print_lvl) tstart(outfile);
  ip_cwk_add(":MVO");

  psio_init(); psio_ipv1_config();
  free(extra_args);
}


void title(void)
{
  if (params.print_lvl) {
    fprintf(outfile, "\n");
    fprintf(outfile,"\t**************************************************\n");
    fprintf(outfile,"\t*     MVO: Obtain Modified Virtual Orbitals      *\n");
    fprintf(outfile,"\t*                                                *\n");
    fprintf(outfile,"\t*               C. David Sherrill                *\n");
    fprintf(outfile,"\t*                  March  2001                   *\n");
    fprintf(outfile,"\t**************************************************\n");
    fprintf(outfile, "\n");
  }
}

void init_ioff(void)
{
  int i;
  ioff = (int *) malloc(IOFF_MAX * sizeof(int));
  if(ioff == NULL) {
    fprintf(stderr, "(transqt): error malloc'ing ioff array\n");
    exit(0);
  }
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) {
    ioff[i] = ioff[i-1] + i;
  }
}

void exit_io(void)
{
  psio_done();
  if (params.print_lvl) tstop(outfile);
  psi_stop(infile,outfile,psi_file_prefix);
}


void get_parameters(void)
{
  int errcod;

  errcod = ip_string("WFN", &(params.wfn),0);
  if(errcod == IPE_KEY_NOT_FOUND) {
    params.wfn = (char *) malloc(sizeof(char)*5);
    strcpy(params.wfn, "CCSD");
  }

  params.h_fzc_file = PSIF_OEI;
  errcod = ip_data("FZC_FILE","%d",&(params.h_fzc_file),0);

  params.print_mos = 0;
  errcod = ip_boolean("PRINT_MOS", &(params.print_mos),0);

  params.print_lvl = 1;
  errcod = ip_data("PRINT", "%d", &(params.print_lvl),0);

  params.oei_erase = 0;
  errcod = ip_boolean("OEI_ERASE",&(params.oei_erase),0); 

  params.fzc = 1;
  errcod = ip_boolean("FZC",&(params.fzc),0);

  /* remove restricted docc from RAS 1 ? */
  params.del_restr_docc = 1;
  errcod = ip_boolean("DELETE_RESTR_DOCC",&(params.del_restr_docc),0);

  params.mp2nos = 0;
  errcod = ip_boolean("MP2NOS",&(params.mp2nos),0);

  params.unos = 0;
  errcod = ip_boolean("UNOS",&(params.unos),0);

  if (strcmp(params.wfn, "CI")==0 || strcmp(params.wfn, "DETCAS")==0 ||
      strcmp(params.wfn, "DETCI")==0) {
    params.ras_type = 1;
  }
  else {
    params.ras_type = 0;
  }

  params.fzc_fock_coeff = 1.0;
  errcod = ip_data("FZC_FOCK_COEFF", "%lf", &(params.fzc_fock_coeff),0);

  params.fock_coeff = 0.0;
  errcod = ip_data("FOCK_COEFF", "%lf", &(params.fock_coeff),0);

  params.ivo = 0;
  errcod = ip_boolean("IVO", &(params.ivo), 0);

  params.canonical = 0;
  errcod = ip_boolean("CANONICAL", &(params.canonical), 0);
  return;

}


void print_parameters(void)
{

  if (params.print_lvl) {
    fprintf(outfile,"\n");
    fprintf(outfile,"\tInput Parameters:\n");
    fprintf(outfile,"\t-----------------\n");
    fprintf(outfile,"\tWavefunction           =  %s\n", params.wfn);
    fprintf(outfile,"\tPrint MOs              =  %s\n", 
                                (params.print_mos ? "Yes": "No"));
    fprintf(outfile,"\tErase OEI file         =  %s\n", 
                                (params.oei_erase ? "Yes": "No"));
    fprintf(outfile,"\tFrozen core            =  %s\n", 
                                (params.fzc ? "Yes": "No"));
    fprintf(outfile,"\tDelete Restricted Docc =  %s\n",
                                (params.del_restr_docc ? "Yes" : "No"));
    fprintf(outfile,"\tIVO                    =  %s\n", 
                                (params.ivo ? "Yes": "No"));
    fprintf(outfile,"\tCanonical              =  %s\n",
		                  (params.canonical ? "Yes": "No"));
    fprintf(outfile,"\tFrozen Core OEI file   =  %d\n", 
                                params.h_fzc_file);
    fprintf(outfile,"\tfzc_fock_coeff         =  %lf\n", 
                                params.fzc_fock_coeff);
    fprintf(outfile,"\tfock_coeff             =  %lf\n", 
                                params.fock_coeff);
    fprintf(outfile,"\tPrint Level            =  %d\n", params.print_lvl);
    fflush(outfile);
  }
  
  return;
}

void get_fzc_operator(void)
{

  moinfo.fzc_operator = (double *) init_array(moinfo.fzc_op_size);
  /*
  iwl_rdone(params.h_fzc_file, moinfo.fzc_operator, &(moinfo.efzc),
            ioff, moinfo.nmo, 0, 0, params.oei_erase,
            (params.print_lvl>4), outfile);
  */
  iwl_rdone(params.h_fzc_file, PSIF_MO_FZC, moinfo.fzc_operator, 
            moinfo.fzc_op_size, params.oei_erase,
            (params.print_lvl>4), outfile);
  
}



void get_moinfo(void)
{
  int i,j,k,h,errcod,size,row,col,p,q,offset,first_offset,last_offset,warned;
  int *tmpi, nopen;
  double **tmpmat, **so2ao;
  double N; /* number of valence electrons */

  chkpt_init(PSIO_OPEN_OLD);
  moinfo.nmo = chkpt_rd_nmo();
  moinfo.nso = chkpt_rd_nso();
  moinfo.nao = chkpt_rd_nao();
  moinfo.nirreps = chkpt_rd_nirreps();
  moinfo.iopen = chkpt_rd_iopen();
  moinfo.labels = chkpt_rd_irr_labs();
  moinfo.sopi   = chkpt_rd_sopi();
  moinfo.orbspi = chkpt_rd_orbspi();
  moinfo.clsdpi = chkpt_rd_clsdpi();
  moinfo.openpi = chkpt_rd_openpi();

  moinfo.fzc_op_size = (moinfo.nmo * (moinfo.nmo + 1)) / 2;

  moinfo.sosym = init_int_array(moinfo.nso);
  for (i=0,k=0; i<moinfo.nirreps; i++) {
    for (j=0; j<moinfo.sopi[i]; j++,k++) {
      moinfo.sosym[k] = i;
    }
  }
  
  moinfo.orbsym = init_int_array(moinfo.nmo);
  for (i=0,k=0; i<moinfo.nirreps; i++) {
    for (j=0; j<moinfo.orbspi[i]; j++,k++) {
      moinfo.orbsym[k] = i;
    }
  }

  moinfo.frdocc = init_int_array(moinfo.nirreps);
  moinfo.fruocc = init_int_array(moinfo.nirreps);
  errcod = ip_int_array("FROZEN_DOCC", moinfo.frdocc, moinfo.nirreps);
  errcod = ip_int_array("FROZEN_UOCC", moinfo.fruocc, moinfo.nirreps);

  moinfo.rstrdocc = init_int_array(moinfo.nirreps);
  moinfo.rstruocc = init_int_array(moinfo.nirreps);
  errcod = ip_int_array("RESTRICTED_DOCC",moinfo.rstrdocc,moinfo.nirreps);
  errcod = ip_int_array("RESTRICTED_UOCC",moinfo.rstruocc,moinfo.nirreps);
                                                                                
  /*
  if (params.treat_cor_as_fzc) {
    for (i=0; i<moinfo.nirreps; i++) {
      moinfo.frdocc[i] += moinfo.rstrdocc[i];
      moinfo.rstrdocc[i] = 0;
    }
  }
  */

  if (!params.fzc) {
    for (i=0; i<moinfo.nirreps; i++) {
      moinfo.frdocc[i] = 0;
    }
  }

  moinfo.nfzc = 0;
  moinfo.nfzv = 0;
  for (i=0; i<moinfo.nirreps; i++) {
    moinfo.nfzc += moinfo.frdocc[i];
    moinfo.nfzv += moinfo.fruocc[i];
  }

  moinfo.ndocc = 0;
  tmpi = init_int_array(moinfo.nirreps);
  errcod = ip_int_array("DOCC", tmpi, moinfo.nirreps);
  if (errcod == IPE_OK) {
    for (i=0,warned=0; i<moinfo.nirreps; i++) {
      if (tmpi[i] != moinfo.clsdpi[i] && !warned) {
        fprintf(outfile, "\tWarning: DOCC doesn't match PSIF_CHKPT\n");
        warned = 1;
      }
      moinfo.clsdpi[i] = tmpi[i];
      moinfo.ndocc += tmpi[i];
    }
  }
  else {
    for (i=0; i<moinfo.nirreps; i++) {
      moinfo.ndocc += moinfo.clsdpi[i];
    }
  }

  moinfo.nsocc = 0;
  errcod = ip_int_array("SOCC", tmpi, moinfo.nirreps);
  if (errcod == IPE_OK) {
    for (i=0,warned=0; i<moinfo.nirreps; i++) {
      if (tmpi[i] != moinfo.openpi[i] && !warned) {
        fprintf(outfile, "\tWarning: SOCC doesn't match PSIF_CHKPT\n");
        warned = 1;
      }
      moinfo.openpi[i] = tmpi[i];
      moinfo.nsocc += tmpi[i];
    }
  }

  moinfo.virtpi = init_int_array(moinfo.nirreps);
  for(i=0; i < moinfo.nirreps; i++) {
    moinfo.virtpi[i] = moinfo.orbspi[i]-moinfo.clsdpi[i]-moinfo.openpi[i];
  }

  if (params.print_lvl) {
    fprintf(outfile,"\n\tCheckpoint file parameters:\n");
    fprintf(outfile,"\t------------------\n");
    fprintf(outfile,"\tNumber of irreps = %d\n",moinfo.nirreps);
    fprintf(outfile,"\tNumber of SOs    = %d\n",moinfo.nso);
    fprintf(outfile,"\tNumber of MOs    = %d\n\n",moinfo.nmo);
    fprintf(outfile,
        "\tLabel\t# SOs\t# MOs\t# FZDC\t# DOCC\t# SOCC\t# VIRT\t# FZVR\n");
    fprintf(outfile,
        "\t-----\t-----\t-----\t------\t------\t------\t------\t------\n");
    for(i=0; i < moinfo.nirreps; i++) {
        fprintf(outfile,
           "\t %s\t   %d\t   %d\t    %d\t    %d\t    %d\t    %d\t    %d\n",
           moinfo.labels[i],moinfo.sopi[i],moinfo.orbspi[i],moinfo.frdocc[i],
           moinfo.clsdpi[i],moinfo.openpi[i],moinfo.virtpi[i],
           moinfo.fruocc[i]);
    }
    fflush(outfile);
  }


  /*
     Construct first and last index arrays for SOs: this defines the first
     absolute orbital index and last absolute orbital
     index for each irrep.  When there are no orbitals for an irrep, the
     value is -1 for first[] and -2 for last[].  Note that there must be
     basis functions in the first irrep (i.e. totally symmetric) for this to 
     work.
  */
  moinfo.first_so = init_int_array(moinfo.nirreps);
  moinfo.last_so = init_int_array(moinfo.nirreps);
  for(h=0; h < moinfo.nirreps; h++) {
    moinfo.first_so[h] = -1;
    moinfo.last_so[h] = -2;
  }
  first_offset = 0;
  last_offset = moinfo.sopi[0] - 1; 
  moinfo.first_so[0] = first_offset;
  moinfo.last_so[0] = last_offset;
  for(h=1; h < moinfo.nirreps; h++) {
    first_offset += moinfo.sopi[h-1];
    last_offset += moinfo.sopi[h];
    if(moinfo.sopi[h]) {
        moinfo.first_so[h] = first_offset;
        moinfo.last_so[h] = last_offset;
    }
  }
  
  /*
     Construct first and last index arrays: this defines the first
     absolute orbital index (Pitzer ordering) and last absolute orbital
     index for each irrep.  When there are no orbitals for an irrep, the
     value is -1 for first[] and -2 for last[].  Note that there must be
     orbitals in the first irrep (i.e. totally symmetric) for this to work.
  */
  moinfo.first = init_int_array(moinfo.nirreps);
  moinfo.last = init_int_array(moinfo.nirreps);
  for(h=0; h < moinfo.nirreps; h++) {
    moinfo.first[h] = -1;
    moinfo.last[h] = -2;
  }
  first_offset = 0;
  last_offset = moinfo.orbspi[0] - 1; 
  moinfo.first[0] = first_offset;
  moinfo.last[0] = last_offset;
  for(h=1; h < moinfo.nirreps; h++) {
    first_offset += moinfo.orbspi[h-1];
    last_offset += moinfo.orbspi[h];
    if(moinfo.orbspi[h]) {
      moinfo.first[h] = first_offset;
      moinfo.last[h] = last_offset;
    }
  }
  /*
     Construct first and last active index arrays: this defines the first
     absolute orbital index (Pitzer ordering) and last absolute orbital
     index for each irrep, excluding frozen orbitals.  When there are no
     orbitals for an irrep, the value is -1 for first[] and -2 for last[].
     Note that there must be orbitals in the first irrep (i.e. totally
     symmetric) for this to work.  
  */
  moinfo.fstact = init_int_array(moinfo.nirreps);
  moinfo.lstact = init_int_array(moinfo.nirreps);
  for(h=0; h < moinfo.nirreps; h++) {
    moinfo.fstact[h] = -1;
    moinfo.lstact[h] = -2;
  }
  first_offset = moinfo.frdocc[0];
  last_offset = moinfo.orbspi[0] - moinfo.fruocc[0] - 1; 
  moinfo.fstact[0] = first_offset;
  moinfo.lstact[0] = last_offset;
  for(h=1; h < moinfo.nirreps; h++) {
    first_offset += moinfo.orbspi[h-1]+moinfo.frdocc[h]-moinfo.frdocc[h-1];
    last_offset += moinfo.orbspi[h] - moinfo.fruocc[h] + moinfo.fruocc[h-1];
    if(moinfo.orbspi[h]) {
      moinfo.fstact[h] = first_offset;
      moinfo.lstact[h] = last_offset;
    }
  }

  /* Now define active[] such that frozen orbitals are taken into account */
  moinfo.active = init_int_array(moinfo.nirreps);
  for(h=0; h < moinfo.nirreps; h++) {
    moinfo.active[h] = moinfo.orbspi[h]-moinfo.frdocc[h]-moinfo.fruocc[h];
  }

  chkpt_close();

  /* in case IVO's have been asked for */
  if (params.ivo) {

    /* count the number of valence electrons */
    for (h=0; h < moinfo.nirreps; h++) {
      N += 2.0 * (double) (moinfo.clsdpi[h] - moinfo.frdocc[h]);
      N += (double) moinfo.openpi[h];
    }  

    /* use the number of valence electrons to compute the mixing ratio */
    params.fock_coeff = (N-1)/N;

    /*
    fprintf(outfile, "inside IVO routine now, N=%f, fock=%f\n", N,
      params.fock_coeff);
    */

    if (params.fock_coeff < 0.0) params.fock_coeff = 0.0;
    params.fzc_fock_coeff = 1.0 - params.fock_coeff;
  }

}




void get_reorder_array(void)
{
  int i, errcod;
  int j, k, l, fzv_offset;

  /* the following will only be set nonzero if it is a CI related WFN */
  moinfo.ras_opi = init_int_matrix(MAX_RAS_SPACES,moinfo.nirreps);
  moinfo.order = init_int_array(moinfo.nmo);

  if (strcmp(params.wfn, "CI") == 0 || strcmp(params.wfn, "DETCI") == 0
       || strcmp(params.wfn, "QDPT") == 0 
       || strcmp(params.wfn, "OOCCD") == 0 
       || strcmp(params.wfn, "DETCAS") == 0) {
    
    moinfo.ras_opi = init_int_matrix(MAX_RAS_SPACES,moinfo.nirreps); 
    
    if (!ras_set2(moinfo.nirreps, moinfo.nmo, params.fzc, 
         params.del_restr_docc, moinfo.orbspi, moinfo.clsdpi, moinfo.openpi,
         moinfo.frdocc, moinfo.fruocc, moinfo.rstrdocc, moinfo.rstruocc,
         moinfo.ras_opi, moinfo.order, params.ras_type, 0)) {
      fprintf(outfile, "Error in ras_set2().  Aborting.\n");
      exit(1);
    } 
  } 
  else { /* default (CC, MP2, other) */
    reorder_qt(moinfo.clsdpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc,
	       moinfo.order, moinfo.orbspi, moinfo.nirreps);
  }
  
  /* construct an array to map the other direction, i.e., from correlated */
  /* to Pitzer order */
  moinfo.corr2pitz = init_int_array(moinfo.nmo);
  for (i=0; i<moinfo.nmo; i++) {
    j = moinfo.order[i];
    moinfo.corr2pitz[j] = i;
  }

}


void cleanup(void)
{
  free(moinfo.fzc_operator);
  free_int_matrix(moinfo.ras_opi);
}

void get_mvos(void)
{

  int h, nirreps, nvir, ntri, offset, row, col;
  int i, j, ij, iabs, jabs, icorr, jcorr, ijcorr, nocc;
  int ndv, ndoccv, *ndvph, k, colv;
  int errcod;
  double *FCvv, *FC, **evecs, *evals, **Cvv, **Cvvp, **Cnew;
  double *eig_unsrt, **scf_vector;

  FC = moinfo.fzc_operator;
  nirreps = moinfo.nirreps;

  chkpt_init(PSIO_OPEN_OLD);
  scf_vector = chkpt_rd_scf();
  eig_unsrt = chkpt_rd_evals();
  chkpt_close();

  /* read in which of the doccs are to be treated as virtuals */
  ndv = 0;

  if (ip_exist("DOCC_VIRT",0)) {
     errcod = ip_count("DOCC_VIRT",&ndv,0);
     params.docc_virt = init_int_array(ndv);
     for(k=0; k < ndv; k++) {
       params.docc_virt[k] = 0;
       errcod = ip_data("DOCC_VIRT","%d",(&(params.docc_virt[k])),1,k);
     }
  }

  /* finished reading in which doccs are to be treated as virtuals */

  /* form the vv block of FC for each irrep h */
  for (h=0,offset=0; h < nirreps; h++) {

  /* find out which doccs are to be treated as virtuals */
    ndoccv = 0;
    ndvph = init_int_array(moinfo.clsdpi[h]);
    for (k=0;k<moinfo.clsdpi[h];k++) {
      ndvph[k] = 0;
      }
    if (ip_exist("DOCC_VIRT",0)) {
      for (k=0;k<(ndv-1);k=(k+2)) {
        if ( (params.docc_virt[k]) == (h+1)) {
          ndoccv++;
          ndvph[(params.docc_virt[k+1])-1] = ndoccv;
        }
      }
    }
  /* ndoccv is the number of doccs to be treated as virts for irrep h */
  /* ndvph is a 1-d array of length "number of doccs for irrep h".    */
  /* it is a string of integers, 0 meaning "leave this docc alone",   */
  /* 1 at position j means that docc number j is the first docc-virt, */
  /* 2 at position j means that docc number j is the second docc-virt */
  /* etc.                                                             */

    nvir = moinfo.virtpi[h] + ndoccv;
    nocc = moinfo.clsdpi[h] + moinfo.openpi[h] - ndoccv;

    if (nvir < 1) break;

    if (params.print_lvl) {
      fprintf(outfile, "Working on irrep %d (%d vir, %d occ)\n", 
              h, nvir, nocc);
    }

    ntri = (nvir * (nvir+1))/2;
    FCvv = init_array(ntri);
    evals = init_array(ntri);
    evecs = block_matrix(nvir,nvir);
    Cvv = block_matrix(moinfo.sopi[h],nvir);
    Cvvp = block_matrix(moinfo.sopi[h],nvir);
    Cnew = block_matrix(moinfo.sopi[h],moinfo.orbspi[h]);

    for (i=0,ij=0; i<nvir; i++) {
      /* convert this MO index into a Correlated MO index */
      if (i>=ndoccv) {
        iabs = i + nocc + offset;
      }
      else {
        for (k=0; k<(nocc + ndoccv); k++) {
          if ( (i+1) == ndvph[k] )
             iabs = k + offset;
        }
      }
      icorr = moinfo.order[iabs];
      for (j=0; j<=i; j++,ij++) {
        if (j>=ndoccv) {
          jabs = j + nocc + offset;
        }
        else {
          for (k=0; k<(nocc + ndoccv); k++) {
            if ( (j+1) == ndvph[k] )
               jabs = k + offset; 
          }
        }
        jcorr = moinfo.order[jabs];  
        ijcorr = INDEX(icorr,jcorr);
        FCvv[ij] = params.fzc_fock_coeff * FC[ijcorr] ;
        if (i==j) FCvv[ij] += params.fock_coeff * eig_unsrt[iabs];
      }
    }

    if (params.print_lvl) {
      fprintf(outfile, "Old FCvv matrix for irrep %d\n", h);
      print_array(FCvv, nvir, outfile);
    }

    /* diagonalize the vv block of FC */
    rsp(nvir, nvir, ntri, FCvv, evals, 1, evecs, 1E-12);

    if (params.print_lvl) {
      fprintf(outfile, "\nDiagonalization matrix for FCvv irrep %d\n", h);
      print_mat(evecs, nvir, nvir, outfile);

      fprintf(outfile, "\nEigenvalues of FCvv for irrep %d\n", h);
      for (i=0; i<nvir; i++) {
        fprintf(outfile, "%lf ", evals[i]);
      }
      fprintf(outfile, "\n");
    }

    /* load up the vv part of SCF matrix C */
    for (i=moinfo.first_so[h],row=0; i<= moinfo.last_so[h]; i++,row++) {
      colv = 0;
      for (j=moinfo.first_so[h],k=0; j < (moinfo.first_so[h]+nocc+ndoccv); j++,k++) {
        if (ndvph[k] != 0) {
           Cvv[row][colv] = scf_vector[i][j];
           colv += 1;
        }
      }
      for (j=moinfo.first_so[h]+nocc+ndoccv,col=ndoccv; j<=moinfo.last_so[h]; j++,col++) {
        Cvv[row][col] = scf_vector[i][j];
      }
    } 

    if (params.print_lvl) {
      fprintf(outfile, "Old Cvv matrix for irrep %d\n", h);
      print_mat(Cvv, moinfo.sopi[h], nvir, outfile);
    }
     
    /* multiply the FCvv eigenvectors by Cvv to get new Cvv */
    mmult(Cvv, 0, evecs, 0, Cvvp, 0, moinfo.sopi[h], nvir, nvir, 0);
 
    if (params.print_lvl) {
      fprintf(outfile, "New Cvv matrix for irrep %d\n", h);
      print_mat(Cvvp, moinfo.sopi[h], nvir, outfile);
    }

    /* Load up the new C matrix and write it out */

    /* load up the occupied orbitals */
    for (i=moinfo.first_so[h],row=0; i<= moinfo.last_so[h]; i++,row++) {
      for (j=moinfo.first_so[h],col=0; j< moinfo.first_so[h]+nocc+ndoccv; j++,col++) {
        if (ndvph[col] != 0) {
          Cnew[row][col] = Cvvp[row][(ndvph[col]-1)];
        }
        else {
          Cnew[row][col] = scf_vector[i][j];
        }
      }
    }

    /* load up the new virtuals */
    for (i=0; i<moinfo.sopi[h]; i++) {
      for (j=ndoccv; j<nvir; j++) {
        Cnew[i][j+nocc] = Cvvp[i][j];
      }
    }

    if (params.print_lvl) {
      fprintf(outfile, "New C matrix for irrep %d\n", h);
      print_mat(Cnew, moinfo.sopi[h], moinfo.orbspi[h], outfile);
    }

    chkpt_init(PSIO_OPEN_OLD);
    chkpt_wt_scf_irrep(Cnew, h);
    chkpt_close();

    offset += moinfo.orbspi[h];
    free(FCvv);
    free(evals);
    free_block(evecs);
  }

  free(eig_unsrt);
  free_block(scf_vector);
}

}} // end namespace psi::mvo

