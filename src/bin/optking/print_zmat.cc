/*! \file
    \ingroup OPTKING
    \brief PRINT_ZMAT() printout z-matrix
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <physconst.h>
#include <libchkpt/chkpt.h>

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"

#include <masses.h>

namespace psi { namespace optking {

/* recompute values for z-matrix and write back to chkpt */
void compute_zmat(cartesians &carts, int *unique_zvars) {
  int i, a, cnt;
  int nallatom, natom, *to_nodummy, *atom_dummy;
  struct internals *simples;
  struct z_entry *zmat;
  char sym[30];
  double *opt_fgeom;

  nallatom = optinfo.nallatom;
  chkpt_init(PSIO_OPEN_OLD);
  zmat = chkpt_rd_zmat();
  chkpt_close();

/*
  for (i=0;i<nallatom;++i) {
    fprintf(outfile,"%d %d %d %d\n",  i, zmat[i].bond_atom, zmat[i].angle_atom, zmat[i].tors_atom);
    fprintf(outfile,"%20.10lf %20.10lf %20.10lf\n", zmat[i].bond_val, zmat[i].angle_val, zmat[i].tors_val);
    fprintf(outfile,"%s %s %s \n", zmat[i].bond_label, zmat[i].angle_label, zmat[i].tors_label);
   }
*/

  /* determine and save the unique variables */
  cnt = -1;
  for (i=1; i<nallatom; ++i) {
    ++cnt;
    if (zmat[i].bond_label[0] != '\0') {
      unique_zvars[cnt] = 1;
      strcpy(sym, zmat[i].bond_label);
      for (a=0;a<i;++a) {
        if (strcmp(sym, zmat[a].bond_label) == 0) {
          unique_zvars[cnt] = 0;
          break;
        }
      }
    }
    else unique_zvars[cnt] = 0;
    if (i>1) {
      ++cnt;
      if (zmat[i].angle_label[0] != '\0') {
        unique_zvars[cnt] = 1;
        strcpy(sym, zmat[i].angle_label);
        for (a=0;a<i;++a) {
          if (strcmp(sym, zmat[a].angle_label) == 0) {
            unique_zvars[cnt] = 0;
            break;
          }
        }
      }
      else unique_zvars[cnt] = 0;
    }
    if (i>2) {
      ++cnt;
      if (zmat[i].tors_label[0] != '\0') {
        unique_zvars[cnt] = 1;
        strcpy(sym, zmat[i].tors_label);
        for (a=0;a<i;++a) {
          if (strcmp(sym, zmat[a].tors_label) == 0) {
            unique_zvars[cnt] = 0;
            break;
          }
        }
      }
      else unique_zvars[cnt] = 0;
    }
  }

  int *nints;
  nints = (int *) malloc(6*sizeof(int));
  nints[0] = nallatom-1; /* stre */
  nints[1] = nallatom-2; /* bend */
  nints[2] = nallatom-3; /* tors */
  nints[3] = 0; /* oop */
  nints[4] = 0; /* linb */
  nints[5] = 0; /* fragment */
  for (i=0; i<6; ++i)
    if (nints[i] < 0) nints[i] = 0;
  internals zints(nints);
  /* compute the value of the unique variables */
  zints.stre.set_num(nints[0]);
  zints.bend.set_num(nints[1]);
  zints.tors.set_num(nints[2]);
  zints.out.set_num(0);
  zints.linb.set_num(0);
  zints.frag.set_num(0);

  cnt = 0;
  for (i=0;i<nallatom;++i) {
    if (i>0) {
      zints.stre.set_id(i-1,++cnt);
      zints.stre.set_A(i-1,i);
      zints.stre.set_B(i-1,zmat[i].bond_atom-1);
    }
    if (i>1) {
      zints.bend.set_id(i-2,++cnt);
      zints.bend.set_A(i-2,i);
      zints.bend.set_B(i-2,zmat[i].bond_atom-1);
      zints.bend.set_C(i-2,zmat[i].angle_atom-1);
    }
    if (i>2) {
      zints.tors.set_id(i-3,++cnt);
      zints.tors.set_A(i-3,i);
      zints.tors.set_B(i-3,zmat[i].bond_atom-1);
      zints.tors.set_C(i-3,zmat[i].angle_atom-1);
      zints.tors.set_D(i-3,zmat[i].tors_atom-1);
    }
  }
  // compute value of the zmatrix coordinates
  opt_fgeom = carts.get_fcoord();
  zints.compute_internals(nallatom, opt_fgeom);

  // insert computed values into zmatrix object
  for (i=0;i<nallatom;++i) {
    if (i>0) {
      zmat[i].bond_val = zints.stre.get_val(i-1);
    }
    if (i>1) {
      zmat[i].angle_val = zints.bend.get_val(i-2);
    }
    if (i>2) {
      zmat[i].tors_val = zints.tors.get_val(i-3);
    }
  }

  // write recomputed z-matrix to PSIF_CHKPT
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_zmat(zmat);
  chkpt_close();
}

/*** ZMAT_TO_INTCO if simples are not already there and zmat_simples
* is not turned off, then generate simple internals from zmatrix ***/
void print_zmat(FILE *outfile, int *unique_zvars) {
  int i, a, b, c, d, cnt = 0;
  int nallatom, natom, *to_nodummy, *atom_dummy;
  char **felement;
  char buf[2], *sym;
  double *zvals;
  struct z_entry *zmat;
  const char *X = "X";

  nallatom = optinfo.nallatom;
  natom = optinfo.natom;
  atom_dummy = optinfo.atom_dummy;
  to_nodummy = optinfo.to_nodummy;

  chkpt_init(PSIO_OPEN_OLD);
  zmat = chkpt_rd_zmat();
  felement = chkpt_rd_felement();
  zvals = chkpt_rd_zvals();
  chkpt_close();

  fprintf(outfile,"  zmat = ( \n");
  for (i=0; i<nallatom; ++i) {
    if (atom_dummy[i]) sym = X;
    else sym = atomic_labels[to_nodummy[i]];
    fprintf(outfile,"    ( %s ", sym);
    if (i > 0) {
      fprintf(outfile," %d", zmat[i].bond_atom);
      if (zmat[i].bond_label[0] != '\0') {
        fprintf(outfile," %s", zmat[i].bond_label);
        if (zmat[i].bond_opt) fprintf(outfile,"$");
      }
      else fprintf(outfile," %10.5lf", zmat[i].bond_val);
    }
    if (i > 1) {
      fprintf(outfile," %d", zmat[i].angle_atom);
      if (zmat[i].angle_label[0] != '\0') {
        fprintf(outfile," %s", zmat[i].angle_label);
        if (zmat[i].angle_opt) fprintf(outfile,"$");
      }
      else fprintf(outfile," %10.5lf", zmat[i].angle_val);
    }
    if (i > 2) {
      fprintf(outfile," %d", zmat[i].tors_atom);
      if (zmat[i].tors_label[0] != '\0') {
        fprintf(outfile," %s", zmat[i].tors_label);
        if (zmat[i].tors_opt) fprintf(outfile,"$");
      }
      else fprintf(outfile," %10.5lf", zmat[i].tors_val);
    }
    fprintf(outfile,")\n");
  }
  fprintf(outfile,"  )\n");
  int cnt_vars = 0;
  for (i=0;i<MAX_ZVARS;++i)
    if (unique_zvars[i]) ++cnt_vars;

  if (cnt_vars != 0) {
    cnt = -1;
    fprintf(outfile,"  zvars = ( \n");
    for (i=0; i<nallatom; ++i) {
      if (i > 0) {
        if (unique_zvars[++cnt])
          fprintf(outfile,"    ( %s %10.5lf )\n",
              zmat[i].bond_label, zmat[i].bond_val);
      }
      if (i > 1) {
        if (unique_zvars[++cnt])
          fprintf(outfile,"    ( %s %10.5lf )\n",
              zmat[i].angle_label, zmat[i].angle_val);
      }
      if (i > 2) {
        if (unique_zvars[++cnt])
          fprintf(outfile,"    ( %s %10.5lf )\n",
              zmat[i].tors_label, zmat[i].tors_val);
      }
    }
    fprintf(outfile,"  )\n");
  }
  return;
}

}} /* namespace psi::optking */

