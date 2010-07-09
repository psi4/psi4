/*! \file
    \ingroup OPTKING
    \brief GET_SYMINFO: gets symmetry info from chkpt and char_table.c
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

#include <libipv1/ip_lib.h>
#include <libchkpt/chkpt.h>

namespace psi { //namespace optking {

void get_syminfo(const simples_class & simples) {
  int a, b, c, d, aa, bb, cc, dd, i, j, sign, natom;
  int id, sub_index, sub_index2, ops,linval;
  int nallatom;
  Intco_type intco_type;

  nallatom = optinfo.nallatom;
  natom = optinfo.natom;

  chkpt_init(PSIO_OPEN_OLD);
  if (optinfo.mode == MODE_DISP_LOAD) {
    fprintf(outfile,"Reading symmetry information from root area of checkpoint.\n");
    chkpt_reset_prefix();
    chkpt_commit_prefix();
  }

  if (ip_exist("SYMMETRY",0))
    ip_string("SYMMETRY",&syminfo.symmetry,0);
  else {
    syminfo.symmetry = chkpt_rd_sym_label();
  }
  syminfo.nirreps = chkpt_rd_nirreps();
  syminfo.ict = chkpt_rd_ict();
  syminfo.irrep_lbls = chkpt_rd_irr_labs();
  chkpt_close();

  //for (j=0; j<syminfo.nirreps; ++j) {
    //for(i=0; i<natom; ++i)
      //fprintf(outfile,"%d", syminfo.ict[j][i]);
    //fprintf(outfile,"\n");
  //}

  // abort if group label is greater than 3 characters
  j = strlen(syminfo.symmetry);
  if (j > 3) abort();
  strncpy(ptgrp,syminfo.symmetry,4);
  for ( ;j<3;++j)
    ptgrp[j] = ' ';
  ptgrp[3] = '\0';
  for (i=0;i<3;++i) 
    ptgrp[i] = toupper(ptgrp[i]);

  syminfo.clean_irrep_lbls = (char **) malloc(syminfo.nirreps * sizeof(char *));
  // new char*[nirreps];
  for (i=0;i<syminfo.nirreps;++i) {
    //  syminfo.clean_irrep_lbls[i] = new char[4];
    syminfo.clean_irrep_lbls[i] = (char *) malloc(4*sizeof(char));
    // new char[4];
    strcpy(syminfo.clean_irrep_lbls[i], syminfo.irrep_lbls[i]);
    for (j=0;j<3;++j) {
      if (syminfo.clean_irrep_lbls[i][j] == '\"')
        syminfo.clean_irrep_lbls[i][j] = 'P' ;
      if (syminfo.clean_irrep_lbls[i][j] == '\'')
        syminfo.clean_irrep_lbls[i][j] = 'p' ;
    }
  }
  syminfo.ct = get_char_table(ptgrp);
  syminfo.op_lbls = get_symm_ops(ptgrp);
  syminfo.ict_ops = init_int_matrix(simples.get_num(),syminfo.nirreps);
  syminfo.ict_ops_sign = init_int_matrix(simples.get_num(),syminfo.nirreps);
  ops_in_class = get_ops_in_class(ptgrp, syminfo.nirreps);

  // make dummy atoms transform into themselves
  syminfo.fict = init_int_matrix(syminfo.nirreps, optinfo.nallatom);
  for (i=0; i<syminfo.nirreps; ++i) 
    for (j=0; j<optinfo.nallatom; ++j) 
      syminfo.fict[i][j] = j;

  for (i=0; i<syminfo.nirreps; ++i) 
    for (j=0; j<natom; ++j)  {
      a = optinfo.to_dummy[j];
      syminfo.fict[i][a] = optinfo.to_dummy[syminfo.ict[i][j]-1];
    }
  // for consistency with ict, start atom numbering at 1
  for (i=0; i<syminfo.nirreps; ++i) 
    for (j=0; j<optinfo.natom; ++j) 
      syminfo.fict[i][j] += 1;

  for (i=0;i<simples.get_num();++i)
    for (j=0;j<syminfo.nirreps;++j)
      syminfo.ict_ops_sign[i][j] = 1;

  // Generate simple internal coordinate transformation matrix
  // it would seem that one could only use fict instead of ict and
  // include dymmy atoms, only is symmetry sets of dummy atoms were generated
  for (i=0;i<simples.get_num();++i) {
    id = simples.index_to_id(i);
    simples.locate_id(id,&intco_type,&sub_index,&sub_index2);
    if (intco_type == STRE) {
      a = simples.get_atom(STRE, sub_index, 0);
      b = simples.get_atom(STRE, sub_index, 1);
      for (ops=0;ops < syminfo.nirreps;++ops) {
        aa = syminfo.ict[ops][a]-1;
        bb = syminfo.ict[ops][b]-1;
        syminfo.ict_ops[i][ops] = simples.get_id_from_atoms_stre(aa,bb);
      }
    }
    if (intco_type == BEND) {
      a = simples.get_atom(BEND, sub_index, 0);
      b = simples.get_atom(BEND, sub_index, 1);
      c = simples.get_atom(BEND, sub_index, 2);
      for (ops=0;ops < syminfo.nirreps;++ops) {
        aa = syminfo.ict[ops][a]-1;
        bb = syminfo.ict[ops][b]-1;
        cc = syminfo.ict[ops][c]-1;
        syminfo.ict_ops[i][ops] = simples.get_id_from_atoms_bend(aa,bb,cc);
      }
    }
    if (intco_type == TORS) {
      a = simples.get_atom(TORS, sub_index, 0);
      b = simples.get_atom(TORS, sub_index, 1);
      c = simples.get_atom(TORS, sub_index, 2);
      d = simples.get_atom(TORS, sub_index, 3);
//fprintf(outfile,"tors %d %d %d %d\n", a, b, c, d);
      for (ops=0;ops < syminfo.nirreps;++ops) {
        aa = syminfo.ict[ops][a]-1;
        bb = syminfo.ict[ops][b]-1;
        cc = syminfo.ict[ops][c]-1;
        dd = syminfo.ict[ops][d]-1;
//fprintf(outfile,"transformed tors %d %d %d %d\n", aa, bb, cc, dd);
        syminfo.ict_ops[i][ops] = simples.get_id_from_atoms_tors(aa,bb,cc,dd);
        if ( ('S' == syminfo.op_lbls[ops][0]) || ('I' == syminfo.op_lbls[ops][0]) )
          syminfo.ict_ops_sign[i][ops] = -1;
      }
    }
    if (intco_type == OUT) {
      a = simples.get_atom(OUT, sub_index, 0);
      b = simples.get_atom(OUT, sub_index, 1);
      c = simples.get_atom(OUT, sub_index, 2);
      d = simples.get_atom(OUT, sub_index, 3);
      for (ops=0;ops < syminfo.nirreps;++ops) {
        aa = syminfo.ict[ops][a]-1;
        bb = syminfo.ict[ops][b]-1;
        cc = syminfo.ict[ops][c]-1;
        dd = syminfo.ict[ops][d]-1;
        syminfo.ict_ops[i][ops] = simples.get_id_from_atoms_out(aa,bb,cc,dd,&sign);
        if ( ('S' == syminfo.op_lbls[ops][0]) || ('I' == syminfo.op_lbls[ops][0]) )
          sign *= -1;
        syminfo.ict_ops_sign[i][ops] = sign;
      }
    }
    // this probably doesn't work anyway but here it is
    if (intco_type == LINB) {
      a = simples.get_atom(LINB, sub_index, 0);
      b = simples.get_atom(LINB, sub_index, 1);
      c = simples.get_atom(LINB, sub_index, 2);
      linval = simples.linb_get_linval(sub_index);
      for (ops=0;ops < syminfo.nirreps;++ops) {
        aa = syminfo.ict[ops][a]-1;
        bb = syminfo.ict[ops][b]-1;
        cc = syminfo.ict[ops][c]-1;
        syminfo.ict_ops[i][ops] = simples.get_id_from_atoms_linb(aa,bb,cc,linval);
      }
    }
  }


  if (optinfo.print_symmetry) {
    fprintf(outfile,"\n+++ Symmetry Information +++\n");
    fprintf(outfile,"The ICT table from chkpt:\n");
    for(i=0;i<syminfo.nirreps;++i) {
      for(j=0;j<natom;++j)
        fprintf(outfile,"%3d",syminfo.ict[i][j]);
      fprintf(outfile,"\n");
    }

    fprintf(outfile,"The FICT table from chkpt:\n");
    for(i=0;i<syminfo.nirreps;++i) {
      for(j=0;j<nallatom;++j)
        fprintf(outfile,"%3d",syminfo.fict[i][j]);
      fprintf(outfile,"\n");
    }

    fprintf(outfile,"Clean irrep labels:");
    for (j=0;j<syminfo.nirreps;++j)
      fprintf(outfile,"%5s ", syminfo.clean_irrep_lbls[j]);
    fprintf(outfile,"\n");

    fprintf(outfile,"\nCharacter table from char_table.c and symmetry: %s\n", ptgrp);
    fprintf(outfile,"      ");
    for (i=0;i<syminfo.nirreps;++i) fprintf(outfile,"%5s",syminfo.op_lbls[i]);
    fprintf(outfile,"\n");

    for (i=0;i<syminfo.nirreps;++i) {
      for (j=0;j<syminfo.nirreps;++j) {
        if (j == 0) fprintf(outfile,"%5s ", syminfo.irrep_lbls[i]);
        fprintf(outfile,"%5d",syminfo.ct[i][j]);
      }
      fprintf(outfile,"\n");
    }
    fprintf(outfile,"Internal coordinate transformation matrix.\n");
    for (i=0;i<syminfo.nirreps;++i) fprintf(outfile,"%5s",syminfo.op_lbls[i]);
    fprintf(outfile,"\n");
    for (id=0;id<simples.get_num();++id) {
      for (ops=0;ops < syminfo.nirreps;++ops)
        fprintf(outfile,"%5d",syminfo.ict_ops[id][ops]);
      fprintf(outfile,"\n");
    }
    fprintf(outfile,
        "Internal transformation sign matrix to fix torsions and out-of-planes.\n");
    for (i=0;i<syminfo.nirreps;++i) fprintf(outfile,"%5s",syminfo.op_lbls[i]);
    fprintf(outfile,"\n");
    for (id=0;id<simples.get_num();++id) {
      for (ops=0;ops < syminfo.nirreps;++ops)
        fprintf(outfile,"%5d",syminfo.ict_ops_sign[id][ops]);
      fprintf(outfile,"\n");
    }
  }
  return;
}

}//} /* namespace psi::optking */

