/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here 
*/

/* #define DEBUG */

#define EXTERN

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libmints/mints.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libchkpt/chkpt.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <physconst.h>
#include "structs.h"
#include "globals.h"
#include "civect.h"

namespace psi { namespace detci {

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define TOL 1E-14

void orbsfile_rd_blk(int targetfile, int root, int irrep, double **orbs_vector);
void orbsfile_wt_blk(int targetfile, int root, int irrep, double **orbs_vector);
void ave(int targetfile);
void opdm_block(struct stringwr **alplist, struct stringwr **betlist,
		double **onepdm_a, double **onepdm_b, double **CJ, double **CI, int Ja_list, 
		int Jb_list, int Jnas, int Jnbs, int Ia_list, int Ib_list, 
		int Inas, int Inbs);
void opdm_ke(double **onepdm);
// void get_mo_dipmom_ints(double **mux_mo, double **muy_mo, double **muz_mo);
// void get_dipmom_nuc(double *mu_x_n, double *mu_y_n, double *mu_z_n);


/*
** Computes the one-particle density matrix for all n roots.  If
** transdens is set, then will compute the transition densities
** from Iroot to Jroot where Jroot will run over all roots
*/ 
void opdm(struct stringwr **alplist, struct stringwr **betlist, 
          int transdens, int dipmom,
          int Inroots, int Iroot, int Inunits, int Ifirstunit, 
	  int Jnroots, int Jroot, int Jnunits, int Jfirstunit, 
	  int targetfile, int writeflag, int printflag)
{

  CIvect Ivec, Jvec;
  int i, j, k, l, klast, roots;
  int maxrows, maxcols;
  unsigned long bufsz;
  int max_orb_per_irrep;
  double **transp_tmp = NULL;
  double **transp_tmp2 = NULL;
  double *buffer1, *buffer2, **onepdm, **onepdm_a, **onepdm_b;
  int i_ci, j_ci, irrep, mo_offset, so_offset, orb_length=0, opdm_length=0;
  double *opdm_eigval, **opdm_eigvec, **opdm_blk, **scfvec;
  int Iblock, Iblock2, Ibuf, Iac, Ibc, Inas, Inbs, Iairr;
  int Jblock, Jblock2, Jbuf, Jac, Jbc, Jnas, Jnbs, Jairr;
  int do_Jblock, do_Jblock2;
  int populated_orbs;
  double **tmp_mat, **opdmso;
  double overlap, max_overlap;
  char opdm_key[80]; /* libpsio TOC entry name for OPDM for each root */
  // double **mux_mo, **muy_mo, **muz_mo;
  // double mu_x, mu_y, mu_z, mu_tot;
  double mux_n, muy_n, muz_n; /* nuclear parts of dipole moments */

  if (!transdens) Iroot = 0;
  if (transdens) 
    Jroot = 1;
  else
    Jroot = 0; 
  if (Jroot > Jnroots) return;
  
  Ivec.set(CIblks.vectlen, CIblks.num_blocks, Parameters.icore, Parameters.Ms0,
           CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
           CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
           CalcInfo.nirreps, AlphaG->subgr_per_irrep, Inroots, Inunits,
           Ifirstunit, CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);

  Jvec.set(CIblks.vectlen, CIblks.num_blocks, Parameters.icore, Parameters.Ms0,
           CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
           CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
           CalcInfo.nirreps, AlphaG->subgr_per_irrep, Jnroots, Jnunits,
           Jfirstunit, CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);

  populated_orbs = CalcInfo.num_ci_orbs + CalcInfo.num_fzc_orbs;
  for (irrep=0; irrep<CalcInfo.nirreps; irrep++) {
     opdm_length += (CalcInfo.orbs_per_irr[irrep] - CalcInfo.frozen_uocc[irrep])
                 * (CalcInfo.orbs_per_irr[irrep] - CalcInfo.frozen_uocc[irrep]);
     orb_length += (CalcInfo.so_per_irr[irrep]*CalcInfo.orbs_per_irr[irrep]);
     }

  /* find biggest blocksize */
  for (irrep=0,max_orb_per_irrep=0; irrep<CalcInfo.nirreps; irrep++) {
    if (CalcInfo.orbs_per_irr[irrep] > max_orb_per_irrep)
      max_orb_per_irrep = CalcInfo.so_per_irr[irrep];
  }

  opdm_eigvec = block_matrix(max_orb_per_irrep, max_orb_per_irrep);
  opdm_eigval = init_array(max_orb_per_irrep);
  opdm_blk = block_matrix(max_orb_per_irrep, max_orb_per_irrep);
  opdmso = block_matrix(CalcInfo.nso, CalcInfo.nso);
  tmp_mat = block_matrix(max_orb_per_irrep, max_orb_per_irrep);

  /* this index stuff is probably irrelevant now ... CDS 6/03 */
  /*
  Parameters.opdm_idxmat =
    init_int_matrix(Parameters.num_roots+2, CalcInfo.nirreps);
  Parameters.orbs_idxmat = 
    init_int_matrix(Parameters.num_roots+2, CalcInfo.nirreps);
  for (l=0; l<=(Parameters.num_roots+1); l++) {
     Parameters.opdm_idxmat[l][0] = l * opdm_length * sizeof(double); 
     Parameters.orbs_idxmat[l][0] = l * orb_length * sizeof(double);
     for (irrep=1; irrep<CalcInfo.nirreps; irrep++) {
        Parameters.orbs_idxmat[l][irrep] =
          Parameters.orbs_idxmat[l][irrep-1]+
          CalcInfo.so_per_irr[irrep-1]*CalcInfo.orbs_per_irr[irrep-1]
          *sizeof(double);  
        Parameters.opdm_idxmat[l][irrep] =
          Parameters.opdm_idxmat[l][irrep-1] +
          (CalcInfo.orbs_per_irr[irrep-1]-CalcInfo.frozen_uocc[irrep]) *
          (CalcInfo.orbs_per_irr[irrep-1]-CalcInfo.frozen_uocc[irrep]) *
          sizeof(double);
        } 
     } 
  */

  buffer1 = Ivec.buf_malloc();
  buffer2 = Jvec.buf_malloc();
  Ivec.buf_lock(buffer1);
  Jvec.buf_lock(buffer2);
  onepdm = block_matrix(populated_orbs, populated_orbs); 
  onepdm_a = block_matrix(populated_orbs, populated_orbs); 
  onepdm_b = block_matrix(populated_orbs, populated_orbs); 

  if ((Ivec.icore==2 && Ivec.Ms0 && CalcInfo.ref_sym != 0) || 
      (Ivec.icore==0 && Ivec.Ms0)) {
    for (i=0, maxrows=0, maxcols=0; i<Ivec.num_blocks; i++) {
      if (Ivec.Ia_size[i] > maxrows) maxrows = Ivec.Ia_size[i];
      if (Ivec.Ib_size[i] > maxcols) maxcols = Ivec.Ib_size[i];
    }
    if (maxcols > maxrows) maxrows = maxcols;
    transp_tmp = (double **) malloc (maxrows * sizeof(double *));
    transp_tmp2 = (double **) malloc (maxrows * sizeof(double *));
    if (transp_tmp == NULL || transp_tmp2 == NULL) {
      printf("(opdm): Trouble with malloc'ing transp_tmp\n");
    }
    bufsz = Ivec.get_max_blk_size();
    transp_tmp[0] = init_array(bufsz);
    transp_tmp2[0] = init_array(bufsz);
    if (transp_tmp[0] == NULL || transp_tmp2[0] == NULL) {
      printf("(opdm): Trouble with malloc'ing transp_tmp[0]\n");
    }
  }
 
  if (writeflag) 
    psio_open(targetfile, PSIO_OPEN_OLD);

  /* if we're trying to follow a root, figure out which one here */
  
  if (Parameters.follow_vec_num > 0 && !transdens) {
    max_overlap = 0.0;
    for (i=0,j=0; i<Parameters.num_roots; i++) {

      overlap = Ivec.compute_follow_overlap(i, 
        Parameters.follow_vec_num, Parameters.follow_vec_coef,
        Parameters.follow_vec_Iac, Parameters.follow_vec_Iaridx, 
        Parameters.follow_vec_Ibc, Parameters.follow_vec_Ibridx);

      if (overlap > max_overlap) {
        max_overlap = overlap;
        j = i;
      } 
    }


    Parameters.root = j;
    fprintf(outfile, "DETCI following root %d (overlap %6.4lf)\n", 
      Parameters.root+1,max_overlap);

    for (i=0; i<Parameters.average_num; i++) {
      Parameters.average_states[i] = 0;
      Parameters.average_weights[i] = 0.0;
    }
    Parameters.average_num = 1;
    Parameters.average_states[0] = Parameters.root;
    Parameters.average_weights[0] = 1.0;
    
    /* correct "the" energy in the checkpoint file */
    chkpt_init(PSIO_OPEN_OLD);
    sprintf(opdm_key,"Root %2d energy",Parameters.root);
    overlap = chkpt_rd_e_labeled(opdm_key);
    chkpt_wt_etot(overlap);
    Process::environment.globals["CURRENT ENERGY"] = overlap;
    Process::environment.globals["CI TOTAL ENERGY"] = overlap;
    // eref is wrong for open-shells so replace it with escf until
    // I fix it, CDS 11/5/11
    Process::environment.globals["CI CORRELATION ENERGY"] = overlap - 
      CalcInfo.escf;
    chkpt_close();
  
  }
    
  /* if getting (transition) moments, need to read in the AO dip mom ints */
  /*
  if (dipmom) {
    mux_mo = block_matrix(CalcInfo.nmo, CalcInfo.nmo);
    muy_mo = block_matrix(CalcInfo.nmo, CalcInfo.nmo);
    muz_mo = block_matrix(CalcInfo.nmo, CalcInfo.nmo);
    get_mo_dipmom_ints(mux_mo, muy_mo, muz_mo);
    get_dipmom_nuc(&mux_n, &muy_n, &muz_n);
  } 
  */

  for (; Jroot<Parameters.num_roots; Jroot++) {
   
    zero_mat(onepdm_a, populated_orbs, populated_orbs); 
    zero_mat(onepdm_b, populated_orbs, populated_orbs); 

    if (!transdens) {
      for (i=0; i<CalcInfo.num_fzc_orbs; i++) {
        onepdm_a[i][i] = 1.0;
        onepdm_b[i][i] = 1.0;
      }
    }

    if (Parameters.icore == 0) {
 
      for (Ibuf=0; Ibuf<Ivec.buf_per_vect; Ibuf++) {
        Ivec.read(Iroot, Ibuf);
        Iblock = Ivec.buf2blk[Ibuf];
        Iac = Ivec.Ia_code[Iblock];
        Ibc = Ivec.Ib_code[Iblock];
        Inas = Ivec.Ia_size[Iblock];
        Inbs = Ivec.Ib_size[Iblock];
       
        for (Jbuf=0; Jbuf<Jvec.buf_per_vect; Jbuf++) {
          do_Jblock=0; do_Jblock2=0;
          Jblock = Jvec.buf2blk[Jbuf];
          Jblock2 = -1;
          Jac = Jvec.Ia_code[Jblock];
          Jbc = Jvec.Ib_code[Jblock];
          if (Jvec.Ms0) Jblock2 = Jvec.decode[Jbc][Jac];
          Jnas = Jvec.Ia_size[Jblock];
          Jnbs = Jvec.Ib_size[Jblock];
          if (s1_contrib[Iblock][Jblock] || s2_contrib[Iblock][Jblock]) 
            do_Jblock = 1;
          if (Jvec.buf_offdiag[Jbuf] && (s1_contrib[Iblock][Jblock2] ||
                                         s2_contrib[Iblock][Jblock2]))
            do_Jblock2 = 1;
          if (!do_Jblock && !do_Jblock2) continue;
	 
          Jvec.read(Jroot, Jbuf);
	 
          if (do_Jblock) {
            opdm_block(alplist, betlist, onepdm_a, onepdm_b, Jvec.blocks[Jblock], 
                       Ivec.blocks[Iblock], Jac, Jbc, Jnas,
                       Jnbs, Iac, Ibc, Inas, Inbs);
            }
	 
          if (do_Jblock2) {
            Jvec.transp_block(Jblock, transp_tmp);
            opdm_block(alplist, betlist, onepdm_a, onepdm_b, transp_tmp,
                       Ivec.blocks[Iblock], Jbc, Jac, Jnbs,
                       Jnas, Iac, Ibc, Inas, Inbs);
          }
	 
        } /* end loop over Jbuf */
       
        if (Ivec.buf_offdiag[Ibuf]) { /* need to get contrib of transpose */
          Iblock2 = Ivec.decode[Ibc][Iac];
          Iac = Ivec.Ia_code[Iblock2];
          Ibc = Ivec.Ib_code[Iblock2];
          Inas = Ivec.Ia_size[Iblock2];
          Inbs = Ivec.Ib_size[Iblock2];
       
          Ivec.transp_block(Iblock, transp_tmp2);

          for (Jbuf=0; Jbuf<Jvec.buf_per_vect; Jbuf++) {
            do_Jblock=0; do_Jblock2=0;
            Jblock = Jvec.buf2blk[Jbuf];
            Jblock2 = -1;
            Jac = Jvec.Ia_code[Jblock];
            Jbc = Jvec.Ib_code[Jblock];
            if (Jvec.Ms0) Jblock2 = Jvec.decode[Jbc][Jac];
            Jnas = Jvec.Ia_size[Jblock];
            Jnbs = Jvec.Ib_size[Jblock];
            if (s1_contrib[Iblock2][Jblock] || s2_contrib[Iblock2][Jblock]) 
              do_Jblock = 1;
            if (Jvec.buf_offdiag[Jbuf] && (s1_contrib[Iblock2][Jblock2] ||
                                           s2_contrib[Iblock2][Jblock2]))
              do_Jblock2 = 1;
            if (!do_Jblock && !do_Jblock2) continue;
	   
            Jvec.read(Jroot, Jbuf);
	 
            if (do_Jblock) {
              opdm_block(alplist, betlist, onepdm_a, onepdm_b, Jvec.blocks[Jblock], 
                         transp_tmp2, Jac, Jbc, Jnas,
                         Jnbs, Iac, Ibc, Inas, Inbs);
            }
	   
            if (do_Jblock2) {
              Jvec.transp_block(Jblock, transp_tmp);
              opdm_block(alplist, betlist, onepdm_a, onepdm_b, transp_tmp,
                         transp_tmp2, Jbc, Jac, Jnbs,
                         Jnas, Iac, Ibc, Inas, Inbs);
            }
          } /* end loop over Jbuf */
        } /* end loop over Ibuf transpose */
      } /* end loop over Ibuf */
    } /* end icore==0 */

    else if (Parameters.icore==1) { /* whole vectors in-core */
      Ivec.read(Iroot, 0);
      Jvec.read(Jroot, 0);
      for (Iblock=0; Iblock<Ivec.num_blocks; Iblock++) {
        Iac = Ivec.Ia_code[Iblock];
        Ibc = Ivec.Ib_code[Iblock];
        Inas = Ivec.Ia_size[Iblock];
        Inbs = Ivec.Ib_size[Iblock];
        if (Inas==0 || Inbs==0) continue;
        for (Jblock=0; Jblock<Jvec.num_blocks; Jblock++) {
          Jac = Jvec.Ia_code[Jblock];
          Jbc = Jvec.Ib_code[Jblock];
          Jnas = Jvec.Ia_size[Jblock];
          Jnbs = Jvec.Ib_size[Jblock];
          if (s1_contrib[Iblock][Jblock] || s2_contrib[Iblock][Jblock])
            opdm_block(alplist, betlist, onepdm_a, onepdm_b, Jvec.blocks[Jblock],
                       Ivec.blocks[Iblock], Jac, Jbc, Jnas,
                       Jnbs, Iac, Ibc, Inas, Inbs);
        }
      } /* end loop over Iblock */
    } /* end icore==1 */

    else if (Parameters.icore==2) { /* icore==2 */
      for (Ibuf=0; Ibuf<Ivec.buf_per_vect; Ibuf++) {
        Ivec.read(Iroot, Ibuf);
        Iairr = Ivec.buf2blk[Ibuf];

        for (Jbuf=0; Jbuf<Jvec.buf_per_vect; Jbuf++) {
          Jvec.read(Jroot, Jbuf);
          Jairr = Jvec.buf2blk[Jbuf];
	
        for (Iblock=Ivec.first_ablk[Iairr]; Iblock<=Ivec.last_ablk[Iairr];
             Iblock++) {
          Iac = Ivec.Ia_code[Iblock];
          Ibc = Ivec.Ib_code[Iblock];
          Inas = Ivec.Ia_size[Iblock];
          Inbs = Ivec.Ib_size[Iblock];
	   
          for (Jblock=Jvec.first_ablk[Jairr]; Jblock<=Jvec.last_ablk[Jairr];
               Jblock++) {
            Jac = Jvec.Ia_code[Jblock];
            Jbc = Jvec.Ib_code[Jblock];
            Jnas = Jvec.Ia_size[Jblock];
            Jnbs = Jvec.Ib_size[Jblock];
	   
            if (s1_contrib[Iblock][Jblock] || s2_contrib[Iblock][Jblock])
              opdm_block(alplist, betlist, onepdm_a, onepdm_b, Jvec.blocks[Jblock],
                         Ivec.blocks[Iblock], Jac, Jbc, Jnas,
                         Jnbs, Iac, Ibc, Inas, Inbs);

            if (Jvec.buf_offdiag[Jbuf]) {
              Jblock2 = Jvec.decode[Jbc][Jac];
              if (s1_contrib[Iblock][Jblock2] ||
                  s2_contrib[Iblock][Jblock2]) {
              Jvec.transp_block(Jblock, transp_tmp);
                opdm_block(alplist, betlist, onepdm_a, onepdm_b, transp_tmp,
                  Ivec.blocks[Iblock], Jbc, Jac,
                  Jnbs, Jnas, Iac, Ibc, Inas, Inbs);
	      }
	    }

          } /* end loop over Jblock */

          if (Ivec.buf_offdiag[Ibuf]) {
            Iblock2 = Ivec.decode[Ibc][Iac];
            Ivec.transp_block(Iblock, transp_tmp2);
            Iac = Ivec.Ia_code[Iblock2];
            Ibc = Ivec.Ib_code[Iblock2];
            Inas = Ivec.Ia_size[Iblock2];
            Inbs = Ivec.Ib_size[Iblock2];
	   
            for (Jblock=Jvec.first_ablk[Jairr]; Jblock<=Jvec.last_ablk[Jairr];
              Jblock++) {
              Jac = Jvec.Ia_code[Jblock];
              Jbc = Jvec.Ib_code[Jblock];
              Jnas = Jvec.Ia_size[Jblock];
              Jnbs = Jvec.Ib_size[Jblock];
	   
              if (s1_contrib[Iblock2][Jblock] || s2_contrib[Iblock2][Jblock])
                opdm_block(alplist, betlist, onepdm_a, onepdm_b, Jvec.blocks[Jblock],
                           transp_tmp2, Jac, Jbc, Jnas, Jnbs, Iac, Ibc, 
                           Inas, Inbs);

                if (Jvec.buf_offdiag[Jbuf]) {
                  Jblock2 = Jvec.decode[Jbc][Jac];
                  if (s1_contrib[Iblock][Jblock2] || 
                    s2_contrib[Iblock][Jblock2]) {
                    Jvec.transp_block(Jblock, transp_tmp);
                    opdm_block(alplist, betlist, onepdm_a, onepdm_b, transp_tmp,
                      transp_tmp2, Jbc, Jac, Jnbs, Jnas, Iac, Ibc, Inas, Inbs);
                  }
	        }

	      } /* end loop over Jblock */
            } /* end Ivec offdiag */

          } /* end loop over Iblock */
        } /* end loop over Jbuf */
      } /* end loop over Ibuf */
    } /* end icore==2 */

    else {
      printf("opdm: unrecognized core option!\n");
      return;
    }

    // form total density as a sum of alpha and beta components
    add_mat(onepdm_a,onepdm_b,onepdm,populated_orbs,populated_orbs);
    
    /* write and/or print the opdm */
    if (printflag) {
      fprintf(outfile, "\n\nOne-particle ");
      if (transdens)
        fprintf(outfile, "transition ");
      fprintf(outfile, "density matrix MO basis for root %d\n", Jroot+1);
      print_mat(onepdm, populated_orbs, populated_orbs, outfile);
      fprintf(outfile, "\n");
    }


    if (writeflag) {
      sprintf(opdm_key,"MO-basis %s Root %d", transdens ? "TDM" : "OPDM",Jroot);
      psio_write_entry(targetfile, opdm_key, (char *) onepdm[0], 
        populated_orbs * populated_orbs * sizeof(double));
      if (Parameters.print_lvl) 
        fprintf(outfile, "\nWrote MO-basis %s %d to disk\n", 
          transdens ? "TDM" : "OPDM", Jroot+1);

      sprintf(opdm_key,"MO-basis Alpha %s Root %d", transdens ? "TDM" : "OPDM",Jroot);
      psio_write_entry(targetfile, opdm_key, (char *) onepdm_a[0], 
        populated_orbs * populated_orbs * sizeof(double));
      if (Parameters.print_lvl) 
        fprintf(outfile, "\nWrote MO-basis Alpha %s %d to disk\n", 
          transdens ? "TDM" : "OPDM", Jroot+1);

      sprintf(opdm_key,"MO-basis Beta %s Root %d", transdens ? "TDM" : "OPDM",Jroot);
      psio_write_entry(targetfile, opdm_key, (char *) onepdm_b[0], 
        populated_orbs * populated_orbs * sizeof(double));
      if (Parameters.print_lvl) 
        fprintf(outfile, "\nWrote MO-basis Beta %s %d to disk\n", 
          transdens ? "TDM" : "OPDM", Jroot+1);

      /* write it without the "Root n" part if it's the desired root      */
      /* plain old "MO-basis OPDM" is what is searched by the rest of PSI */
      if (Jroot==Parameters.root) {
        sprintf(opdm_key,"MO-basis %s", transdens ? "TDM" : "OPDM");
        psio_write_entry(targetfile, opdm_key, (char *) onepdm[0],
          populated_orbs * populated_orbs * sizeof(double));
        if (Parameters.print_lvl) 
          fprintf(outfile, "Wrote MO-basis %s to disk\n", 
            transdens ? "TDM" : "OPDM");
        
        // print alpha and beta densities for the target root also
        sprintf(opdm_key,"MO-basis Alpha %s", transdens ? "TDM" : "OPDM");
        psio_write_entry(targetfile, opdm_key, (char *) onepdm_a[0],
          populated_orbs * populated_orbs * sizeof(double));
        if (Parameters.print_lvl) 
          fprintf(outfile, "Wrote MO-basis Alpha %s to disk\n", 
            transdens ? "TDM" : "OPDM");
        sprintf(opdm_key,"MO-basis Beta %s", transdens ? "TDM" : "OPDM");
        psio_write_entry(targetfile, opdm_key, (char *) onepdm_b[0],
          populated_orbs * populated_orbs * sizeof(double));
        if (Parameters.print_lvl) 
          fprintf(outfile, "Wrote MO-basis Beta %s to disk\n", 
            transdens ? "TDM" : "OPDM");

      }
      fprintf(outfile, "\n");
    }

    /* Get the kinetic energy if requested */
    if (Parameters.opdm_ke && !transdens) {
      opdm_ke(onepdm);
    }

    /* get the (transition) dipole moment */
    /*
    if (dipmom) {
      mu_x = 0.0; mu_y = 0.0; mu_z = 0.0; 
      // should I be including nuclear contributions to TM's
      for (i=0; i<populated_orbs; i++) {
        for (j=0; j<populated_orbs; j++) {
          mu_x += mux_mo[i][j] * onepdm[i][j];
          mu_y += muy_mo[i][j] * onepdm[i][j];
          mu_z += muz_mo[i][j] * onepdm[i][j];
        }
      }
      fprintf(outfile, "%sipole moment root %d \n", 
        transdens ? "\nTransition d" : "\nD", Jroot+1); 
      if (!transdens) {
        fprintf(outfile, "Nuclear:    %9.5lf x, %9.5lf y, %9.5lf z au\n",
          mux_n, muy_n, muz_n);
        fprintf(outfile, "            %9.5lf x, %9.5lf y, %9.5lf z D\n",
          mux_n*_dipmom_au2debye,muy_n*_dipmom_au2debye,muz_n*_dipmom_au2debye);
      }
      fprintf(outfile, "Electronic: %9.5lf x, %9.5lf y, %9.5lf z au\n",
        mu_x, mu_y, mu_z);
      fprintf(outfile, "            %9.5lf x, %9.5lf y, %9.5lf z D\n",
        mu_x*_dipmom_au2debye, mu_y*_dipmom_au2debye, mu_z*_dipmom_au2debye);
      if (!transdens) {
        mu_x += mux_n; mu_y += muy_n; mu_z += muz_n; 
        fprintf(outfile, "Total:      %9.5lf x, %9.5lf y, %9.5lf z au\n",
          mu_x, mu_y, mu_z);
        fprintf(outfile, "            %9.5lf x, %9.5lf y, %9.5lf z D\n",
          mu_x*_dipmom_au2debye, mu_y*_dipmom_au2debye, 
          mu_z*_dipmom_au2debye);
      }
      mu_tot = sqrt(mu_x * mu_x + mu_y * mu_y + mu_z * mu_z);
      fprintf(outfile, "|mu|  =     %9.5lf au %9.5lf D\n", mu_tot,
        mu_tot * _dipmom_au2debye);
       
      if (Jroot + 1 == Parameters.num_roots) fprintf(outfile, "\n");
    }
    */

    /* Call OEProp here for each root opdm */
    boost::shared_ptr<OEProp> oe(new OEProp());
    boost::shared_ptr<Wavefunction> wfn = 
      Process::environment.reference_wavefunction(); 
    boost::shared_ptr<Matrix> Ca = wfn->Ca(); 
    std::stringstream ss;
    ss << "CI " << (transdens ? "TDM" : "OPDM");
    if (transdens) {
      ss << " Root " << (Iroot+1) << " -> Root " << (Jroot+1); 
    } 
    else {
      ss << " Root " << (Iroot+1); 
    }

    std::stringstream ss_a;
    ss_a << ss.str() << " alpha";

    SharedMatrix opdm_a(new Matrix(ss_a.str(), Ca->colspi(), Ca->colspi())); 
    int mo_offset = 0;
    for (int h = 0; h < Ca->nirrep(); h++) {
      int nmo = CalcInfo.orbs_per_irr[h];
      int nfv = CalcInfo.frozen_uocc[h];
      int nmor = nmo - nfv;
      //int nmo = Ca->colspi()[h];
      //int nmor = nmo - ref->frzvpi()[h];
      if (!nmo || !nmor) continue;
      double** opdmap = opdm_a->pointer(h);
        
      for (int i=0; i<CalcInfo.orbs_per_irr[h]- CalcInfo.frozen_uocc[h]; i++) {
        for (int j=0; j<CalcInfo.orbs_per_irr[h]-
          CalcInfo.frozen_uocc[h]; j++) {
          int i_ci = CalcInfo.reorder[i+mo_offset];
          int j_ci = CalcInfo.reorder[j+mo_offset]; 
          opdmap[i][j] = onepdm_a[i_ci][j_ci];
        } 
      }
      mo_offset += CalcInfo.orbs_per_irr[h];
    }
    oe->set_Da_mo(opdm_a);

    if (Parameters.ref == "ROHF") {
      std::stringstream ss_b;
      ss_b << ss.str() << " beta";
      SharedMatrix opdm_b(new Matrix(ss_b.str(), Ca->colspi(), Ca->colspi())); 
      mo_offset = 0;
      for (int h = 0; h < Ca->nirrep(); h++) {
        int nmo = CalcInfo.orbs_per_irr[h];
        int nfv = CalcInfo.frozen_uocc[h];
        int nmor = nmo - nfv;
        //int nmo = Ca->colspi()[h];
        //int nmor = nmo - ref->frzvpi()[h];
        if (!nmo || !nmor) continue;
        double** opdmbp = opdm_b->pointer(h);
            
        for (int i=0; i<CalcInfo.orbs_per_irr[h]-CalcInfo.frozen_uocc[h]; i++) {
          for (int j=0; j<CalcInfo.orbs_per_irr[h]-
            CalcInfo.frozen_uocc[h]; j++) {
            int i_ci = CalcInfo.reorder[i+mo_offset];
            int j_ci = CalcInfo.reorder[j+mo_offset]; 
            opdmbp[i][j] = onepdm_b[i_ci][j_ci];
          } 
        }
        mo_offset += CalcInfo.orbs_per_irr[h];
      }
      oe->set_Db_mo(opdm_b);
    }

    std::stringstream oeprop_label;
    if (transdens) {
      oeprop_label << "CI ROOT " << (Iroot+1) << " -> ROOT " << (Jroot+1); 
    } 
    else {
      oeprop_label << "CI ROOT " << (Iroot+1); 
    }
    oe->set_title(oeprop_label.str());
    if (!transdens) {
        oe->add("DIPOLE");
        oe->add("MULLIKEN_CHARGES");
        oe->add("NO_OCCUPATIONS");
        if (Parameters.print_lvl > 1) { 
            oe->add("QUADRUPOLE");
        }
    } 
    else {
        oe->add("TRANSITION_DIPOLE");
        if (Parameters.print_lvl > 1) { 
            oe->add("TRANSITION_QUADRUPOLE");
        }
    }
    
    fprintf(outfile, "  ==> Properties %s <==\n", ss.str().c_str());
    oe->compute();

    // std::pair<SharedMatrix,SharedVector> nos = oe->Na_mo();

    // if this is the "special" root, then copy over OEProp 
    // Process::environment variables from the current root into
    // more general locations
    if (Iroot == Parameters.root) {
      std::stringstream ss2;
      ss2 << oeprop_label.str() << " DIPOLE X"; 
      Process::environment.globals["CI DIPOLE X"] = 
        Process::environment.globals[ss2.str()]; 
      ss2.str(std::string());
      ss2 << oeprop_label.str() << " DIPOLE Y"; 
      Process::environment.globals["CI DIPOLE Y"] = 
        Process::environment.globals[ss2.str()]; 
      ss2.str(std::string());
      ss2 << oeprop_label.str() << " DIPOLE Z"; 
      Process::environment.globals["CI DIPOLE Z"] = 
        Process::environment.globals[ss2.str()]; 
      if (Parameters.print_lvl > 1) { 
         ss2.str(std::string());
         ss2 << oeprop_label.str() << " QUADRUPOLE XX"; 
         Process::environment.globals["CI QUADRUPOLE XX"] = 
           Process::environment.globals[ss2.str()]; 
         ss2.str(std::string());
         ss2 << oeprop_label.str() << " QUADRUPOLE YY"; 
         Process::environment.globals["CI QUADRUPOLE YY"] = 
           Process::environment.globals[ss2.str()]; 
         ss2.str(std::string());
         ss2 << oeprop_label.str() << " QUADRUPOLE ZZ"; 
         Process::environment.globals["CI QUADRUPOLE ZZ"] = 
           Process::environment.globals[ss2.str()]; 
         ss2.str(std::string());
         ss2 << oeprop_label.str() << " QUADRUPOLE XY"; 
         Process::environment.globals["CI QUADRUPOLE XY"] = 
           Process::environment.globals[ss2.str()]; 
         ss2.str(std::string());
         ss2 << oeprop_label.str() << " QUADRUPOLE XZ"; 
         Process::environment.globals["CI QUADRUPOLE XZ"] = 
           Process::environment.globals[ss2.str()]; 
         ss2.str(std::string());
         ss2 << oeprop_label.str() << " QUADRUPOLE YZ"; 
         Process::environment.globals["CI QUADRUPOLE YZ"] = 
           Process::environment.globals[ss2.str()]; 
      }
    }

    fflush(outfile);
    if (!transdens) Iroot++;
  } /* end loop over num_roots Jroot */  

  if (writeflag) {
    sprintf(opdm_key,"Num MO-basis %s", transdens ? "TDM" : "OPDM");
    i = Parameters.num_roots; /* num max index, not the number for TDM
                                 b/c we don't write transition 0->0 */
    psio_write_entry(targetfile, opdm_key, (char *) &i, sizeof(int));
    psio_close(targetfile, 1);
  }

  if (transp_tmp != NULL) free_block(transp_tmp);
  if (transp_tmp2 != NULL) free_block(transp_tmp2);
  Ivec.buf_unlock();
  Jvec.buf_unlock();
  free(buffer1);
  free(buffer2);

  /*
  if (dipmom) {
    free_block(mux_mo);
    free_block(muy_mo);
    free_block(muz_mo);
  }
  */

  if (transdens) {
    fflush(outfile);
    free_block(opdm_blk);
    return;
  }

  /* Average the opdm's */
  /* if (Parameters.opdm_diag) rfile(targetfile); */
  if (Parameters.opdm_ave) {
    psio_open(targetfile, PSIO_OPEN_OLD); 
    ave(targetfile); 
    psio_close(targetfile, 1);
  }

  /* get CI Natural Orbitals */
  if (Parameters.opdm_diag) {

    psio_open(targetfile, PSIO_OPEN_OLD);
    chkpt_init(PSIO_OPEN_OLD);

    /* reorder opdm from ci to pitzer and diagonalize each 
    ** symmetry blk
    */
    if (Parameters.opdm_ave) {
      klast = 1;
    }
    else {
      klast = Parameters.num_roots;
    }

    /* loop over roots or averaged opdm */
    for(k=0; k<klast; k++) {
      if (Parameters.opdm_ave && Parameters.print_lvl > 1) {
        fprintf(outfile,"\n\n\t\t\tCI Natural Orbitals for the Averaged\n");
        fprintf(outfile,"\t\t\tOPDM of %d Roots in terms of Molecular"
                 " Orbitals\n\n",k); 
      }
      else if (Parameters.print_lvl > 1) {
        fprintf(outfile,
             "\n\t\t\tCI Natural Orbitals in terms of Molecular Orbitals\n\n");
        fprintf(outfile,"\t\t\t Root %d\n\n",k+1);
        fflush(outfile);
      }

      mo_offset = 0;

      if (!Parameters.opdm_ave)
        sprintf(opdm_key, "MO-basis OPDM Root %d", k);
      else 
        sprintf(opdm_key, "MO-basis OPDM Ave");

      psio_read_entry(targetfile, opdm_key, (char *) onepdm[0],
        populated_orbs * populated_orbs * sizeof(double));

      for (irrep=0; irrep<CalcInfo.nirreps; irrep++) {
        if (CalcInfo.orbs_per_irr[irrep] == 0) continue; 
        for (i=0; i<CalcInfo.orbs_per_irr[irrep]-
                    CalcInfo.frozen_uocc[irrep]; i++) {
          for (j=0; j<CalcInfo.orbs_per_irr[irrep]-
                    CalcInfo.frozen_uocc[irrep]; j++) {
            i_ci = CalcInfo.reorder[i+mo_offset];
            j_ci = CalcInfo.reorder[j+mo_offset]; 
            opdm_blk[i][j] = onepdm[i_ci][j_ci];
          } 
        }
 
        /* Writing SCF vector to OPDM file for safe keeping 
         * because you might overwrite it with NO's for earlier root 
         * and we're only storing one irrep in core.  Only need to do once.
         */
        if (k==0) {

          scfvec = chkpt_rd_scf_irrep(irrep);

            #ifdef DEBUG
            fprintf(outfile,"Cvec for k==0, read in from chkpt original\n");
            fprintf(outfile," %s Block \n", CalcInfo.labels[irrep]);
            print_mat(scfvec, CalcInfo.orbs_per_irr[irrep],
                      CalcInfo.orbs_per_irr[irrep], outfile);
            #endif

          sprintf(opdm_key, "Old SCF Matrix Irrep %d", irrep);
          psio_write_entry(targetfile, opdm_key, (char *) scfvec[0],
            CalcInfo.orbs_per_irr[irrep] * CalcInfo.orbs_per_irr[irrep] *
            sizeof(double));
	  free_block(scfvec);
        }

        zero_mat(opdm_eigvec, max_orb_per_irrep, max_orb_per_irrep);

        if (CalcInfo.orbs_per_irr[irrep]-CalcInfo.frozen_uocc[irrep] > 0) {
          sq_rsp(CalcInfo.orbs_per_irr[irrep]-CalcInfo.frozen_uocc[irrep],
                 CalcInfo.orbs_per_irr[irrep]-CalcInfo.frozen_uocc[irrep],
                 opdm_blk, opdm_eigval, 1, opdm_eigvec, TOL); 
          }
        for (i=CalcInfo.orbs_per_irr[irrep]-CalcInfo.frozen_uocc[irrep]; 
             i<CalcInfo.orbs_per_irr[irrep]; i++) {
           opdm_eigvec[i][i] = 1.0;
           opdm_eigval[i] = 0.0;
           }
        eigsort(opdm_eigval, opdm_eigvec, -(CalcInfo.orbs_per_irr[irrep]));



// FAE
// eigsort sometimes will swap the order of orbitals, for example
// in frozen core computations the focc may be mixed with docc and 
// change the final result

//Loop over "populated"
  for (i=0;i<CalcInfo.orbs_per_irr[irrep]-CalcInfo.frozen_uocc[irrep];i++){
    max_overlap = 0;
    int m       = 0;
     for (j=i;j<CalcInfo.orbs_per_irr[irrep]-CalcInfo.frozen_uocc[irrep];j++){
       overlap = opdm_eigvec[i][j] * opdm_eigvec[i][j];
       if(overlap > max_overlap){
         m = j;
         max_overlap = overlap;
       }
     }
     for (j=0;j<CalcInfo.orbs_per_irr[irrep];j++){
         double temporary  = opdm_eigvec[j][i];
         opdm_eigvec[j][i] = opdm_eigvec[j][m];
         opdm_eigvec[j][m] = temporary;
     }
     double temporary = opdm_eigval[i];
     opdm_eigval[i] = opdm_eigval[m];
     opdm_eigval[m] = temporary;

  }
  // End FAE changes, May 3 2007





        if (Parameters.print_lvl > 0) {
          if (irrep==0) {
            if (Parameters.opdm_ave) { 
              fprintf(outfile, "\n Averaged CI Natural Orbitals in terms "
                "of Molecular Orbitals\n\n");
              }
            else fprintf(outfile, "\n CI Natural Orbitals in terms of "
                   "Molecular Orbitals: Root %d\n\n", k+1);
          }
          fprintf(outfile,"\n %s Block (MO basis)\n", CalcInfo.labels[irrep]);
          eivout(opdm_eigvec, opdm_eigval, CalcInfo.orbs_per_irr[irrep],
                 CalcInfo.orbs_per_irr[irrep], outfile);
        }

        /* Write them if we need to */
        if (Parameters.opdm_wrtnos && (k==Parameters.opdm_orbs_root)) {
          if (irrep==0) {
            if (!Parameters.opdm_ave) {
              fprintf(outfile,"\n Writing CI Natural Orbitals for root %d"
                      " to checkpoint in terms of Symmetry Orbitals\n\n",k+1);
            }
          }
	  
	  scfvec = block_matrix(CalcInfo.orbs_per_irr[irrep], 
                       CalcInfo.orbs_per_irr[irrep]);
          sprintf(opdm_key, "Old SCF Matrix Irrep %d", irrep);
          psio_read_entry(targetfile, opdm_key, (char *) scfvec[0],
            CalcInfo.orbs_per_irr[irrep] * CalcInfo.orbs_per_irr[irrep] *
            sizeof(double));

          #ifdef DEBUG
          fprintf(outfile,"\nCvec read for MO to SO trans\n\n");
          fprintf(outfile," %s Block \n", CalcInfo.labels[irrep]);
          print_mat(scfvec, CalcInfo.orbs_per_irr[irrep],
                    CalcInfo.orbs_per_irr[irrep], outfile);
          fprintf(outfile,"\nOpdm_eigvec before MO to SO trans\n\n");
          fprintf(outfile," %s Block \n", CalcInfo.labels[irrep]);
          print_mat(opdm_eigvec, CalcInfo.orbs_per_irr[irrep],
                    CalcInfo.orbs_per_irr[irrep], outfile); 
          #endif
          mmult(scfvec, 0, opdm_eigvec, 0, opdm_blk, 0,
                CalcInfo.so_per_irr[irrep], CalcInfo.orbs_per_irr[irrep],
                CalcInfo.orbs_per_irr[irrep], 0); 
          free_block(scfvec);
          if (Parameters.print_lvl > 0) {  // FAE
            fprintf(outfile," %s Block (SO basis)\n", CalcInfo.labels[irrep]);
            print_mat(opdm_blk, CalcInfo.so_per_irr[irrep],
                      CalcInfo.orbs_per_irr[irrep], outfile);
          }
          chkpt_wt_scf_irrep(opdm_blk, irrep);
          fprintf(outfile, "\n Warning: Natural Orbitals for the ");
	  if (Parameters.opdm_ave)
            fprintf(outfile, "Averaged OPDM ");
          else
            fprintf(outfile, "Root %d OPDM ", k);
          fprintf(outfile, "have been written to the checkpoint file!\n\n"); 
        } /* end code to write the NO's to disk */
        mo_offset += CalcInfo.orbs_per_irr[irrep];
      } /* end loop over irreps */
    } /* end loop over roots */

    free_block(onepdm);
    free_block(onepdm_a);
    free_block(onepdm_b);
    free_block(opdm_eigvec);
    free(opdm_eigval); 
    chkpt_close();
    psio_close(targetfile, 1);
  } /* CINOS completed */


  fflush(outfile);
  free_block(opdm_blk);

}

void opdm_block(struct stringwr **alplist, struct stringwr **betlist,
		double **onepdm_a, double **onepdm_b, double **CJ, double **CI, int Ja_list, 
		int Jb_list, int Jnas, int Jnbs, int Ia_list, int Ib_list, 
		int Inas, int Inbs)
{
  int Ia_idx, Ib_idx, Ja_idx, Jb_idx, Ja_ex, Jb_ex, Jbcnt, Jacnt; 
  struct stringwr *Jb, *Ja;
  signed char *Jbsgn, *Jasgn;
  unsigned int *Jbridx, *Jaridx;
  double C1, C2, Ib_sgn, Ia_sgn;
  int i, j, oij, nfzc, *Jboij, *Jaoij;
 
  nfzc = CalcInfo.num_fzc_orbs;

  /* loop over Ia in Ia_list */
  if (Ia_list == Ja_list) {
    for (Ia_idx=0; Ia_idx<Inas; Ia_idx++) {
      for (Jb=betlist[Jb_list], Jb_idx=0; Jb_idx<Jnbs; Jb_idx++, Jb++) {
	C1 = CJ[Ia_idx][Jb_idx];

	/* loop over excitations E^b_{ij} from |B(J_b)> */
	Jbcnt = Jb->cnt[Ib_list];
	Jbridx = Jb->ridx[Ib_list];
	Jbsgn = Jb->sgn[Ib_list];
	Jboij = Jb->oij[Ib_list];
	for (Jb_ex=0; Jb_ex < Jbcnt; Jb_ex++) {
	  oij = *Jboij++;
	  Ib_idx = *Jbridx++;
	  Ib_sgn = (double) *Jbsgn++;
	  C2 = CI[Ia_idx][Ib_idx];
          i = oij/CalcInfo.num_ci_orbs + nfzc;
          j = oij%CalcInfo.num_ci_orbs + nfzc;
	  onepdm_b[i][j] += C1 * C2 * Ib_sgn;
	}
      }
    }
  }

  /* loop over Ib in Ib_list */
  if (Ib_list == Jb_list) {
    for (Ib_idx=0; Ib_idx<Inbs; Ib_idx++) {
      for (Ja=alplist[Ja_list], Ja_idx=0; Ja_idx<Jnas; Ja_idx++, Ja++) {
	C1 = CJ[Ja_idx][Ib_idx];
	
	/* loop over excitations */
	Jacnt = Ja->cnt[Ia_list];
	Jaridx = Ja->ridx[Ia_list];
	Jasgn = Ja->sgn[Ia_list];
	Jaoij = Ja->oij[Ia_list];
	for (Ja_ex=0; Ja_ex < Jacnt; Ja_ex++) {
	  oij = *Jaoij++;
	  Ia_idx = *Jaridx++;
	  Ia_sgn = (double) *Jasgn++;
	  C2 = CI[Ia_idx][Ib_idx];
          i = oij/CalcInfo.num_ci_orbs + nfzc;
          j = oij%CalcInfo.num_ci_orbs + nfzc;
	  onepdm_a[i][j] += C1 * C2 * Ia_sgn;
	}
      }
    }
  }
}


/*
** ave()
** 
** Parameters:
**  targetfile = file number to obtain matrices from 
**
** Modified 7/13/04 to use new state average parameters --- CDS
**
*/
void ave(int targetfile)
{
  int root, root_count, i, j, populated_orbs;
  double **tmp_mat1, **tmp_mat2;
  char opdm_key[80];
  
  const char spinlabels[][10] = { "", "Alpha ", "Beta " };
  const int nspincases = 3;

  populated_orbs = CalcInfo.nmo-CalcInfo.num_fzv_orbs;
  tmp_mat1 = block_matrix(populated_orbs, populated_orbs);  
  tmp_mat2 = block_matrix(populated_orbs, populated_orbs);
  
  for(int spincase=0; spincase<nspincases; ++spincase) {
    
    zero_mat(tmp_mat1, populated_orbs, populated_orbs);

    for(root_count=0; root_count<Parameters.average_num; root_count++) {

      root = Parameters.average_states[root_count];

      sprintf(opdm_key, "MO-basis %sOPDM Root %d", spinlabels[spincase], root);

      psio_read_entry(targetfile, opdm_key, (char *) tmp_mat2[0],
                      populated_orbs * populated_orbs * sizeof(double));

      if (Parameters.opdm_print) {
        fprintf(outfile,"\n\n\t\t%sOPDM for Root %d",spinlabels[spincase],root+1);
        print_mat(tmp_mat2, populated_orbs, populated_orbs, outfile);
      }

      for (i=0; i<populated_orbs; i++)
        for (j=0; j<populated_orbs; j++)
          tmp_mat1[i][j] += Parameters.average_weights[root]*tmp_mat2[i][j];

    }

    sprintf(opdm_key, "MO-basis %sOPDM", spinlabels[spincase]);
    psio_write_entry(targetfile, opdm_key, (char *) tmp_mat1[0],
                     populated_orbs*populated_orbs*sizeof(double));
    
    if (Parameters.print_lvl> 0 || Parameters.opdm_print) {
      fprintf(outfile,
              "\n\t\t\t Averaged %sOPDM's for %d Roots written to opdm_file \n\n",
              spinlabels[spincase], Parameters.average_num);
    }
    if (Parameters.opdm_print) {
      print_mat(tmp_mat1, populated_orbs, populated_orbs, outfile);
    }
    
  } // loop over spincases

  free_block(tmp_mat1);
  free_block(tmp_mat2);

}


/*
** opdm_ke
**
** Compute the kinetic energy contribution from the correlated part of the
** one-particle density matrix.  For Daniel Crawford
*/
void opdm_ke(double **onepdm)
{
  int errcod;
  int src_T_file, mo_offset, so_offset, irrep, i, j, i_ci, j_ci, i2, j2, ij;
  int maxorbs;
  int noeints;
  double *T, **scfmat, **opdm_blk, **tmp_mat, ke, kei;

  ke = kei = 0.0;

  /* read in the kinetic energy integrals */
  noeints = CalcInfo.nso*(CalcInfo.nso+1)/2;
  src_T_file = PSIF_OEI;

  T = init_array(noeints);
  if (Parameters.print_lvl>2) 
    fprintf(outfile, "Kinetic energy integrals (SO basis):\n");
  errcod = iwl_rdone(src_T_file,PSIF_SO_T,T,noeints,0,
                     (Parameters.print_lvl>2),outfile);
  if (!errcod) {
    printf("(detci): Error reading kinetic energy ints\n");
    exit(1);
  }

  /* find biggest blocksize */
  for (irrep=0,maxorbs=0; irrep<CalcInfo.nirreps; irrep++) {
    if (CalcInfo.orbs_per_irr[irrep] > maxorbs)
      maxorbs = CalcInfo.so_per_irr[irrep];
  }
  opdm_blk = block_matrix(maxorbs, maxorbs);
  tmp_mat = block_matrix(maxorbs, maxorbs);

  /* transform the onepdm into SO form, one irrep at a time */
  so_offset = mo_offset = 0;
  fprintf(outfile,"Correlation Piece of OPDM in SO basis\n");
  for (irrep=0; irrep<CalcInfo.nirreps; irrep++) {
    if (CalcInfo.orbs_per_irr[irrep] == 0) continue;
    for (i=0; i<CalcInfo.orbs_per_irr[irrep]; i++) {
      for (j=0; j<CalcInfo.orbs_per_irr[irrep]; j++) {
        i_ci = CalcInfo.reorder[i+mo_offset];
        j_ci = CalcInfo.reorder[j+mo_offset];
        opdm_blk[i][j] = onepdm[i_ci][j_ci];
      }
    }
    /* need to subtract out reference piece, assume single det ref */
    for (i=0; i<CalcInfo.docc[irrep]; i++) {
      opdm_blk[i][i] -= 2.0;
    }
    /* keep counting from i to take out socc part */
    for (j=0; j<CalcInfo.socc[irrep]; j++,i++) {
      opdm_blk[i][i] -= 1.0;
    }

    
    fprintf(outfile, "Irrep %d\n", irrep);
    fprintf(outfile, "MO basis, Pitzer order\n");
    print_mat(opdm_blk,CalcInfo.so_per_irr[irrep],
              CalcInfo.so_per_irr[irrep],outfile);

    /* transform back to SO basis */
    scfmat = chkpt_rd_scf_irrep(irrep);
    mmult(opdm_blk,0,scfmat,1,tmp_mat,0,CalcInfo.orbs_per_irr[irrep],
          CalcInfo.orbs_per_irr[irrep],CalcInfo.so_per_irr[irrep],0);
    mmult(scfmat,0,tmp_mat,0,opdm_blk,0,CalcInfo.so_per_irr[irrep],
          CalcInfo.orbs_per_irr[irrep],CalcInfo.so_per_irr[irrep],0); 
    
    fprintf(outfile, "SO basis, Pitzer order\n");
    print_mat(opdm_blk,CalcInfo.so_per_irr[irrep],
              CalcInfo.so_per_irr[irrep],outfile);

    /* get kinetic energy contribution */
    kei = 0.0;
    for (i=0,i2=so_offset; i<CalcInfo.so_per_irr[irrep]; i++,i2++) {
      for (j=0,j2=so_offset; j<CalcInfo.so_per_irr[irrep]; j++,j2++) {
        ij = ioff[MAX0(i2,j2)] + MIN0(i2,j2);
        kei += opdm_blk[i][j] * T[ij];
      }
    }

    fprintf(outfile,"Contribution to correlation kinetic energy = %15.10lf\n", 
            kei);

    ke += kei;
    mo_offset += CalcInfo.orbs_per_irr[irrep];
    so_offset += CalcInfo.so_per_irr[irrep];

  } /* end loop over irreps */

  fprintf(outfile, "\nTotal correlation kinetic energy = %15.10lf\n", ke);

}

/*
**
** Note: This function is deprecated (and will no longer work in PSI4).
** 
** Its functionality is now provided directly by PSI4
**
** This function will read in the dipole moment integrals from AO basis
** off disk (obtained from cints --oeprop) and transform them to the MO
** basis (CI ordering) for subsequent contraction with the density or
** transition density matrices.  Some code adapted from ccdensity/dipole.c
*/
/*
void get_mo_dipmom_ints(double **MUX_MO, double **MUY_MO, double **MUZ_MO)
{
  int nao, nso, nmo, noei;
  int i, j, I, stat;
  double *mu_x_ints, *mu_y_ints, *mu_z_ints;
  double **usotao, **scf_pitzer, **scf_mo;
  double **MUX_AO, **MUY_AO, **MUZ_AO;
  double **MUX_SO, **MUY_SO, **MUZ_SO;
  double **X;

  // Run dip mom ints if needed 
  stat = 0;
  psio_open(PSIF_OEI, PSIO_OPEN_OLD);
  if (psio_tocscan(PSIF_OEI, PSIF_AO_MX) == NULL) 
    stat = 1; // not on disk yet 
  psio_close(PSIF_OEI, 1);  

  if (stat && system("cints --oeprop")) {
    fprintf(outfile, "DETCI (get_mo_dipmom_ints): Can't run cints --oeprop\n");
    exit(1);
  }

  chkpt_init(PSIO_OPEN_OLD);
  nso = chkpt_rd_nso();
  nao = chkpt_rd_nao();
  nmo = chkpt_rd_nmo();
  usotao = chkpt_rd_usotao();
  scf_pitzer = chkpt_rd_scf();                                                         
  chkpt_close();

  // reorder SCF eigenvectors to MO (correlated CI) ordering
  scf_mo = block_matrix(nso, nmo);
  for (i=0; i<nmo; i++) {
    I = CalcInfo.reorder[i]; // Pitzer -> MO ordering
    for (j=0; j<nso; j++) {
      scf_mo[j][I] = scf_pitzer[j][i];
    }
  }
  free_block(scf_pitzer);

  // Read in dipole moment integrals in the AO basis 
  noei = nao * (nao + 1)/2;
                                                                                
  mu_x_ints = init_array(noei);
  stat = iwl_rdone(PSIF_OEI,PSIF_AO_MX,mu_x_ints,noei,0,0,outfile);
  mu_y_ints = init_array(noei);
  stat = iwl_rdone(PSIF_OEI,PSIF_AO_MY,mu_y_ints,noei,0,0,outfile);
  mu_z_ints = init_array(noei);
  stat = iwl_rdone(PSIF_OEI,PSIF_AO_MZ,mu_z_ints,noei,0,0,outfile);
                                                                                
  MUX_AO = block_matrix(nao,nao);
  MUY_AO = block_matrix(nao,nao);
  MUZ_AO = block_matrix(nao,nao);

  MUX_SO = block_matrix(nso,nso);
  MUY_SO = block_matrix(nso,nso);
  MUZ_SO = block_matrix(nso,nso);

  for(i=0; i<nao; i++) {
    for(j=0; j<nao; j++) {
      MUX_AO[i][j] = mu_x_ints[INDEX(i,j)];
      MUY_AO[i][j] = mu_y_ints[INDEX(i,j)];
      MUZ_AO[i][j] = mu_z_ints[INDEX(i,j)];
    }
  }

  // Transform the AO dipole integrals to the SO basis 
  X = block_matrix(nso,nao); // just a temporary matrix
                                                                                
  C_DGEMM('n','n',nso,nao,nao,1.0,&(usotao[0][0]),nao,&(MUX_AO[0][0]),nao,
          0,&(X[0][0]),nao);
  C_DGEMM('n','t',nso,nso,nao,1.0,&(X[0][0]),nao,&(usotao[0][0]),nao,
          0,&(MUX_SO[0][0]),nso);
                                                                                
  C_DGEMM('n','n',nso,nao,nao,1.0,&(usotao[0][0]),nao,&(MUY_AO[0][0]),nao,
          0,&(X[0][0]),nao);
  C_DGEMM('n','t',nso,nso,nao,1.0,&(X[0][0]),nao,&(usotao[0][0]),nao,
          0,&(MUY_SO[0][0]),nso);
                                                                                
  C_DGEMM('n','n',nso,nao,nao,1.0,&(usotao[0][0]),nao,&(MUZ_AO[0][0]),nao,
          0,&(X[0][0]),nao);
  C_DGEMM('n','t',nso,nso,nao,1.0,&(X[0][0]),nao,&(usotao[0][0]),nao,
          0,&(MUZ_SO[0][0]),nso);
                                                                                
  free(mu_x_ints); free(mu_y_ints); free(mu_z_ints);
  free_block(MUX_AO); free_block(MUY_AO); free_block(MUZ_AO);
  free_block(X);


  // Transform the SO dipole integrals to the MO basis 
                                                                                
  X = block_matrix(nmo,nso); // just a temporary matrix
                                                                                
  C_DGEMM('t','n',nmo,nso,nso,1.0,&(scf_mo[0][0]),nmo,&(MUX_SO[0][0]),nso,
          0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1.0,&(X[0][0]),nso,&(scf_mo[0][0]),nmo,
          0,&(MUX_MO[0][0]),nmo);
                                                                                
  C_DGEMM('t','n',nmo,nso,nso,1.0,&(scf_mo[0][0]),nmo,&(MUY_SO[0][0]),nso,
          0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1.0,&(X[0][0]),nso,&(scf_mo[0][0]),nmo,
          0,&(MUY_MO[0][0]),nmo);

  C_DGEMM('t','n',nmo,nso,nso,1.0,&(scf_mo[0][0]),nmo,&(MUZ_SO[0][0]),nso,
          0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1.0,&(X[0][0]),nso,&(scf_mo[0][0]),nmo,
          0,&(MUZ_MO[0][0]),nmo);

  free_block(scf_mo);
  free_block(MUX_SO); free_block(MUY_SO); free_block(MUZ_SO);
  free_block(X);

  return;
}
*/

/*
** Note: This function is deprecated (and may no longer work in PSI4).
** 
** Its functionality is now provided directly by PSI4
** get the nuclear part of the dipole moment 
*/
/*
void get_dipmom_nuc(double *mu_x_n, double *mu_y_n, double *mu_z_n)
{
  int i, natom;
  double *zvals, **geom;
  double x, y, z;

  chkpt_init(PSIO_OPEN_OLD);
  natom = chkpt_rd_natom();
  zvals = chkpt_rd_zvals();
  geom = chkpt_rd_geom();
  chkpt_close();

  *mu_x_n = 0.0; *mu_y_n = 0.0; *mu_z_n = 0.0;
  for(i=0;i<natom;i++) {
    *mu_x_n += zvals[i]*geom[i][0];
    *mu_y_n += zvals[i]*geom[i][1];
    *mu_z_n += zvals[i]*geom[i][2];
  }

  free(zvals);
  free_block(geom);

}  
*/

}} // namespace psi::detci

