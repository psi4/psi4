/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here 
*/

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
#include "civect.h"
#include "ciwave.h"

namespace psi { namespace detci {

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define TOL 1E-14



void CIWavefunction::form_opdm(void)
{

  std::vector<std::vector<SharedMatrix> > opdm_list;
  /* don't need Parameters_->root since it writes all opdm's */
  if (Parameters_->transdens) {
    opdm_list = opdm(0, Parameters_->num_roots, Parameters_->first_d_tmp_unit, 
                     Parameters_->first_d_tmp_unit, true);
    for (int i=0; i<Parameters_->num_roots-1; i++){
        opdm_map_[opdm_list[i][0]->name()] = opdm_list[i][0];
        opdm_map_[opdm_list[i][1]->name()] = opdm_list[i][1];
        opdm_map_[opdm_list[i][2]->name()] = opdm_list[i][2];
    }
  }
  if (Parameters_->opdm) {
    opdm_list = opdm(0, Parameters_->num_roots, Parameters_->first_d_tmp_unit, 
                     Parameters_->first_d_tmp_unit, false);
    for (int i=0; i<Parameters_->num_roots; i++){
        opdm_map_[opdm_list[i][0]->name()] = opdm_list[i][0];
        opdm_map_[opdm_list[i][1]->name()] = opdm_list[i][1];
        opdm_map_[opdm_list[i][2]->name()] = opdm_list[i][2];
    }

    // Figure out which OPDM should be current
    if (Parameters_->opdm_ave){
        Dimension act_dim = get_dimension("ACT");
        opdm_a_ = SharedMatrix(new Matrix("MO-basis Alpha OPDM", nirrep_, act_dim, act_dim));
        opdm_b_ = SharedMatrix(new Matrix("MO-basis Beta OPDM", nirrep_, act_dim, act_dim));
        opdm_   = SharedMatrix(new Matrix("MO-basis OPDM", nirrep_, act_dim, act_dim));
    
        for(int i=0; i<Parameters_->average_num; i++) {
            int croot = Parameters_->average_states[i];
            double weight = Parameters_->average_weights[i];
            opdm_a_->axpy(weight, opdm_list[croot][0]);
            opdm_b_->axpy(weight, opdm_list[croot][1]);
            opdm_->axpy(weight, opdm_list[croot][2]);
        }
    }
    else{
        int croot = Parameters_->root;
        opdm_a_ = opdm_list[croot][0]->clone();
        opdm_b_ = opdm_list[croot][1]->clone();
        opdm_ = opdm_list[croot][2]->clone();
    }
    Da_ = opdm_add_inactive(opdm_a_, 1.0, true);
    Db_ = opdm_add_inactive(opdm_b_, 1.0, true);
  }

}
/*
** Resizes an act by act OPDM matrix to include the inactive porition.
** The diagonal of the inactive portion is set to value. If virt is true
** the return OPDM with of size nmo x nmo. 
*/
SharedMatrix CIWavefunction::opdm_add_inactive(SharedMatrix opdm, double value, bool virt)
{

    Dimension dim_drc = get_dimension("DRC");
    Dimension dim_act = get_dimension("ACT");
    Dimension dim_ia = dim_drc + dim_act;
    Dimension dim_ret;
    
    if (virt){
        dim_ret = dim_ia + get_dimension("DRV");
    }
    else{
        dim_ret = dim_ia;
    }

    
    SharedMatrix ret(new Matrix(opdm->name(), dim_ret, dim_ret));

    for (int h=0; h<nirrep_; h++){
        if (!dim_ia[h]) continue;
        
        double** retp = ret->pointer(h);
        double** opdmp = opdm->pointer(h);

        // Diagonal inact part
        for (int i=0; i<dim_drc[h]; i++){
            retp[i][i] = value;
        }

        // Non-diagonal act part
        for (int t=0; t<dim_act[h]; t++){
            for (int u=0; u<dim_act[h]; u++){
                retp[dim_drc[h] + t][dim_drc[h] + u] = opdmp[t][u];
            }
        }

    }
    return ret;

}

/*
** Computes the one-particle density matrix for all nroots starting at root_start.
** If transden is true, opdm will then compute transition matrices for Iroot = root_start
** and Jroot from Iroot to Iroot + Nroots.
** \gamma_{tu} = < Iroot | \hat{E}^{tu} | Jroot >
*/ 
std::vector<std::vector<SharedMatrix> > CIWavefunction::opdm(int root_start, int nroots, int Ifile,
                                                             int Jfile, bool transden)
{

  //CIvect Ivec, Jvec;
  int i, j, k, l, klast, roots;
  int maxrows, maxcols;
  unsigned long bufsz;
  int max_orb_per_irrep;
  double **transp_tmp = NULL;
  double **transp_tmp2 = NULL;
  double *buffer1, *buffer2, **onepdm, **onepdm_a, **onepdm_b;
  int i_ci, j_ci, irrep, mo_offset, so_offset, orb_length=0, opdm_length=0;
  int Iblock, Iblock2, Ibuf, Iac, Ibc, Inas, Inbs, Iairr;
  int Jblock, Jblock2, Jbuf, Jac, Jbc, Jnas, Jnbs, Jairr;
  int do_Jblock, do_Jblock2;
  int populated_orbs;
  double overlap, max_overlap;
  char opdm_key[80]; /* libpsio TOC entry name for OPDM for each root */

  int Iroot = root_start;
  int Jroot;
  if (transden) Jroot = root_start + 1;
  else Jroot = root_start;

  CIvect Ivec(Parameters_->icore, nroots, 1, Ifile, CIblks_, CalcInfo_, Parameters_,
              H0block_, false);
  Ivec.init_io_files(true);

  CIvect Jvec(Parameters_->icore, nroots, 1, Jfile, CIblks_, CalcInfo_, Parameters_,
              H0block_, false);
  Jvec.init_io_files(true);

  std::vector<std::vector<SharedMatrix> > opdm_list;


  buffer1 = Ivec.buf_malloc();
  buffer2 = Jvec.buf_malloc();
  Ivec.buf_lock(buffer1);
  Jvec.buf_lock(buffer2);

  if ((Ivec.icore_==2 && Ivec.Ms0_ && CalcInfo_->ref_sym != 0) || 
      (Ivec.icore_==0 && Ivec.Ms0_)) {
    for (i=0, maxrows=0, maxcols=0; i<Ivec.num_blocks_; i++) {
      if (Ivec.Ia_size_[i] > maxrows) maxrows = Ivec.Ia_size_[i];
      if (Ivec.Ib_size_[i] > maxcols) maxcols = Ivec.Ib_size_[i];
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
 

  /* if we're trying to follow a root, figure out which one here */
  //if (Parameters_->follow_vec_num > 0 && !transdens) {
  //  max_overlap = 0.0;
  //  for (i=0,j=0; i<Parameters_->num_roots; i++) {

  //    overlap = Ivec.compute_follow_overlap(i, 
  //      Parameters_->follow_vec_num, Parameters_->follow_vec_coef,
  //      Parameters_->follow_vec_Iac, Parameters_->follow_vec_Iaridx, 
  //      Parameters_->follow_vec_Ibc, Parameters_->follow_vec_Ibridx);

  //    if (overlap > max_overlap) {
  //      max_overlap = overlap;
  //      j = i;
  //    } 
  //  }

  //  Parameters_->root = j;
  //  outfile->Printf( "DETCI following root %d (overlap %6.4lf)\n", 
  //    Parameters_->root+1,max_overlap);

  //  for (i=0; i<Parameters_->average_num; i++) {
  //    Parameters_->average_states[i] = 0;
  //    Parameters_->average_weights[i] = 0.0;
  //  }
  //  Parameters_->average_num = 1;
  //  Parameters_->average_states[0] = Parameters_->root;
  //  Parameters_->average_weights[0] = 1.0;
  //  
  //  /* correct "the" energy in the checkpoint file */
  //  Process::environment.globals["CURRENT ENERGY"] = overlap;
  //  Process::environment.globals["CI TOTAL ENERGY"] = overlap;

  //  // eref is wrong for open-shells so replace it with escf until
  //  // I fix it, CDS 11/5/11
  //  Process::environment.globals["CI CORRELATION ENERGY"] = overlap - CalcInfo_->escf;
  //  Process::environment.globals["CURRENT CORRELATION ENERGY"] = overlap - CalcInfo_->escf;
  //  Process::environment.globals["CURRENT REFERENCE ENERGY"] = CalcInfo_->escf;
  //
  //}
  int nci = CalcInfo_->num_ci_orbs;
  Dimension act_dim = get_dimension("ACT");
  SharedMatrix scratch_a(new Matrix("OPDM A Scratch", nci, nci));
  SharedMatrix scratch_b(new Matrix("OPDM B Scratch", nci, nci));
  double** scratch_ap = scratch_a->pointer();
  double** scratch_bp = scratch_b->pointer();

  for (; Jroot<nroots; Jroot++) {
  
    scratch_a->zero();
    scratch_b->zero();
     

    if (Parameters_->icore == 0) {
 
      for (Ibuf=0; Ibuf<Ivec.buf_per_vect_; Ibuf++) {
        Ivec.read(Iroot, Ibuf);
        Iblock = Ivec.buf2blk_[Ibuf];
        Iac = Ivec.Ia_code_[Iblock];
        Ibc = Ivec.Ib_code_[Iblock];
        Inas = Ivec.Ia_size_[Iblock];
        Inbs = Ivec.Ib_size_[Iblock];
       
        for (Jbuf=0; Jbuf<Jvec.buf_per_vect_; Jbuf++) {
          do_Jblock=0; do_Jblock2=0;
          Jblock = Jvec.buf2blk_[Jbuf];
          Jblock2 = -1;
          Jac = Jvec.Ia_code_[Jblock];
          Jbc = Jvec.Ib_code_[Jblock];
          if (Jvec.Ms0_) Jblock2 = Jvec.decode_[Jbc][Jac];
          Jnas = Jvec.Ia_size_[Jblock];
          Jnbs = Jvec.Ib_size_[Jblock];
          if (s1_contrib_[Iblock][Jblock] || s2_contrib_[Iblock][Jblock]) 
            do_Jblock = 1;
          if (Jvec.buf_offdiag_[Jbuf] && (s1_contrib_[Iblock][Jblock2] ||
                                         s2_contrib_[Iblock][Jblock2]))
            do_Jblock2 = 1;
          if (!do_Jblock && !do_Jblock2) continue;
	 
          Jvec.read(Jroot, Jbuf);
	 
          if (do_Jblock) {
            opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, Jvec.blocks_[Jblock], 
                       Ivec.blocks_[Iblock], Jac, Jbc, Jnas,
                       Jnbs, Iac, Ibc, Inas, Inbs);
            }
	 
          if (do_Jblock2) {
            Jvec.transp_block(Jblock, transp_tmp);
            opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, transp_tmp,
                       Ivec.blocks_[Iblock], Jbc, Jac, Jnbs,
                       Jnas, Iac, Ibc, Inas, Inbs);
          }
	 
        } /* end loop over Jbuf */
       
        if (Ivec.buf_offdiag_[Ibuf]) { /* need to get contrib of transpose */
          Iblock2 = Ivec.decode_[Ibc][Iac];
          Iac = Ivec.Ia_code_[Iblock2];
          Ibc = Ivec.Ib_code_[Iblock2];
          Inas = Ivec.Ia_size_[Iblock2];
          Inbs = Ivec.Ib_size_[Iblock2];
       
          Ivec.transp_block(Iblock, transp_tmp2);

          for (Jbuf=0; Jbuf<Jvec.buf_per_vect_; Jbuf++) {
            do_Jblock=0; do_Jblock2=0;
            Jblock = Jvec.buf2blk_[Jbuf];
            Jblock2 = -1;
            Jac = Jvec.Ia_code_[Jblock];
            Jbc = Jvec.Ib_code_[Jblock];
            if (Jvec.Ms0_) Jblock2 = Jvec.decode_[Jbc][Jac];
            Jnas = Jvec.Ia_size_[Jblock];
            Jnbs = Jvec.Ib_size_[Jblock];
            if (s1_contrib_[Iblock2][Jblock] || s2_contrib_[Iblock2][Jblock]) 
              do_Jblock = 1;
            if (Jvec.buf_offdiag_[Jbuf] && (s1_contrib_[Iblock2][Jblock2] ||
                                           s2_contrib_[Iblock2][Jblock2]))
              do_Jblock2 = 1;
            if (!do_Jblock && !do_Jblock2) continue;
	   
            Jvec.read(Jroot, Jbuf);
	 
            if (do_Jblock) {
              opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, Jvec.blocks_[Jblock], 
                         transp_tmp2, Jac, Jbc, Jnas,
                         Jnbs, Iac, Ibc, Inas, Inbs);
            }
	   
            if (do_Jblock2) {
              Jvec.transp_block(Jblock, transp_tmp);
              opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, transp_tmp,
                         transp_tmp2, Jbc, Jac, Jnbs,
                         Jnas, Iac, Ibc, Inas, Inbs);
            }
          } /* end loop over Jbuf */
        } /* end loop over Ibuf transpose */
      } /* end loop over Ibuf */
    } /* end icore==0 */

    else if (Parameters_->icore==1) { /* whole vectors in-core */
      Ivec.read(Iroot, 0);
      Jvec.read(Jroot, 0);
      for (Iblock=0; Iblock<Ivec.num_blocks_; Iblock++) {
        Iac = Ivec.Ia_code_[Iblock];
        Ibc = Ivec.Ib_code_[Iblock];
        Inas = Ivec.Ia_size_[Iblock];
        Inbs = Ivec.Ib_size_[Iblock];
        if (Inas==0 || Inbs==0) continue;
        for (Jblock=0; Jblock<Jvec.num_blocks_; Jblock++) {
          Jac = Jvec.Ia_code_[Jblock];
          Jbc = Jvec.Ib_code_[Jblock];
          Jnas = Jvec.Ia_size_[Jblock];
          Jnbs = Jvec.Ib_size_[Jblock];
          if (s1_contrib_[Iblock][Jblock] || s2_contrib_[Iblock][Jblock])
            opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, Jvec.blocks_[Jblock],
                       Ivec.blocks_[Iblock], Jac, Jbc, Jnas,
                       Jnbs, Iac, Ibc, Inas, Inbs);
        }
      } /* end loop over Iblock */
    } /* end icore==1 */

    else if (Parameters_->icore==2) { /* icore==2 */
      for (Ibuf=0; Ibuf<Ivec.buf_per_vect_; Ibuf++) {
        Ivec.read(Iroot, Ibuf);
        Iairr = Ivec.buf2blk_[Ibuf];

        for (Jbuf=0; Jbuf<Jvec.buf_per_vect_; Jbuf++) {
          Jvec.read(Jroot, Jbuf);
          Jairr = Jvec.buf2blk_[Jbuf];
	
        for (Iblock=Ivec.first_ablk_[Iairr]; Iblock<=Ivec.last_ablk_[Iairr];
             Iblock++) {
          Iac = Ivec.Ia_code_[Iblock];
          Ibc = Ivec.Ib_code_[Iblock];
          Inas = Ivec.Ia_size_[Iblock];
          Inbs = Ivec.Ib_size_[Iblock];
	   
          for (Jblock=Jvec.first_ablk_[Jairr]; Jblock<=Jvec.last_ablk_[Jairr];
               Jblock++) {
            Jac = Jvec.Ia_code_[Jblock];
            Jbc = Jvec.Ib_code_[Jblock];
            Jnas = Jvec.Ia_size_[Jblock];
            Jnbs = Jvec.Ib_size_[Jblock];
	   
            if (s1_contrib_[Iblock][Jblock] || s2_contrib_[Iblock][Jblock])
              opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, Jvec.blocks_[Jblock],
                         Ivec.blocks_[Iblock], Jac, Jbc, Jnas,
                         Jnbs, Iac, Ibc, Inas, Inbs);

            if (Jvec.buf_offdiag_[Jbuf]) {
              Jblock2 = Jvec.decode_[Jbc][Jac];
              if (s1_contrib_[Iblock][Jblock2] ||
                  s2_contrib_[Iblock][Jblock2]) {
              Jvec.transp_block(Jblock, transp_tmp);
                opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, transp_tmp,
                  Ivec.blocks_[Iblock], Jbc, Jac,
                  Jnbs, Jnas, Iac, Ibc, Inas, Inbs);
	      }
	    }

          } /* end loop over Jblock */

          if (Ivec.buf_offdiag_[Ibuf]) {
            Iblock2 = Ivec.decode_[Ibc][Iac];
            Ivec.transp_block(Iblock, transp_tmp2);
            Iac = Ivec.Ia_code_[Iblock2];
            Ibc = Ivec.Ib_code_[Iblock2];
            Inas = Ivec.Ia_size_[Iblock2];
            Inbs = Ivec.Ib_size_[Iblock2];
	   
            for (Jblock=Jvec.first_ablk_[Jairr]; Jblock<=Jvec.last_ablk_[Jairr];
              Jblock++) {
              Jac = Jvec.Ia_code_[Jblock];
              Jbc = Jvec.Ib_code_[Jblock];
              Jnas = Jvec.Ia_size_[Jblock];
              Jnbs = Jvec.Ib_size_[Jblock];
	   
              if (s1_contrib_[Iblock2][Jblock] || s2_contrib_[Iblock2][Jblock])
                opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, Jvec.blocks_[Jblock],
                           transp_tmp2, Jac, Jbc, Jnas, Jnbs, Iac, Ibc, 
                           Inas, Inbs);

                if (Jvec.buf_offdiag_[Jbuf]) {
                  Jblock2 = Jvec.decode_[Jbc][Jac];
                  if (s1_contrib_[Iblock][Jblock2] || 
                    s2_contrib_[Iblock][Jblock2]) {
                    Jvec.transp_block(Jblock, transp_tmp);
                    opdm_block(alplist_, betlist_, scratch_ap, scratch_bp, transp_tmp,
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
      throw PSIEXCEPTION("CIWavefunction::opdm: unrecognized core option!\n");
    }

    std::stringstream opdm_name;
    opdm_name << "MO-basis Alpha OPDM <" << Iroot+1 << "| Etu |" << Jroot+1 << ">"; 
    SharedMatrix new_OPDM_a(new Matrix(opdm_name.str(), nirrep_, act_dim, act_dim));

    opdm_name.str(std::string());
    opdm_name << "MO-basis Beta OPDM <" << Iroot+1 << "| Etu |" << Jroot+1 << ">"; 
    SharedMatrix new_OPDM_b(new Matrix(opdm_name.str(), nirrep_, act_dim, act_dim));

    opdm_name.str(std::string());
    opdm_name << "MO-basis OPDM <" << Iroot+1 << "| Etu |" << Jroot+1 << ">"; 
    SharedMatrix new_OPDM(new Matrix(opdm_name.str(), nirrep_, act_dim, act_dim));

    int offset = 0;
    for (int h=0; h<nirrep_; h++){
      if (!CalcInfo_->ci_orbs[h]) continue;

      double* opdm_a = new_OPDM_a->pointer(h)[0];
      double* opdm_b = new_OPDM_b->pointer(h)[0];
      double* opdm  = new_OPDM->pointer(h)[0];

      for (int i=0, target=0; i<CalcInfo_->ci_orbs[h]; i++){
        int ni = CalcInfo_->act_reorder[i + offset];
        for (int j=0; j<CalcInfo_->ci_orbs[h]; j++){
          int nj = CalcInfo_->act_reorder[j + offset];

          opdm_a[target] = scratch_ap[ni][nj];
          opdm_b[target] = scratch_bp[ni][nj];
          opdm[target++] = scratch_ap[ni][nj] + scratch_bp[ni][nj]; 
        }
      }
      offset += CalcInfo_->ci_orbs[h];
    }

    //new_OPDM_a->print();
    //new_OPDM->print();

    std::vector<SharedMatrix> opdm_root_vec;
    opdm_root_vec.push_back(new_OPDM_a);
    opdm_root_vec.push_back(new_OPDM_b);
    opdm_root_vec.push_back(new_OPDM);

    opdm_list.push_back(opdm_root_vec); 

    
    if (!transden) Iroot++;
  } /* end loop over num_roots Jroot */  


  if (transp_tmp != NULL) free_block(transp_tmp);
  if (transp_tmp2 != NULL) free_block(transp_tmp2);
  Ivec.buf_unlock();
  Jvec.buf_unlock();
  free(buffer1);
  free(buffer2);

  scratch_a.reset();
  scratch_b.reset();

  opdm_called_ = true;
  return opdm_list;
}

void CIWavefunction::opdm_block(struct stringwr **alplist, struct stringwr **betlist,
		double **onepdm_a, double **onepdm_b, double **CJ, double **CI, int Ja_list, 
		int Jb_list, int Jnas, int Jnbs, int Ia_list, int Ib_list, 
		int Inas, int Inbs)
{
  int Ia_idx, Ib_idx, Ja_idx, Jb_idx, Ja_ex, Jb_ex, Jbcnt, Jacnt; 
  struct stringwr *Jb, *Ja;
  signed char *Jbsgn, *Jasgn;
  unsigned int *Jbridx, *Jaridx;
  double C1, C2, Ib_sgn, Ia_sgn;
  int i, j, oij, ndrc, *Jboij, *Jaoij;
 
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
          i = oij/CalcInfo_->num_ci_orbs;
          j = oij%CalcInfo_->num_ci_orbs;
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
          i = oij/CalcInfo_->num_ci_orbs;
          j = oij%CalcInfo_->num_ci_orbs;
	  onepdm_a[i][j] += C1 * C2 * Ia_sgn;
	}
      }
    }
  }
}

/*
** Here we compute all properties for opdms in opdm_map_ from form_opdm 
*/
void CIWavefunction::opdm_properties()
{
    if (!opdm_called_){
        throw PSIEXCEPTION("CIWavefunction::opdm_properties: Called before form_opdm, no opdms to run properties on!"); 
    }
    outfile->Printf("Figuring out root list\n");

    // Figure out which opdms we needs
    std::vector<std::pair<int, int> > root_list;
    for (int i=0; i<Parameters_->num_roots; i++){
        root_list.push_back(std::pair<int, int>(i, i));
    }
    if (Parameters_->transdens) {
        for (int i=1; i<Parameters_->num_roots; i++){
            root_list.push_back(std::pair<int, int>(0, i));
        }
    }

    boost::shared_ptr<OEProp> oe(new OEProp());
    oe->set_Ca(get_orbitals("ALL"));
    SharedMatrix opdm_a;
    SharedMatrix opdm_b;
    std::stringstream opdm_name, oeprop_label;
    double inactive_value = 0.0;

    std::stringstream ss;

    // Loop over roots
    for (int i=0; i<root_list.size(); i++){
        int Iroot = root_list[i].first;
        int Jroot = root_list[i].second;

        if (Iroot == Jroot) inactive_value = 1.0;
        else inactive_value = 0.0;

        opdm_name.str(std::string());
        opdm_name << "MO-basis Alpha OPDM <" << Iroot+1 << "| Etu |" << Jroot+1 << ">"; 
        opdm_a = opdm_add_inactive(opdm_map_[opdm_name.str()], inactive_value, true);
        oe->set_Da_mo(opdm_a);

        if (Parameters_->ref == "ROHF") {
            opdm_name.str(std::string());
            opdm_name << "MO-basis Beta OPDM <" << Iroot+1 << "| Etu |" << Jroot+1 << ">"; 
            opdm_b = opdm_add_inactive(opdm_map_[opdm_name.str()], inactive_value, true);
            oe->set_Db_mo(opdm_b);
        }

        oe->clear();
        oeprop_label.str(std::string());
        if (Iroot == Jroot) {
            oe->add("DIPOLE");
            oe->add("MULLIKEN_CHARGES");
            oe->add("NO_OCCUPATIONS");
            oe->add("QUADRUPOLE");
            oeprop_label << "CI ROOT " << (Iroot+1); 
        } 
        else {
            oe->add("TRANSITION_DIPOLE");
            oe->add("TRANSITION_QUADRUPOLE");
            oeprop_label << "CI ROOT " << (Iroot+1) << " -> ROOT " << (Jroot+1); 
        }
        oe->set_title(oeprop_label.str());
        
        
        outfile->Printf( "\n  ==> Properties %s <==\n", oeprop_label.str().c_str());
        oe->compute();

        // if this is the "special" root, then copy over OEProp 
        // Process::environment variables from the current root into
        // more general locations

        if ((Iroot == Parameters_->root) && (Jroot == Parameters_->root)) {
          std::stringstream ss2;
          ss2 << oeprop_label.str() << " DIPOLE X"; 
          Process::environment.globals["CI DIPOLE X"] = Process::environment.globals[ss2.str()]; 

          ss2.str(std::string());
          ss2 << oeprop_label.str() << " DIPOLE Y"; 
          Process::environment.globals["CI DIPOLE Y"] = Process::environment.globals[ss2.str()]; 

          ss2.str(std::string());
          ss2 << oeprop_label.str() << " DIPOLE Z"; 
          Process::environment.globals["CI DIPOLE Z"] = Process::environment.globals[ss2.str()]; 

          ss2.str(std::string());
          ss2 << oeprop_label.str() << " QUADRUPOLE XX"; 
          Process::environment.globals["CI QUADRUPOLE XX"] = Process::environment.globals[ss2.str()]; 

          ss2.str(std::string());
          ss2 << oeprop_label.str() << " QUADRUPOLE YY"; 
          Process::environment.globals["CI QUADRUPOLE YY"] = Process::environment.globals[ss2.str()]; 

          ss2.str(std::string());
          ss2 << oeprop_label.str() << " QUADRUPOLE ZZ"; 
          Process::environment.globals["CI QUADRUPOLE ZZ"] = Process::environment.globals[ss2.str()]; 

          ss2.str(std::string());
          ss2 << oeprop_label.str() << " QUADRUPOLE XY"; 
          Process::environment.globals["CI QUADRUPOLE XY"] = Process::environment.globals[ss2.str()]; 

          ss2.str(std::string());
          ss2 << oeprop_label.str() << " QUADRUPOLE XZ"; 
          Process::environment.globals["CI QUADRUPOLE XZ"] = Process::environment.globals[ss2.str()]; 

          ss2.str(std::string());
          ss2 << oeprop_label.str() << " QUADRUPOLE YZ"; 
          Process::environment.globals["CI QUADRUPOLE YZ"] = Process::environment.globals[ss2.str()]; 
        }

    } // End loop over roots
}

/*
void CIWavefunction::ci_nat_orbs()
{
    double *opdm_eigval, **opdm_eigvec, **opdm_blk, **scfvec;
    for (int irrep=0, max_orb_per_irrep=0; irrep<CalcInfo_->nirreps; irrep++) {
      if (CalcInfo_->orbs_per_irr[irrep] > max_orb_per_irrep)
        max_orb_per_irrep = CalcInfo_->so_per_irr[irrep];
    }

    opdm_eigvec = block_matrix(max_orb_per_irrep, max_orb_per_irrep);
    opdm_eigval = init_array(max_orb_per_irrep);
    opdm_blk = block_matrix(max_orb_per_irrep, max_orb_per_irrep);

    psio_open(targetfile, PSIO_OPEN_OLD);
    chkpt_init(PSIO_OPEN_OLD);

    // reorder opdm from ci to pitzer and diagonalize each 
    if (Parameters_->opdm_ave) {
      klast = 1;
    }
    else {
      klast = Parameters_->num_roots;
    }

    // loop over roots or averaged opdm
    for(k=0; k<klast; k++) {
      if (Parameters_->opdm_ave && Parameters_->print_lvl > 1) {
        outfile->Printf("\n\n\t\t\tCI Natural Orbitals for the Averaged\n");
        outfile->Printf("\t\t\tOPDM of %d Roots in terms of Molecular"
                 " Orbitals\n\n",k); 
      }
      else if (Parameters_->print_lvl > 1) {
        outfile->Printf(
             "\n\t\t\tCI Natural Orbitals in terms of Molecular Orbitals\n\n");
        outfile->Printf("\t\t\t Root %d\n\n",k+1);
        
      }

      mo_offset = 0;

      if (!Parameters_->opdm_ave)
        sprintf(opdm_key, "MO-basis OPDM Root %d", k);
      else 
        sprintf(opdm_key, "MO-basis OPDM Ave");

      psio_read_entry(targetfile, opdm_key, (char *) onepdm[0],
        populated_orbs * populated_orbs * sizeof(double));

      for (irrep=0; irrep<CalcInfo_->nirreps; irrep++) {
        if (CalcInfo_->orbs_per_irr[irrep] == 0) continue; 
        for (i=0; i<CalcInfo_->orbs_per_irr[irrep]-
                    CalcInfo_->dropped_uocc[irrep]; i++) {
          for (j=0; j<CalcInfo_->orbs_per_irr[irrep]-
                    CalcInfo_->dropped_uocc[irrep]; j++) {
            i_ci = CalcInfo_->reorder[i+mo_offset];
            j_ci = CalcInfo_->reorder[j+mo_offset]; 
            opdm_blk[i][j] = onepdm[i_ci][j_ci];
          } 
        }
 
        // Writing SCF vector to OPDM file for safe keeping 
        // because you might overwrite it with NO's for earlier root 
        // and we're only storing one irrep in core.  Only need to do once.
        
        if (k==0) {

          scfvec = chkpt_rd_scf_irrep(irrep);

            if (Parameters_->print_lvl > 3) {
              outfile->Printf("Cvec for k==0, read in from chkpt original\n");
              outfile->Printf(" %s Block \n", CalcInfo_->labels[irrep]);
              print_mat(scfvec, CalcInfo_->orbs_per_irr[irrep],
                        CalcInfo_->orbs_per_irr[irrep], "outfile");
            }

          sprintf(opdm_key, "Old SCF Matrix Irrep %d", irrep);
          psio_write_entry(targetfile, opdm_key, (char *) scfvec[0],
            CalcInfo_->orbs_per_irr[irrep] * CalcInfo_->orbs_per_irr[irrep] *
            sizeof(double));
	      free_block(scfvec);
        }

        zero_mat(opdm_eigvec, max_orb_per_irrep, max_orb_per_irrep);

        if (CalcInfo_->orbs_per_irr[irrep]-CalcInfo_->dropped_uocc[irrep] > 0) {
          sq_rsp(CalcInfo_->orbs_per_irr[irrep]-CalcInfo_->dropped_uocc[irrep],
                 CalcInfo_->orbs_per_irr[irrep]-CalcInfo_->dropped_uocc[irrep],
                 opdm_blk, opdm_eigval, 1, opdm_eigvec, TOL); 
          }
        for (i=CalcInfo_->orbs_per_irr[irrep]-CalcInfo_->dropped_uocc[irrep]; 
             i<CalcInfo_->orbs_per_irr[irrep]; i++) {
           opdm_eigvec[i][i] = 1.0;
           opdm_eigval[i] = 0.0;
           }
        eigsort(opdm_eigval, opdm_eigvec, -(CalcInfo_->orbs_per_irr[irrep]));


        // FAE
        // eigsort sometimes will swap the order of orbitals, for example
        // in frozen core computations the focc may be mixed with docc and 
        // change the final result
        
        //Loop over "populated"
          for (i=0;i<CalcInfo_->orbs_per_irr[irrep]-CalcInfo_->dropped_uocc[irrep];i++){
            max_overlap = 0;
            int m       = 0;
             for (j=i;j<CalcInfo_->orbs_per_irr[irrep]-CalcInfo_->dropped_uocc[irrep];j++){
               overlap = opdm_eigvec[i][j] * opdm_eigvec[i][j];
               if(overlap > max_overlap){
                 m = j;
                 max_overlap = overlap;
               }
             }
             for (j=0;j<CalcInfo_->orbs_per_irr[irrep];j++){
                 double temporary  = opdm_eigvec[j][i];
                 opdm_eigvec[j][i] = opdm_eigvec[j][m];
                 opdm_eigvec[j][m] = temporary;
             }
             double temporary = opdm_eigval[i];
             opdm_eigval[i] = opdm_eigval[m];
             opdm_eigval[m] = temporary;
        
          }
          // End FAE changes, May 3 2007


        if (Parameters_->print_lvl > 0) {
          if (irrep==0) {
            if (Parameters_->opdm_ave) { 
              outfile->Printf( "\n Averaged CI Natural Orbitals in terms "
                "of Molecular Orbitals\n\n");
              }
            else outfile->Printf( "\n CI Natural Orbitals in terms of "
                   "Molecular Orbitals: Root %d\n\n", k+1);
          }
          outfile->Printf("\n %s Block (MO basis)\n", CalcInfo_->labels[irrep]);
          eivout(opdm_eigvec, opdm_eigval, CalcInfo_->orbs_per_irr[irrep],
                 CalcInfo_->orbs_per_irr[irrep], "outfile");
        }

        if (Parameters_->opdm_wrtnos && (k==Parameters_->opdm_orbs_root)) {
          if (irrep==0) {
            if (!Parameters_->opdm_ave) {
              outfile->Printf("\n Writing CI Natural Orbitals for root %d"
                      " to checkpoint in terms of Symmetry Orbitals\n\n",k+1);
            }
          }
	  
	      scfvec = block_matrix(CalcInfo_->orbs_per_irr[irrep], 
                                CalcInfo_->orbs_per_irr[irrep]);
          sprintf(opdm_key, "Old SCF Matrix Irrep %d", irrep);
          psio_read_entry(targetfile, opdm_key, (char *) scfvec[0],
            CalcInfo_->orbs_per_irr[irrep] * CalcInfo_->orbs_per_irr[irrep] *
            sizeof(double));

          if (Parameters_->print_lvl > 3) {
            outfile->Printf("\nCvec read for MO to SO trans\n\n");
            outfile->Printf(" %s Block \n", CalcInfo_->labels[irrep]);
            print_mat(scfvec, CalcInfo_->orbs_per_irr[irrep],
                      CalcInfo_->orbs_per_irr[irrep], "outfile");
            outfile->Printf("\nOpdm_eigvec before MO to SO trans\n\n");
            outfile->Printf(" %s Block \n", CalcInfo_->labels[irrep]);
            print_mat(opdm_eigvec, CalcInfo_->orbs_per_irr[irrep],
                    CalcInfo_->orbs_per_irr[irrep], "outfile");
          }
          mmult(scfvec, 0, opdm_eigvec, 0, opdm_blk, 0,
                CalcInfo_->so_per_irr[irrep], CalcInfo_->orbs_per_irr[irrep],
                CalcInfo_->orbs_per_irr[irrep], 0); 
          free_block(scfvec);
          if (Parameters_->print_lvl > 0) {  // FAE
            outfile->Printf(" %s Block (SO basis)\n", CalcInfo_->labels[irrep]);
            print_mat(opdm_blk, CalcInfo_->so_per_irr[irrep],
                      CalcInfo_->orbs_per_irr[irrep], "outfile");
          }
          chkpt_wt_scf_irrep(opdm_blk, irrep);
          outfile->Printf( "\n Warning: Natural Orbitals for the ");
	  if (Parameters_->opdm_ave)
            outfile->Printf( "Averaged OPDM ");
          else
            outfile->Printf( "Root %d OPDM ", k);
          outfile->Printf( "have been written to the checkpoint file!\n\n"); 
        } // end code to write the NO's to disk
        mo_offset += CalcInfo_->orbs_per_irr[irrep];
      } // end loop over irreps
    } // end loop over roots 

    free_block(opdm_eigvec);
    free(opdm_eigval); 
    chkpt_close();
    psio_close(targetfile, 1);
}
*/


/*
** opdm_ave()
** 
** Parameters:
**  targetfile = file number to obtain matrices from 
**
** Modified 7/13/04 to use new state average parameters --- CDS
**
*/
//void CIWavefunction::opdm_ave(int targetfile)
//{
//  int root, root_count, i, j, populated_orbs;
//  double **tmp_mat1, **tmp_mat2;
//  char opdm_key[80];
//  
//  const char spinlabels[][10] = { "", "Alpha ", "Beta " };
//  const int nspincases = 3;
//
//  populated_orbs = CalcInfo_->nmo-CalcInfo_->num_drv_orbs;
//  tmp_mat1 = block_matrix(populated_orbs, populated_orbs);  
//  tmp_mat2 = block_matrix(populated_orbs, populated_orbs);
//
//  
//  for(int spincase=0; spincase<nspincases; ++spincase) {
//    
//    zero_mat(tmp_mat1, populated_orbs, populated_orbs);
//
//    for(root_count=0; root_count<Parameters_->average_num; root_count++) {
//
//      root = Parameters_->average_states[root_count];
//
//      sprintf(opdm_key, "MO-basis %sOPDM Root %d", spinlabels[spincase], root);
//
//      psio_read_entry(targetfile, opdm_key, (char *) tmp_mat2[0],
//                      populated_orbs * populated_orbs * sizeof(double));
//
//      if (Parameters_->opdm_print) {
//        outfile->Printf("\n\n\t\t%sOPDM for Root %d",spinlabels[spincase],root+1);
//        print_mat(tmp_mat2, populated_orbs, populated_orbs, "outfile");
//      }
//
//      for (i=0; i<populated_orbs; i++)
//        for (j=0; j<populated_orbs; j++)
//          tmp_mat1[i][j] += Parameters_->average_weights[root]*tmp_mat2[i][j];
//
//    }
//
//    sprintf(opdm_key, "MO-basis %sOPDM", spinlabels[spincase]);
//    psio_write_entry(targetfile, opdm_key, (char *) tmp_mat1[0],
//                     populated_orbs*populated_orbs*sizeof(double));
//    
//    if (Parameters_->print_lvl> 0 || Parameters_->opdm_print) {
//      outfile->Printf(
//              "\n\t\t\t Averaged %sOPDM's for %d Roots written to opdm_file \n\n",
//              spinlabels[spincase], Parameters_->average_num);
//    }
//    if (Parameters_->opdm_print) {
//      print_mat(tmp_mat1, populated_orbs, populated_orbs, "outfile");
//    }
//    
//  } // loop over spincases
//
//  free_block(tmp_mat1);
//  free_block(tmp_mat2);
//
//}


/*
** opdm_ke
**
** Compute the kinetic energy contribution from the correlated part of the
** one-particle density matrix.  For Daniel Crawford
*/
void CIWavefunction::opdm_ke(double **onepdm)
{
  int errcod;
  int src_T_file, mo_offset, so_offset, irrep, i, j, i_ci, j_ci, i2, j2, ij;
  int maxorbs;
  int noeints;
  double *T, **scfmat, **opdm_blk, **tmp_mat, ke, kei;

  ke = kei = 0.0;

  /* read in the kinetic energy integrals */
  noeints = CalcInfo_->nso*(CalcInfo_->nso+1)/2;
  src_T_file = PSIF_OEI;

  T = init_array(noeints);
  if (Parameters_->print_lvl>2) 
    outfile->Printf( "Kinetic energy integrals (SO basis):\n");
  errcod = iwl_rdone(src_T_file,PSIF_SO_T,T,noeints,0,
                     (Parameters_->print_lvl>2),"outfile");
  if (!errcod) {
    throw PsiException("(detci): Error reading kinetic energy ints",__FILE__,__LINE__);
  }

  /* find biggest blocksize */
  for (irrep=0,maxorbs=0; irrep<CalcInfo_->nirreps; irrep++) {
    if (CalcInfo_->orbs_per_irr[irrep] > maxorbs)
      maxorbs = CalcInfo_->so_per_irr[irrep];
  }
  opdm_blk = block_matrix(maxorbs, maxorbs);
  tmp_mat = block_matrix(maxorbs, maxorbs);

  /* transform the onepdm into SO form, one irrep at a time */
  so_offset = mo_offset = 0;
  outfile->Printf("Correlation Piece of OPDM in SO basis\n");
  for (irrep=0; irrep<CalcInfo_->nirreps; irrep++) {
    if (CalcInfo_->orbs_per_irr[irrep] == 0) continue;
    for (i=0; i<CalcInfo_->orbs_per_irr[irrep]; i++) {
      for (j=0; j<CalcInfo_->orbs_per_irr[irrep]; j++) {
        i_ci = CalcInfo_->reorder[i+mo_offset];
        j_ci = CalcInfo_->reorder[j+mo_offset];
        opdm_blk[i][j] = onepdm[i_ci][j_ci];
      }
    }
    /* need to subtract out reference piece, assume single det ref */
    for (i=0; i<CalcInfo_->docc[irrep]; i++) {
      opdm_blk[i][i] -= 2.0;
    }
    /* keep counting from i to take out socc part */
    for (j=0; j<CalcInfo_->socc[irrep]; j++,i++) {
      opdm_blk[i][i] -= 1.0;
    }

    
    outfile->Printf( "Irrep %d\n", irrep);
    outfile->Printf( "MO basis, Pitzer order\n");
    print_mat(opdm_blk,CalcInfo_->so_per_irr[irrep],
              CalcInfo_->so_per_irr[irrep],"outfile");

    /* transform back to SO basis */
    scfmat = chkpt_rd_scf_irrep(irrep);
    mmult(opdm_blk,0,scfmat,1,tmp_mat,0,CalcInfo_->orbs_per_irr[irrep],
          CalcInfo_->orbs_per_irr[irrep],CalcInfo_->so_per_irr[irrep],0);
    mmult(scfmat,0,tmp_mat,0,opdm_blk,0,CalcInfo_->so_per_irr[irrep],
          CalcInfo_->orbs_per_irr[irrep],CalcInfo_->so_per_irr[irrep],0); 
    
    outfile->Printf( "SO basis, Pitzer order\n");
    print_mat(opdm_blk,CalcInfo_->so_per_irr[irrep],
              CalcInfo_->so_per_irr[irrep],"outfile");

    /* get kinetic energy contribution */
    kei = 0.0;
    for (i=0,i2=so_offset; i<CalcInfo_->so_per_irr[irrep]; i++,i2++) {
      for (j=0,j2=so_offset; j<CalcInfo_->so_per_irr[irrep]; j++,j2++) {
        ij = ioff[MAX0(i2,j2)] + MIN0(i2,j2);
        kei += opdm_blk[i][j] * T[ij];
      }
    }

    outfile->Printf("Contribution to correlation kinetic energy = %15.10lf\n", 
            kei);

    ke += kei;
    mo_offset += CalcInfo_->orbs_per_irr[irrep];
    so_offset += CalcInfo_->so_per_irr[irrep];

  } /* end loop over irreps */

  outfile->Printf( "\nTotal correlation kinetic energy = %15.10lf\n", ke);

}


}} // namespace psi::detci
