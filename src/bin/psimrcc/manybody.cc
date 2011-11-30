/**
 *  @file manybody.cc
 *  @ingroup (PSIMRCC)
 *  @brief The base class for all the many-body methods
*/
#include <iostream>

#include <algorithm>
#include <cmath>
#include <functional>
#include <utility>
#include <vector>

#include <liboptions/liboptions.h>
#include <libutil/libutil.h>
#include <libmoinfo/libmoinfo.h>

#include "algebra_interface.h"
#include "blas.h"
#include "index.h"
#include "debugging.h"
#include "manybody.h"
#include "matrix.h"
#include "sort.h"

namespace psi{
    extern FILE *outfile;
    namespace psimrcc{
    extern MOInfo *moinfo;
    extern MemoryManager* memory_manager;


using namespace std;

/**
 * Allocate the effective Hamiltonian matrices and eigenvectors
 * @todo wrap the current operations in an init() function
 */
CCManyBody::CCManyBody(Options &options):
        options_(options)
{
  // Allocate memory for the eigenvector and the effective Hamiltonian
  allocate1(double,zeroth_order_eigenvector,moinfo->get_nrefs());
  allocate1(double,right_eigenvector,moinfo->get_nrefs());
  allocate1(double,left_eigenvector,moinfo->get_nrefs());
  allocate2(double,Heff,moinfo->get_nrefs(),moinfo->get_nrefs());
  allocate2(double,Heff_mrpt2,moinfo->get_nrefs(),moinfo->get_nrefs());

  pert_cbs = false;
  pert_cbs_coupling = false;
  huge  = 1.0e100;
  norm_amps = 0.0;
  delta_t1_amps = 0.0;
  delta_t2_amps = 0.0;
  d3_ooo = d3_ooO = d3_oOO = d3_OOO =  d3_vvv = d3_vvV = d3_vVV = d3_VVV = NULL;
}

/**
 * Deallocate the effective Hamiltonian matrices and eigenvectors
 * @todo wrap the current operations in an cleanup() function
 */
CCManyBody::~CCManyBody()
{
  release1(zeroth_order_eigenvector);
  release1(left_eigenvector);
  release1(right_eigenvector);
  release2(Heff);
  release2(Heff_mrpt2);

  if(triples_type>ccsd)
    deallocate_triples_denominators();
}

/**
 * Creates a CCSort object and stores the address in the global pointer sorter
 */
void CCManyBody::generate_integrals()
{
  Timer timer;
  DEBUGGING(1,
    fprintf(outfile,"\n\tvoid CCManyBody::generate_integrals()");
    fflush(outfile);
  )
  // CCSort reads the one and two electron integrals
  // and creates the Fock matrices
  sorter   = new CCSort(out_of_core_sort);
//   blas->show_storage();
  blas->compute_storage_strategy();
//   blas->show_storage();

  DEBUGGING(1,
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  )
}

void CCManyBody::generate_triples_denominators()
{
  generate_d3_ijk(d3_ooo,true,true,true);
  generate_d3_ijk(d3_ooO,true,true,false);
  generate_d3_ijk(d3_oOO,true,false,false);
  generate_d3_ijk(d3_OOO,false,false,false);
  generate_d3_abc(d3_vvv,true,true,true);
  generate_d3_abc(d3_vvV,true,true,false);
  generate_d3_abc(d3_vVV,true,false,false);
  generate_d3_abc(d3_VVV,false,false,false);
}

void CCManyBody::generate_d3_ijk(double***& d3,bool alpha_i,bool alpha_j,bool alpha_k)
{
  allocate2(double*,d3,moinfo->get_nunique(),moinfo->get_nirreps());
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int reference = moinfo->get_ref_number(ref,UniqueRefs);

    // N.B. Never introduce Matrices/Vectors with O or V in the name before you compute the Fock matrix elements
    std::vector<int> aocc     = moinfo->get_aocc(reference,AllRefs);
    std::vector<int> bocc     = moinfo->get_bocc(reference,AllRefs);

    // Build the is_ arrays for reference ref
    bool* is_aocc = new bool[moinfo->get_nocc()];
    bool* is_bocc = new bool[moinfo->get_nocc()];
    for(int i=0;i<moinfo->get_nocc();i++){
      is_aocc[i]=false;
      is_bocc[i]=false;
    }
    for(size_t i=0;i<aocc.size();i++) is_aocc[aocc[i]]=true;
    for(size_t i=0;i<bocc.size();i++) is_bocc[bocc[i]]=true;

    // Read the Fock matrices
    CCMatTmp f_oo_Matrix = blas->get_MatTmp("fock[oo]",reference,none);
    CCMatTmp f_OO_Matrix = blas->get_MatTmp("fock[OO]",reference,none);

    CCMatrix* f_ii_Matrix;
    CCMatrix* f_jj_Matrix;
    CCMatrix* f_kk_Matrix;

    if(alpha_i)
      f_ii_Matrix = f_oo_Matrix.get_CCMatrix();
    else
      f_ii_Matrix = f_OO_Matrix.get_CCMatrix();

    if(alpha_j)
      f_jj_Matrix = f_oo_Matrix.get_CCMatrix();
    else
      f_jj_Matrix = f_OO_Matrix.get_CCMatrix();

    if(alpha_k)
      f_kk_Matrix = f_oo_Matrix.get_CCMatrix();
    else
      f_kk_Matrix = f_OO_Matrix.get_CCMatrix();

    CCIndex* ooo_indexing = blas->get_index("[ooo]");
    short**  ooo_tuples   = ooo_indexing->get_tuples();

    for(int h=0;h<moinfo->get_nirreps();h++){
      size_t ooo_offset  = ooo_indexing->get_first(h);
      allocate1(double,d3[ref][h],ooo_indexing->get_pairpi(h));
      for(size_t ijk = 0;ijk<ooo_indexing->get_pairpi(h);ijk++){
        short i = ooo_tuples[ooo_offset + ijk][0];
        short j = ooo_tuples[ooo_offset + ijk][1];
        short k = ooo_tuples[ooo_offset + ijk][2];

        bool external = true;
        if((alpha_i && !is_aocc[i]) || (!alpha_i && !is_bocc[i]))
          external = false;
        if((alpha_j && !is_aocc[j]) || (!alpha_j && !is_bocc[j]))
          external = false;
        if((alpha_k && !is_aocc[k]) || (!alpha_k && !is_bocc[k]))
          external = false;

        if(external)
          d3[ref][h][ijk]=f_ii_Matrix->get_two_address_element(i,i)
                    +f_jj_Matrix->get_two_address_element(j,j)
                   +f_kk_Matrix->get_two_address_element(k,k);
        else
          d3[ref][h][ijk]=huge;
      } // End loop over h,ijk
    }
    delete[] is_aocc;
    delete[] is_bocc;
  }
}

void CCManyBody::generate_d3_abc(double***& d3,bool alpha_a,bool alpha_b,bool alpha_c)
{
  allocate2(double*,d3,moinfo->get_nunique(),moinfo->get_nirreps());
  // Loop over references
  for(int ref=0;ref<moinfo->get_nunique();ref++){
    int reference = moinfo->get_ref_number(ref,UniqueRefs);

    // N.B. Never introduce Matrices/Vectors with O or V in the name before you compute the Fock matrix elements
    std::vector<int> avir     = moinfo->get_avir(reference,AllRefs);
    std::vector<int> bvir     = moinfo->get_bvir(reference,AllRefs);

    // Build the is_ arrays for reference ref
    bool* is_avir = new bool[moinfo->get_nvir()];
    bool* is_bvir = new bool[moinfo->get_nvir()];
    for(int i=0;i<moinfo->get_nvir();i++){
      is_avir[i]=false;
      is_bvir[i]=false;
    }
    for(size_t i=0;i<avir.size();i++) is_avir[avir[i]]=true;
    for(size_t i=0;i<bvir.size();i++) is_bvir[bvir[i]]=true;

    // Read the Fock matrices
    CCMatTmp f_vv_Matrix = blas->get_MatTmp("fock[vv]",reference,none);
    CCMatTmp f_VV_Matrix = blas->get_MatTmp("fock[VV]",reference,none);

    CCMatrix* f_aa_Matrix;
    CCMatrix* f_bb_Matrix;
    CCMatrix* f_cc_Matrix;


    if(alpha_a)
      f_aa_Matrix = f_vv_Matrix.get_CCMatrix();
    else
      f_aa_Matrix = f_VV_Matrix.get_CCMatrix();

    if(alpha_b)
      f_bb_Matrix = f_vv_Matrix.get_CCMatrix();
    else
      f_bb_Matrix = f_VV_Matrix.get_CCMatrix();

    if(alpha_c)
      f_cc_Matrix = f_vv_Matrix.get_CCMatrix();
    else
      f_cc_Matrix = f_VV_Matrix.get_CCMatrix();

    CCIndex* vvv_indexing = blas->get_index("[vvv]");
    short**  vvv_tuples   = vvv_indexing->get_tuples();

    for(int h=0;h<moinfo->get_nirreps();h++){
      size_t vvv_offset  = vvv_indexing->get_first(h);
      allocate1(double,d3[ref][h],vvv_indexing->get_pairpi(h));
      for(size_t abc = 0;abc<vvv_indexing->get_pairpi(h);abc++){
        short a = vvv_tuples[vvv_offset + abc][0];
        short b = vvv_tuples[vvv_offset + abc][1];
        short c = vvv_tuples[vvv_offset + abc][2];

        bool external = true;
        if((alpha_a && !is_avir[a]) || (!alpha_a && !is_bvir[a]))
          external = false;
        if((alpha_b && !is_avir[b]) || (!alpha_b && !is_bvir[b]))
          external = false;
        if((alpha_c && !is_avir[c]) || (!alpha_c && !is_bvir[c]))
          external = false;

        if(external)
          d3[ref][h][abc]=f_aa_Matrix->get_two_address_element(a,a)
                     +f_bb_Matrix->get_two_address_element(b,b)
                     +f_cc_Matrix->get_two_address_element(c,c);
        else
          d3[ref][h][abc]=-huge;
      } // End lvvp over h,ijk
    }
    delete[] is_avir;
    delete[] is_bvir;
  }
}

void CCManyBody::deallocate_triples_denominators()
{
  for(int ref=0;ref<moinfo->get_nunique();ref++)
    for(int h=0;h<moinfo->get_nirreps();h++){
      release1(d3_ooo[ref][h]);
      release1(d3_ooO[ref][h]);
      release1(d3_oOO[ref][h]);
      release1(d3_OOO[ref][h]);

      release1(d3_vvv[ref][h]);
      release1(d3_vvV[ref][h]);
      release1(d3_vVV[ref][h]);
      release1(d3_VVV[ref][h]);
    }

  release2(d3_ooo);
  release2(d3_ooO);
  release2(d3_oOO);
  release2(d3_OOO);
  release2(d3_vvv);
  release2(d3_vvV);
  release2(d3_vVV);
  release2(d3_VVV);
}

/**
 * Computes the energy for each unique reference determinant
 */
void CCManyBody::compute_reference_energy()
{
  Timer timer;
  DEBUGGING(3,
    fprintf(outfile,"\n\tvoid CCManyBody::compute_reference_energy()");
    fflush(outfile);
  )

  // Compute the zeroth-order energy for the unique references
  for(int n=0;n<moinfo->get_nunique();n++){
    int unique_n = moinfo->get_ref_number(n,UniqueRefs);
    double ref_energy=moinfo->get_nuclear_energy()+moinfo->get_fzcore_energy();
    // Grab reference n and the list of occupied orbitals
    std::vector<int> aocc     = moinfo->get_aocc(n,UniqueRefs);
    std::vector<int> bocc     = moinfo->get_bocc(n,UniqueRefs);

    // Read these matrices
    CCMatTmp f_oo_Matrix = blas->get_MatTmp("fock[o][o]",unique_n,none);
    CCMatTmp f_OO_Matrix = blas->get_MatTmp("fock[O][O]",unique_n,none);
    CCMatTmp V_oooo_Matrix = blas->get_MatTmp("<[oo]:[oo]>",none);
    CCMatTmp V_oOoO_Matrix = blas->get_MatTmp("<[oo]|[oo]>",none);

    for(size_t i=0;i<aocc.size();i++)
      ref_energy+=f_oo_Matrix->get_two_address_element(aocc[i],aocc[i]);
    for(size_t i=0;i<bocc.size();i++)
      ref_energy+=f_OO_Matrix->get_two_address_element(bocc[i],bocc[i]);

    for(size_t i=0;i<aocc.size();i++)
      for(size_t j=0;j<aocc.size();j++)
        ref_energy -= 0.5 * V_oooo_Matrix->get_four_address_element(aocc[i],aocc[j],aocc[i],aocc[j]);
    for(size_t i=0;i<bocc.size();i++)
      for(size_t j=0;j<bocc.size();j++)
        ref_energy -= 0.5 * V_oooo_Matrix->get_four_address_element(bocc[i],bocc[j],bocc[i],bocc[j]);
    for(size_t i=0;i<aocc.size();i++)
      for(size_t j=0;j<bocc.size();j++)
        ref_energy -= V_oOoO_Matrix->get_four_address_element(aocc[i],bocc[j],aocc[i],bocc[j]);
    // Write the energy to the ERef
    CCMatTmp ERef_Matrix = blas->get_MatTmp("ERef",unique_n,none);
    ERef_Matrix->set_scalar(ref_energy);
  }


  DEBUGGING(3,
    blas->print("ERef{u}");
    fprintf(outfile," done. Timing %20.6f s",timer.get());
    fflush(outfile);
  )
}

void CCManyBody::print_method(const char* text)
{
  fprintf(outfile,"\n");
  fprintf(outfile,"\n  ==============================================================================");
  fprintf(outfile,"\n  %s",text);
  fprintf(outfile,"\n  ==============================================================================");
  fprintf(outfile,"\n");
  fflush(outfile);
}

void CCManyBody::print_eigensystem(int ndets, double** Heff,double*& eigenvector)
{
  if(ndets < 8){
    fprintf(outfile,"\n\n  Heff Matrix\n");
    for(int i=0;i<ndets;i++){
      fprintf(outfile,"\n  ");
      for(int j=0;j<ndets;j++)
        fprintf(outfile," %22.15f",Heff[i][j]);
    }
  }

  std::vector<std::pair<double,int> > eigenvector_index_pair;
  for(int i = 0; i < ndets; ++i){
    eigenvector_index_pair.push_back(make_pair(eigenvector[i]*eigenvector[i],i));
  }
  sort(eigenvector_index_pair.begin(),eigenvector_index_pair.end(),greater<pair<double,int> >());
  int max_size_list = std::min(10,static_cast<int>(eigenvector_index_pair.size()));
  fprintf(outfile,"\n\n  Most important determinants in the wave function");
  fprintf(outfile,"\n\n  determinant  eigenvector   eigenvector^2\n");
  for(int i = 0; i < max_size_list; ++i){
    fprintf(outfile,"\n  %11d   %9.6f    %9.6f  %s",eigenvector_index_pair[i].second
                                             ,eigenvector[eigenvector_index_pair[i].second]
                                             ,eigenvector_index_pair[i].first
                                             ,moinfo->get_determinant_label(eigenvector_index_pair[i].second).c_str());
  }
}

/**
 * This function computes \f$ E = \mathbf{c}^{\dagger} \mathbf{H} \mathbf{c} \f$
 * @param ndets size of the \f$ \mathbf{c} \f$ vector
 * @param H the \f$ \mathbf{H} \f$ matrix stored as a double**
 * @param c the \f$ \mathbf{c} \f$ vector stored as a double*
 * @return \f$ E \f$
 */
double CCManyBody::c_H_c(int ndets, double** H,double*& c)
{
    double energy = 0.0;
    for(int i=0;i<ndets;i++)
      for(int j=0;j<ndets;j++)
        energy += c[i]*H[i][j]*c[j];
    return(energy);
}

/**
 * This function computes the left and right eigenvalues of a generic real matrix
 * @param root selects the root for which the left-eigenvector must be saved
 * @param ndets size of the matrix
 * @param Heff the \f$ \mathbf{H}^{\mathrm{eff}} \f$ matrix stored as a double**
 * @param eigenvector the \f$ \mathbf{c} \f$ left-eigenvector stored as a double*  * @param initial a bool used to enable root following. initial = true allows you to select a root while initial = false follows the root that has the largest overlap with the previous eigenvector
 * @return
 */
double CCManyBody::diagonalize_Heff(int root,int ndets, double** Heff,double*& right_eigenvector,double*& left_eigenvector, bool initial)
{
  double      energy;
  double*     real;
  double*     imaginary;
  double*     work;
  double**    left;
  double**    right;
  double**    H;

  int lwork = 6*ndets*ndets;
  allocate1(double,work,lwork);
  allocate1(double,real,ndets);
  allocate1(double,imaginary,ndets);

  allocate2(double,H,ndets,ndets);
  allocate2(double,left,ndets,ndets);
  allocate2(double,right,ndets,ndets);

  for(int i=0;i<ndets;i++)
    for(int j=0;j<ndets;j++)
      H[j][i] = Heff[i][j];

  int info;

  F_DGEEV("V","V",&ndets, &(H[0][0]), &ndets, &(real[0]), &(imaginary[0]),
      &(left[0][0]), &ndets, &(right[0][0]), &ndets, &(work[0]), &lwork, &info);

  sort_eigensystem(ndets,real,imaginary,left,right);

  if(initial){
    if(ndets < 8){
      fprintf(outfile,"\n\n  Heff Matrix\n");
      for(int i=0;i<ndets;i++){
        fprintf(outfile,"\n  ");
        for(int j=0;j<ndets;j++)
          fprintf(outfile," %22.12f",Heff[i][j]);
      }

      fprintf(outfile,"\n\n  Left Matrix\n");
      for(int i=0;i<ndets;i++){
        fprintf(outfile,"\n  ");
        for(int j=0;j<ndets;j++)
          fprintf(outfile," %22.12f",left[j][i]);
      }

      fprintf(outfile,"\n\n  Right Matrix\n");
      for(int i=0;i<ndets;i++){
        fprintf(outfile,"\n  ");
        for(int j=0;j<ndets;j++)
          fprintf(outfile," %22.12f",right[j][i]);
      }

      fprintf(outfile,"\n\n  Real                  Imaginary\n");
      for(int i=0;i<ndets;i++)
        fprintf(outfile,"\n  %22.12f   %22.12f",real[i],imaginary[i]);
      fprintf(outfile,"\n");
    }else{
      fprintf(outfile,"\n\n  There are too many determinants to print the eigensystem");
    }
    fprintf(outfile,"\n\n  The eigenvalue for root %d is %.12f (%.12f)",root,real[root],imaginary[root]);
  }

  // Select the eigenvector to follow
  if(initial){
    for(int k=0;k<ndets;k++){
      zeroth_order_eigenvector[k] = right[root][k];
      right_eigenvector[k]        = right[root][k];
      left_eigenvector[k]         =  left[root][k];
    }
    energy = real[root];
    // Eliminate the triplet solution if required
    if((options_.get_bool("LOCK_SINGLET")==1)&&(ndets==4)){
      if((fabs(right_eigenvector[0])<5.0e-2)&& (fabs(right_eigenvector[3])<5.0e-2) && ((right_eigenvector[1]/right_eigenvector[2])<-0.5)){
        fprintf(outfile,"\n\tSelecting root %d since original root is a triplet\n",root+1);
        root++;
        for(int k=0;k<ndets;k++){
          right_eigenvector[k] = right[root][k];
          left_eigenvector[k]  =  left[root][k];
        }
        energy = real[root];
      }
    }
  }
  else // find vector with maximum overlap
  {
    int select_vect=0;
    double max_overlap=0.0;
    double overlap=0.0;
    for(int i=0;i<ndets;i++){
      overlap=0.0;
      for(int m=0;m<ndets;m++)
        overlap += zeroth_order_eigenvector[m] * right[i][m];
      overlap=sqrt(overlap*overlap);
      if(overlap>max_overlap){
        select_vect=i;
        max_overlap=overlap;
      }
    }
    for(int m=0;m<ndets;m++){
      right_eigenvector[m] = right[select_vect][m];
      left_eigenvector[m]  =  left[select_vect][m];
    }
    energy = real[select_vect];
  }

  // Normalize the left-eigenvector to <L|R> = 1
  double lnorm = 0.0;
  for(int m = 0; m < ndets; m++){
    lnorm += right_eigenvector[m] * left_eigenvector[m];
  }

  for(int m = 0; m < ndets; m++){
    left_eigenvector[m] = left_eigenvector[m] / lnorm;
  }


  release1(work);
  release1(real);
  release1(imaginary);
  release2(H);
  release2(left);
  release2(right);
  return(energy);
}

void CCManyBody::sort_eigensystem(int ndets,double*& real,double*& imaginary,double**& left,double**& right)
{
  std::vector<std::pair<double, int> > pairs;
  for(int i=0;i<ndets;i++)
    pairs.push_back(make_pair(real[i],i));
  sort(pairs.begin(),pairs.end());

  double*  tempv;
  double** tempm;
  allocate1(double,tempv,ndets);
  allocate2(double,tempm,ndets,ndets);

  for(int i=0;i<ndets;i++) tempv[i] = real[pairs[i].second];
  for(int i=0;i<ndets;i++) real[i]  = tempv[i];

  for(int i=0;i<ndets;i++) tempv[i]     = imaginary[pairs[i].second];
  for(int i=0;i<ndets;i++) imaginary[i] = tempv[i];

  for(int i=0;i<ndets;i++)
    for(int j=0;j<ndets;j++)
      tempm[i][j] = left[pairs[i].second][j];
  for(int i=0;i<ndets;i++)
    for(int j=0;j<ndets;j++)
      left[i][j] = tempm[i][j];

  for(int i=0;i<ndets;i++)
    for(int j=0;j<ndets;j++)
      tempm[i][j] = right[pairs[i].second][j];
  for(int i=0;i<ndets;i++)
    for(int j=0;j<ndets;j++)
      right[i][j] = tempm[i][j];

  release1(tempv);
  release2(tempm);
}

//void CCManyBody::zero_internal_amps()
//{
//  if(options_get_bool("ZERO_INTERNAL_AMPS")){
//    // Zero internal amplitudes for unique reference i
//    for(int i=0;i<moinfo->get_nunique();i++){
//      int unique_i = moinfo->get_ref_number(i,UniqueRefs);
//      // Loop over reference j
//      for(int j=0;j<moinfo->get_ref_size(AllRefs);j++){
//        vector<pair<int,int> >  alpha_internal_excitation = moinfo->get_alpha_internal_excitation(unique_i,j);
//        vector<pair<int,int> >   beta_internal_excitation = moinfo->get_beta_internal_excitation(unique_i,j);
//
//        // Zero alpha-alpha single excitations
//        if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==0)){
//          blas->get_MatTmp("t1[o][v]",unique_i,none)->set_two_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[0].second,
//                                            0.0);
//        }
//
//        // Zero beta-beta single excitations
//        if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==1))
//          blas->get_MatTmp("t1[O][V]",unique_i,none)->set_two_address_element(
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//
//        // Zero (alpha,alpha)->(alpha,alpha) double excitations (all permutations)
//        if((alpha_internal_excitation.size()==2)&&(beta_internal_excitation.size()==0)){
//          blas->get_MatTmp("t2[oo][vv]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[1].first,
//                                            alpha_internal_excitation[0].second,
//                                            alpha_internal_excitation[1].second,
//                                            0.0);
//          blas->get_MatTmp("t2[oo][vv]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[1].first,
//                                            alpha_internal_excitation[1].second,
//                                            alpha_internal_excitation[0].second,
//                                            0.0);
//          blas->get_MatTmp("t2[oo][vv]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[1].first,
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[0].second,
//                                            alpha_internal_excitation[1].second,
//                                            0.0);
//          blas->get_MatTmp("t2[oo][vv]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[1].first,
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[1].second,
//                                            alpha_internal_excitation[0].second,
//                                            0.0);
//        }
//
//        // Zero (alpha,beta)->(alpha,beta) double excitations
//        if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==1)){
//          blas->get_MatTmp("t2[oO][vV]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            beta_internal_excitation[0].first,
//                                            alpha_internal_excitation[0].second,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//        }
//
//        // Zero (beta,beta)->(beta,beta) double excitations (all permutations)
//        if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==2)){
//          blas->get_MatTmp("t2[OO][VV]",unique_i,none)->set_four_address_element(
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[1].first,
//                                            beta_internal_excitation[0].second,
//                                            beta_internal_excitation[1].second,
//                                            0.0);
//          blas->get_MatTmp("t2[OO][VV]",unique_i,none)->set_four_address_element(
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[1].first,
//                                            beta_internal_excitation[1].second,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//          blas->get_MatTmp("t2[OO][VV]",unique_i,none)->set_four_address_element(
//                                            beta_internal_excitation[1].first,
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[0].second,
//                                            beta_internal_excitation[1].second,
//                                            0.0);
//          blas->get_MatTmp("t2[OO][VV]",unique_i,none)->set_four_address_element(
//                                            beta_internal_excitation[1].first,
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[1].second,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//        }
//      }
//    }
//
//    // Print the t-amplitudes
//    DEBUGGING(3,
//      blas->print("t1[o][v]{u}");
//      blas->print("t1[O][V]{u}");
//      blas->print("t2[oo][vv]{u}");
//      blas->print("t2[oO][vV]{u}");
//      blas->print("t2[OO][VV]{u}");
//    )
//  }else{
//    fprintf(outfile,"\n  Warning: the internal amplitudes are not zeroed.\n  This is not proper Mk-MRCC. Size-extensivity might be lost\n");
//  }
//}
//
//
//void CCManyBody::zero_t1_internal_amps()
//{
//  if(options_get_bool("ZERO_INTERNAL_AMPS")){
//    // Zero internal amplitudes for unique reference i
//    for(int i=0;i<moinfo->get_nunique();i++){
//      int unique_i = moinfo->get_ref_number(i,UniqueRefs);
//      // Loop over reference j
//      for(int j=0;j<moinfo->get_ref_size(AllRefs);j++){
//        vector<pair<int,int> >  alpha_internal_excitation = moinfo->get_alpha_internal_excitation(unique_i,j);
//        vector<pair<int,int> >   beta_internal_excitation = moinfo->get_beta_internal_excitation(unique_i,j);
//
//        // Zero alpha-alpha single excitations
//        if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==0))
//          blas->get_MatTmp("t1[o][v]",unique_i,none)->set_two_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[0].second,
//                                            0.0);
//
//        // Zero beta-beta single excitations
//        if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==1))
//          blas->get_MatTmp("t1[O][V]",unique_i,none)->set_two_address_element(
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//      }
//    }
//
//    // Print the t-amplitudes
//    DEBUGGING(3,
//      blas->print("t1[o][v]{u}");
//      blas->print("t1[O][V]{u}");
//    )
//  }else{
//    fprintf(outfile,"\n  Warning: the internal amplitudes are not zeroed.\n  This is not proper Mk-MRCC. Size-extensivity might be lost\n");
//  }
//}
//
//void CCManyBody::zero_internal_delta_amps()
//{
//  if(options_get_bool("ZERO_INTERNAL_AMPS")){
//    // Zero internal amplitudes for unique reference i
//    for(int i=0;i<moinfo->get_nunique();i++){
//      int unique_i = moinfo->get_ref_number(i,UniqueRefs);
//      // Loop over reference j
//      for(int j=0;j<moinfo->get_ref_size(AllRefs);j++){
//        vector<pair<int,int> >  alpha_internal_excitation = moinfo->get_alpha_internal_excitation(unique_i,j);
//        vector<pair<int,int> >   beta_internal_excitation = moinfo->get_beta_internal_excitation(unique_i,j);
//
//        // Zero alpha-alpha single excitations
//        if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==0))
//          blas->get_MatTmp("t1_delta[o][v]",unique_i,none)->set_two_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[0].second,
//                                            0.0);
//
//        // Zero beta-beta single excitations
//        if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==1))
//          blas->get_MatTmp("t1_delta[O][V]",unique_i,none)->set_two_address_element(
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//
//        // Zero (alpha,alpha)->(alpha,alpha) double excitations (all permutations)
//        if((alpha_internal_excitation.size()==2)&&(beta_internal_excitation.size()==0)){
//          blas->get_MatTmp("t2_delta[oo][vv]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[1].first,
//                                            alpha_internal_excitation[0].second,
//                                            alpha_internal_excitation[1].second,
//                                            0.0);
//          blas->get_MatTmp("t2_delta[oo][vv]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[1].first,
//                                            alpha_internal_excitation[1].second,
//                                            alpha_internal_excitation[0].second,
//                                            0.0);
//          blas->get_MatTmp("t2_delta[oo][vv]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[1].first,
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[0].second,
//                                            alpha_internal_excitation[1].second,
//                                            0.0);
//          blas->get_MatTmp("t2_delta[oo][vv]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[1].first,
//                                            alpha_internal_excitation[0].first,
//                                            alpha_internal_excitation[1].second,
//                                            alpha_internal_excitation[0].second,
//                                            0.0);
//        }
//
//        // Zero (alpha,beta)->(alpha,beta) double excitations
//        if((alpha_internal_excitation.size()==1)&&(beta_internal_excitation.size()==1)){
//          blas->get_MatTmp("t2_delta[oO][vV]",unique_i,none)->set_four_address_element(
//                                            alpha_internal_excitation[0].first,
//                                            beta_internal_excitation[0].first,
//                                            alpha_internal_excitation[0].second,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//        }
//
//        // Zero (beta,beta)->(beta,beta) double excitations (all permutations)
//        if((alpha_internal_excitation.size()==0)&&(beta_internal_excitation.size()==2)){
//          blas->get_MatTmp("t2_delta[OO][VV]",unique_i,none)->set_four_address_element(
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[1].first,
//                                            beta_internal_excitation[0].second,
//                                            beta_internal_excitation[1].second,
//                                            0.0);
//          blas->get_MatTmp("t2_delta[OO][VV]",unique_i,none)->set_four_address_element(
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[1].first,
//                                            beta_internal_excitation[1].second,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//          blas->get_MatTmp("t2_delta[OO][VV]",unique_i,none)->set_four_address_element(
//                                            beta_internal_excitation[1].first,
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[0].second,
//                                            beta_internal_excitation[1].second,
//                                            0.0);
//          blas->get_MatTmp("t2_delta[OO][VV]",unique_i,none)->set_four_address_element(
//                                            beta_internal_excitation[1].first,
//                                            beta_internal_excitation[0].first,
//                                            beta_internal_excitation[1].second,
//                                            beta_internal_excitation[0].second,
//                                            0.0);
//        }
//      }
//    }
//
//    // Print the t-amplitudes
//    DEBUGGING(3,
//      blas->print("t1_delta[o][v]{u}");
//      blas->print("t1_delta[O][V]{u}");
//      blas->print("t2_delta[oo][vv]{u}");
//      blas->print("t2_delta[oO][vV]{u}");
//      blas->print("t2_delta[OO][VV]{u}");
//    )
//  }else{
//    fprintf(outfile,"\n  Warning: the internal amplitudes are not zeroed.\n  This is not proper Mk-MRCC. Size-extensivity might be lost\n");
//  }
//}

}} /* End Namespaces */
