/**
 *  @file ccsort_out_of_core.cpp
 *  @ingroup (PSIMRCC)
*/

#include <cstdio>
#include <libmoinfo/libmoinfo.h>

#include "blas.h"
#include "matrix.h"
#include "sort.h"
#include "transform.h"

extern FILE* outfile;

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;
    extern MemoryManager *memory_manager;

using namespace std;

/**
 * Builds the integral matrices on disk using an out-of-core algorithm
 */
void CCSort::build_integrals_out_of_core()
{
  trans->read_oei_from_transqt();
  // One electron contribution to the frozen core energy from each irrep
  efzc = 0.0;
  for(int i = 0; i < nfzc; ++i){
    int ii = frozen_core[i];
    efzc += 2.0 * trans->oei(ii,ii);
  }

  MatrixMap matrix_map= blas->get_MatrixMap();
  MatMapIt mat_it     = matrix_map.begin();
  MatMapIt mat_end    = matrix_map.end();
  int mat_irrep = 0;
  int cycle     = 0;

  size_t ccintegrals_memory = static_cast<size_t>(static_cast<double>(memory_manager->get_FreeMemory())*fraction_of_memory_for_sorting);


  fprintf(outfile,"\n\n  Sorting integrals:");
  fprintf(outfile,"\n    Memory available                       = %14lu bytes",
                  (unsigned long)memory_manager->get_FreeMemory());
  fprintf(outfile,"\n    Memory available for sorting           = %14lu bytes (%.1f%%)",
                  (unsigned long)ccintegrals_memory,fraction_of_memory_for_sorting*100.0);
  fflush(outfile);

  while(mat_it!=mat_end){
    fprintf(outfile,"\n\n    Pass %d:",cycle + 1);
    // Find how many matrices blocks we can store in 95% of the free memory and allocate them
    MatrixBlks to_be_processed;
    setup_out_of_core_list(mat_it,mat_irrep,mat_end,to_be_processed);
    int first_irrep = 0;
    int last_irrep  = 0;
    // The one-particle integrals are added at the beginning to avoid interfering with the
    // way the transformation code handles the process
    form_fock_one_out_of_core(to_be_processed);
    while(last_irrep < moinfo->get_nirreps()){
      last_irrep = trans->read_tei_mo_integrals_block(first_irrep);
      if(cycle==0) frozen_core_energy_out_of_core();
      sort_integrals_out_of_core(first_irrep,last_irrep,to_be_processed);
      trans->free_tei_mo_integrals_block(first_irrep,last_irrep);
      first_irrep = last_irrep;
    }
    // Write to disk and free memory
    dump_integrals_to_disk(to_be_processed);
    cycle++;
  }
}

void CCSort::setup_out_of_core_list(MatMapIt& mat_it,int& mat_irrep,MatMapIt& mat_end,MatrixBlks& to_be_processed)
{
  fprintf(outfile,"\n    Setting up the matrix list:");
  fflush(outfile);
  size_t ccintegrals_memory = static_cast<size_t>(static_cast<double>(memory_manager->get_FreeMemory())*fraction_of_memory_for_sorting);

  int blocks_added = 0;
  bool out_of_memory = false;
  while((mat_it != mat_end) && !out_of_memory){
    if(mat_it->second->is_integral() || mat_it->second->is_fock()){
      CCMatrix* Matrix = mat_it->second;
      while(mat_irrep < moinfo->get_nirreps() && !out_of_memory){
        size_t block_memory = Matrix->get_memorypi2(mat_irrep);
        if(block_memory < ccintegrals_memory){
          to_be_processed.push_back(make_pair(Matrix,mat_irrep));
          // Allocate the matrix, this will also take care of MOInfo::allocated_memory
          Matrix->allocate_block(mat_irrep);
          ccintegrals_memory -= block_memory;
          mat_irrep++;
          blocks_added++;
        }else{
          if(blocks_added == 0){
            fprintf(outfile,"\n    Matrix: %s irrep %d does not fit into memory",Matrix->get_label().c_str(),mat_irrep);
            fprintf(outfile,"\n            memory required = %14lu bytes",(unsigned long)block_memory);
            fflush(outfile);
          }
          out_of_memory = true;
        }
      }
      if(!out_of_memory)
        mat_irrep=0;
    }
    if(!out_of_memory)
      ++mat_it;
  }
  fprintf(outfile," added %d matrices blocks",blocks_added);
  fflush(outfile);
}

void CCSort::sort_integrals_out_of_core(int first_irrep, int last_irrep, MatrixBlks& to_be_processed)
{
  for(MatBlksIt block_it=to_be_processed.begin();block_it!=to_be_processed.end();++block_it){
    form_fock_out_of_core(block_it->first,block_it->second);
    form_two_electron_integrals_out_of_core(block_it->first,block_it->second);
  }
}

void CCSort::form_fock_one_out_of_core(MatrixBlks& to_be_processed)
{
  for(MatBlksIt block_it=to_be_processed.begin();block_it!=to_be_processed.end();++block_it){
    CCMatrix* Matrix = block_it->first;
    if(Matrix->is_fock()){
      int h = block_it->second;
      double*** matrix = Matrix->get_matrix();
      short* pq = new short[2];

      for(size_t i = 0; i < Matrix->get_left_pairpi(h); ++i){
        for(size_t j = 0; j < Matrix->get_right_pairpi(h); ++j){
          // Find p and q from the pairs
          Matrix->get_two_indices_pitzer(pq,h,i,j);
          // Add the h(p,q) contribution
          matrix[h][i][j] = trans->oei(pq[0],pq[1]);
        }
      }
      delete[] pq;
    }
  }
}

void CCSort::form_fock_out_of_core(CCMatrix* Matrix, int h)
{
  if(Matrix->is_fock()){
    string label     = Matrix->get_label();
    double*** matrix = Matrix->get_matrix();
    short* pq = new short[2];
    const intvec& oa2p = moinfo->get_occ_to_mo();

    bool alpha = true;
    if((label.find("O")!=string::npos) || (label.find("V")!=string::npos) || (label.find("A")!=string::npos) || (label.find("F")!=string::npos)) // NB This was missing the last bit, this might be a problem
      alpha = false;

    // N.B. Never introduce Matrices/Vectors with O or V in the name before you compute the Fock matrix elements
    vector<int> aocc = moinfo->get_aocc(Matrix->get_reference(),AllRefs);
    vector<int> bocc = moinfo->get_bocc(Matrix->get_reference(),AllRefs);

    for(size_t i = 0; i < Matrix->get_left_pairpi(h); ++i)
      for(size_t j = 0; j < Matrix->get_right_pairpi(h); ++j){
        // Find p and q from the pairs
        Matrix->get_two_indices_pitzer(pq,h,i,j);
        // Add the core contribution//
        for(int k=0;k<nfzc;k++){
          int kk=frozen_core[k];
          matrix[h][i][j]+=add_fock_two_out_of_core(pq[0],pq[1],kk,true);
          matrix[h][i][j]+=add_fock_two_out_of_core(pq[0],pq[1],kk,false);
        }
        for(size_t k = 0; k < aocc.size(); ++k){
          int kk=oa2p[aocc[k]];
          if(alpha)
            matrix[h][i][j]+=add_fock_two_out_of_core(pq[0],pq[1],kk,true);
          else
            matrix[h][i][j]+=add_fock_two_out_of_core(pq[0],pq[1],kk,false);
        }
        for(size_t k = 0; k < bocc.size(); ++k){
          int kk=oa2p[bocc[k]];
          if(!alpha)
            matrix[h][i][j]+=add_fock_two_out_of_core(pq[0],pq[1],kk,true);
          else
            matrix[h][i][j]+=add_fock_two_out_of_core(pq[0],pq[1],kk,false);
        }
      }
    delete[] pq;
  }
}


void CCSort::form_two_electron_integrals_out_of_core(CCMatrix* Matrix, int h)
{
  if(Matrix->is_integral()){
    short*      pqrs = new short[4];
    double*** matrix = Matrix->get_matrix();
    bool antisymmetric = Matrix->is_antisymmetric();
    if(Matrix->is_chemist()){
      for(size_t i = 0; i < Matrix->get_left_pairpi(h); ++i)
        for(size_t j = 0; j < Matrix->get_right_pairpi(h); ++j){
          Matrix->get_four_indices_pitzer(pqrs,h,i,j);
          // From (pq|rs) = <pr|qs> we define
          // (pq:rs) = <pr:qs> = (pq|rs) - (ps|qr)

          // Add the +<pr|qs> = (pq|rs) contribution
          matrix[h][i][j] += trans->tei_block(pqrs[0],pqrs[1],pqrs[2],pqrs[3]);

          // Add the -<pq|sr> = -(ps|qr) contribution
          if(antisymmetric)
            matrix[h][i][j] -= trans->tei_block(pqrs[0],pqrs[3],pqrs[1],pqrs[2]);
        }
    }else{
      for(size_t i = 0; i < Matrix->get_left_pairpi(h); ++i)
        for(size_t j = 0; j < Matrix->get_right_pairpi(h); ++j){
          Matrix->get_four_indices_pitzer(pqrs,h,i,j);
          // Add the +<pq|rs> = (pr|qs) contribution
          matrix[h][i][j] += trans->tei_block(pqrs[0],pqrs[2],pqrs[1],pqrs[3]);

          // Add the -<pq|sr> = -(ps|qr) contribution
          if(antisymmetric)
            matrix[h][i][j] -= trans->tei_block(pqrs[0],pqrs[3],pqrs[1],pqrs[2]);
        }
    }
    delete[] pqrs;
  }
}

double CCSort::add_fock_two_out_of_core(int p, int q, int k, bool exchange)
{
  // Add the (pq|kk) contribution
  double term = trans->tei_block(p,q,k,k);
  // Add the -(pk|qk) contribution
  if(exchange)
    term -= trans->tei_block(p,k,q,k);
  return(term);
}

void CCSort::frozen_core_energy_out_of_core()
{
  // Two electron contribution to the frozen core energy
  for(int i=0;i<nfzc;i++){
    for(int j=0;j<nfzc;j++){
      int ii=frozen_core[i];
      int jj=frozen_core[j];
      efzc+=2.0*trans->tei_block(ii,ii,jj,jj);
      efzc-=trans->tei_block(ii,jj,ii,jj);
    }
  }
}

/**
 * Dump (write + free memory) the two electron integral matrix blocks contained in the list to_be_processed
 * @param to_be_processed
 */
void CCSort::dump_integrals_to_disk(MatrixBlks& to_be_processed)
{
  for(MatBlksIt block_it=to_be_processed.begin();block_it!=to_be_processed.end();++block_it){
    CCMatrix* Matrix = block_it->first;
    Matrix->dump_block_to_disk(block_it->second);
  }
}

}} /* End Namespaces */
