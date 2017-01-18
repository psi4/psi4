/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#include <limits>
#include <cmath>

 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP
#include "psi4/psifiles.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"

#include "algebra_interface.h"
#include "blas.h"
#include "matrix.h"


namespace psi{

    namespace psimrcc{
    extern MOInfo *moinfo;
    extern MemoryManager* memory_manager;

using namespace std;

vector<pair<string,string> > diis_matrices;
const double diis_singular_tollerance = 1.0e-12;

void CCBLAS::diis_add(string amps, string delta_amps)
{
  vector<string> amps_names = moinfo->get_matrix_names(amps);
  vector<string> delta_amps_names = moinfo->get_matrix_names(delta_amps);
  for(size_t n=0;n<amps_names.size();n++){
    diis_matrices.push_back(make_pair(amps_names[n],delta_amps_names[n]));
  }
}

void CCBLAS::diis_save_t_amps(int cycle)
{
  if(options_.get_int("DIIS_MAX_VECS") != 0){
    int diis_step = cycle % options_.get_int("DIIS_MAX_VECS");
    for(vector<pair<string,string> >::iterator it=diis_matrices.begin();it!=diis_matrices.end();++it){
      for(int h=0;h<moinfo->get_nirreps();h++){
        CCMatIrTmp Amps = get_MatIrTmp(it->first,h,none);
        double** matrix = Amps->get_matrix()[h];
        size_t   block_sizepi = Amps->get_block_sizepi(h);
        if(block_sizepi>0){
          char data_label[80];
          sprintf(data_label,"%s_%s_%d_%d",(it->first).c_str(),"DIIS",h,diis_step);
          _default_psio_lib_->write_entry(PSIF_PSIMRCC_INTEGRALS,data_label,(char*)&(matrix[0][0]),block_sizepi*sizeof(double));
        }
      }
    }
  }
}

void CCBLAS::diis(int cycle, double delta, DiisType diis_type)
{
  if(options_.get_int("DIIS_MAX_VECS") != 0){
    int diis_step = cycle % options_.get_int("DIIS_MAX_VECS");

    for(vector<pair<string,string> >::iterator it=diis_matrices.begin();it!=diis_matrices.end();++it){
      if(it->second.find("t3_delta")==string::npos){
        for(int h=0;h<moinfo->get_nirreps();h++){
          CCMatIrTmp DeltaAmps = get_MatIrTmp(it->second,h,none);
          double** matrix = DeltaAmps->get_matrix()[h];
          size_t   block_sizepi = DeltaAmps->get_block_sizepi(h);
          if(block_sizepi>0){
            char data_label[80];
            sprintf(data_label,"%s_%s_%d_%d",(it->second).c_str(),"DIIS",h,diis_step);
            _default_psio_lib_->write_entry(PSIF_PSIMRCC_INTEGRALS,data_label,(char*)&(matrix[0][0]),block_sizepi*sizeof(double));
          }
        }
      }
    }
    outfile->Printf("   S");


    // Decide if we are doing a DIIS extrapolation in this cycle
    bool do_diis_extrapolation = false;
    if(diis_type == DiisEachCycle){
      if(cycle >= options_.get_int("DIIS_MAX_VECS") + options_.get_int("DIIS_START"))
        do_diis_extrapolation = true;
    }else if(diis_type == DiisCC){
      if(diis_step == options_.get_int("DIIS_MAX_VECS")-1)
        do_diis_extrapolation = true;
    }

    // Do a DIIS step
    if(do_diis_extrapolation){
      double** diis_B;
      double*  diis_A;
      allocate1(double,diis_A,options_.get_int("DIIS_MAX_VECS")+1);
      allocate2(double,diis_B,options_.get_int("DIIS_MAX_VECS")+1,options_.get_int("DIIS_MAX_VECS")+1);
      bool singularities_found = false;
      for(vector<pair<string,string> >::iterator it=diis_matrices.begin();it!=diis_matrices.end();++it){
        // Zero A and B
        for(int i=0;i<options_.get_int("DIIS_MAX_VECS");i++){
          diis_A[i]=0.0;
          diis_B[i][options_.get_int("DIIS_MAX_VECS")]=diis_B[options_.get_int("DIIS_MAX_VECS")][i]=-1.0;
          for(int j=0;j<options_.get_int("DIIS_MAX_VECS");j++)
            diis_B[i][j]=0.0;
        }
        diis_B[options_.get_int("DIIS_MAX_VECS")][options_.get_int("DIIS_MAX_VECS")]=0.0;
        diis_A[options_.get_int("DIIS_MAX_VECS")]=-1.0;

        // Build B
        for(int h=0;h<moinfo->get_nirreps();h++){
          CCMatIrTmp Amps       = get_MatIrTmp(it->first,h,none);
          size_t   block_sizepi = Amps->get_block_sizepi(h);
          if(block_sizepi>0){
            double*  i_matrix;
            double*  j_matrix;
            allocate1(double,i_matrix,block_sizepi);
            allocate1(double,j_matrix,block_sizepi);

            // Build the diis_B matrix
            for(int i=0;i<options_.get_int("DIIS_MAX_VECS");i++){
              // Load vector i irrep h
              char i_data_label[80];
              sprintf(i_data_label,"%s_%s_%d_%d",(it->second).c_str(),"DIIS",h,i);
              _default_psio_lib_->read_entry(PSIF_PSIMRCC_INTEGRALS,i_data_label,(char*)&(i_matrix[0]),block_sizepi*sizeof(double));

              for(int j=i;j<options_.get_int("DIIS_MAX_VECS");j++){
                // Load vector j irrep h
                char j_data_label[80];
                sprintf(j_data_label,"%s_%s_%d_%d",(it->second).c_str(),"DIIS",h,j);
                _default_psio_lib_->read_entry(PSIF_PSIMRCC_INTEGRALS,j_data_label,(char*)&(j_matrix[0]),block_sizepi*sizeof(double));

                int dx = 1;
                int lenght = block_sizepi;
                if( block_sizepi < static_cast<size_t>(numeric_limits<int>::max()) ){
                  diis_B[i][j] += F_DDOT(&lenght,i_matrix,&dx,j_matrix,&dx);
                  diis_B[j][i] = diis_B[i][j];
                }else{
                  throw PSIEXCEPTION("The numeric limits for int was reached for F_DDOT");
                }
              }
            }
            release1(i_matrix);
            release1(j_matrix);
          }
        }

        // Solve B x = A
        int  matrix_size = options_.get_int("DIIS_MAX_VECS") + 1;
        int* IPIV = new int[matrix_size];
        int nrhs = 1;
        int info = 0;
        F_DGESV(&matrix_size, &nrhs, &(diis_B[0][0]),&matrix_size, &(IPIV[0]), &(diis_A[0]),&matrix_size, &info);
        delete[] IPIV;

        // Update T = sum t(i) * A(i);
        if(!info){
          for(int h=0;h<moinfo->get_nirreps();h++){
            CCMatIrTmp Amps       = get_MatIrTmp(it->first,h,none);
            size_t   block_sizepi = Amps->get_block_sizepi(h);
            if(block_sizepi>0){
              // Update the amplitudes
              double*  i_matrix;
              double*  j_matrix;
              allocate1(double,i_matrix,block_sizepi);
              allocate1(double,j_matrix,block_sizepi);
              double* t_matrix = &(Amps->get_matrix()[h][0][0]);
              Amps->zero_matrix_block(h);
              for(int i=0;i<options_.get_int("DIIS_MAX_VECS");i++){
                char i_data_label[80];
                sprintf(i_data_label,"%s_%s_%d_%d",(it->first).c_str(),"DIIS",h,i);
                _default_psio_lib_->read_entry(PSIF_PSIMRCC_INTEGRALS,i_data_label,(char*)&(i_matrix[0]),block_sizepi*sizeof(double));
                for(size_t n=0;n<block_sizepi;n++){
                  t_matrix[n] += diis_A[i]*i_matrix[n];
                }
              }
              release1(i_matrix);
              release1(j_matrix);
            }
          }
        }else{
          singularities_found = true;
        }

      }
      outfile->Printf("/E");
      if(singularities_found)
        outfile->Printf(" (singularities found)");
      release1(diis_A);
      release2(diis_B);
    }
  }
}

}} /* End Namespaces */
