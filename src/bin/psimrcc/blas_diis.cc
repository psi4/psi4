#include <limits>
#include <cmath>

#include <boost/shared_ptr.hpp>
#include <psifiles.h>
#include <liboptions/liboptions.h>
#include <libmoinfo/libmoinfo.h>
#include <libutil/libutil.h>
#include <libpsio/psio.hpp>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

#include "algebra_interface.h"
#include "blas.h"
#include "matrix.h"


namespace psi{
    extern FILE *outfile;
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
    fprintf(outfile,"   S");


    // Decide if we are doing a DIIS extrapolation in this cycle
    bool do_diis_extrapolation = false;
    if(diis_type == DiisEachCycle){
      if(cycle >= options_.get_int("DIIS_MAX_VECS") + options_.get_int("START_DIIS"))
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
      fprintf(outfile,"/E");
      if(singularities_found)
        fprintf(outfile," (singularities found)");
      release1(diis_A);
      release2(diis_B);
    }
  }
}

}} /* End Namespaces */
