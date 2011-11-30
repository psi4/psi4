#include <iostream>
#include <cmath>

#include <libutil/libutil.h>
#include <libmoinfo/libmoinfo.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.hpp>

#include "algebra_interface.h"
#include "blas.h"
#include "debugging.h"
#include "index.h"
#include "matrix.h"
#include "operation.h"

namespace psi{
    extern FILE *outfile;
    namespace psimrcc{
    extern MOInfo *moinfo;

using namespace std;

void CCOperation::setup_contractions()
{
  Timer PartA;

  CCMatTmp AMatTmp = blas->get_MatTmp(A_Matrix,none);
  check_and_zero_target();

  bool need_sort = false;
  if(reindexing.size()>0)
    need_sort = true;
  CCIndex*  T_left  =A_Matrix->get_left();
  CCIndex*  T_right =A_Matrix->get_right();
  double*** T_matrix=A_Matrix->get_matrix();

  // Determine the target indexing if sorting is needed
  if(need_sort){
    if(operation[0]=='1'){
      T_left=B_Matrix->get_right();
    }  else {
      T_left=B_Matrix->get_left();
    }
    if(operation[2]=='1'){
      T_right=C_Matrix->get_right();
    }  else {
      T_right=C_Matrix->get_left();
    }
    size_t T_matrix_offset = 0;
    T_matrix = new double**[moinfo->get_nirreps()];
    for(int irrep=0;irrep<moinfo->get_nirreps();irrep++){
      T_matrix[irrep] = new double*[T_left->get_pairpi(irrep)];
      size_t block_size = T_left->get_pairpi(irrep)*T_right->get_pairpi(irrep);
        for(size_t i=0;i<T_left->get_pairpi(irrep);i++){
          T_matrix[irrep][i]=&(local_work[T_matrix_offset+i*T_right->get_pairpi(irrep)]);
        }
        T_matrix_offset+=block_size;
    }
    if(T_matrix_offset>0)
      zero_arr(&(local_work[0]),T_matrix_offset);
  }

  PartA_timing += PartA.get();
  Timer PartB;

  double**  A_matrix;
  double**  B_matrix;
  double**  C_matrix;

  for(int h=0;h<moinfo->get_nirreps();h++){
    bool B_on_disk = false;
    bool C_on_disk = false;
    if(!B_Matrix->is_block_allocated(h)){
      if(B_Matrix->is_integral()){
        B_on_disk = true;
      }else{
        B_Matrix->load_irrep(h);
      }
    }
    if(!C_Matrix->is_block_allocated(h)){
      if(C_Matrix->is_integral()){
        C_on_disk = true;
      }else{
        C_Matrix->load_irrep(h);
      }
    }
    if(B_on_disk && C_on_disk)
//      print_developing(outfile,"Both on disk matrix multiply",__FILE__,__LINE__);

    //////////////////////////////////////////////////////////
    // Case I. A,B,C are in core. Perform direct a BLAS call
    // in one pass.
    //////////////////////////////////////////////////////////
    if(!B_on_disk && !C_on_disk){
      size_t offset = 0;
      A_matrix = T_matrix[h];
      B_matrix = B_Matrix->get_matrix()[h];
      C_matrix = C_Matrix->get_matrix()[h];
      size_t rows_A = T_left->get_pairpi(h);
      size_t cols_A = T_right->get_pairpi(h);
      size_t rows_B = B_Matrix->get_left_pairpi(h);
      size_t cols_B = B_Matrix->get_right_pairpi(h);
      size_t rows_C = C_Matrix->get_left_pairpi(h);
      size_t cols_C = C_Matrix->get_right_pairpi(h);
      // Now call BLAS
      // Start a timer
      Timer timer;
      // Do the job
      contract_in_core(A_matrix,B_matrix,C_matrix,B_on_disk,C_on_disk,rows_A,rows_B,rows_C,cols_A,cols_B,cols_C,offset);
      // Store the timing in moinfo
      moinfo->add_dgemm_timing(timer.get());
    }

    //////////////////////////////////////////////////////////
    // Case II. A,C are in core. B is on disk. Perform several
    // BLAS calls.
    //////////////////////////////////////////////////////////
    if(B_on_disk && !C_on_disk){
      // Assign pointers to in core matrices
      double** A_matrix = T_matrix[h];
      double** B_matrix = new double*[1];
      B_matrix[0] = &out_of_core_buffer[0];
      double** C_matrix = C_Matrix->get_matrix()[h];

      int strip = 0;
      size_t offset = 0;
      bool done = false;
      while(!done){
        size_t strip_length = B_Matrix->read_strip_from_disk(h,strip,out_of_core_buffer);
        if(strip_length == 0){
          done = true;
        }else{
          // Compute the size of the sub-matrix B
          size_t rows_A = T_left->get_pairpi(h);
          size_t cols_A = T_right->get_pairpi(h);
          size_t rows_B = strip_length;
          size_t cols_B = B_Matrix->get_right_pairpi(h);
          size_t rows_C = C_Matrix->get_left_pairpi(h);
          size_t cols_C = C_Matrix->get_right_pairpi(h);

          // Now call BLAS
          // Start a timer
          Timer timer;
          // Do the job
          contract_in_core(A_matrix,B_matrix,C_matrix,B_on_disk,C_on_disk,rows_A,rows_B,rows_C,cols_A,cols_B,cols_C,offset);
          // Store the timing in moinfo
          moinfo->add_dgemm_timing(timer.get());
          offset += strip_length;
        }
        strip++;
      }
      delete[] B_matrix;
    }
    //////////////////////////////////////////////////////////
    // Case III. A,B are in core. C is on disk. Perform several
    // BLAS calls.
    //////////////////////////////////////////////////////////
    if(!B_on_disk && C_on_disk){
      // Assign pointers to in core matrices
      double** A_matrix = T_matrix[h];
      double** B_matrix = B_Matrix->get_matrix()[h];
      double** C_matrix = new double*[1];
      C_matrix[0] = &out_of_core_buffer[0];

      int strip = 0;
      size_t offset = 0;
      bool done = false;
      while(!done){
        size_t strip_length = C_Matrix->read_strip_from_disk(h,strip,out_of_core_buffer);
        if(strip_length == 0){
          done = true;
        }else{
          // Compute the size of the sub-matrix B
          size_t rows_A = T_left->get_pairpi(h);
          size_t cols_A = T_right->get_pairpi(h);
          size_t rows_B = B_Matrix->get_left_pairpi(h);
          size_t cols_B = B_Matrix->get_right_pairpi(h);
          size_t rows_C = strip_length;
          size_t cols_C = C_Matrix->get_right_pairpi(h);

          // Now call BLAS
          // Start a timer
          Timer timer;
          // Do the job
          contract_in_core(A_matrix,B_matrix,C_matrix,B_on_disk,C_on_disk,rows_A,rows_B,rows_C,cols_A,cols_B,cols_C,offset);
          // Store the timing in moinfo
          moinfo->add_dgemm_timing(timer.get());
          offset += strip_length;
        }
        strip++;
      }
      delete[] C_matrix;
    }
  }  // end of for loop over irreps

  PartB_timing += PartB.get();
  Timer PartC;
  if(need_sort){
    sort(T_left,T_right,T_matrix,1.0);
    for(int h=0;h<moinfo->get_nirreps();h++)
//       if(T_left->get_pairpi(h)*T_right->get_pairpi(h)>0)
        delete[] T_matrix[h];
    delete[] T_matrix;
  }
  PartC_timing += PartC.get();
}

void CCOperation::contract_in_core(double** A_matrix,double** B_matrix,double** C_matrix,bool B_on_disk,bool C_on_disk,int rows_A,int rows_B,int rows_C,int cols_A,int cols_B,int cols_C,int offset)
{
  // CASE A, worst case scenario
  //
  //            -----------------------------
  //           |         ---------------     |
  //           |        |               |    |
  //   A[A_l][A_r] = B[B_i][B_c] 1@1 C[C_c][C_i]
  //      |                  |
  //       ------------------
  //  Solution:
  //  1) Transpose B,C
  //  2) Perform contraction of the indices
  //  3) Transpose the result and store in A ???
  if(operation=="1@1"){
    double beta = 1.0;
    int m = rows_A;
    int n = cols_A;
    int k;
    if(m*n == 0) return;
    if(!B_on_disk && !C_on_disk){
      k = rows_B;
    }
    if( B_on_disk && !C_on_disk){
      k = rows_B;
    }
    if(!B_on_disk &&  C_on_disk){
      k = rows_C;
    }
    if (k != 0){
      F_DGEMM("n","t",&n,&m,&k,&factor,&(C_matrix[(B_on_disk ? offset : 0)][0]),&n,&(B_matrix[(C_on_disk ? offset : 0)][0]),&m,&beta,&(A_matrix[0][0]),&n);
    }
    DEBUGGING(5,
      for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
          for(int l=0;l<k;l++)
            A_matrix[i][j]+=factor*B_matrix[(C_on_disk ? l+offset : l)][i]*C_matrix[(B_on_disk ? l+offset : l)][j];
    );
  }
  // CASE B
  //
  //            ------------------------
  //           |         ---------------|----
  //           |        |               |    |
  //   A[A_l][A_r] = B[B_i][B_c] 1@2 C[C_c][C_i]
  //      |                  |
  //       ------------------
  //  Solution:
  //  1) Transpose B
  //  2) Perform contraction of the indices ???
  //  3) Transpose the result and store in A ??
  if(operation=="1@2"){
    double beta = 1.0;
    int m = rows_A;
    int n = rows_C;
    int k;
    if(m*n == 0) return;
    k = rows_B;
    if (k != 0){
      F_DGEMM("t","t",&n,&m,&k,&factor,&(C_matrix[0][(B_on_disk ? offset : 0)]),&cols_C,&(B_matrix[0][0]),&m,&beta,&(A_matrix[0][(C_on_disk ? offset : 0)]),&cols_A);
    }
    DEBUGGING(5,
      for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
          for(int l=0;l<k;l++)
            A_matrix[i][(C_on_disk ? offset+j : j)]+=factor*B_matrix[l][i]*C_matrix[j][(B_on_disk ? offset+l : l)];
    );

  }
  // CASE C,
  //
  //       -------------       ---------
  //      |             |     |         |
  //   A[A_l][A_r] = B[B_i][B_c] 2@1 C[C_c][C_i]
  //           |                             |
  //            -----------------------------
  //  Solution:
  //  Perform contraction of the indices by calling BLAS with swapped C and B
  if(operation=="2@1"){
    int m = rows_B;
    int n = cols_A;
    int k = rows_C;
    if(m*n == 0) return;
    double beta = 1.0;
    if (k != 0){
      F_DGEMM("n","n",&n,&m,&k,&factor,&(C_matrix[0][0]),&n,&(B_matrix[0][(C_on_disk ? offset : 0)]),&cols_B,&beta,&(A_matrix[(B_on_disk ? offset : 0)][0]),&n);
    }
    DEBUGGING(5,
      for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
          for(int l=0;l<k;l++)
            A_matrix[(B_on_disk ? offset+i : i)][j]+=factor*B_matrix[i][(C_on_disk ? l+offset : l)]*C_matrix[l][j];
    );

  }
  // CASE D, best case scenario
  //
  //       -------------       ---------------
  //      |             |     |               |
  //   A[A_l][A_r] = B[B_i][B_c] 2@2 C[C_c][C_i]
  //           |                        |
  //            ------------------------
  //  Solution:
  //  1) Perform contraction of the indices
  //  2) Store the result in A
  if(operation=="2@2"){
    int m = rows_B;
    int n = rows_C;
    int k = cols_C;
    if(m*n == 0) return;
    double beta = 1.0;
    if (k != 0){
      F_DGEMM("t","n",&n,&m,&k,&factor,&(C_matrix[0][0]),&k,&(B_matrix[0][0]),&k,&beta,&(A_matrix[(B_on_disk ? offset : 0)][(C_on_disk ? offset : 0)]),&cols_A);
    }
    DEBUGGING(5,
      for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
          for(int l=0;l<k;l++)
            A_matrix[(B_on_disk ? offset+i : i)][(C_on_disk ? offset+j : j)]+=factor*B_matrix[i][l]*C_matrix[j][l];
    );

  }
}



// Just a copy of the unmodified version
/*
void CCOperation::setup_contractions()
{
  Timer PartA;
  bool need_sort = false;
  if(reindexing.size()>0)
    need_sort = true;
  CCIndex*  T_left  =A_Matrix->get_left();
  CCIndex*  T_right =A_Matrix->get_right();
  double*** T_matrix=A_Matrix->get_matrix();

  // Determine the target indexing if sorting is needed
  if(need_sort){
    if(operation[0]=='1'){
      T_left=B_Matrix->get_right();
    }  else {
      T_left=B_Matrix->get_left();
    }
    if(operation[2]=='1'){
      T_right=C_Matrix->get_right();
    }  else {
      T_right=C_Matrix->get_left();
    }
//      T_matrix = blas->get_sortmap(T_left,T_right,thread_id);
    int T_matrix_offset = 0;
    T_matrix = new double**[moinfo->get_nirreps()];
    for(int irrep=0;irrep<moinfo->get_nirreps();irrep++){
      T_matrix[h] = new double*[T_left->get_pairpi(h)];
      int block_size = T_left->get_pairpi(h)*T_right->get_pairpi(h);
        for(int i=0;i<T_left->get_pairpi(h);i++){
          T_matrix[h][i]=&(thread_work[T_matrix_offset+i*T_right->get_pairpi(h)]);
        }
        T_matrix_offset+=block_size;
    }
    if(T_matrix_offset>0)
      zero_arr(&(thread_work[0]),T_matrix_offset);
  }

  PartA_timing += PartA.get();
  Timer PartB;

  double**  A_matrix;
  double**  B_matrix;
  double**  C_matrix;

  char      label[80];

  bool      B_on_disk;
  int       B_in_core_rows;
  int       B_in_core_cols;

  bool      C_on_disk;
  int       C_in_core_rows;
  int       C_in_core_cols;

  if(!B_Matrix->get_on_disk() && !C_Matrix->get_on_disk()){
    B_on_disk = false;
    C_on_disk = false;
  }
  if(B_Matrix->get_on_disk() && !C_Matrix->get_on_disk()){
    B_on_disk = true;
    C_on_disk = false;
  }
  if(!B_Matrix->get_on_disk() && C_Matrix->get_on_disk()){
    B_on_disk = false;
    C_on_disk = true;
  }
  if(B_Matrix->get_on_disk() && C_Matrix->get_on_disk()){
    print_developing("Both on disk matrix multiply",__FILE__,__LINE__);
  }
  for(int irrep=0;irrep<moinfo->get_nirreps();irrep++){
    bool done       = false;
    int  first_block= -1;
    int  last_block = 0;
    int  rows_B,rows_C,offset;
    while(!done){
      //////////////////////////////////////////////////////////
      // Case I. A,B,C are in core. Perform direct a BLAS call
      // in one pass.
      //////////////////////////////////////////////////////////
      if(!B_on_disk && !C_on_disk){
        done = true;
        A_matrix = T_matrix[h];
        B_matrix = B_Matrix->get_matrix()[h];
        C_matrix = C_Matrix->get_matrix()[h];
        rows_B = B_Matrix->get_left_pairpi(h);
        rows_C = C_Matrix->get_left_pairpi(h);
        offset = 0;
      }
      //////////////////////////////////////////////////////////
      // Case II. A,C are in core. B is on disk. Perform several
      // BLAS calls.
      //////////////////////////////////////////////////////////
      if(B_on_disk && !C_on_disk){
        int total_blocks = B_Matrix->get_blocks_per_irrep(h);
        // Determine how many blocks can be held in core
        double free_memory = memory_manager->get_free_memory();
        first_block = last_block;
        last_block  = first_block;
        for(int i=first_block;i<total_blocks;i++){
          if(free_memory > B_Matrix->get_memory_per_block(irrep,i)){
            free_memory -= B_Matrix->get_memory_per_block(irrep,i);
            last_block++;
          }
        }
        if(last_block == total_blocks)
          done = true;

        // Assign pointers to in core matrices
        A_matrix = T_matrix[h];
        C_matrix = C_Matrix->get_matrix()[h];

        // Compute the size of the sub-matrix B
        B_in_core_cols = B_Matrix->get_right_pairpi(h);
        B_in_core_rows = 0;
        for(int i=first_block;i<last_block;i++)
          B_in_core_rows += B_Matrix->get_rows_per_block(irrep,i);
        allocate2(double,B_matrix,B_in_core_rows,B_in_core_cols);

        // Read the sub-matrix B
        int block_first_row = 0;
        for(int i=first_block;i<last_block;i++){
           int rows_per_block = B_Matrix->get_rows_per_block(irrep,i);
           sprintf(label,"%s_%d_%d",B_Matrix->get_label().c_str(),irrep,i);
           psio_read_entry(PSIF_PSIMRCC_INTEGRALS,label,(char*)&B_matrix[block_first_row][0],rows_per_block*B_in_core_cols*sizeof(double));
           block_first_row += rows_per_block;
        }
        rows_B = B_in_core_rows;
        rows_C = C_Matrix->get_left_pairpi(h);
        offset = B_Matrix->get_first_row_per_block(irrep,first_block);
      }
      //////////////////////////////////////////////////////////
      // Case III. A,B are in core. C is on disk. Perform several
      // BLAS calls.
      //////////////////////////////////////////////////////////
      if(!B_on_disk && C_on_disk){
        int total_blocks = C_Matrix->get_blocks_per_irrep(h);
        // Determine how many blocks can be held in core
        double free_memory = memory_manager->get_free_memory();
        first_block = last_block;
        last_block  = first_block;
        for(int i=first_block;i<total_blocks;i++){
          if(free_memory > C_Matrix->get_memory_per_block(irrep,i)){
            free_memory -= C_Matrix->get_memory_per_block(irrep,i);
            last_block++;
          }
        }
        if(last_block == total_blocks)
          done = true;

        // Assign pointers to in core matrices
        A_matrix = T_matrix[h];
        B_matrix = B_Matrix->get_matrix()[h];

        // Compute the size of the sub-matrix C
        C_in_core_cols = C_Matrix->get_right_pairpi(h);
        C_in_core_rows = 0;
        for(int i=first_block;i<last_block;i++)
          C_in_core_rows += C_Matrix->get_rows_per_block(irrep,i);
        allocate2(double,C_matrix,C_in_core_rows,C_in_core_cols);

        // Read the sub-matrix C
        int block_first_row = 0;
        for(int i=first_block;i<last_block;i++){
           int rows_per_block = C_Matrix->get_rows_per_block(irrep,i);
           sprintf(label,"%s_%d_%d",C_Matrix->get_label().c_str(),irrep,i);
           psio_read_entry(PSIF_PSIMRCC_INTEGRALS,label,(char*)&C_matrix[block_first_row][0],rows_per_block*C_in_core_cols*sizeof(double));
           block_first_row += rows_per_block;
        }
        rows_B = B_Matrix->get_left_pairpi(h);
        rows_C = C_in_core_rows;
        offset = C_Matrix->get_first_row_per_block(irrep,first_block);
      }

      //////////////////////////
      // Now call BLAS
      //////////////////////////
      int rows_A = T_left->get_pairpi(h);
      int cols_A = T_right->get_pairpi(h);
      int cols_B = B_Matrix->get_right_pairpi(h);
      int cols_C = C_Matrix->get_right_pairpi(h);

      // Start a timer
      Timer timer;
      // Do the job
      contract_in_core(A_matrix,B_matrix,C_matrix,B_on_disk,C_on_disk,rows_A,rows_B,rows_C,cols_A,cols_B,cols_C,offset);
      // Store the timing in moinfo
      moinfo->add_dgemm_timing(timer.get());

      /////////////////////////////////////////
      // Free the memory allocated for B or C
      /////////////////////////////////////////
      if(B_on_disk)
        free_dmatrix(B_matrix,B_in_core_rows,B_in_core_cols);
      if(C_on_disk)
        free_dmatrix(C_matrix,C_in_core_rows,C_in_core_cols);

    }  // end of while(!done)
  }  // end of for loop over irreps

  PartB_timing += PartB.get();
  Timer PartC;
  if(need_sort){
    sort(T_left,T_right,T_matrix,1.0);
    for(int irrep=0;irrep<moinfo->get_nirreps();irrep++)
//       if(T_left->get_pairpi(h)*T_right->get_pairpi(h)>0)
        delete[] T_matrix[h];
    delete[] T_matrix;
  }
  PartC_timing += PartC.get();
}

void CCOperation::contract_in_core(double** A_matrix,double** B_matrix,double** C_matrix,bool B_on_disk,bool C_on_disk,int rows_A,int rows_B,int rows_C,int cols_A,int cols_B,int cols_C,int offset)
{
  // CASE A, worst case scenario
  //
  //            -----------------------------
  //           |         ---------------     |
  //           |        |               |    |
  //   A[A_l][A_r] = B[B_i][B_c] 1@1 C[C_c][C_i]
  //      |                  |
  //       ------------------
  //  Solution:
  //  1) Transpose B,C
  //  2) Perform contraction of the indices
  //  3) Transpose the result and store in A ???
  if(operation=="1@1"){
    double beta = 1.0;
    int m = rows_A;
    int n = cols_A;
    int k;
    if(m*n == 0) return;
    if(!B_on_disk && !C_on_disk){
      k = rows_B;
    }
    if( B_on_disk && !C_on_disk){
      k = rows_B;
    }
    if(!B_on_disk &&  C_on_disk){
      k = rows_C;
    }
    DEBUGGING(5,
      for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
          for(int l=0;l<k;l++)
            A_matrix[i][j]+=factor*B_matrix[(C_on_disk ? l+offset : l)][i]*C_matrix[(B_on_disk ? l+offset : l)][j];
    }
    else if (k != 0){
      F_DGEMM("n","t",&n,&m,&k,&factor,&(C_matrix[(B_on_disk ? offset : 0)][0]),&n,&(B_matrix[(C_on_disk ? offset : 0)][0]),&m,&beta,&(A_matrix[0][0]),&n);
    }
  }
  // CASE B
  //
  //            ------------------------
  //           |         ---------------|----
  //           |        |               |    |
  //   A[A_l][A_r] = B[B_i][B_c] 1@2 C[C_c][C_i]
  //      |                  |
  //       ------------------
  //  Solution:
  //  1) Transpose B
  //  2) Perform contraction of the indices ???
  //  3) Transpose the result and store in A ??
  if(operation=="1@2"){
    double beta = 1.0;
    int m = rows_A;
    int n = rows_C;
    int k;
    if(m*n == 0) return;
    k = rows_B;
    DEBUGGING(5,
      for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
          for(int l=0;l<k;l++)
            A_matrix[i][(C_on_disk ? offset+j : j)]+=factor*B_matrix[l][i]*C_matrix[j][(B_on_disk ? offset+l : l)];
    }
    else if (k != 0){
      F_DGEMM("t","t",&n,&m,&k,&factor,&(C_matrix[0][(B_on_disk ? offset : 0)]),&cols_C,&(B_matrix[0][0]),&m,&beta,&(A_matrix[0][(C_on_disk ? offset : 0)]),&cols_A);
    }
  }
  // CASE C,
  //
  //       -------------       ---------
  //      |             |     |         |
  //   A[A_l][A_r] = B[B_i][B_c] 2@1 C[C_c][C_i]
  //           |                             |
  //            -----------------------------
  //  Solution:
  //  Perform contraction of the indices by calling BLAS with swapped C and B
  if(operation=="2@1"){
    int m = rows_B;
    int n = cols_A;
    int k = rows_C;
    if(m*n == 0) return;
    double beta = 1.0;
    DEBUGGING(5,
      for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
          for(int l=0;l<k;l++)
            A_matrix[(B_on_disk ? offset+i : i)][j]+=factor*B_matrix[i][(C_on_disk ? l+offset : l)]*C_matrix[l][j];
    }
    else if (k != 0){
      F_DGEMM("n","n",&n,&m,&k,&factor,&(C_matrix[0][0]),&n,&(B_matrix[0][(C_on_disk ? offset : 0)]),&cols_B,&beta,&(A_matrix[(B_on_disk ? offset : 0)][0]),&n);
    }
  }
  // CASE D, best case scenario
  //
  //       -------------       ---------------
  //      |             |     |               |
  //   A[A_l][A_r] = B[B_i][B_c] 2@2 C[C_c][C_i]
  //           |                        |
  //            ------------------------
  //  Solution:
  //  1) Perform contraction of the indices
  //  2) Store the result in A
  if(operation=="2@2"){
    int m = rows_B;
    int n = rows_C;
    int k = cols_C;
    if(m*n == 0) return;
    double beta = 1.0;
    DEBUGGING(5,
      for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
          for(int l=0;l<k;l++)
            A_matrix[(B_on_disk ? offset+i : i)][(C_on_disk ? offset+j : j)]+=factor*B_matrix[i][l]*C_matrix[j][l];
    }
    else if (k != 0){
      F_DGEMM("t","n",&n,&m,&k,&factor,&(C_matrix[0][0]),&k,&(B_matrix[0][0]),&k,&beta,&(A_matrix[(B_on_disk ? offset : 0)][(C_on_disk ? offset : 0)]),&cols_A);
    }
  }
}
*/


/**
 * This routine will setup the contraction and call library blas routines.
 * When there is no resorting required, this routine can process one block
 * at the time, otherwise matrix A has to be loaded in memory, after the
 * contraction is done.
 */
/*
void CCOperation::setup_contractions()
{
  Timer PartA;
  CCMatTmp AMatTmp = blas->get_MatTmp(A_Matrix,none);
  check_and_zero_target(); // TODO : this is a temporary solution!!!
  CCMatTmp BMatTmp = blas->get_MatTmp(B_Matrix,none);
  CCMatTmp CMatTmp = blas->get_MatTmp(C_Matrix,none);

  bool need_sort = false;
  if(reindexing.size()>0)
    need_sort = true;
  CCIndex*  T_left  =A_Matrix->get_left();
  CCIndex*  T_right =A_Matrix->get_right();
  double*** T_matrix=A_Matrix->get_matrix();

  // Determine the target indexing if sorting is needed
  if(need_sort){
    if(operation[0]=='1'){
      T_left=B_Matrix->get_right();
    }  else {
      T_left=B_Matrix->get_left();
    }
    if(operation[2]=='1'){
      T_right=C_Matrix->get_right();
    }  else {
      T_right=C_Matrix->get_left();
    }
//      T_matrix = blas->get_sortmap(T_left,T_right,thread_id);
    int T_matrix_offset = 0;
    T_matrix = new double**[moinfo->get_nirreps()];
    for(int irrep=0;irrep<moinfo->get_nirreps();irrep++){
      T_matrix[h] = new double*[T_left->get_pairpi(h)];
      int block_size = T_left->get_pairpi(h)*T_right->get_pairpi(h);
        for(int i=0;i<T_left->get_pairpi(h);i++){
          T_matrix[h][i]=&(local_work[T_matrix_offset+i*T_right->get_pairpi(h)]);
        }
        T_matrix_offset+=block_size;
    }
    if(T_matrix_offset>0)
      zero_arr(&(local_work[0]),T_matrix_offset);
  }

  PartA_timing += PartA.get();
  Timer PartB;

//   if(need_sort)
//     A_Matrix->allocate_memory()
//   bool      A_is_allocated = A_Matrix->is_allocated();
//   bool      B_is_allocated = B_Matrix->is_allocated();
//   bool      C_is_allocated = C_Matrix->is_allocated();

  char      label[80];
  ///////////////////////////////////
  // Multiply one irrep at time
  ///////////////////////////////////
  // TODO use tricks, i.e. when there is no reindexing you can do one irrep at time
  for(int irrep=0;irrep<moinfo->get_nirreps();irrep++){
    double** A_matrix = T_matrix[h];
    double** B_matrix = B_Matrix->get_matrix()[h];
    double** C_matrix = C_Matrix->get_matrix()[h];
    int rows_A        = T_left->get_pairpi(h);
    int cols_A        = T_right->get_pairpi(h);
    int rows_B        = B_Matrix->get_left_pairpi(h);
    int cols_B        = B_Matrix->get_right_pairpi(h);
    int cols_C        = C_Matrix->get_right_pairpi(h);
    int rows_C        = C_Matrix->get_left_pairpi(h);

    //////////////////////////
    // Now call BLAS
    //////////////////////////

    // Start a timer
    Timer timer;
    // Do the job
    contract_in_core(A_matrix,B_matrix,C_matrix,rows_A,rows_B,rows_C,cols_A,cols_B,cols_C);
    // Store the timing in moinfo
    moinfo->add_dgemm_timing(timer.get());


  }  // end of for loop over irreps

  PartB_timing += PartB.get();
  Timer PartC;
  if(need_sort){
    sort(T_left,T_right,T_matrix,1.0);
    for(int irrep=0;irrep<moinfo->get_nirreps();irrep++)
//       if(T_left->get_pairpi(h)*T_right->get_pairpi(h)>0)
        delete[] T_matrix[h];
    delete[] T_matrix;
  }
  PartC_timing += PartC.get();
}

void CCOperation::contract_in_core(double** A_matrix,double** B_matrix,double** C_matrix,int rows_A,int rows_B,int rows_C,int cols_A,int cols_B,int cols_C)
{
  // CASE A, worst case scenario
  //
  //            -----------------------------
  //           |         ---------------     |
  //           |        |               |    |
  //   A[A_l][A_r] = B[B_i][B_c] 1@1 C[C_c][C_i]
  //      |                  |
  //       ------------------
  //  Solution:
  //  1) Transpose B,C
  //  2) Perform contraction of the indices
  //  3) Transpose the result and store in A ???
  if(operation=="1@1"){
    double beta = 1.0;
    int m = rows_A;
    int n = cols_A;
    int k = rows_B;
    if(m*n == 0) return;
    DEBUGGING(5,
      for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
          for(int l=0;l<k;l++)
            A_matrix[i][j]+=factor*B_matrix[l][i]*C_matrix[l][j];
    }
    else if (k != 0){
      F_DGEMM("n","t",&n,&m,&k,&factor,&(C_matrix[0][0]),&n,&(B_matrix[0][0]),&m,&beta,&(A_matrix[0][0]),&n);
    }
  }
  // CASE B
  //
  //            ------------------------
  //           |         ---------------|----
  //           |        |               |    |
  //   A[A_l][A_r] = B[B_i][B_c] 1@2 C[C_c][C_i]
  //      |                  |
  //       ------------------
  //  Solution:
  //  1) Transpose B
  //  2) Perform contraction of the indices ???
  //  3) Transpose the result and store in A ??
  if(operation=="1@2"){
    double beta = 1.0;
    int m = rows_A;
    int n = rows_C;
    int k = rows_B;
    if(m*n == 0) return;
    DEBUGGING(5,
      for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
          for(int l=0;l<k;l++)
            A_matrix[i][j]+=factor*B_matrix[l][i]*C_matrix[j][l];
    }
    else if (k != 0){
      F_DGEMM("t","t",&n,&m,&k,&factor,&(C_matrix[0][0]),&cols_C,&(B_matrix[0][0]),&m,&beta,&(A_matrix[0][0]),&cols_A);
    }
  }
  // CASE C,
  //
  //       -------------       ---------
  //      |             |     |         |
  //   A[A_l][A_r] = B[B_i][B_c] 2@1 C[C_c][C_i]
  //           |                             |
  //            -----------------------------
  //  Solution:
  //  Perform contraction of the indices by calling BLAS with swapped C and B
  if(operation=="2@1"){
    int m = rows_B;
    int n = cols_A;
    int k = rows_C;
    if(m*n == 0) return;
    double beta = 1.0;
    DEBUGGING(5,
      for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
          for(int l=0;l<k;l++)
            A_matrix[i][j]+=factor*B_matrix[i][l]*C_matrix[l][j];
    }
    else if (k != 0){
      F_DGEMM("n","n",&n,&m,&k,&factor,&(C_matrix[0][0]),&n,&(B_matrix[0][0]),&cols_B,&beta,&(A_matrix[0][0]),&n);
    }
  }
  // CASE D, best case scenario
  //
  //       -------------       ---------------
  //      |             |     |               |
  //   A[A_l][A_r] = B[B_i][B_c] 2@2 C[C_c][C_i]
  //           |                        |
  //            ------------------------
  //  Solution:
  //  1) Perform contraction of the indices
  //  2) Store the result in A
  if(operation=="2@2"){
    int m = rows_B;
    int n = rows_C;
    int k = cols_C;
    if(m*n == 0) return;
    double beta = 1.0;
    DEBUGGING(5,
      for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
          for(int l=0;l<k;l++)
            A_matrix[i][j]+=factor*B_matrix[i][l]*C_matrix[j][l];
    }
    else if (k != 0){
      F_DGEMM("t","n",&n,&m,&k,&factor,&(C_matrix[0][0]),&k,&(B_matrix[0][0]),&k,&beta,&(A_matrix[0][0]),&cols_A);
    }
  }
}

*/


}} /* End Namespaces */
