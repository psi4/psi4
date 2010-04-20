#if 0
/***************************************************************************
 *   Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *   frank@ccc.uga.edu
 *   SR/MRCC Code
 ***************************************************************************/

#include <libmoinfo/libmoinfo.h>
#include "blas.h"
#include <libutil/libutil.h>
#include <pthread.h>

namespace psi{ namespace psimrcc{


using namespace std;

extern FILE *infile, *outfile;



void CCBLAS::Operation::reindex(double constant)
{
  double*** A_matrix = A_Matrix->get_matrix();
  double*** B_matrix = B_Matrix->get_matrix();

  // TBM : make this more general
  // Determine wether any of the indices needs to be expanded
  // Example A[vv][vv] <- B[vv][v>v]
  bool expand_left  = false;
  bool expand_right = false;
  if(B_Matrix->get_left()->get_greater_than() && (!A_Matrix->get_left()->get_greater_than()))
    expand_left = true;
  if(B_Matrix->get_right()->get_greater_than() && (!A_Matrix->get_right()->get_greater_than()))
    expand_right = true;

  if(reindexing.size()==2){
    // Setup reindexing_array
    // This assumes that the reindexing starts from 1 !!! This can cost you an headache
    short* reindexing_array = new short[2];
    if(reindexing.size()==0)
      reindexing="12";
    for(int i = 0; i< reindexing.size(); i++)
      reindexing_array[i]=(int)ToDouble(reindexing.substr(i,1))-1;

    if((!expand_left) && (!expand_right)){
      // Do not expand
      short* pq   = new short[2];
      for(int n=0;n<moinfo->get_nirreps();n++){
        for(int i = 0;i<A_Matrix->get_left_pairpi(n);i++){
          for(int j = 0;j<A_Matrix->get_right_pairpi(n);j++){
            // Get the pq indices  [v][o],[vo] - > v, o
            A_Matrix->get_two_indices(pq,n,i,j);
            A_matrix[n][i][j]+=constant*B_Matrix->get_two_address_element(pq[reindexing_array[0]],pq[reindexing_array[1]]);
          }
        }
      }
      delete[] pq;
    }else{
      print_developing("Reindexing algorithm with 2-index expansion",__FILE__,__LINE__);
    }
  }else{
    /* ==== Use the PThreads library to parallelize the 4-index resorting ==== */
    if(reindexing.size()==0)
      reindexing="1234";
    data = new data_t[moinfo->get_num_threads()];
    for(int i=0; i<moinfo->get_num_threads(); i++){
      data[i].num_threads     = moinfo->get_num_threads();
      data[i].A_Matrix        = A_Matrix;
      data[i].B_Matrix        = B_Matrix;
      data[i].constant        = constant;
      data[i].num_indices     = reindexing.size();
      data[i].reindexing_array = new short[reindexing.size()];
      for(int j = 0; j< reindexing.size(); j++)
        data[i].reindexing_array[j]=(int)ToDouble(reindexing.substr(j,1))-1;
        // This assumes that the reindexing starts from 1 !!! This can cost you an headache
    }

    // Allocate memory for the threads
    pthread_t *thread_id = new pthread_t[moinfo->get_num_threads()];
    // Now we create data arrays and start the threads
    for(int i=0;i<moinfo->get_num_threads();i++){
      if(!expand_left && !expand_right)
        pthread_create(&(thread_id[i]),NULL,reindex_thread,(void*)i);
      else if(!expand_left && expand_right)
        pthread_create(&(thread_id[i]),NULL,reindex_thread_right_expand,(void*)i);
      else if(expand_left){
        print_developing("Reindexing algorithm with 4-index left expansion",__FILE__,__LINE__);
      }
    }

    // Now we wait for all threads to finish before exiting
    for(int i=0;i<moinfo->get_num_threads();i++){
      pthread_join(thread_id[i],NULL);
      delete[] data[i].reindexing_array;
    }
    delete[] data;
    delete[] thread_id;
  }
}

void *CCBLAS::Operation::reindex_thread(void * my_id)
{
  int thread_id           = (int) my_id;
  short *reindexing_array = data[thread_id].reindexing_array;
  CCMatrix* A_Matrix       = data[thread_id].A_Matrix;
  CCMatrix* B_Matrix       = data[thread_id].B_Matrix;
  double constant         = data[thread_id].constant;
  int n_threads           = data[thread_id].num_threads;
  int num_indices         = data[thread_id].num_indices;

  CCIndex* A_left   = A_Matrix->get_left();
  CCIndex* A_right  = A_Matrix->get_right();
  CCIndex* B_left   = B_Matrix->get_left();
  CCIndex* B_right  = B_Matrix->get_right();
  double*** A_matrix= A_Matrix->get_matrix();
  double*** B_matrix= B_Matrix->get_matrix();

  int A_left_nelements = A_left->get_nelements();
  int B_left_nelements = B_left->get_nelements();

  int*    A_left_pairpi   = A_left->get_pairpi_thread(thread_id);
  int*    A_right_pairpi  = A_right->get_pairpi_thread(thread_id);
  int*    A_left_first    = A_left->get_first_thread(thread_id);
  int*    A_right_first   = A_right->get_first_thread(thread_id);
  short** A_left_tuples   = A_left->get_tuples_thread(thread_id);
  short** A_right_tuples  = A_right->get_tuples_thread(thread_id);

  int*    B_left_first    = B_left->get_first_thread(thread_id);
  int*    B_right_first   = B_right->get_first_thread(thread_id);

  int b_irrep,b_left,b_right;
  if(num_indices==4){
    short* pqrs = new short[4];
    // A[x][xxx] <- B[x][xxx]
    if((A_left_nelements == 1) && (B_left_nelements == 1)){
      int*    B_left_one_index_to_irrep    = B_left->get_one_index_to_irrep_thread(thread_id);
      int*    B_left_one_index_to_tuple    = B_left->get_one_index_to_tuple_thread(thread_id);
      int***  B_right_three_index_to_tuple = B_right->get_three_index_to_tuple_thread(thread_id);
      for(int n=0;n<moinfo->get_nirreps();n++){
        for(int i = 0;i<A_left_pairpi[n];i++){
          //Check that this should be handled by this thread
          pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
          if (thread_id == i%n_threads)
            for(int j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[1]=A_right_tuples[j+A_right_first[n]][0];
              pqrs[2]=A_right_tuples[j+A_right_first[n]][1];
              pqrs[3]=A_right_tuples[j+A_right_first[n]][2];
              b_irrep = B_left_one_index_to_irrep[pqrs[reindexing_array[0]]];
              b_left  = B_left_one_index_to_tuple[pqrs[reindexing_array[0]]];
              b_right = B_right_three_index_to_tuple[pqrs[reindexing_array[1]]]
                                                    [pqrs[reindexing_array[2]]]
                                                    [pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*B_matrix[b_irrep][b_left][b_right];
            }
        }
      }
    }

    // A[xx][xx] <- B[x][xxx]
    if((A_left_nelements == 2) && (B_left_nelements == 1)){
      int*    B_left_one_index_to_irrep    = B_left->get_one_index_to_irrep_thread(thread_id);
      int*    B_left_one_index_to_tuple    = B_left->get_one_index_to_tuple_thread(thread_id);
      int***  B_right_three_index_to_tuple = B_right->get_three_index_to_tuple_thread(thread_id);
      for(int n=0;n<moinfo->get_nirreps();n++){
        for(int i = 0;i<A_left_pairpi[n];i++){
          //Check that this should be handled by this thread
          pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
          pqrs[1]=A_left_tuples[i+A_left_first[n]][1];
          if (thread_id == i%n_threads)
            for(int j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[2]=A_right_tuples[j+A_right_first[n]][0];
              pqrs[3]=A_right_tuples[j+A_right_first[n]][1];
              b_irrep = B_left_one_index_to_irrep[pqrs[reindexing_array[0]]];
              b_left  = B_left_one_index_to_tuple[pqrs[reindexing_array[0]]];
              b_right = B_right_three_index_to_tuple[pqrs[reindexing_array[1]]]
                                                    [pqrs[reindexing_array[2]]]
                                                    [pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*B_matrix[b_irrep][b_left][b_right];
            }
        }
      }
    }

    // A[xxx][x] <- B[x][xxx]
    if((A_left_nelements == 3) && (B_left_nelements == 1)){
      int*    B_left_one_index_to_irrep    = B_left->get_one_index_to_irrep_thread(thread_id);
      int*    B_left_one_index_to_tuple    = B_left->get_one_index_to_tuple_thread(thread_id);
      int***  B_right_three_index_to_tuple = B_right->get_three_index_to_tuple_thread(thread_id);
      for(int n=0;n<moinfo->get_nirreps();n++){
        for(int i = 0;i<A_left_pairpi[n];i++){
          //Check that this should be handled by this thread
          pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
          pqrs[1]=A_left_tuples[i+A_left_first[n]][1];
          pqrs[2]=A_left_tuples[i+A_left_first[n]][2];
          if (thread_id == i%n_threads)
            for(int j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[3]=A_right_tuples[j+A_right_first[n]][0];
              b_irrep = B_left_one_index_to_irrep[pqrs[reindexing_array[0]]];
              b_left  = B_left_one_index_to_tuple[pqrs[reindexing_array[0]]];
              b_right = B_right_three_index_to_tuple[pqrs[reindexing_array[1]]]
                                                    [pqrs[reindexing_array[2]]]
                                                    [pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*B_matrix[b_irrep][b_left][b_right];
            }
        }
      }
    }

    // A[x][xxx] <- B[xx][xx]
    if((A_left_nelements == 1) && (B_left_nelements == 2)){
      int** B_left_two_index_to_irrep    = B_left->get_two_index_to_irrep_thread(thread_id);
      int** B_left_two_index_to_tuple    = B_left->get_two_index_to_tuple_thread(thread_id);
      int** B_right_two_index_to_tuple   = B_right->get_two_index_to_tuple_thread(thread_id);
      for(int n=0;n<moinfo->get_nirreps();n++){
        for(int i = 0;i<A_left_pairpi[n];i++){
          //Check that this should be handled by this thread
          pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
          if (thread_id == i%n_threads)
            for(int j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[1]=A_right_tuples[j+A_right_first[n]][0];
              pqrs[2]=A_right_tuples[j+A_right_first[n]][1];
              pqrs[3]=A_right_tuples[j+A_right_first[n]][2];
              b_irrep = B_left_two_index_to_irrep[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              b_left  = B_left_two_index_to_tuple[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              b_right = B_right_two_index_to_tuple[pqrs[reindexing_array[2]]][pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*B_matrix[b_irrep][b_left][b_right];
            }
        }
      }
    }

    // A[xx][xx] <- B[xx][xx]
    if((A_left_nelements == 2) && (B_left_nelements == 2)){
      int** B_left_two_index_to_irrep    = B_left->get_two_index_to_irrep_thread(thread_id);
      int** B_left_two_index_to_tuple    = B_left->get_two_index_to_tuple_thread(thread_id);
      int** B_right_two_index_to_tuple   = B_right->get_two_index_to_tuple_thread(thread_id);
      for(int n=0;n<moinfo->get_nirreps();n++){
        for(int i = 0;i<A_left_pairpi[n];i++){
          //Check that this should be handled by this thread
          pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
          pqrs[1]=A_left_tuples[i+A_left_first[n]][1];
          if (thread_id == i%n_threads)
            for(int j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[2]=A_right_tuples[j+A_right_first[n]][0];
              pqrs[3]=A_right_tuples[j+A_right_first[n]][1];
              b_irrep = B_left_two_index_to_irrep[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              b_left  = B_left_two_index_to_tuple[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              b_right = B_right_two_index_to_tuple[pqrs[reindexing_array[2]]][pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*B_matrix[b_irrep][b_left][b_right];
            }
        }
      }
    }

    // A[xxx][x] <- B[xx][xx]
    if((A_left_nelements == 3) && (B_left_nelements == 2)){
      int** B_left_two_index_to_irrep    = B_left->get_two_index_to_irrep_thread(thread_id);
      int** B_left_two_index_to_tuple    = B_left->get_two_index_to_tuple_thread(thread_id);
      int** B_right_two_index_to_tuple   = B_right->get_two_index_to_tuple_thread(thread_id);
      for(int n=0;n<moinfo->get_nirreps();n++){
        for(int i = 0;i<A_left_pairpi[n];i++){
          //Check that this should be handled by this thread
          pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
          pqrs[1]=A_left_tuples[i+A_left_first[n]][1];
          pqrs[2]=A_left_tuples[i+A_left_first[n]][2];
          if (thread_id == i%n_threads)
            for(int j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[3]=A_right_tuples[j+A_right_first[n]][0];
              b_irrep = B_left_two_index_to_irrep[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              b_left  = B_left_two_index_to_tuple[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              b_right = B_right_two_index_to_tuple[pqrs[reindexing_array[2]]][pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*B_matrix[b_irrep][b_left][b_right];
            }
        }
      }
    }

    // A[x][xxx] <- B[xxx][x]
    if((A_left_nelements == 1) && (B_left_nelements == 3)){
      int*** B_left_three_index_to_tuple  = B_left->get_three_index_to_tuple_thread(thread_id);
      int*   B_right_one_index_to_irrep   = B_right->get_one_index_to_irrep_thread(thread_id);
      int*   B_right_one_index_to_tuple   = B_right->get_one_index_to_tuple_thread(thread_id);
      for(int n=0;n<moinfo->get_nirreps();n++){
        for(int i = 0;i<A_left_pairpi[n];i++){
          //Check that this should be handled by this thread
          pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
          if (thread_id == i%n_threads)
            for(int j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[1]=A_right_tuples[j+A_right_first[n]][0];
              pqrs[2]=A_right_tuples[j+A_right_first[n]][1];
              pqrs[3]=A_right_tuples[j+A_right_first[n]][2];
              b_irrep = B_right_one_index_to_irrep[pqrs[reindexing_array[3]]];
              b_left  = B_left_three_index_to_tuple[pqrs[reindexing_array[0]]]
                                                   [pqrs[reindexing_array[1]]]
                                                   [pqrs[reindexing_array[2]]];
              b_right = B_right_one_index_to_tuple[pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*B_matrix[b_irrep][b_left][b_right];
            }
        }
      }
    }

    // A[xx][xx] <- B[xxx][x]
    if((A_left_nelements == 2) && (B_left_nelements == 3)){
      int*** B_left_three_index_to_tuple  = B_left->get_three_index_to_tuple_thread(thread_id);
      int*   B_right_one_index_to_irrep   = B_right->get_one_index_to_irrep_thread(thread_id);
      int*   B_right_one_index_to_tuple   = B_right->get_one_index_to_tuple_thread(thread_id);
      for(int n=0;n<moinfo->get_nirreps();n++){
        for(int i = 0;i<A_left_pairpi[n];i++){
          //Check that this should be handled by this thread
          pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
          pqrs[1]=A_left_tuples[i+A_left_first[n]][1];
          if (thread_id == i%n_threads)
            for(int j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[2]=A_right_tuples[j+A_right_first[n]][0];
              pqrs[3]=A_right_tuples[j+A_right_first[n]][1];
              b_irrep = B_right_one_index_to_irrep[pqrs[reindexing_array[3]]];
              b_left  = B_left_three_index_to_tuple[pqrs[reindexing_array[0]]]
                                                   [pqrs[reindexing_array[1]]]
                                                   [pqrs[reindexing_array[2]]];
              b_right = B_right_one_index_to_tuple[pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*B_matrix[b_irrep][b_left][b_right];
            }
        }
      }
    }

    // A[xxx][x] <- B[xxx][x]
    if((A_left_nelements == 3) && (B_left_nelements == 3)){
      int*** B_left_three_index_to_tuple  = B_left->get_three_index_to_tuple_thread(thread_id);
      int*   B_right_one_index_to_irrep   = B_right->get_one_index_to_irrep_thread(thread_id);
      int*   B_right_one_index_to_tuple   = B_right->get_one_index_to_tuple_thread(thread_id);
      for(int n=0;n<moinfo->get_nirreps();n++){
        for(int i = 0;i<A_left_pairpi[n];i++){
          //Check that this should be handled by this thread
          pqrs[0]=A_left_tuples[i+A_left_first[n]][0];
          pqrs[1]=A_left_tuples[i+A_left_first[n]][1];
          pqrs[2]=A_left_tuples[i+A_left_first[n]][2];
          if (thread_id == i%n_threads)
            for(int j = 0;j<A_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[3]=A_right_tuples[j+A_right_first[n]][0];
              b_irrep = B_right_one_index_to_irrep[pqrs[reindexing_array[3]]];
              b_left  = B_left_three_index_to_tuple[pqrs[reindexing_array[0]]]
                                                   [pqrs[reindexing_array[1]]]
                                                   [pqrs[reindexing_array[2]]];
              b_right = B_right_one_index_to_tuple[pqrs[reindexing_array[3]]];
              A_matrix[n][i][j]+=constant*B_matrix[b_irrep][b_left][b_right];
            }
        }
      }
    }

    delete[] pqrs;
  }
  return NULL;
}



void *CCBLAS::Operation::reindex_thread_right_expand(void * my_id)
{
  int thread_id           = (int) my_id;
  short *reindexing_array = data[thread_id].reindexing_array;
  CCMatrix* A_Matrix      = data[thread_id].A_Matrix;
  CCMatrix* B_Matrix      = data[thread_id].B_Matrix;
  double constant         = data[thread_id].constant;
  int n_threads           = data[thread_id].num_threads;
  int num_indices         = data[thread_id].num_indices;

  CCIndex* A_left   = A_Matrix->get_left();
  CCIndex* A_right  = A_Matrix->get_right();
  CCIndex* B_left   = B_Matrix->get_left();
  CCIndex* B_right  = B_Matrix->get_right();
  double*** A_matrix= A_Matrix->get_matrix();
  double*** B_matrix= B_Matrix->get_matrix();

  int A_left_nelements = A_left->get_nelements();
  int B_left_nelements = B_left->get_nelements();


  // Matrix A data
  int*    A_left_first    = A_left->get_first_thread(thread_id);
  int*    A_right_first   = A_right->get_first_thread(thread_id);

  // Matrix B data
  int*    B_left_pairpi   = B_left->get_pairpi_thread(thread_id);
  int*    B_right_pairpi  = B_right->get_pairpi_thread(thread_id);
  int*    B_left_first    = B_left->get_first_thread(thread_id);
  int*    B_right_first   = B_right->get_first_thread(thread_id);
  short** B_left_tuples   = B_left->get_tuples_thread(thread_id);
  short** B_right_tuples  = B_right->get_tuples_thread(thread_id);

  int     a_irrep,a_left,a_right;
  short   swap_temp;
  double  sym_constant = constant * B_Matrix->get_symmetry(); 

  if(num_indices==4){
    short* pqrs = new short[4];
    // A[x][xxx] <- B[x][xxx]
    if((A_left_nelements == 1) && (B_left_nelements == 1)){
      print_developing("A[x][xxx] <- B[x][xxx] with expansion",__FILE__,__LINE__);
    }

    // A[xx][xx] <- B[x][xxx]
    if((A_left_nelements == 2) && (B_left_nelements == 1)){
      print_developing("A[xx][xx] <- B[x][xxx] with expansion",__FILE__,__LINE__);
    }

    // A[xxx][x] <- B[x][xxx]
    if((A_left_nelements == 3) && (B_left_nelements == 1)){
      print_developing("A[xxx][x] <- B[x][xxx] with expansion",__FILE__,__LINE__);
    }

    // A[x][xxx] <- B[xx][xx]
    if((A_left_nelements == 1) && (B_left_nelements == 2)){
      print_developing("A[x][xxx] <- B[xx][xx] with expansion",__FILE__,__LINE__);
    }

    // A[xx][xx] <- B[xx][xx]
    if((A_left_nelements == 2) && (B_left_nelements == 2)){
      int** A_left_two_index_to_irrep    = A_left->get_two_index_to_irrep_thread(thread_id);
      int** A_left_two_index_to_tuple    = A_left->get_two_index_to_tuple_thread(thread_id);
      int** A_right_two_index_to_tuple   = A_right->get_two_index_to_tuple_thread(thread_id);
      for(int n=0;n<moinfo->get_nirreps();n++){
        for(int i = 0;i<B_left_pairpi[n];i++){
          //Check that this should be handled by this thread
          pqrs[0]=B_left_tuples[i+B_left_first[n]][0];
          pqrs[1]=B_left_tuples[i+B_left_first[n]][1];
          if (thread_id == i%n_threads)
            for(int j = 0;j<B_right_pairpi[n];j++){
              // Get the pqrs indices
              pqrs[2]=B_right_tuples[j+B_right_first[n]][0];
              pqrs[3]=B_right_tuples[j+B_right_first[n]][1];
              a_irrep = A_left_two_index_to_irrep[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              a_left  = A_left_two_index_to_tuple[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              a_right = A_right_two_index_to_tuple[pqrs[reindexing_array[2]]][pqrs[reindexing_array[3]]];
              A_matrix[a_irrep][a_left][a_right]+=constant*B_matrix[n][i][j];
              // Add the expandend term with the correct symmetry
              swap_temp =pqrs[2];
              pqrs[2] = pqrs[3];
              pqrs[3] = swap_temp;
              a_irrep = A_left_two_index_to_irrep[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              a_left  = A_left_two_index_to_tuple[pqrs[reindexing_array[0]]][pqrs[reindexing_array[1]]];
              a_right = A_right_two_index_to_tuple[pqrs[reindexing_array[2]]][pqrs[reindexing_array[3]]];
              A_matrix[a_irrep][a_left][a_right]+=sym_constant*B_matrix[n][i][j];
            }
        }
      }
    }

    // A[xxx][x] <- B[xx][xx]
    if((A_left_nelements == 3) && (B_left_nelements == 2)){
      print_developing("A[xxx][x] <- B[xx][xx] with expansion",__FILE__,__LINE__);
    }

    // A[x][xxx] <- B[xxx][x]
    if((A_left_nelements == 1) && (B_left_nelements == 3)){
      print_developing("A[x][xxx] <- B[xxx][x] with expansion",__FILE__,__LINE__);
    }

    // A[xx][xx] <- B[xxx][x]
    if((A_left_nelements == 2) && (B_left_nelements == 3)){
      print_developing("A[xx][xx] <- B[xxx][x] with expansion",__FILE__,__LINE__);
    }

    // A[xxx][x] <- B[xxx][x]
    if((A_left_nelements == 3) && (B_left_nelements == 3)){
      print_developing("A[xxx][x] <- B[xxx][x] with expansion",__FILE__,__LINE__);
    }

    delete[] pqrs;
  }
  return NULL;
}

}} /* End Namespaces */

#endif
