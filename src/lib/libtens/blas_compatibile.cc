#include <cstdlib>
#include <cstdio>

#include <libutil/libutil.h>

#include "blas.h"
#include "index.h"
#include "matrix.h"

extern FILE *outfile;

namespace psi{ namespace psimrcc{

using namespace std;

bool CCOperation::compatible_dot()
{
  //
  //  Here we are doing the following check:
  //            --------------
  //           |              |
  //   A = B[B_l][B_r] . C[C_l][C_r]
  //               |               |
  //                ---------------
  //
  bool same = false;

  int A_left   = A_Matrix->get_left()->get_ntuples();
  int A_right  = A_Matrix->get_right()->get_ntuples();
  int B_left   = B_Matrix->get_left()->get_ntuples();
  int B_right  = B_Matrix->get_right()->get_ntuples();
  int C_left   = C_Matrix->get_left()->get_ntuples();
  int C_right  = C_Matrix->get_right()->get_ntuples();

  if((A_left==1) && (B_left==C_left) && (A_right==1) && (B_right==C_right))
    same = true;
  if(!same){
    fprintf(outfile,"\n\nSolve couldn't perform the operation ");
    print_operation();
    fflush(outfile);
    exit(1);
  }
  return(same);
}

bool CCOperation::compatible_element_by_element()
{
  //
  //  Here we are doing the following check:
  //            ---------------------------
  //           |             |             |
  //   A[A_l][A_r] = B[B_l][B_r] / C[C_l][C_r]
  //      |             |             |
  //       ---------------------------
  //
  bool same = false;

  int A_left   = A_Matrix->get_left()->get_ntuples();
  int A_right  = A_Matrix->get_right()->get_ntuples();
  int B_left   = B_Matrix->get_left()->get_ntuples();
  int B_right  = B_Matrix->get_right()->get_ntuples();

  // First case: We are comparing A and B
  if(C_Matrix==NULL){
    if((A_left==B_left) && (A_right==B_right))
      same = true;
  }else{
  // Second case: We are comparing A, B, and C
    int C_left   = C_Matrix->get_left()->get_ntuples();
    int C_right  = C_Matrix->get_right()->get_ntuples();
    if((A_left==B_left) && (B_left==C_left) && (A_right==B_right) && (B_right==C_right))
      same = true;
    if((B_left!=C_left) || (B_right!=C_right)){
      fprintf(outfile,"\n\nSolve couldn't perform the operation ");
      print_operation();
      fflush(outfile);
      exit(1);
    }
  }
  return(same);
}

bool CCOperation::compatible_contract()
{
  //
  //  Here we are doing the following check:
  //            -----------------------------
  //           |                             |
  //   A[A_l][A_r] = B[B_i][B_c] 2@1 C[C_c][C_i]
  //      |             |    |          |
  //       -------------      ----------
  //
  bool same = false;

  int A_left   = A_Matrix->get_left()->get_ntuples();
  int A_right  = A_Matrix->get_right()->get_ntuples();
  int B_contracted;
  int B_index;
  int C_contracted;
  int C_index;

  if(operation[0]=='1'){
    B_contracted   = B_Matrix->get_left()->get_ntuples();
    B_index        = B_Matrix->get_right()->get_ntuples();
  } else{
    B_contracted   = B_Matrix->get_right()->get_ntuples();
    B_index        = B_Matrix->get_left()->get_ntuples();
  }
  if(operation[2]=='1'){
    C_contracted   = C_Matrix->get_left()->get_ntuples();
    C_index        = C_Matrix->get_right()->get_ntuples();
  } else{
    C_contracted   = C_Matrix->get_right()->get_ntuples();
    C_index        = C_Matrix->get_left()->get_ntuples();
  }

  if((B_contracted==C_contracted) && (A_left==B_index) && (A_right==C_index))
    same = true;
  if(B_contracted!=C_contracted){
    fprintf(outfile,"\n\nSolve couldn't perform the operation ");
    print_operation();
    fflush(outfile);
    exit(1);
  }
  return(same);
}

}} /* End Namespaces */
