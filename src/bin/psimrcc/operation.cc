#include <libutil/libutil.h>
#include <cstdio>

#include "blas.h"
#include "debugging.h"
#include "matrix.h"

namespace psi{
    extern FILE *outfile;
    namespace psimrcc{

double* CCOperation::local_work = NULL;
double* CCOperation::out_of_core_buffer = NULL;
double CCOperation::zero_timing=0.0;
double CCOperation::numerical_timing=0.0;
double CCOperation::contract_timing=0.0;
double CCOperation::tensor_timing=0.0;
double CCOperation::dot_timing=0.0;
double CCOperation::plus_timing=0.0;
double CCOperation::product_timing=0.0;
double CCOperation::division_timing=0.0;
double CCOperation::sort_timing=0.0;
double CCOperation::PartA_timing=0.0;
double CCOperation::PartB_timing=0.0;
double CCOperation::PartC_timing=0.0;

CCOperation::CCOperation(double in_factor,std::string in_assignment,
            std::string in_reindexing,std::string in_operation,
            CCMatrix* in_A_Matrix, CCMatrix* in_B_Matrix, CCMatrix* in_C_Matrix,double* work,double* buffer)
: factor(in_factor), assignment(in_assignment), reindexing(in_reindexing),operation(in_operation),
A_Matrix(in_A_Matrix),B_Matrix(in_B_Matrix),C_Matrix(in_C_Matrix)
{
  local_work = work;
  out_of_core_buffer = buffer;
}

CCOperation::~CCOperation()
{
}

void CCOperation::print()
{
  if(reindexing.size())
    fprintf(outfile,"\n\tReindexing = %s",reindexing.c_str());
  fprintf(outfile,"\n\tNumericalFactor = %lf",factor);
  fprintf(outfile,"\tAssigment = %s",assignment.c_str());
  fprintf(outfile,"\tOperation = %s",operation.c_str());
  fprintf(outfile,"\n\tA = %s",A_Matrix->get_label().c_str());
  if(B_Matrix!=NULL)
    fprintf(outfile,"\tB = %s",B_Matrix->get_label().c_str());
  if(C_Matrix!=NULL)
    fprintf(outfile,"\tC = %s",C_Matrix->get_label().c_str());
}

void CCOperation::print_operation()
{
  fprintf(outfile,"%s",A_Matrix->get_label().c_str());
  fprintf(outfile," %s",assignment.c_str());
  if(reindexing.size())
    fprintf(outfile," %s",reindexing.c_str());
  fprintf(outfile," %lf",factor);
  if(B_Matrix!=NULL)
    fprintf(outfile," %s",B_Matrix->get_label().c_str());
  fprintf(outfile," %s",operation.c_str());
  if(C_Matrix!=NULL)
    fprintf(outfile," %s",C_Matrix->get_label().c_str());
}

void CCOperation::print_timing()
{
  DEBUGGING(1,
  fprintf(outfile,"\n-----------------------------------------");
  fprintf(outfile,"\nzero_timing             = %f",zero_timing);
  fprintf(outfile,"\nnumerical_timing        = %f",numerical_timing);
  fprintf(outfile,"\ncontract_timing         = %f",contract_timing);
  fprintf(outfile,"\ntensor_timing           = %f",tensor_timing);
  fprintf(outfile,"\ndot_timing              = %f",dot_timing);
  fprintf(outfile,"\nplus_timing             = %f",plus_timing);
  fprintf(outfile,"\nproduct_timing          = %f",product_timing);
  fprintf(outfile,"\ndivision_timing         = %f",division_timing);
  fprintf(outfile,"\nsort_timing             = %f",sort_timing);
  fprintf(outfile,"\nPartA_timing            = %f",PartA_timing);
  fprintf(outfile,"\nPartB_timing            = %f",PartB_timing);
  fprintf(outfile,"\nPartC_timing            = %f",PartC_timing);
  fprintf(outfile,"\n-----------------------------------------\n");
  );
}

}} /* End Namespaces */
