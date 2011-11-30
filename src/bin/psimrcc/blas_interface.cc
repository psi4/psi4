#include <cstdlib>

#include <libmoinfo/libmoinfo.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libutil/libutil.h>

#include "blas.h"
#include "index.h"
#include "debugging.h"
#include "matrix.h"

namespace psi{
    extern FILE *outfile;
    namespace psimrcc{
    extern MOInfo *moinfo;
    extern MemoryManager *memory_manager;

using namespace std;

void CCBLAS::add_index(const char* cstr)
{
  // Make sure that the element that we are adding is not present
  string str(cstr);
  to_lower(str);
  if(indices.find(str)==indices.end()){
    indices.insert(make_pair(str,new CCIndex(str)));
  }
}

void CCBLAS::add_Matrix(const char* cstr)
{
  string str(cstr);
  vector<string> names = moinfo->get_matrix_names(str);
  for(size_t n = 0; n < names.size(); ++n)
    add_Matrix_ref(names[n]);
}

void CCBLAS::add_Matrix(string str)
{
  vector<string> names = moinfo->get_matrix_names(str);
  for(size_t n = 0; n < names.size(); ++n)
    add_Matrix_ref(names[n]);
}

void CCBLAS::add_Matrix_ref(std::string& str)
{
  // Make sure that the element that we are adding is not present
  MatrixMap::iterator iter = matrices.find(str);
  if(iter==matrices.end()){
    CCIndex* index_pointer[2];
    // Default: assume the [] indexing
    index_pointer[0]=get_index("[]");
    index_pointer[1]=get_index("[]");
    vector<string> index_string_vec = split_indices(str);
    for(size_t i = 0; i < index_string_vec.size(); ++i)
      index_pointer[i]=get_index(index_string_vec[i]);
    matrices.insert(make_pair(str,new CCMatrix(str,index_pointer[0],index_pointer[1])));
  }
}

CCIndex* CCBLAS::get_index(const char* cstr)
{
  string str(cstr);
  to_lower(str);
  // Make sure that the element that we are retrieving is present
  IndexMap::iterator iter = indices.find(str);
  if(iter!=indices.end()){
    return(indices[str]);
  }
  throw PSIEXCEPTION("\nCCBLAS::get_index() couldn't find index " + str);
  return(NULL);
}

CCIndex* CCBLAS::get_index(string& str)
{
  to_lower(str);
  // Make sure that the element that we are retrieving is present
  IndexMap::iterator iter = indices.find(str);
  if(iter!=indices.end()){
    return(indices[str]);
  }
  throw PSIEXCEPTION("\nCCBLAS::get_index() couldn't find index " + str);
  return(NULL);
}

CCMatTmp CCBLAS::get_MatTmp(std::string str, int reference, DiskOpt disk_option)
{
  append_reference(str,reference);
  load(get_Matrix(str));
  return(CCMatTmp(get_Matrix(str),disk_option));
}

CCMatTmp CCBLAS::get_MatTmp(std::string str, DiskOpt disk_option)
{
  load(get_Matrix(str));
  return(CCMatTmp(get_Matrix(str),disk_option));
}

CCMatTmp CCBLAS::get_MatTmp(CCMatrix* Matrix, DiskOpt disk_option)
{
  load(Matrix);
  return(CCMatTmp(Matrix,disk_option));
}

CCMatIrTmp CCBLAS::get_MatIrTmp(std::string str, int reference, int irrep,  DiskOpt disk_option)
{
  append_reference(str,reference);
  load_irrep(get_Matrix(str),irrep);
  return(CCMatIrTmp(get_Matrix(str),irrep,disk_option));
}

CCMatIrTmp CCBLAS::get_MatIrTmp(std::string str, int irrep, DiskOpt disk_option)
{
  load_irrep(get_Matrix(str),irrep);
  return(CCMatIrTmp(get_Matrix(str),irrep,disk_option));
}

CCMatIrTmp CCBLAS::get_MatIrTmp(CCMatrix* Matrix, int irrep, DiskOpt disk_option)
{
  load_irrep(Matrix,irrep);
  return(CCMatIrTmp(Matrix,irrep,disk_option));
}

CCMatrix* CCBLAS::get_Matrix(const char* cstr, int reference)
{
  string str(cstr);
  append_reference(str,reference);
  return(get_Matrix(str));
}

CCMatrix* CCBLAS::get_Matrix(const char* cstr)
{
  string str(cstr);
  return(get_Matrix(str));
}

CCMatrix* CCBLAS::get_Matrix(string& str)
{
  // Make sure that the element that we are retrieving is present
  MatrixMap::iterator iter = matrices.find(str);
  if(iter!=matrices.end())
    return(matrices[str]);
  throw PSIEXCEPTION("\nCCBLAS::get_matrix() couldn't find matrix " + str);
  return(NULL);
}

CCMatrix* CCBLAS::get_Matrix(string& str, string& expression)
{
  // Make sure that the element that we are retrieving is present
  MatrixMap::iterator iter = matrices.find(str);
  if(iter!=matrices.end()){
    return(matrices[str]);
  }
  throw PSIEXCEPTION("\n\nCCBLAS::parse() couldn't find the matrix " + str + " in the CCMatrix list\n\nwhile parsing the string:\n\t " + expression + "\n\n");
  return NULL;
}

void CCBLAS::set_scalar(const char* cstr,int reference,double value)
{
  string str(cstr);
  set_scalar(str,reference,value);
}

void CCBLAS::set_scalar(string& str,int reference,double value)
{
  string matrix_str = add_reference(str,reference);
  // Make sure that the element that we are retrieving is present
  MatrixMap::iterator iter = matrices.find(matrix_str);
  if(iter!=matrices.end()){
    load(iter->second);
    iter->second->set_scalar(value);
    return;
  }
  throw PSIEXCEPTION("\nCCBLAS::set_scalar() couldn't find matrix " + matrix_str);
}

double CCBLAS::get_scalar(const char* cstr,int reference)
{
  string str(cstr);
  return(get_scalar(str,reference));
}

double CCBLAS::get_scalar(string& str,int reference)
{
  string matrix_str(str);
  append_reference(matrix_str,reference);
  // Make sure that the element that we are retrieving is present
  MatrixMap::iterator iter = matrices.find(matrix_str);
  if(iter!=matrices.end()){
    load(iter->second);
    return(iter->second->get_scalar());
  }
  throw PSIEXCEPTION("\nCCBLAS::get_scalar() couldn't find matrix " + matrix_str);
  return (0.0);
}

double CCBLAS::get_scalar(string str)
{
  // Make sure that the element that we are retrieving is present
  MatrixMap::iterator iter = matrices.find(str);
  if(iter!=matrices.end()){
    load(iter->second);
    return(iter->second->get_scalar());
  }
  throw PSIEXCEPTION("\nCCBLAS::get_scalar() couldn't find matrix " + str);
  return (0.0);
}

void CCBLAS::load(CCMatrix* Matrix)
{
  if(Matrix->is_allocated()){
    DEBUGGING(2,
      fprintf(outfile,"\nCCBLAS::load(%s): matrix is in core.",Matrix->get_label().c_str());
    );
  }else{
    DEBUGGING(2,
      fprintf(outfile,"\nCCBLAS::load(%s): matrix is not in core. Loading it :[",Matrix->get_label().c_str());
    );
    // Do we have enough memory to fit the entire matrix in core?
    size_t memory_required = Matrix->get_memory2();
    make_space(memory_required);
    Matrix->load();
    DEBUGGING(2,
      fprintf(outfile,"\n] <- done.");
    );
  }
}

void CCBLAS::load_irrep(CCMatrix* Matrix,int h)
{
  if(Matrix->is_block_allocated(h)){
    DEBUGGING(2,
      fprintf(outfile,"\nCCBLAS::load_irrep(%s,%d): matrix block is in core.",Matrix->get_label().c_str(),h);
    )
  }else{
    DEBUGGING(2,
      fprintf(outfile,"\nCCBLAS::load_irrep(%s,%d): matrix block is not in core. Loading it : [",Matrix->get_label().c_str(),h);
    )
    // Do we have enough memory to fit the entire matrix in core?
    size_t memory_required = Matrix->get_memorypi2(h);
    make_space(memory_required);
    Matrix->load_irrep(h);
    DEBUGGING(2,
      fprintf(outfile,"\n] <- done.");
    )
  }
}

void CCBLAS::make_space(size_t memory_required)
{
  if(memory_required < memory_manager->get_FreeMemory())
    return;
  else{
    fprintf(outfile,"\nCCBLAS::make_space() not implemented yet!!!");
    // Attempt #1
  }
}

}} /* End Namespaces */
