#include <libmoinfo/libmoinfo.h>
#include <libutil/libutil.h>
#include <cstdio>

#include "blas.h"
#include "debugging.h"
#include "matrix.h"


namespace psi{
    extern FILE *outfile;
    namespace psimrcc{
    extern MOInfo *moinfo;

extern MemoryManager* memory_manager;

using namespace std;

void CCBLAS::zero(const char* cstr)
{
  string str(cstr);
  // To zero diagonals of things like "Fae[v][v]{u}"
  vector<string> names = moinfo->get_matrix_names(str);
  for(size_t n = 0; n < names.size(); ++n){
    CCMatrix* Matrix = get_Matrix(names[n]);
    Matrix->zero_matrix();
    DEBUGGING(5,
      fprintf(outfile,"\n...setting %s to zero",names[n].c_str());
    );
  }
}

void CCBLAS::zero_right_four_diagonal(const char* cstr)
{
  string str(cstr);
  // To zero diagonals of things like "Fae[v][v]{u}"
  vector<string> names = moinfo->get_matrix_names(str);
  for(size_t n = 0; n < names.size(); ++n){
    CCMatrix* Matrix = get_Matrix(names[n]);
    Matrix->zero_right_four_diagonal();
    DEBUGGING(5,
      fprintf(outfile,"\n...setting the right diagonal terms of %s to zero",names[n].c_str());
    );
  }
}

void CCBLAS::zero_non_doubly_occupied(const char* cstr)
{
  string str(cstr);
  // To zero non-doubly occupied MOs of things like "Fae[v][v]{u}"
  vector<string> names = moinfo->get_matrix_names(str);
  for(size_t n = 0; n < names.size(); ++n){
    CCMatrix* Matrix = get_Matrix(names[n]);
    Matrix->zero_non_doubly_occupied();
    DEBUGGING(5,
      fprintf(outfile,"\n...setting the right diagonal terms of %s to zero",names[n].c_str());
    );
  }
}

void CCBLAS::zero_non_external(const char* cstr)
{
  string str(cstr);
  // To zero non-external MOs of things like "Fae[v][v]{u}"
  vector<string> names = moinfo->get_matrix_names(str);
  for(size_t n = 0; n < names.size(); ++n){
    CCMatrix* Matrix = get_Matrix(names[n]);
    Matrix->zero_non_external();
    DEBUGGING(5,
      fprintf(outfile,"\n...setting the right diagonal terms of %s to zero",names[n].c_str());
    );
  }
}

void CCBLAS::scale(const char* cstr,int reference,double value)
{
  string str(cstr);
  scale(str,reference,value);
}

void CCBLAS::scale(string& str,int reference,double value)
{
  string matrix_str = add_reference(str,reference);
  // Make sure that the element that we are retrieving is present
  MatrixMap::iterator iter = matrices.find(matrix_str);
  if(iter!=matrices.end()){
    load(iter->second);
    iter->second->scale(value);
    return;
  }
  throw PSIEXCEPTION("\nCCBLAS::scale() couldn't find matrix " + matrix_str);
}

void CCBLAS::reduce_spaces(const char* out,const char* in)
{
  string  in_str(in);
  string out_str(out);
  // To zero diagonals of things like "Fae[v][v]{u}"
  vector<string>  in_names = moinfo->get_matrix_names(in_str);
  vector<string> out_names = moinfo->get_matrix_names(out_str);
  if(in_names.size()!=out_names.size())
    throw PSIEXCEPTION("CCBLAS::map_spaces, number of references mismatch");
  for(size_t n = 0; n < in_names.size(); ++n){
    CCMatrix*  in_Matrix = get_Matrix(in_names[n]);
    CCMatrix* out_Matrix = get_Matrix(out_names[n]);
    process_reduce_spaces(out_Matrix,in_Matrix);
  }
}

void CCBLAS::process_reduce_spaces(CCMatrix* out_Matrix,CCMatrix* in_Matrix)
{
  double*** out_matrix = out_Matrix->get_matrix();
  const intvec&  act_to_occ = moinfo->get_actv_to_occ();
  const intvec&  act_to_vir = moinfo->get_actv_to_vir();

  string& out_index_label = out_Matrix->get_index_label();
  string&  in_index_label =  in_Matrix->get_index_label();

  int index_label_size = out_index_label.size();

  int** map;
  allocate2(int,map,index_label_size,moinfo->get_nmo());

  for(int k=0;k<index_label_size;k++){
    if(out_index_label[k]=='a' && in_index_label[k]=='o'){
      for(int l=0;l<moinfo->get_nactv();l++){
        map[k][l] = act_to_occ[l];
      }
    }else if(out_index_label[k]=='a' && in_index_label[k]=='v'){
      for(int l=0;l<moinfo->get_nactv();l++){
        map[k][l] = act_to_vir[l];
      }
    } else {
      for(int l=0;l<moinfo->get_nmo();l++){
        map[k][l] = l;
      }
    }
  }

  if(index_label_size==2){
    short* pq = new short[2];
    for(int h=0;h<moinfo->get_nirreps();h++){
      for(size_t i = 0;i < out_Matrix->get_left_pairpi(h); ++i){
        for(size_t j = 0; j < out_Matrix->get_right_pairpi(h); ++j){
          out_Matrix->get_two_indices(pq,h,i,j);
          out_matrix[h][i][j] = in_Matrix->get_two_address_element(map[0][pq[0]],map[1][pq[1]]);
        }
      }
    }
    delete[] pq;
  }else if(index_label_size==4){
    short* pqrs = new short[4];
    for(int h=0;h<moinfo->get_nirreps();h++){
      for(size_t i = 0; i < out_Matrix->get_left_pairpi(h); ++i){
        for(size_t j = 0; j < out_Matrix->get_right_pairpi(h); ++j){
          out_Matrix->get_four_indices(pqrs,h,i,j);
          out_matrix[h][i][j] = in_Matrix->get_four_address_element(map[0][pqrs[0]],map[1][pqrs[1]],map[2][pqrs[2]],map[3][pqrs[3]]);
        }
      }
    }
    delete[] pqrs;
  }
  release2(map);
}

void CCBLAS::expand_spaces(const char* out,const char* in)
{
  string  in_str(in);
  string out_str(out);

  vector<string>  in_names = moinfo->get_matrix_names(in_str);
  vector<string> out_names = moinfo->get_matrix_names(out_str);
  if(in_names.size()!=out_names.size())
    throw PSIEXCEPTION("CCBLAS::map_spaces, number of references mismatch");
  for(size_t n = 0; n < in_names.size(); ++n){
    CCMatrix*  in_Matrix = get_Matrix(in_names[n]);
    CCMatrix* out_Matrix = get_Matrix(out_names[n]);
    process_expand_spaces(out_Matrix,in_Matrix);
  }
}

void CCBLAS::process_expand_spaces(CCMatrix* out_Matrix,CCMatrix* in_Matrix)
{
  double*** out_matrix = out_Matrix->get_matrix();
  const intvec&    act_to_occ = moinfo->get_actv_to_occ();
  const intvec&    act_to_vir = moinfo->get_actv_to_vir();

  string& out_index_label = out_Matrix->get_index_label();
  string&  in_index_label =  in_Matrix->get_index_label();

  int index_label_size = out_index_label.size();

  int** map;
  allocate2(int,map,index_label_size,moinfo->get_nmo());

  for(int k=0;k<index_label_size;k++){
    if(out_index_label[k]=='a' && in_index_label[k]=='o'){
      for(int l=0;l<moinfo->get_nactv();l++){
        map[k][l] = act_to_occ[l];
      }
    }else if(out_index_label[k]=='a' && in_index_label[k]=='v'){
      for(int l=0;l<moinfo->get_nactv();l++){
        map[k][l] = act_to_vir[l];
      }
    } else {
      for(int l=0;l<moinfo->get_nmo();l++){
        map[k][l] = l;
      }
    }
  }

  if(index_label_size==2){
    short* pq = new short[2];

    for(int h=0;h<moinfo->get_nirreps();h++){
      for(size_t i = 0; i < out_Matrix->get_left_pairpi(h); ++i){
        for(size_t j = 0; j < out_Matrix->get_right_pairpi(h); ++j){
          out_Matrix->get_two_indices(pq,h,i,j);
          in_Matrix->set_two_address_element(map[0][pq[0]],
                                             map[1][pq[1]],
                                             out_matrix[h][i][j]);
        }
      }
    }

    delete[] pq;
  }else if(index_label_size==4){
    short* pqrs = new short[4];
    for(int h=0;h<moinfo->get_nirreps();h++){
      for(size_t i = 0; i < out_Matrix->get_left_pairpi(h); ++i){
        for(size_t j = 0; j < out_Matrix->get_right_pairpi(h); ++j){
          out_Matrix->get_four_indices(pqrs,h,i,j);
          in_Matrix->set_four_address_element(map[0][pqrs[0]],
                                              map[1][pqrs[1]],
                                              map[2][pqrs[2]],
                                              map[3][pqrs[3]],
                                              out_matrix[h][i][j]);
        }
      }
    }
    delete[] pqrs;
  }
  release2(map);
}

}} /* End Namespaces */
