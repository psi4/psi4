#ifndef _psi_src_lib_libmemtrix_block_matrix_h_
#define _psi_src_lib_libmemtrix_block_matrix_h_

#include "block_vector.h"
#include "matrix_base.h"
#include <string>

typedef std::vector<int> vecint;

extern FILE* outfile;

namespace psi{ namespace mcscf{

class BlockMatrix
{
public:
  BlockMatrix();
  BlockMatrix(std::string label, int nirreps, size_t*& rows_size, size_t*& cols_size);
  BlockMatrix(std::string label, int nirreps, int*& rows_size, int*& cols_size);
  BlockMatrix(std::string label, int nirreps, vecint& rows_size, vecint& cols_size);
  ~BlockMatrix();

  // Inlines
  void        set(int h, int i, int j, double value) {matrix_base_[h]->set(i,j,value);}
  void        add(int h, int i, int j, double value) {matrix_base_[h]->add(i,j,value);}
  double      get(int h, int i, int j)               {return(matrix_base_[h]->get(i,j));}
  size_t      get_rows(int h)                        {return(matrix_base_[h]->get_rows());}
  size_t      get_cols(int h)                        {return(matrix_base_[h]->get_cols());}
  size_t      get_abs_row(int h,int i)               {return(rows_offset_[h] + i);}
  size_t      get_abs_col(int h,int i)               {return(cols_offset_[h] + i);}

  // Overloaded operators
  BlockMatrix& operator=(BlockMatrix& rhs);
  BlockMatrix& operator+=(const BlockMatrix& rhs);
  BlockMatrix& operator-=(const BlockMatrix& rhs);
  friend double dot(BlockMatrix* A,BlockMatrix* B);

  void        print();
  void        zero();
  void        zero_diagonal();
  void        scale(double factor);
  void        transpose();

  void        multiply(bool transpose_A, bool transpose_B, BlockMatrix* A, BlockMatrix* B);
  void        diagonalize(BlockMatrix* eigenvectors, BlockVector* eigenvalues);
  MatrixBase* getMatrixBase(int h) {return(matrix_base_[h]);}


  // Reference counting related
  unsigned int ref ()  const { return ref_;}   // Number of references
  void add_reference      () { ref_++;}
  bool subtract_reference () { if (--ref_ == 0){ delete this; return true;} return false; }
  // Reference count
  unsigned int ref_;

private:
  // Matrix label and pointer
  std::string label_;
  MatrixBase** matrix_base_;

  // Block sizes etc.
  size_t*  rows_size_;
  size_t*  cols_size_;
  size_t*  rows_offset_;
  size_t*  cols_offset_;
  int      nirreps_;

  void startup(std::string label, int nirreps, size_t*& rows_size, size_t*& cols_size );
  void startup(std::string label, int nirreps, int*& rows_size, int*& cols_size);
  void startup(std::string label, int nirreps, vecint& rows_size, vecint& cols_size);
  void cleanup();
};

}}

#endif // _psi_src_lib_libmemtrix_block_matrix_h_
