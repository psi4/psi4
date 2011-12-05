#ifndef _psi_src_lib_libmemtrix_sblock_matrix_h_
#define _psi_src_lib_libmemtrix_sblock_matrix_h_

#include <string>
#include <vector>

#include "block_matrix.h"
#include "sblock_vector.h"

namespace psi{ namespace mcscf{

// Smart version of BlockMatrix
class SBlockMatrix
{
public:
  SBlockMatrix();
  SBlockMatrix(std::string label, int nirreps, size_t*& rows_size, size_t*& cols_size);
  SBlockMatrix(std::string label, int nirreps, int*& rows_size, int*& cols_size);
  SBlockMatrix(std::string label, int nirreps, vecint& rows_size, vecint& cols_size);
  ~SBlockMatrix() { if(block_matrix_)
                        if(block_matrix_->subtract_reference())
                            block_matrix_ = 0;}

  // Manual allocation
  void allocate(std::string label, int nirreps, size_t*& rows_size, size_t*& cols_size);
  void allocate(std::string label, int nirreps, int*& rows_size, int*& cols_size);
  void allocate(std::string label, int nirreps, vecint& rows_size, vecint& cols_size);

  void subtract_reference(){
      if(block_matrix_){
          if(block_matrix_->subtract_reference())
              block_matrix_ = 0;
      }
  }

  // Copy constructor and assignment operator
  SBlockMatrix             (SBlockMatrix& src);
  SBlockMatrix& operator=  (SBlockMatrix& src);
  SBlockMatrix& operator+= (SBlockMatrix& src);
  SBlockMatrix& operator-= (SBlockMatrix& src);

  // Allow access to the implementation object functions
  const BlockMatrix* operator-> () const {return block_matrix_;}
  BlockMatrix*       operator-> ()       {return block_matrix_;}

  // Copy the implementation object
//   void copy(SBlockMatrix& src);

  // Access the implementation object
  BlockMatrix* getBlockMatrix() {return block_matrix_;}

  //Operations
  void multiply(bool transpose_A, bool transpose_B, SBlockMatrix& A, SBlockMatrix& B);
  void diagonalize(SBlockMatrix& eigenmatrix,SBlockVector& eigenvalues);
  void transpose() {block_matrix_->transpose();}
  void scale(double factor) {block_matrix_->scale(factor);}
  friend double dot(SBlockMatrix& A,SBlockMatrix& B);

  // Checking functions
  bool is_allocated() {return (block_matrix_);}
  void check(const char* cstr);

private:
  SBlockMatrix(BlockMatrix* block_matrix);

  BlockMatrix* block_matrix_;
};

}}

#endif // _psi_src_lib_libmemtrix_sblock_matrix_h_

/*
  SBlockMatrix(std::string label_, int nirreps_, int*& block_size_);
  ~SBlockMatrix();
  void print();
  int       get_nirreps()                          {return(nirreps);}
  int       get_nirreps() const                    {return(nirreps);}
  double*** get_matrix()                           {return(matrix);}
  double**  get_block(int h)                       {return(matrix[h]);}
  const double**  get_block(int h) const           {return((const double**)matrix[h]);}
  double    get(int h, int i, int j)               {return(matrix[h][i][j]);}
  void      set(int h, int i, int j, double value) {matrix[h][i][j]  = value;}
  void      add(int h, int i, int j, double value) {matrix[h][i][j] += value;}
  int       get_block_size(int h)                  {return(block_size[h]);}
  int       get_block_size(int h) const            {return(block_size[h]);}
  void      diagonalize(SBlockMatrix* eigenvectors, double* eigenvalues);
  void      DGEMM(bool transpose_A, bool transpose_B, SBlockMatrix* A, SBlockMatrix* B);
  void      minus(SBlockMatrix* B);
  void      zero();
  void      scale(double factor);
  void      transpose();
  friend double  operator^(const SBlockMatrix& rhs,const SBlockMatrix& lhs);
  SBlockMatrix& operator=(const SBlockMatrix& rhs);
  SBlockMatrix& operator+=(const SBlockMatrix& lhs);
private:
  // Matrix label and pointer
  std::string label;
  double***   matrix;

  // Block sizes etc.
  int   nirreps;
  int*  block_size;
  int*  block_offset;

  // Private functions
  void  sort_eigensystem(int n,double*& eigenvalues,double**& eigenvectors);
  void  cleanup();
*/
