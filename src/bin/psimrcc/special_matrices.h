#ifndef _psi_src_bin_psimrcc_special_matrices_h_
#define _psi_src_bin_psimrcc_special_matrices_h_

#include "manybody.h"
#include "index.h"
#include "index_iterator.h"

namespace psi{ namespace psimrcc{

class MatrixBase
{
public:
  // Constructor and destructor
  MatrixBase(size_t nrows_,size_t ncols);
  ~MatrixBase();
  size_t   get_nrows() {return nrows;}
  size_t   get_ncols() {return ncols;}
  void     add(size_t row,size_t col,double value) {matrix[row][col] += value;}
  void     set(size_t row,size_t col,double value) {matrix[row][col] = value;}
  double   get(size_t row,size_t col) {return matrix[row][col];}
  double** get_matrix() {return matrix;}
  void     zero();
  void     print();
  void     add(MatrixBase* A, double alpha, double beta);
  void     multiply(MatrixBase* A, MatrixBase* B, double alpha, double beta);
  void     contract(MatrixBase* A, MatrixBase* B, double const alpha, double const beta);
  double   norm();
private:
  size_t   nrows;
  size_t   ncols;
  double** matrix;
};


class BlockMatrix
{
public:
  // Constructor and destructor
  BlockMatrix(int nirreps, std::vector<size_t>& rows_size_, std::vector<size_t>& cols_size_,int sym);
  ~BlockMatrix();

  void print();

  void     add(BlockMatrix* A, double alpha, double beta);
  void     add_acb(double z,BlockMatrix* A, CCIndex* pqr_index,CCIndex* p_index,CCIndex* qr_index,double a);
  void     add_cab(double z,BlockMatrix* A, CCIndex* pqr_index,CCIndex* p_index,CCIndex* qr_index,double a);

  void     add(int h,size_t row,size_t col,double value) {blocks[h]->add(row,col,value);}
  void     set(int h,size_t row,size_t col,double value) {blocks[h]->set(row,col,value);}
  double   get(int h,size_t row,size_t col) {return blocks[h]->get(row,col);}
  double** get_matrix(int h) {return blocks[h]->get_matrix();}
  MatrixBase* get_matrixbase(int h) {return blocks[h];}
  void     multiply(BlockMatrix* A, BlockMatrix* B,double alpha, double beta);
  void     contract(BlockMatrix* A, BlockMatrix* B, double alpha, double beta);
  void     cyclical_permutation_1_2(BlockMatrix* A, CCIndex* pqr_index,CCIndex* p_index,CCIndex* qr_index);
  void     a_b_permutation_1_2(BlockMatrix* A, CCIndex* pqr_index,CCIndex* p_index,CCIndex* qr_index);
  void     add_c_ab_permutation_1_2(BlockMatrix* A, CCIndex* pqr_index,CCIndex* p_index,CCIndex* qr_index);
  void     add_permutation_1_2(double z,BlockMatrix* A, CCIndex* pqr_index,CCIndex* p_index,CCIndex* qr_index,
      double a,double b,double c,double d,double e,double f);
  void     a_b_permutation(CCIndex* pqr_index,CCIndex* p_index,CCIndex* qr_index);
//  void     add_a_b_permutation_1_2(BlockMatrix* A, CCIndex* pqr_index,CCIndex* p_index,CCIndex* qr_index);
  void     zero();
  double   norm();
private:
  MatrixBase** blocks;
  // Block sizes etc.
  std::vector<size_t> rows_size;
  std::vector<size_t> cols_size;
  std::vector<size_t> rows_offset;
  std::vector<size_t> cols_offset;
  int      nirreps;
  int      sym;
};


class IndexMatrix
{
  typedef std::pair<size_t,int>          IMIndex;
  typedef std::map<IMIndex,BlockMatrix*> BMMap;
public:
  // Constructor and destructor
  IndexMatrix();
  ~IndexMatrix();

  void add_block_matrix(size_t index,int ref,BlockMatrix* block_matrix);
  BlockMatrix* get_block_matrix(size_t index,int ref = 0);
  void print();
private:
  BMMap  matrices;
};

void multiply(BlockMatrix* A, BlockMatrix* B);

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_special_matrices_h_

