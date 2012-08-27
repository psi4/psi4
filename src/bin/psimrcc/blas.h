#ifndef _psi_src_bin_psimrcc_ccblas_h
#define _psi_src_bin_psimrcc_ccblas_h

/*! \file
    \ingroup (PSIMRCC)
    \brief   A class to perform contractions
*/

#include <deque>
#include <string>
#include <vector>
#include <map>
#include <utility>

#include "liboptions/liboptions.h"
#include "matrixtmp.h"
#include "operation.h"

#include "index_types.h"
#include "matrix_types.h"
#include "types.h"

namespace psi{ namespace psimrcc{

class CCIndex;
class CCMatrix;

enum DiisType {DiisEachCycle,DiisCC};

/**
	@author Francesco A. Evangelista and Andrew C. Simmonett <frank@ccc.uga.edu>
*/
class CCBLAS{
public:
  typedef std::vector<std::string>            strvec;
  typedef std::vector<int>                    intvec;
  typedef std::vector<std::pair<int,int> >    intpairvec;
  typedef std::deque<CCOperation>             OpDeque;


  CCBLAS(Options &options);
  ~CCBLAS();

  Options &options_;
  // Add routines
  void       add_Matrix(const char* cstr);
  void       add_Matrix(std::string str);
  void       add_index(const char* cstr);
  // Solve and sort
  void       solve(const char* cstr);
  void       solve(std::string str);
  void       solve_zero_two_diagonal(const char* cstr);
  void       zero_right_four_diagonal(const char* cstr);
  void       zero_left_four_diagonal(const char* cstr);
  void       zero_non_doubly_occupied(const char* cstr);
  void       zero_non_external(const char* cstr);
  void       zero(const char* cstr);
  void       scale(const char* cstr,int reference,double value);
  void       scale(std::string& str,int reference,double value);
  void       reduce_spaces(const char* out,const char* in);
  void       expand_spaces(const char* out,const char* in);
  void       append(const char* cstr);
  void       append(std::string str);
  void       append_zero_two_diagonal(const char* cstr);
  void       compute();
  int        compute_storage_strategy();
  void       show_storage();
  // DIIS
  void       diis_add(std::string amps, std::string delta_amps);
  void       diis_save_t_amps(int cycle);
  void       diis(int cycle, double delta, DiisType diis_type);
  // Printing
  void       print(const char* cstr);
  void       print_ref(std::string& str);
  void       print_memory();
  // Safe get and set
  CCIndex*   get_index(const char* cstr);
  CCIndex*   get_index(std::string& str);
  CCMatTmp   get_MatTmp(std::string str, int reference, DiskOpt disk_option);
  CCMatTmp   get_MatTmp(std::string str, DiskOpt disk_option);
  CCMatTmp   get_MatTmp(CCMatrix* Matrix, DiskOpt disk_option);
  CCMatIrTmp get_MatIrTmp(std::string str, int reference, int irrep,  DiskOpt disk_option);
  CCMatIrTmp get_MatIrTmp(std::string str, int irrep, DiskOpt disk_option);
  CCMatIrTmp get_MatIrTmp(CCMatrix* Matrix, int irrep, DiskOpt disk_option);

  double     get_scalar(std::string str);
  double     get_scalar(const char* cstr,int reference);
  double     get_scalar(std::string& str,int reference);
  void       set_scalar(const char* cstr,int reference,double value);
  void       set_scalar(std::string& str,int reference,double value);

  //DIIS
  std::vector<std::pair<std::string,std::string> > diis_matrices;

  // These have to be improved
  MatrixMap& get_MatrixMap() {return(matrices);}
private:
  bool       full_in_core;
  size_t     work_size;
  size_t     buffer_size;
  MatrixMap  matrices;
  IndexMap   indices;
  OpDeque    operations;
  ArrayVec   work;
  ArrayVec   buffer;
  MatCnt     matrices_in_deque;
  MatCnt     matrices_in_deque_target;
  MatCnt     matrices_in_deque_source;
  SortMap    sortmap;
private:

  IndexMap&  get_IndexMap()  {return(indices);}
  CCMatrix*  get_Matrix(std::string& str);
  CCMatrix*  get_Matrix(const char* cstr);
  CCMatrix*  get_Matrix(const char* cstr, int reference);
  CCMatrix*  get_Matrix(std::string& str,std::string& expression); // Prints a clear error message
  double*    get_work(int n)   {return(work[n]);}
//   double***  get_sortmap(CCIndex* T_left,CCIndex* T_right,int thread);

  void       allocate_matrices_in_core();
  void       load(CCMatrix* Matrix);
  void       load_irrep(CCMatrix* Matrix,int h);
  ///////////////////////////////////////////////////////////////////////////////
  // Class private functions
  ///////////////////////////////////////////////////////////////////////////////
  void       add_Matrix_ref(std::string& str);
  void       add_indices();
  void       add_matrix_ref(std::string& str);
  void       solve_ref(std::string& str);
  int        parse(std::string& str);
  void       process_operations();
  void       process_reduce_spaces(CCMatrix* out_Matrix,CCMatrix* in_Matrix);
  void       process_expand_spaces(CCMatrix* out_Matrix,CCMatrix* in_Matrix);
  bool       get_factor(const std::string& str,double& factor);
  // General routines

  void       make_space(size_t memory_required);
  // Low level memory routines
  void       init();
  void       cleanup();
  void       allocate_work();
  void       allocate_buffer();
  void       free_sortmap();
  void       free_work();
  void       free_indices();
  void       free_matrices();
  void       free_buffer();
};

extern CCBLAS *blas;

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_ccblas_h
