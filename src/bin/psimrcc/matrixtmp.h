#ifndef _psi_src_bin_psimrcc_ccmatrixtmp_h
#define _psi_src_bin_psimrcc_ccmatrixtmp_h
/***************************************************************************
 *   Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *   frank@ccc.uga.edu
 *   SR/MRCC Code
 ***************************************************************************/


namespace psi{ namespace psimrcc{

class CCMatrix;

enum DiskOpt {none,dump,release};

/**
	@author Francesco A. Evangelista and Andrew C. Simmonett <frank@ccc.uga.edu>
*/
class CCMatTmp{
public:
  CCMatTmp(CCMatrix* Matrix, DiskOpt disk_option = none);
  CCMatrix* operator->()   const {return(Matrix_);}
  CCMatrix* get_CCMatrix() const {return(Matrix_);}
  ~CCMatTmp();
private:
  DiskOpt     disk_option_;
  CCMatrix*   Matrix_;
};

/**
	@author Francesco A. Evangelista and Andrew C. Simmonett <frank@ccc.uga.edu>
*/
class CCMatIrTmp{
public:
  CCMatIrTmp(CCMatrix* Matrix, int irrep, DiskOpt disk_option = none);
  CCMatrix* operator->()   const {return(Matrix_);}
  CCMatrix* get_CCMatrix() const {return(Matrix_);}
  ~CCMatIrTmp();
private:
  DiskOpt     disk_option_;
  int         irrep_;
  CCMatrix*   Matrix_;
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_ccmatrixtmp_h
