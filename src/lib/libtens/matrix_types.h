#ifndef _psi_src_bin_psimrcc_matrix_types_h_
#define _psi_src_bin_psimrcc_matrix_types_h_

typedef std::map<psi::psimrcc::CCMatrix*,int>                           MatCnt;
typedef std::map<std::string,psi::psimrcc::CCMatrix*>                   MatrixMap;
typedef std::map<std::string,psi::psimrcc::CCMatrix*>::iterator         MatMapIt;
typedef std::vector<std::pair<psi::psimrcc::CCMatrix*,int> >            MatrixBlks;
typedef std::vector<std::pair<psi::psimrcc::CCMatrix*,int> >::iterator  MatBlksIt;

#endif // _psi_src_bin_psimrcc_matrix_types_h_
