#ifndef _psi_src_bin_cints_DFT_init_close_shell_info_h
#define _psi_src_bin_cints_DFT_init_close_shell_info_h

/*! \file init_close_shell_info.h
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include"data_structs.h"
namespace psi { namespace cints {
  struct close_shell_info_s init_close_shell_info(void);
  void free_close_shell_info(struct close_shell_info_s close_shell_info);
  void print_close_shell_info(struct close_shell_info_s close);
}}
#endif
