#ifndef _psi_src_bin_cints_Tools_iwl_tebuf_h
#define _psi_src_bin_cints_Tools_iwl_tebuf_h

/*! \file iwl_tebuf.h
    \ingroup CINTS
*/namespace psi { namespace CINTS {

void iwl_buf_wrt_struct_nocut(struct iwlbuf *Buf, struct tebuf *Tebuf, int size);
void iwl_buf_wrt_struct(struct iwlbuf *Buf, struct tebuf *Tebuf, int size, double cutoff);

}}
#endif
