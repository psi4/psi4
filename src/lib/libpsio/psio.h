#ifndef PSIO_H
#define PSIO_H

#include <cstdio>
#include <libpsio/config.h>

namespace psi {

  int psio_init(void);
  int psio_ipv1_config(void);
  int psio_state(void);
  int psio_done(void);
  void psio_error(unsigned int unit, unsigned int errval);
  int psio_open(unsigned int unit, int status);
  int psio_close(unsigned int unit, int keep);

  unsigned int psio_get_numvols_default(void);
  int psio_get_volpath_default(unsigned int volume, char **path);
  int psio_get_filename_default(char **name);
  psio_address psio_get_address(psio_address start, ULI shift);
  psio_address psio_get_global_address(psio_address entry_start,
                                       psio_address rel_address);
  int psio_volseek(psio_vol *vol, ULI page, ULI offset, ULI numvols);
  ULI psio_get_length(psio_address sadd, psio_address eadd);
  psio_address psio_get_entry_end(unsigned int unit, const char *key);

  int psio_tocwrite(unsigned int unit);
  int psio_tocread(unsigned int unit);
  void psio_tocprint(unsigned int unit, FILE *output);
  psio_tocentry *psio_tocscan(unsigned int unit, const char *key);
  bool psio_tocentry_exists(unsigned int unit, const char *key);
  psio_tocentry *psio_toclast(unsigned int unit);
  int psio_tocclean(unsigned int unit, const char *key);

  int psio_write(unsigned int unit, const char *key, char *buffer, ULI size,
                 psio_address sadd, psio_address *eadd);
  int psio_read(unsigned int unit, const char *key, char *buffer, ULI size,
                psio_address sadd, psio_address *eadd);
  int psio_write_entry(unsigned int unit, const char *key, char *buffer, ULI size);
  int psio_read_entry(unsigned int unit, const char *key, char *buffer, ULI size);
  int psio_write_block(unsigned int unit, const char *key, char *buffer, ULI blksiz,
                       ULI start_blk, ULI end_blk);
  int psio_read_block(unsigned int unit, const char *key, char *buffer, ULI blksiz,
                      ULI start_blk, ULI end_blk);
  int psio_rw(unsigned int unit, char *buffer, psio_address address, ULI size,
              int wrt);
  int psio_zero_disk(unsigned int unit, const char *key, ULI rows, ULI cols);

  int psio_open_check(unsigned int unit);

  unsigned long int psio_rd_toclen(unsigned int unit);
  void psio_wt_toclen(unsigned int unit, ULI toclen);

  int psio_set_filescfg_kwd(const char* kwdgrp, const char* kwd, int unit,
                            const char* kwdval);
  const char* psio_get_filescfg_kwd(const char* kwdgrp, const char* kwd,
                                    int unit);

}

#endif    /* #ifndef PSIO_H */

