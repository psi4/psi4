#ifndef _psi3_libiwl_config_h_
#define _psi3_libiwl_config_h_

#include <libpsio/config.h>

typedef short int Label;
typedef double Value;

struct iwlbuf {
  int itap;                   /* tape number for input file */
  psio_address bufpos;        /* current page/offset */
  int ints_per_buf;           /* integrals per buffer */
  int bufszc;                 /* buffer size in characters (bytes) */
  double cutoff;              /* cutoff value for writing */
  int lastbuf;                /* is this the last IWL buffer? 1=yes,0=no */
  int inbuf;                  /* how many ints in current buffer? */
  int idx;                    /* index of integral in current buffer */
  Label *labels;              /* pointer to where integral values begin */
  Value *values;              /* integral values */
};

#define IWL_KEY_BUF "IWL Buffers"
#define IWL_KEY_ONEL "IWL One-electron matrix elements"

#define IWL_INTS_PER_BUF 2980

#endif // _psi3_libiwl_config_h_
