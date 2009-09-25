/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

void classify(int p, int q, int r, int s, double value,
	      struct iwlbuf *ABuf, struct iwlbuf *BBuf,
	      struct iwlbuf *CBuf, struct iwlbuf *DBuf,
	      struct iwlbuf *EBuf, struct iwlbuf *F1Buf, 
	      struct iwlbuf *F2Buf);

void distribute_rhf(int filenum, int first_tmp, double tolerance, int keep_input)
{
  struct iwlbuf InBuf;
  struct iwlbuf ABuf, BBuf, CBuf, DBuf, EBuf, F1Buf, F2Buf;
  int lastbuf;
  Value *valptr;
  Label *lblptr;
  int idx, p, q, r, s;
  double value;

  iwl_buf_init(&InBuf, filenum, tolerance, 1, 1);
  iwl_buf_init(&ABuf, first_tmp, tolerance, 0, 0);
  if(params.make_abcd) iwl_buf_init(&BBuf, first_tmp+1, tolerance, 0, 0);
  iwl_buf_init(&CBuf, first_tmp+2, tolerance, 0, 0);
  iwl_buf_init(&DBuf, first_tmp+3, tolerance, 0, 0);
  iwl_buf_init(&EBuf, first_tmp+4, tolerance, 0, 0);
  iwl_buf_init(&F1Buf, first_tmp+5, tolerance, 0, 0);
  if(params.make_aibc) iwl_buf_init(&F2Buf, first_tmp+6, tolerance, 0, 0);

  /* Run through the buffer that's already available */
  lblptr = InBuf.labels;
  valptr = InBuf.values;
  lastbuf = InBuf.lastbuf;

  for (idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
    p = (int) lblptr[idx++];
    q = (int) lblptr[idx++];
    r = (int) lblptr[idx++];
    s = (int) lblptr[idx++];

    value = (double) valptr[InBuf.idx];

    /* Check integral into each class */
    classify(p,q,r,s,value,&ABuf,&BBuf,&CBuf,&DBuf,&EBuf,&F1Buf, &F2Buf);

    /*    fprintf(outfile, "(%d %d|%d %d) = %20.10lf\n", p, q, r, s, value);  */

  } /* end loop through current buffer */

  /* Now run through the rest of the buffers in the file */
  while (!lastbuf) {
    iwl_buf_fetch(&InBuf);
    lastbuf = InBuf.lastbuf;

    for (idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
      p = (int) lblptr[idx++];
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];

      value = (double) valptr[InBuf.idx];

      /* Check integral into each class */
      classify(p,q,r,s,value,&ABuf,&BBuf,&CBuf,&DBuf,&EBuf,&F1Buf, &F2Buf);

      /*      fprintf(outfile, "(%d %d|%d %d) = %20.10lf\n", p, q, r, s, value); */

    } /* end loop through current buffer */
  } /* end loop over reading buffers */


  iwl_buf_close(&InBuf, keep_input);
  iwl_buf_flush(&ABuf, 1);
  if(params.make_abcd) iwl_buf_flush(&BBuf, 1);
  iwl_buf_flush(&CBuf, 1);
  iwl_buf_flush(&DBuf, 1);
  iwl_buf_flush(&EBuf, 1);
  iwl_buf_flush(&F1Buf, 1);
  if(params.make_aibc) iwl_buf_flush(&F2Buf, 1);
  iwl_buf_close(&ABuf, 1);
  if(params.make_abcd) iwl_buf_close(&BBuf, 1);
  iwl_buf_close(&CBuf, 1);
  iwl_buf_close(&DBuf, 1);
  iwl_buf_close(&EBuf, 1);
  iwl_buf_close(&F1Buf, 1);
  if(params.make_aibc) iwl_buf_close(&F2Buf, 1);

  fflush(outfile);
}

void classify_uhf(int p, int q, int r, int s, double value, const char *spin,
		  struct iwlbuf *ABuf1, struct iwlbuf *BBuf1,
		  struct iwlbuf *CBuf1, struct iwlbuf *CBuf2, 
		  struct iwlbuf *DBuf1, struct iwlbuf *EBuf1, 
		  struct iwlbuf *EBuf2, struct iwlbuf *FBuf1, 
		  struct iwlbuf *FBuf2, struct iwlbuf *FBuf3,
		  struct iwlbuf *FBuf4);

void distribute_uhf(const char *spin, int filenum, int first_tmp, double tolerance, int keep_input)
{
  struct iwlbuf InBuf;
  struct iwlbuf ABuf1, BBuf1, CBuf1, CBuf2, DBuf1;
  struct iwlbuf EBuf1, EBuf2, FBuf1, FBuf2, FBuf3, FBuf4;
  int lastbuf;
  Value *valptr;
  Label *lblptr;
  int idx, p, q, r, s;
  double value;

  iwl_buf_init(&InBuf, filenum, tolerance, 1, 1);
  iwl_buf_init(&ABuf1, first_tmp, tolerance, 0, 0);
  iwl_buf_init(&BBuf1, first_tmp+1, tolerance, 0, 0);
  iwl_buf_init(&CBuf1, first_tmp+2, tolerance, 0, 0);
  iwl_buf_init(&CBuf2, first_tmp+3, tolerance, 0, 0);
  iwl_buf_init(&DBuf1, first_tmp+4, tolerance, 0, 0);
  iwl_buf_init(&EBuf1, first_tmp+5, tolerance, 0, 0);
  iwl_buf_init(&EBuf2, first_tmp+6, tolerance, 0, 0);
  iwl_buf_init(&FBuf1, first_tmp+7, tolerance, 0, 0);
  iwl_buf_init(&FBuf2, first_tmp+8, tolerance, 0, 0);
  iwl_buf_init(&FBuf3, first_tmp+9, tolerance, 0, 0);
  iwl_buf_init(&FBuf4, first_tmp+10, tolerance, 0, 0);

  /* Run through the buffer that's already available */
  lblptr = InBuf.labels;
  valptr = InBuf.values;
  lastbuf = InBuf.lastbuf;

  for (idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
    p = (int) lblptr[idx++];
    q = (int) lblptr[idx++];
    r = (int) lblptr[idx++];
    s = (int) lblptr[idx++];

    value = (double) valptr[InBuf.idx];

    /* Check integral into each class */
    classify_uhf(p,q,r,s,value,spin,&ABuf1,&BBuf1,&CBuf1,&CBuf2,
		 &DBuf1,&EBuf1,&EBuf2,&FBuf1,&FBuf2,&FBuf3,&FBuf4);

    /*    fprintf(outfile, "(%d %d|%d %d) = %20.10lf\n", p, q, r, s, value);  */

  } /* end loop through current buffer */

  /* Now run through the rest of the buffers in the file */
  while (!lastbuf) {
    iwl_buf_fetch(&InBuf);
    lastbuf = InBuf.lastbuf;

    for (idx=4*InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
      p = (int) lblptr[idx++];
      q = (int) lblptr[idx++];
      r = (int) lblptr[idx++];
      s = (int) lblptr[idx++];

      value = (double) valptr[InBuf.idx];

      /* Check integral into each class */
      classify_uhf(p,q,r,s,value,spin,&ABuf1,&BBuf1,&CBuf1,&CBuf2,
		   &DBuf1,&EBuf1,&EBuf2,&FBuf1,&FBuf2,&FBuf3,&FBuf4);

      /*      fprintf(outfile, "(%d %d|%d %d) = %20.10lf\n", p, q, r, s, value); */

    } /* end loop through current buffer */
  } /* end loop over reading buffers */


  iwl_buf_close(&InBuf, keep_input);
  iwl_buf_flush(&ABuf1, 1);
  iwl_buf_flush(&BBuf1, 1);
  iwl_buf_flush(&CBuf1, 1);
  iwl_buf_flush(&CBuf2, 1);
  iwl_buf_flush(&DBuf1, 1);
  iwl_buf_flush(&EBuf1, 1);
  iwl_buf_flush(&EBuf2, 1);
  iwl_buf_flush(&FBuf1, 1);
  iwl_buf_flush(&FBuf2, 1);
  iwl_buf_flush(&FBuf3, 1);
  iwl_buf_flush(&FBuf4, 1);
  iwl_buf_close(&ABuf1, 1);
  iwl_buf_close(&BBuf1, 1);
  iwl_buf_close(&CBuf1, 1);
  iwl_buf_close(&CBuf2, 1);
  iwl_buf_close(&DBuf1, 1);
  iwl_buf_close(&EBuf1, 1);
  iwl_buf_close(&EBuf2, 1);
  iwl_buf_close(&FBuf1, 1);
  iwl_buf_close(&FBuf2, 1);
  iwl_buf_close(&FBuf3, 1);
  iwl_buf_close(&FBuf4, 1);

  fflush(outfile);
}

}} // namespace psi::ccsort
