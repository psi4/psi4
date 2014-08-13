/* begin_generated_IBM_copyright_prolog                             */
/*                                                                  */
/* --------------------------------------------------------------   */
/* (C)Copyright IBM Corp. 2007, 2008                                */
/* IBM BSD License.                                                 */
/* --------------------------------------------------------------   */
/* end_generated_IBM_copyright_prolog                               */
/********************************************************************/

/* $Id: bgmldefs.h 6904 2010-09-16 18:55:57Z manoj $ */

#ifndef _bgmldefs_h
#define _bgmldefs_h

void wait_callback(void *clientdata);
void BGML_Wait(unsigned *clientdata);
void BGML_WaitProc(int proc);
void BGML_WaitAll();

void bgml_init_locks (void * local_memlock_table);
void bgml_lockmem (void *start, void *end, int proc);
void bgml_unlockmem (int proc);

typedef void (*BGML_Barrier)  (unsigned);
extern BGML_Barrier bgml_barrier;

#endif
