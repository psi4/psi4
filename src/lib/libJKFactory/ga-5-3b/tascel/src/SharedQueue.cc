#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <cstdio>
#include <cstring>

#include <armci.h>

#include "Comm.h"
#include "massert.h"
#include "SharedQueue.h"

using namespace std;
using namespace tascel;
using namespace tascel::comm;


SharedQueue::SharedQueue(int _tsk_size, int _max_ntsks)
  : max_ntsks(_max_ntsks)
  , tsk_size(_tsk_size)
  , q(NULL)
  , sq_state(NULL)
  , head(NULL)
  , tail(NULL)
  , size(NULL)
  , dirty(NULL)
  , td() {
  q = new char *[nproc()];
  massert(q);
  ARMCI_Malloc((void **)q, max_ntsks * tsk_size);
  massertl(q[me()], error1);
  ARMCI_Create_mutexes(1);
  sq_state = new sq_state_t*[nproc()];
  massertl(sq_state, error1);
  ARMCI_Malloc((void **)sq_state, nproc()*sizeof(sq_state_t));
  massertl(sq_state[me()], error2);
  head = &sq_state[me()]->head;
  tail = &sq_state[me()]->tail;
  size = &sq_state[me()]->size;
  dirty = &sq_state[me()]->dirty;
  *head = *tail = *size = *dirty = 0;
  //printf("(cons) %d: *head=%p *tail=%p *size=%p\n", me(), head, tail, size);
  return;

error2:
  printf("Error2 in SharedQueue constructor\n");
  delete [] sq_state;
error1:
  printf("Error1 in SharedQueue constructor\n");
  delete [] q;
error:
  printf("Error in SharedQueue constructor\n");
  throw TSL_ERR;
}

SharedQueue::~SharedQueue() {
  barrier();
  ARMCI_Free(sq_state[me()]);
  ARMCI_Destroy_mutexes();
  ARMCI_Free(q[me()]);
  delete [] sq_state;
  delete [] q;
}

bool
SharedQueue::empty() const {
  /*SK: speculative implementation. It could become empty while this
    operation is ongoing. However, the converse is not true, since
    only the host process call empty() or can add tasks (thus
    incrementing size) */
  massert(*size >= 0);
  return (*size == 0);

error:
  throw TSL_ERR;
}

bool
SharedQueue::getTask(void *dscr, int dlen) {
  massert(dscr);
  massert(dlen == tsk_size);
  ARMCI_Lock(0, me());
  if (*size == 0) {
    ARMCI_Unlock(0, me());
    return false;
  }
  *head = (*head - 1 + max_ntsks) % max_ntsks;
  memcpy(dscr, &q[me()][*head*tsk_size], tsk_size);
  *size -= 1;
  ARMCI_Unlock(0, me());
  return true;
  //  return false;

error:
  throw TSL_ERR;
}

void
SharedQueue::addTask(void *dscr, int dlen) {
  massert(dscr);
  massert(dlen == tsk_size);
  ARMCI_Lock(0, me());
  //   printf("addTask. %d: *head=%p *tail=%p *size=%p\n", me(), head, tail, size);
  massertl(*size != max_ntsks, error2);
  memcpy(&q[me()][*head*tsk_size], dscr, tsk_size);
  *head  = (*head + 1) % max_ntsks;
  *size += 1;
  ARMCI_Unlock(0, me());
  return;

error2:
  ARMCI_Unlock(0, me());
error:
  throw TSL_ERR;
}

bool
SharedQueue::steal(int proc) {
  char *buf;
  massert(*size == 0);
  buf = new char[tsk_size];
  massert(buf);

  //FIXME: might need a lock+unlock to ensure both size and dirty are
  //updated atomically by any thief (say through RDMA)
  td.progress(*dirty == 0);
  *dirty = 0;
  if (td.hasTerminated()) {
    return false;
  }

  ARMCI_Lock(0, proc);
  sq_state_t sq;
  ARMCI_Get(sq_state[proc], &sq, sizeof(sq), proc);
  //   printf("%d: Locked %d for steal. size=%d\n", me(), proc, sq.size);
  if (sq.size == 0) {
    //nothing to steal
    ARMCI_Unlock(0, proc);
    return false;
  }
  sq.size = sq.size - 1;
  sq.dirty = 1; //mark dirty for termination detection
  sq.head = (sq.head - 1 + max_ntsks) % max_ntsks;
  ARMCI_Get(&q[proc][sq.head*tsk_size], buf, tsk_size, proc);
  ARMCI_Put(&sq, sq_state[proc], sizeof(sq), proc);
  //FIXME: does unlock imply an ARMCI_Fence. If not, add one
  ARMCI_Unlock(0, proc);

  addTask(buf, tsk_size);
  delete [] buf;
  return true;

error:
  throw TSL_ERR;
}

bool
SharedQueue::hasTerminated() {
  return td.hasTerminated();
}

void
SharedQueue::td_progress() {
  td.progress(*dirty == 0);
  *dirty = 0;
}

