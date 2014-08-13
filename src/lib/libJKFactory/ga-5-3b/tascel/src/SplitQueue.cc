#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <algorithm>
#include <cstdio>
#include <cstring>

#include <armci.h>

#include "Comm.h"
#include "massert.h"
#include "SplitQueue.h"

using namespace std;
using namespace tascel;
using namespace tascel::comm;


SplitQueue::SplitQueue(int _tsk_size, int _max_ntsks)
  : max_ntsks(_max_ntsks)
  , tsk_size(_tsk_size)
  , q(NULL)
  , sq_state(NULL)
  , head(NULL)
  , tail(NULL)
  , size_local(NULL)
  , size_shared(NULL)
  , split(NULL)
  , dirty(NULL)
  , size_local_val(0)
  , head_val(0)
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
  tail = &sq_state[me()]->tail;
  split = &sq_state[me()]->split;
  head = &head_val;
  size_local = &size_local_val;
  size_shared = &sq_state[me()]->size_shared;
  dirty = &sq_state[me()]->dirty;
  *head = *tail = *split = *size_local = *size_shared = *dirty = 0;
  //printf("(cons) %d: *head=%p *tail=%p *size=%p\n", me(), head, tail, size);
  return;

error2:
  printf("Error2 in SplitQueue constructor\n");
  delete [] sq_state;
error1:
  printf("Error1 in SplitQueue constructor\n");
  delete [] q;
error:
  printf("Error in SplitQueue constructor\n");
  /*FIXME: cannot throw exception in constructor!*/
}

SplitQueue::~SplitQueue() {
  barrier();
  ARMCI_Free(sq_state[me()]);
  ARMCI_Destroy_mutexes();
  ARMCI_Free(q[me()]);
  delete [] sq_state;
  delete [] q;
}

bool
SplitQueue::empty() const {
  /*SK: speculative implementation. It could become empty while this
    operation is ongoing. However, the converse is not true, since
    only the host process call empty() or can add tasks (thus
    incrementing size) */
  massert(*size_local >= 0 && *size_shared >= 0);
  if (*size_local == 0 && *size_shared == 0) {
    massert(*split == *head);
    massert(*split == *tail);
  }

  return (*size_local == 0 && *size_shared == 0);

error:
  throw TSL_ERR;
}

bool
SplitQueue::full() const {
  if (*head == *tail && (*split != *head || *split != *tail)) {
    return true;
  }
  return false;
}

/**
 * If local portion has tasks, get one. If not, lock and try to find
 * work in the shared portion. If no work in shared portion as well,
 * return empty handed.
 */
bool
SplitQueue::getTask(void *dscr, int dlen) {
  massert(dscr);
  massert(dlen == tsk_size);
  if (*size_local > 0) {
    if (*size_shared == 0) {
      releaseToShared(*size_local / 2);
    }
    *head = (*head - 1 + max_ntsks) % max_ntsks;
    memcpy(dscr, &q[me()][*head*tsk_size], tsk_size);
    *size_local -= 1;
    return true;
  }

  if (acquireFromShared(1)) {
    //SK: note that this breaks with a push-based model. We assert
    //that when this condition is true, we will not recurse to get
    //work from shared again.
    massert(getTask(dscr, dlen) == true);
    return true;
  }
  return false;

error:
  throw TSL_ERR;
}

int
SplitQueue::acquireFromShared(int _nacquire) {
  int nacquire;
  if (_nacquire == 0) {
    return 0;
  }
  ARMCI_Lock(0, me()); //SK: is it needed?
  if (*size_shared == 0) {
    ARMCI_Unlock(0, me());
    return 0;
  }
  massert(*split != *tail);
  massert(*head == *split);
  nacquire = min(_nacquire, *size_shared);
  *split = (*split - nacquire + max_ntsks) % max_ntsks;
  massert((*split <= *head) || (*split > *tail));
  *size_shared -= nacquire;
  *size_local += nacquire;
  ARMCI_Unlock(0, me());
  return nacquire;

error:
  throw TSL_ERR;
}

void
SplitQueue::releaseToShared(int ndonate) {
  massert(ndonate <= *size_local);
  if (ndonate == 0) {
    return;
  }
  ARMCI_Lock(0, me());
  *split = (*split + ndonate) % max_ntsks;
  massert((*split <= *head) || (*split > *tail));
  //mem-fence
  *size_shared += ndonate;
  *size_local -= ndonate;
  ARMCI_Unlock(0, me());
  return;

error:
  throw TSL_ERR;
}

/**
 * Check is there is space. If not, wait. If there is, add at head.
 */
void
SplitQueue::addTask(void *dscr, int dlen) {
  massert(dscr);
  massert(dlen == tsk_size);
  while (full()) {}
  memcpy(&q[me()][*head*tsk_size], dscr, tsk_size);
  *head  = (*head + 1) % max_ntsks;
  *size_local += 1;
  if (*size_shared == 0) {
    releaseToShared(*size_local / 2);
  }
  return;

error:
  throw TSL_ERR;
}

bool
SplitQueue::steal(int proc) {
  char *buf;
  massert(*size_local == 0);
  massert(!full());
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
  if (sq.size_shared == 0) {
    //nothing to steal
    ARMCI_Unlock(0, proc);
    return false;
  }
  sq.size_shared = sq.size_shared - 1;
  sq.dirty = 1; //mark dirty for termination detection
  ARMCI_Get(&q[proc][sq.tail*tsk_size], buf, tsk_size, proc);
  sq.tail = (sq.tail + 1) % max_ntsks;
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
SplitQueue::hasTerminated() {
  return td.hasTerminated();
}

void SplitQueue::td_progress() {  
  td.progress(*dirty == 0);
  *dirty = 0;
}
