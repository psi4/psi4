/**
 * This is an optimized implementation of SplitQueue, based on the
 * algorithm presented in:
 * "Scalable work stealing", SC'09.
 *
 * This version includes support for tasks with data movement.
 *
 */

#include "SplitQueueOpt.h"
#include <armci.h>
#include "Comm.h"
#include <massert.h>
#include <cstring>
#include <cstdio>
#include <algorithm>

using namespace std;
using namespace tascel;
using namespace comm;

SplitQueueOpt::SplitQueueOpt(int _tsk_size, int _max_ntsks)
  : tsk_size(_tsk_size), max_ntsks(_max_ntsks) {
  q = new char *[nproc()];
  massert(q);
  ARMCI_Malloc((void **)q, max_ntsks * tsk_size);
  massertl(q[me()], error1);
  ARMCI_Create_mutexes(1);
  sq_state = new sq_state_t*[nproc()];
  massertl(sq_state, error1);
  ARMCI_Malloc((void **)sq_state, nproc()*sizeof(sq_state_t));
  massertl(sq_state[me()], error2);
  tail = &sq_state[me()]->ps.tail;
  split = &sq_state[me()]->ps.split;
  size_shared = &sq_state[me()]->ps.size_shared;
  dirty = &sq_state[me()]->ps.dirty;
  vtail2 = &sq_state[me()]->vtail2;
  head = &head_val;
  size_local = &size_local_val;
  vtail1 = &vtail1_val;

  *head = *tail = *vtail1 = *vtail2 = *split = *size_local = *size_shared = *dirty = 0;
  ARMCI_Barrier(); //for the mem fence. we don't want steals before this is init-ed
  //   printf("(cons) %d: *head=%p *tail=%p *size=%p\n", me(), head, tail, size);
  return;

error2:
  printf("Error2 in SplitQueueOpt constructor\n");
  delete [] sq_state;
error1:
  printf("Error1 in SplitQueueOpt constructor\n");
  delete [] q;
error:
  printf("Error in SplitQueueOpt constructor\n");
  /*FIXME: cannot throw exception in constructor!*/
}

SplitQueueOpt::~SplitQueueOpt() {
  barrier();
  ARMCI_Free(sq_state[me()]);
  ARMCI_Destroy_mutexes();
  ARMCI_Free(q[me()]);
  delete [] sq_state;
  delete [] q;
}

bool
SplitQueueOpt::empty() const {
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
SplitQueueOpt::full() const {
  if (*head == *tail && (*split != *head || *split != *tail)) {
    return true;
  }
  return false;
}

bool
SplitQueueOpt::spaceAvailable() const {
  if (*head == *vtail1 && (*split != *head || *split != *vtail1)) {
    return false;
  }
  return true;
}

/**
 * If local portion has tasks, get one. If not, lock and try to find
 * work in the shared portion. If no work in shared portion as well,
 * return empty handed.
 */
bool
SplitQueueOpt::getTask(void *dscr, int dlen) {
  massert(dscr);
  massert(dlen == tsk_size);
  if (*size_local > 0) {
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
SplitQueueOpt::acquireFromShared(int _nacquire) {
  int nacquire;
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
SplitQueueOpt::releaseToShared(int ndonate) {
  massert(ndonate <= *size_local);
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

void
SplitQueueOpt::updateVtail1() {
  if (*vtail2 == *tail) {
    *vtail1 = *tail;
  }
}


/**
 * Check is there is space. If not, wait. If there is, add at head.
 */
void
SplitQueueOpt::addTask(void *dscr, int dlen) {
  massert(dscr);
  massert(dlen == tsk_size);

  while (!spaceAvailable()) {
    updateVtail1();
  }
  memcpy(&q[me()][*head*tsk_size], dscr, tsk_size);
  *head  = (*head + 1) % max_ntsks;
  *size_local += 1;
  if (*size_shared == 0) {
    releaseToShared(1);
  }
  return;

error:
  throw TSL_ERR;
}

bool
SplitQueueOpt::steal(int proc) {
  char *buf;
  int newtail, nshift, oldtail;
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
  if (sq.ps.size_shared == 0) {
    //nothing to steal
    ARMCI_Unlock(0, proc);
    return false;
  }
  sq.ps.size_shared = sq.ps.size_shared - 1;
  sq.ps.dirty = 1; //mark dirty for termination detection
  newtail = (sq.ps.tail + 1) % max_ntsks;
  oldtail = sq.ps.tail;
  nshift = newtail - sq.ps.tail;
  sq.ps.tail = newtail;
  ARMCI_Put(&sq.ps, sq_state[proc], sizeof(sq), proc);
  //FIXME: does unlock imply an ARMCI_Fence. If not, add one
  ARMCI_Unlock(0, proc);

  ARMCI_Get(&q[proc][oldtail*tsk_size], buf, tsk_size, proc);
  int dummy;
  ARMCI_Rmw(ARMCI_FETCH_AND_ADD, &dummy, &sq_state[proc]->vtail2, nshift, proc);

  addTask(buf, tsk_size);
  delete [] buf;
  return true;

error:
  throw TSL_ERR;
}

bool
SplitQueueOpt::hasTerminated() {
  return td.hasTerminated();
}
