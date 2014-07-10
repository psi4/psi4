#include "DenseArray.h"
#include "massert.h"
#include "ga.h"

#include <cstdio>
#include <cassert>
#include <algorithm>

using namespace std;
using namespace tascel;

int
DenseArray::getProc(const void *idx, int idxlen) {
  int lo[ndim], hi[ndim], ld[ndim], proc;
  massert(idxlen == ndim * sizeof(int));

  computeBounds(idx, idxlen, lo, hi, ld);
#define OWNER_LOWER_LEFT 1
#if OWNER_LOWER_LEFT
  // Oversimplification. The owner of the block is considered to be the proc
  // which owns the lowest index.
  proc = NGA_Locate(ga, lo);
  massert(proc >= 0);
#else
  // The owner is the proc which owns most of the data, if more than one proc
  // is found to own the data.
  {
    const int nnodes = GA_Nnodes();
    // not optimal since nnodes could be huge
    int nown, map[nnodes*2*ndim], procs[nnodes];
    nown = NGA_Locate_region(ga, lo, hi, map, procs);
    if (1 == nown) {
      proc = procs[0];
    } else {
      long biggest_size = 0;
      for (int i=0; i<nown; ++i) {
        long current_size = 1;
        for (int j=0; j<ndim; ++j) {
          current_size *= map[i*2*ndim+ndim+j] - map[i*2*ndim+j] + 1;
          assert(current_size > 0);
        }
        if (current_size > biggest_size) {
          biggest_size = current_size;
          proc = procs[i];
        }
      }
      assert(biggest_size > 0);
    }
  }
#endif

  return proc;

error:
  throw TSL_ERR;
}

void
DenseArray::computeLo(const void *idx, int idxlen, int *lo) {
  int *iidx = (int *)idx;
  massert(idxlen = ndim * sizeof(int));
  for (int i = 0; i < ndim; i++) {
    lo[i] = iidx[i] * block[i];
  }
  return;

error:
  throw TSL_ERR;
}

void
DenseArray::computeHi(const void *idx, int idxlen, int *hi) {
  int *iidx = (int *)idx;
  massert(idxlen = ndim * sizeof(int));
  for (int i = 0; i < ndim; i++) {
    hi[i] = min((iidx[i] + 1) * block[i], dims[i]) - 1;
  }
  return;

error:
  throw TSL_ERR;
}

void
DenseArray::computeLd(const void *idx, int idxlen, int *ld) {
  int lo[ndim], hi[ndim];
  massert(idxlen == ndim * sizeof(int));

  computeLo(idx, idxlen, lo);
  computeHi(idx, idxlen, hi);

  for (int i = 0; i < ndim; i++) {
    ld[i] = hi[i] - lo[i] + 1;
  }
  return;

error:
  throw TSL_ERR;
}

void
DenseArray::computeBounds(const void *idx, int idxlen, int *lo, int *hi, int *ld) {
  int *iidx = (int *)idx;
  massert(idxlen = ndim * sizeof(int));
  for (int i = 0; i < ndim; i++) {
    lo[i] = iidx[i] * block[i];
    hi[i] = min((iidx[i] + 1) * block[i], dims[i]) - 1;
    ld[i] = hi[i] - lo[i] + 1;
  }
  return;

error:
  throw TSL_ERR;
}

DenseArray::DenseArray(int _ga, const int *_block, int len)
  : ga(_ga),
    ndim(GA_Ndim(_ga)) {
  int _dims[ndim], tmp, type;
  //ndim = GA_Ndim(ga);
  assert(len == ndim);
  NGA_Inquire(ga, &type, &tmp, _dims);
  assert(tmp == ndim);

  for (int i = 0; i < ndim; i++) {
    block.push_back(_block[i]);
    dims.push_back(_dims[i]);
  }
}

int
DenseArray::getSize(const void *idx, int idxlen) {
  int lo[ndim], hi[ndim], ld[ndim], sz;
  massert(idxlen == ndim * sizeof(int));

  computeBounds(idx, idxlen, lo, hi, ld);

  sz = sizeof(double);
  for (int i = 0; i < ndim; i++) {
    sz *= hi[i] - lo[i] + 1;
  }
  return sz;

error:
  throw TSL_ERR;
}

void
DenseArray::get(const void *idx, int idxlen, void *buf, int buflen) {
  int lo[ndim], hi[ndim], ld[ndim];
  massert(idxlen == ndim * sizeof(int));

  computeBounds(idx, idxlen, lo, hi, ld);
  massert(getSize(idx, idxlen) == buflen);

  NGA_Get(ga, lo, hi, buf, ld);
  return;

error:
  throw TSL_ERR;
}

void
DenseArray::put(const void *idx, int idxlen, const void *buf, int buflen) {
  int lo[ndim], hi[ndim], ld[ndim];
  massert(idxlen == ndim * sizeof(int));

  computeBounds(idx, idxlen, lo, hi, ld);
  massert(getSize(idx, idxlen) == buflen);

  NGA_Put(ga, lo, hi, (void *)buf, ld);
  return;

error:
  throw TSL_ERR;
}

void
DenseArray::add(const void *idx, int idxlen, const void *buf, int buflen) {
  int lo[ndim], hi[ndim], ld[ndim];
  double alpha = 1.0;
  massert(idxlen == ndim * sizeof(int));

  computeBounds(idx, idxlen, lo, hi, ld);
  massert(getSize(idx, idxlen) == buflen);

  NGA_Acc(ga, lo, hi, (void *)buf, ld, &alpha);
  return;

error:
  throw TSL_ERR;
}

DenseArray::nbh_t
DenseArray::nbGet(const void *idx, int idxlen, void *buf, int buflen, nbh_t nbh) {
  int lo[ndim], hi[ndim], ld[ndim];
  ga_nbhdl_t *rnbh;
  massert(idxlen == ndim * sizeof(int));

  rnbh = new ga_nbhdl_t();
  massert(rnbh!=NULL);

  computeBounds(idx, idxlen, lo, hi, ld);
  massert(getSize(idx, idxlen) == buflen);

  NGA_NbGet(ga, lo, hi, buf, ld, rnbh);
  return rnbh;

error:
  throw TSL_ERR;
}

DenseArray::nbh_t
DenseArray::nbPut(const void *idx, int idxlen, const void *buf, int buflen, nbh_t nbh) {
  int lo[ndim], hi[ndim], ld[ndim];
  ga_nbhdl_t *rnbh;
  massert(idxlen == ndim * sizeof(int));

  rnbh = new ga_nbhdl_t();
  massert(rnbh!=NULL);

  computeBounds(idx, idxlen, lo, hi, ld);
  massert(getSize(idx, idxlen) == buflen);

  NGA_NbPut(ga, lo, hi, (void *)buf, ld, rnbh);
  return rnbh;

error:
  throw TSL_ERR;
}

DenseArray::nbh_t
DenseArray::nbAdd(const void *idx, int idxlen, const void *buf, int buflen, nbh_t nbh) {
  int lo[ndim], hi[ndim], ld[ndim];
  ga_nbhdl_t *rnbh;
  double alpha = 1.0;
  massert(idxlen == ndim * sizeof(int));

  rnbh = new ga_nbhdl_t();
  massert(rnbh!=NULL);

  computeBounds(idx, idxlen, lo, hi, ld);
  massert(getSize(idx, idxlen) == buflen);

  NGA_NbAcc(ga, lo, hi, (void *)buf, ld, &alpha, rnbh);
  return rnbh;

error:
  throw TSL_ERR;
}

void
DenseArray::waitHandle(nbh_t nbh) {
  ga_nbhdl_t *rnbh = (ga_nbhdl_t*)nbh;
  NGA_NbWait(rnbh);
  delete rnbh;
}

