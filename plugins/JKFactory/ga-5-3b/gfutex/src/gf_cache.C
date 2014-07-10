#include <assert.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include <mpi.h>

#include <ga.h>
#include <matypes.h>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/tick_count.h>

#include "gfutex.h"
#include "gfutex_internal.h"
#include "gf_cache_internal.h"

#undef CACHE_STATS

namespace globalFutures {
  using namespace globalFutures_implementation;

  typedef tbb::enumerable_thread_specific<double> TSum;

  static GACache cache;
  static GACacheNBOps nbops;

#ifdef CACHE_STATS
  static GACacheStats cacheStats;
#endif

  template <typename T> void accumulate(T *cbuf, const T *buf, size_t nelems);

  static bool is_slice_local(int g_a, int lo[], int hi[]);
  static size_t type_size(int type);
  static size_t slice_size(int g_a, int lo[], int hi[]);
  static bool get_accessor(int g_a, int lo[], int hi[], int ld[],
			   GAInstanceCacheMap::const_accessor &acc);
  static void get_accessor(int g_a, int lo[], int hi[], int ld[],
			   GAInstanceCacheMap::accessor &acc);
  static bool try_read_cache(int g_a, int lo[], int hi[], int ld[], void *buf);
  static void local_accumulate(int g_a, size_t sz, void *cbuf, const void *buf,
			       const void *alpha);
  static void typed_accumulate(int type, void *cbuf, const void *buf, size_t nelems,
			       const void *alpha);
  static void create_one(int type, void *alpha);
  static void dump_cache(std::ofstream &ofs, const int line);

  std::ofstream ofs;

  void GF_CachedGet(int g_a, int lo[], int hi[], void *buf, int ld[])
  {
#if 0
    if (!ofs.is_open()) {
      std::ostringstream ostr;

      ostr << "diags.dat." << me;
      ofs.open(ostr.str().c_str(), std::ios_base::out | std::ios_base::app);
    }
#endif

    if (is_slice_local(g_a, lo, hi)) {
      GF_Get(g_a, lo, hi, buf, ld);
      return;
    }

    // Optimistic cache search, assume the slice is there...
    bool incache = try_read_cache(g_a, lo, hi, ld, buf);

    if (incache)
      return;

    GAInstanceCacheMap::accessor inst_mod_acc;

    get_accessor(g_a, lo, hi, ld, inst_mod_acc);

    char *&dat = inst_mod_acc->second.getData();
    size_t &sz = inst_mod_acc->second.getSize();

    size_t ssz = slice_size(g_a, lo, hi);

    if (!dat) {
      dat = new char[ssz];
      sz = ssz;

      GF_Get(g_a, lo, hi, buf, ld);

      memcpy(dat, buf, ssz);
    }
    else {
      assert(dat && (sz == ssz));

      memcpy(buf, dat, sz);
    }
  } // GF_CachedGet

  void GF_CachedAcc(int g_a, int lo[], int hi[], void *buf, int ld[], void *alpha)
  {
    if (is_slice_local(g_a, lo, hi)) {
      GF_Acc(g_a, lo, hi, buf, ld, alpha);
      return;
    }

    // Optimistic cache search, assume the slice is there...
    GAInstanceCacheMap::accessor inst_mod_acc;

    get_accessor(g_a, lo, hi, ld, inst_mod_acc);

    char *&dat = inst_mod_acc->second.getData();
    size_t &sz = inst_mod_acc->second.getSize();

    if (!dat) {
      size_t ssz = slice_size(g_a, lo, hi);

      dat = new char[ssz];
      sz = ssz;

      memset(dat, 0, ssz);
    }
    else
      assert(dat && (sz == slice_size(g_a, lo, hi)));

    local_accumulate(g_a, sz, dat, buf, alpha);
  } // GF_CachedAcc

  void GF_CachedNbGet(int g_a, int lo[], int hi[], void *buf, int ld[], ga_nbhdl_t *nbhandle)
  {
    try {
      bool is_local = is_slice_local(g_a, lo, hi);
      bool incache = try_read_cache(g_a, lo, hi, ld, buf);

      GACacheNBOps::const_accessor cacc;

      bool found = nbops.find(cacc, nbhandle);

      if (found) {
	// ofs << "Non-blocking operation found in cache: " << __LINE__ << std::flush << std::endl;
	throw GFutException("Non-blocking operation should not be in the cache!");
      }

      GASliceNBInfo nbinfo(g_a, lo, hi, buf, ld, is_local, incache);

      GACacheNBOps::accessor acc;

      bool ins = nbops.insert(acc, nbhandle);

      if (!ins) {
	// ofs << "Non-blocking operation found in cache: " << __LINE__ << std::flush << std::endl;
	throw GFutException("Non-blocking operation should not be in the cache!");
      }

      // ofs << "Inserted nbop in cache: " << nbhandle << ", is_local: " << is_local <<
      // ", incache: " << incache << std::flush << std::endl;

      acc->second = nbinfo;

      if (!incache || is_local)
	GF_NbGet(g_a, lo, hi, buf, ld, nbhandle);
    }
    catch (const GFutException &excp) {
      std::cerr << excp.what() << std::endl;
      exit(EXIT_FAILURE);
    }
  } // GF_CachedNbGet

  void GF_CachedNbWait(ga_nbhdl_t *nbhandle)
  {
    try {
      GACacheNBOps::const_accessor acc;

      bool found = nbops.find(acc, nbhandle);

      if (!found) {
	// ofs << "Problems with nbop: " << nbhandle << std::flush << std::endl;
	throw GFutException("Non-blocking operation should be in the cache!");
      }

      bool is_local = acc->second.get_local();
      bool incache = acc->second.get_incache();

      if (is_local || !incache)
	GF_NbWait(nbhandle);

      if (!is_local && !incache) {
	size_t ssz = slice_size(acc->second.get_g_a(), acc->second.get_lo(),
				acc->second.get_hi());

	GAInstanceCacheMap::accessor inst_acc;

	get_accessor(acc->second.get_g_a(), acc->second.get_lo(), acc->second.get_hi(),
		     acc->second.get_ld(), inst_acc);

	char *&dat = inst_acc->second.getData();
	size_t &sz = inst_acc->second.getSize();

	if (!dat) {
	  dat = new char[ssz];
	  sz = ssz;

	  memcpy(dat, acc->second.get_buf(), ssz);
	}
	else
	  assert(dat && (sz == slice_size(acc->second.get_g_a(), acc->second.get_lo(),
					  acc->second.get_hi())));
      }

      nbops.erase(acc); // Erase record of non-blocking operation
    }
    catch (const GFutException &excp) {
      std::cerr << excp.what() << std::endl;
      exit(EXIT_FAILURE);
    }
  } // GF_CachedNbWait

  void GF_CacheReadOnlyEmpty(int g_a)
  {
    GACache::accessor mod_acc;

    bool found = cache.find(mod_acc, g_a);

    if (found) {
#ifdef CACHE_STATS
      GACacheStats::accessor sacc;
#endif

      mod_acc->second.cmap.clear(); // Erase all cache entries for GA instance g_a

#ifdef CACHE_STATS
      cacheStats.find(sacc, g_a);
      sacc->second.hits = 0U;
      sacc->second.misses = 0U;
#endif
    }
  } // GF_CacheReadOnlyEmpty

  void GF_CacheReadWriteFlush(int g_a)
  {
    GACache::accessor acc;

    bool found = cache.find(acc, g_a);

    if (!found)
      return;

#ifdef CACHE_STATS
    GACacheStats::accessor sacc;
#endif

    size_t sz = acc->second.cmap.size();
    ga_nbhdl_t *hdls = new ga_nbhdl_t[sz];
    int pos = 0;

    for (GAInstanceCacheMap::const_iterator iter = acc->second.cmap.begin();
         iter != acc->second.cmap.end(); iter++) {
      GF_NbPut(g_a, iter->first.get_lo(), iter->first.get_hi(), iter->second.getData(),
	       iter->first.get_ld(), &hdls[pos]);
      pos++;
    }

    for (size_t i = 0U; i < sz; i++)
      GF_NbWait(&hdls[pos]);

    acc->second.cmap.clear();

#ifdef CACHE_STATS
    cacheStats.find(sacc, g_a);
    sacc->second.hits = 0U;
    sacc->second.misses = 0U;
#endif

    delete [] hdls; 
  } // GF_CacheReadWriteFlush

  void GF_CacheAccFlush(int g_a)
  {
    GACache::accessor acc;

    bool found = cache.find(acc, g_a);

    if (!found)
      return;

    int dims[GA_MAX_DIM];
    int ndim, type;

    GF_Inquire(g_a, &type, &ndim, dims);

    char *alpha = new char[type_size(type)];

#ifdef CACHE_STATS
    GACacheStats::accessor sacc;
#endif

    size_t sz = acc->second.cmap.size();
    ga_nbhdl_t *hdls = new ga_nbhdl_t[sz];
    int pos = 0;

    create_one(type, alpha);

    for (GAInstanceCacheMap::const_iterator iter = acc->second.cmap.begin();
         iter != acc->second.cmap.end(); iter++) {
      GF_NbAcc(g_a, iter->first.get_lo(), iter->first.get_hi(), iter->second.getData(),
	       iter->first.get_ld(), alpha, &hdls[pos]);
      pos++;
    }

    for (size_t i = 0U; i < sz; i++)
      GF_NbWait(&hdls[pos]);

    acc->second.cmap.clear();

#ifdef CACHE_STATS
    cacheStats.find(sacc, g_a);
    sacc->second.hits = 0U;
    sacc->second.misses = 0U;
#endif

    delete [] alpha;
    delete [] hdls; 
  } // GF_CacheAccFlush

  unsigned int GF_CacheGetMisses(int g_a)
  {
#ifdef CACHE_STATS
    GACacheStats::const_accessor acc;

    bool found = cacheStats.find(acc, g_a);

    if (found)
      return acc->second.misses;
#endif

    return 0U; 
  } // GF_CacheGetMisses

  unsigned int GF_CacheGetHits(int g_a)
  {
#ifdef CACHE_STATS
    GACacheStats::const_accessor acc;

    bool found = cacheStats.find(acc, g_a);

    if (found)
      return acc->second.hits;
#endif

    return 0U; 
  } // GF_CacheGetHits

  bool is_slice_local(int g_a, int lo[], int hi[])
  {
    int ndims = GF_Ndim(g_a);

    int *map = new int[2 * ndims * nproc];
    int *procs = new int[nproc];
    bool ret = false;

    int np = GF_Locate_region(g_a, lo, hi, map, procs);

    if (np == 1 && procs[0] == me)
      ret = true;

    delete [] map;
    delete [] procs;

    return ret;
  } // is_slice_local

  size_t type_size(int type)
  {
    size_t tsz;

    switch (type) {
    case C_CHAR:
      tsz = sizeof(char);
      break;
    case C_DBL:
      tsz = sizeof(double);
      break;
    case C_DCPL:
      tsz = sizeof(MA_DoubleComplex);
      break;
    case C_FLOAT:
      tsz = sizeof(float);
      break;
    case C_INT:
      tsz = sizeof(int);
      break;
    case C_LDBL:
      tsz = sizeof(MA_LongDouble);
      break;
    case C_LDCPL:
      tsz = sizeof(MA_LongDoubleComplex);
      break;
    case C_LONGLONG:
      tsz = sizeof(long long int);
      break;
    case C_LONG:
      tsz = sizeof(long int);
      break;
    case C_SCPL:
      tsz = sizeof(MA_SingleComplex);
      break;
    default:
      std::cerr << "Proc: " << me << ", GA data type not found! " << type << std::endl;
      exit(EXIT_FAILURE);
    }

    return tsz;
  } // type_size

  size_t slice_size(int g_a, int lo[], int hi[])
  {
    int dims[GA_MAX_DIM];
    int ndim, type;
    size_t sz;

    GF_Inquire(g_a, &type, &ndim, dims);   

    sz = type_size(type);

    for (int i = 0; i < ndim; i++)
      sz *= (hi[i] - lo[i] + 1);

    return sz;
  } // slice_size

  bool get_accessor(int g_a, int lo[], int hi[], int ld[],
		    GAInstanceCacheMap::const_accessor &acc)
  {
    GACache::const_accessor accessor;

    cache.insert(accessor, g_a);

    GASliceKey slkey(g_a, lo, hi, ld);

    bool slice_found = accessor->second.cmap.find(acc, slkey);

    accessor.release();

#ifdef CACHE_STATS
    GACacheStats::accessor sacc;

    cacheStats.insert(sacc, g_a);

    if (slice_found)
      sacc->second.hits++;
    else
      sacc->second.misses++;
#endif

    return slice_found;
  } // get_accessor

  void get_accessor(int g_a, int lo[], int hi[], int ld[], GAInstanceCacheMap::accessor &acc)
  {
    GACache::accessor mod_acc;

    cache.insert(mod_acc, g_a);

    GASliceKey slkey(g_a, lo, hi, ld);

    bool ins = mod_acc->second.cmap.insert(acc, slkey);

    mod_acc.release();

#ifdef CACHE_STATS
    GACacheStats::accessor sacc;

    cacheStats.insert(sacc, g_a);

    if (ins)
      sacc->second.misses++;
    else
      sacc->second.hits++;
#endif
  } // get_accessor

  bool try_read_cache(int g_a, int lo[], int hi[], int ld[], void *buf)
  {
    GAInstanceCacheMap::const_accessor inst_acc;

    bool slice_found = get_accessor(g_a, lo, hi, ld, inst_acc);

    if (slice_found) {
      size_t sz = inst_acc->second.getSize();

      memcpy(buf, inst_acc->second.getData(), sz);
    }

    return slice_found;
  } // try_read_cache

  template <typename T> void accumulate(T *cbuf, const T *buf, size_t nelems, const T *alpha)
  {
    for (size_t i = 0; i < nelems; i++)
      cbuf[i] += (*alpha * buf[i]);
  } // accumulate

  void local_accumulate(int g_a, size_t sz, void *cbuf, const void *buf, const void *alpha)
  {
    int dims[GA_MAX_DIM];
    int ndim, type;
    size_t nelems;

    GF_Inquire(g_a, &type, &ndim, dims);

    nelems = sz / type_size(type);

    typed_accumulate(type, cbuf, buf, nelems, alpha);
  } // local_accumulate

  void typed_accumulate(int type, void *cbuf, const void *buf, size_t nelems,
			const void *alpha)
  {
    switch (type) {
    case C_CHAR:
      accumulate(static_cast<char *>(cbuf), static_cast<const char *>(buf), nelems,
		 static_cast<const char *>(alpha));
      break;
    case C_DBL:
      accumulate(static_cast<double *>(cbuf), static_cast<const double *>(buf), nelems,
		 static_cast<const double *>(alpha));
      break;
    case C_FLOAT:
      accumulate(static_cast<float *>(cbuf), static_cast<const float *>(buf), nelems,
		 static_cast<const float *>(alpha));
      break;
    case C_INT:
      accumulate(static_cast<int *>(cbuf), static_cast<const int *>(buf), nelems,
		 static_cast<const int *>(alpha));
      break;
    case C_LONGLONG:
      accumulate(static_cast<long long int *>(cbuf), static_cast<const long long int *>(buf),
		 nelems, static_cast<const long long int *>(alpha));
      break;
    case C_LONG:
      accumulate(static_cast<long int *>(cbuf), static_cast<const long int *>(buf), nelems,
		 static_cast<const long int *>(alpha));
      break;
    default:
      std::cerr << "Proc: " << me << ", unsupported type for cached accumulate: " << type <<
	std::endl;
      exit(EXIT_FAILURE);
    }
  } // typed_accumulate

  void create_one(int type, void *alpha)
  {
    char one[sizeof(long long int)];

    memset(one, 0, sizeof(one));

    switch (type) {
    case C_CHAR:
      one[0] = 1;
      memcpy(alpha, one, sizeof(char));
      break;
    case C_DBL:
      *(reinterpret_cast<double *>(one)) = 1.0;
      memcpy(alpha, one, sizeof(double));
      break;
    case C_FLOAT:
      *(reinterpret_cast<float *>(one)) = 1.0f;
      memcpy(alpha, one, sizeof(float));
      break;
    case C_INT:
      *(reinterpret_cast<int *>(one)) = 1;
      memcpy(alpha, one, sizeof(int));
      break;
    case C_LONGLONG:
      *(reinterpret_cast<long long int *>(one)) = 1L;
      memcpy(alpha, one, sizeof(long long int));
      break;
    case C_LONG:
      *(reinterpret_cast<long int *>(one)) = 1L;
      memcpy(alpha, one, sizeof(long int));
      break;
    default:
      std::cerr << "Proc: " << me << ", unsupported type for cached accumulate: " << type <<
	std::endl;
      exit(EXIT_FAILURE);
    }
  } // create_one

  void dump_cache(std::ofstream &ofs, const int line)
  {
    for (GACache::const_iterator iter = cache.begin(); iter != cache.end(); iter++) {
      int ndim = GF_Ndim(iter->first);

      ofs << line << " slices for g_a: " << iter->first << ", number: " <<
        iter->second.cmap.size() << std::endl << std::endl;

      for (GAInstanceCacheMap::const_iterator iiter = iter->second.cmap.begin();
	   iiter != iter->second.cmap.end(); iiter++) {
	ofs << "g_a: << " << iiter->first.get_g_a() << ", lo: ";
	for (int i = 0; i < ndim - 1; i++)
	  ofs << iiter->first.get_lo()[i] << " ";

	ofs << iiter->first.get_lo()[ndim - 1] << ", hi: ";

	for (int i = 0; i < ndim - 1; i++)
	  ofs << iiter->first.get_hi()[i] << " ";

	ofs << iiter->first.get_hi()[ndim - 1] << ", sz: " << iiter->second.getSize();
	ofs << ", dat: " << std::hex << static_cast<void *>(iiter->second.getData()) <<
	  std::dec;
	ofs << std::endl;
      }
    }

    ofs << std::flush;
  } // dump_cache
}
