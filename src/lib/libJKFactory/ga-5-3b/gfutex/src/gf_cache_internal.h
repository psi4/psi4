#ifndef __GF_CACHE_INTERNAL_H
#define __GF_CACHE_INTERNAL_H

#include <string.h>

#include <fstream>
#include <iostream>

#include <ga.h>

#include <tbb/concurrent_hash_map.h>

#include "gfutex.h"

namespace globalFutures {
  extern std::ofstream ofs;
}

namespace globalFutures_implementation {

  using namespace tbb;
  using namespace globalFutures;

  class GASliceKey {
  protected:
    int g_a;
    int *lo;
    int *hi;
    int *ld;

  public:
    GASliceKey() : g_a(0), lo(NULL), hi(NULL), ld(NULL) { };
    GASliceKey(const int og_a, const int *olo, const int *ohi, const int *old);
    GASliceKey(const GASliceKey &othis);
    ~GASliceKey();

    int get_g_a() const { return g_a; };
    int *get_lo() const { return lo; };
    int *get_hi() const { return hi; };
    int *get_ld() const { return ld; };

    friend class GASliceKeyCompare;
  };

  class GASliceCache {
  protected:
    char *data;
    size_t sz;

  public:
    GASliceCache() : data(NULL), sz(0U) { };
    GASliceCache(const GASliceCache &othis);
    ~GASliceCache();

    char *&getData() { return data; };
    char *getData() const { return data; };
    size_t getSize() const { return sz; };
    size_t &getSize() { return sz; };
  };

  class GASliceKeyCompare {
  public:
    GASliceKeyCompare() { };
    GASliceKeyCompare(const GASliceKeyCompare &othis) { };
    ~GASliceKeyCompare() { };

    bool equal(const GASliceKey &k, const GASliceKey &j) const;
    size_t hash(const GASliceKey &k) const;
  };

  class GASliceNBInfo {
  protected:
    int g_a;
    int *lo;
    int *hi;
    void *buf;
    int *ld;
    bool local;
    bool incache;

  public:
    GASliceNBInfo() : g_a(0), lo(NULL), hi(NULL), buf(NULL), ld(NULL), local(false),
		      incache(false) { };
    GASliceNBInfo(int og_a, int *olo, int *ohi, void *obuf, int *old, bool olocal,
		  bool oincache);
    GASliceNBInfo(const GASliceNBInfo &othis) :
      g_a(0), lo(NULL), hi(NULL), buf(NULL), ld(NULL), local(false), incache(false)
    { *this = othis; };
    ~GASliceNBInfo();

    GASliceNBInfo &operator = (const GASliceNBInfo &othis);

    int get_g_a() const { return g_a; };
    int *get_lo() const { return lo; };
    int *get_hi() const { return hi; };
    void *get_buf() const { return buf; };
    int *get_ld() const { return ld; };
    bool get_local() const { return local; };
    bool get_incache() const { return incache; };
  };

  typedef concurrent_hash_map<GASliceKey, GASliceCache, GASliceKeyCompare> GAInstanceCacheMap;

  struct GAInstanceCache {
    GAInstanceCacheMap cmap;

    GAInstanceCache() { };
    GAInstanceCache(const GAInstanceCache &othis) : cmap(othis.cmap) { };
  };

  struct GAInstanceStats {
    unsigned int hits;
    unsigned int misses;

    GAInstanceStats() : hits(0U), misses(0U) { };
  };

  typedef concurrent_hash_map<int, GAInstanceCache> GACache;
  typedef concurrent_hash_map<ga_nbhdl_t *, GASliceNBInfo> GACacheNBOps;
  typedef concurrent_hash_map<int, GAInstanceStats> GACacheStats;

  inline GASliceKey::GASliceKey(const int og_a, const int *olo, const int *ohi, const int *old)
    : g_a(og_a)
  {
    int ndim = GF_Ndim(g_a);

    lo = new int[ndim];
    hi = new int[ndim];
    ld = new int[ndim - 1];

    for (int i = 0; i < ndim; i++) {
      lo[i] = olo[i];
      hi[i] = ohi[i];
    }

    for (int i = 0; i < ndim - 1; i++)
      ld[i] = old[i];
  } // GASliceKey

  inline GASliceKey::GASliceKey(const GASliceKey &othis)
    : g_a(othis.g_a)
  {
    int ndim = GF_Ndim(g_a);

    lo = new int[ndim];
    hi = new int[ndim];
    ld = new int[ndim - 1];

    for (int i = 0; i < ndim; i++) {
      lo[i] = othis.lo[i];
      hi[i] = othis.hi[i];
    }

    for (int i = 0; i < ndim - 1; i++)
      ld[i] = othis.ld[i];
  } // GASliceKey

  inline GASliceKey::~GASliceKey()
  {
    if (lo)
      delete [] lo;

    if (hi)
      delete [] hi;

    if (ld)
      delete [] ld;
  } // ~GASliceKey

  inline GASliceCache::GASliceCache(const GASliceCache &othis)
    : data(NULL), sz(othis.sz)
  {
    if (sz > 0) {
      data = new char[sz];

      memcpy(data, othis.data, sz);
    }
  } // GASliceCache

  inline GASliceCache::~GASliceCache()
  {
    if (data)
      delete [] data;
  } // ~GASliceCache

  inline bool GASliceKeyCompare::equal(const GASliceKey &k, const GASliceKey &j) const
  {
    if (k.g_a != j.g_a)
      return false;

    int ndim = GF_Ndim(k.g_a);

    for (int i = 0; i < ndim; i++) {
      if (k.lo[i] != j.lo[i])
	return false;

      if (k.hi[i] != j.hi[i])
	return false;
    }

    return true;
  } // equal

  inline size_t GASliceKeyCompare::hash(const GASliceKey &k) const
  {
    int ndim = GF_Ndim(k.g_a);
    size_t curr = tbb_hasher(k.g_a);

    for (int i = 0; i < ndim; i++)
      curr ^= (tbb_hasher(k.lo[i]) ^ tbb_hasher(k.hi[i]));

    return curr;
  } // hash

  inline GASliceNBInfo::GASliceNBInfo(int og_a, int *olo, int *ohi, void *obuf, int *old,
				      bool olocal, bool oincache)
    : g_a(og_a), buf(obuf), local(olocal), incache(oincache)
  {
    assert(NGA_Valid_handle(g_a));
    int ndim = GF_Ndim(g_a);

    lo = new int[ndim];
    hi = new int[ndim];
    ld = new int[ndim - 1];

    for (int i = 0; i < ndim; i++) {
      lo[i] = olo[i];
      hi[i] = ohi[i];
    }

    for (int i = 0; i < ndim - 1; i++)
      ld[i] = old[i];
  } // GASliceNBInfo

  inline GASliceNBInfo::~GASliceNBInfo()
  {
    if (lo)
      delete [] lo;

    if (hi)
      delete [] hi;

    if (ld)
      delete [] ld;
  } // ~GASliceNBInfo

  inline GASliceNBInfo &GASliceNBInfo::operator = (const GASliceNBInfo &othis)
  {
    g_a = othis.g_a;
    buf = othis.buf;
    local = othis.local;
    incache = othis.incache;

    if (lo)
      delete [] lo;

    if (hi)
      delete [] hi;

    if (ld)
      delete ld;


    lo = NULL;
    hi = NULL;
    ld = NULL;

    if (NGA_Valid_handle(g_a)) {
      int ndim = GF_Ndim(g_a);

      lo = new int[ndim];
      hi = new int[ndim];
      ld = new int[ndim - 1];

      for (int i = 0; i < ndim; i++) {
	lo[i] = othis.lo[i];
	hi[i] = othis.hi[i];
      }

      for (int i = 0; i < ndim - 1; i++)
	ld[i] = othis.ld[i];
    }

    return *this;
  } // operator =
}

#endif // __GF_CACHE_INTERNAL_H
