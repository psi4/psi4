#ifndef __DenseArray_h__
#define __DenseArray_h__

#include <vector>
#include "DataColl.h"

namespace tascel {

  /**
   * DataColl abstraction for a GlobalArray.
   *
   * A limited abstraction for GA data movement to be managed by the
   * tascel runtime. Note that the same GA could be wrapped by many
   * DenseArray objects, each providing a potentially different block
   * view of the same array.
   */
  class DenseArray : public DataColl {
  private:
    const int ga;            /**< GA being wrapped */
    const int ndim;          /**< Number of dimensions of GA */
    std::vector<int> block;  /**< The blocking for this wrapper*/
    std::vector<int> dims;   /**< The dimensions of the GA */
    
    /**
     * Compute low index for use in GA call.
     *
     * Low index for the block pointed to by idx.
     * @param[in]  idx    Pointer to index
     * @param[in]  idxlen Length of index (in bytes)
     * @param[out] lo     Lower bound for the block
     * @return            none
     */
    void computeLo(const void *idx, int idxlen, int *lo);

    /**
     * Compute high index for use in GA call.
     *
     * High index for the block pointed to by idx.
     * @param[in]  idx    Pointer to index
     * @param[in]  idxlen Length of index (in bytes)
     * @param[out] hi     Upper bound for the block
     * @return            none
     */
    void computeHi(const void *idx, int idxlen, int *hi);

    /**
     * Compute stride (ld) for use in GA call.
     *
     * High index for the block pointed to by idx.
     * @param[in]  idx    Pointer to index
     * @param[in]  idxlen Length of index (in bytes)
     * @param[out] ld     ld values for the block
     * @return            none
     */
    void computeLd(const void *idx, int idxlen, int *ld);

    /**
     * Combimes computation of lo, hi, and ld arguments of a GA call. 
     *
     * @param[in]  idx    Pointer to index
     * @param[in]  idxlen Length of index (in bytes)
     * @param[out] lo     lo values for the block
     * @param[out] hi     hi values for the block
     * @param[out] ld     ld values for the block
     * @return            none
     */
    void computeBounds(const void *idx, int idxlen, int *lo, int *hi, int *ld);
    
  public:

    /**
     * Construct a dense array object to wrap a GA
     *
     * A GA itself it not created. The user is still responsible to
     * pass a valid GA, and destroy the GA after this object if
     * destroyed. 
     * @param[in] ga    Handle to ga to be wrapped
     * @param[in] block The size of each block being put of get. 
     * @param[in] len   Number of integers pointed to by block. This
     * should match the number of dimensions of ga.
     */
    DenseArray(int ga, const int *block, int len);

    /**
     * @copybrief   DataColl::getSize(const void*,int)
     * @copydetails DataColl::getSize(const void*,int)
     */
    int getSize(const void *idx, int idxlen);

    /**
     * @copybrief   DataColl::getProc(const void*,int)
     * @copydetails DataColl::getProc(const void*,int)
     */
    int getProc(const void *idx, int idxlen);

    /**
     * @copybrief   DataColl::get(const void*,int,void*,int)
     * @copydetails DataColl::get(const void*,int,void*,int)
     */
    void get(const void *idx, int idxlen, void *buf, int buflen);

    /**
     * @copybrief   DataColl::put(const void*,int,const void*,int)
     * @copydetails DataColl::put(const void*,int,const void*,int)
     */
    void put(const void *idx, int idxlen, const void *buf, int buflen);

    /**
     * @copybrief   DataColl::add(const void*,int,const void*,int)
     * @copydetails DataColl::add(const void*,int,const void*,int)
     */
    void add(const void *idx, int idxlen, const void *buf, int buflen);

    /**
     * @copybrief   DataColl::nbGet(const void*,int,void*,int,nbh_t)
     * @copydetails DataColl::nbGet(const void*,int,void*,int,nbh_t)
     */
    virtual nbh_t nbGet(const void *idx, int idxlen, void *buf, int buflen, nbh_t nbh=0);

    /**
     * @copybrief   DataColl::nbPut(const void*,int,const void*,int,nbh_t)
     * @copydetails DataColl::nbPut(const void*,int,const void*,int,nbh_t)
     */
    virtual nbh_t nbPut(const void *idx, int idxlen, const void *buf, int buflen, nbh_t nbh=0);

    /**
     * @copybrief   DataColl::nbAdd(const void*,int,const void*,int,nbh_t)
     * @copydetails DataColl::nbAdd(const void*,int,const void*,int,nbh_t)
     */
    virtual void *nbAdd(const void *idx, int idxlen, const void *buf, int buflen, nbh_t nbh=0);

    /**
     * @copybrief   DataColl::waitHandle(nbh_t)
     * @copydetails DataColl::waitHandle(nbh_t)
     */
    virtual void waitHandle(nbh_t nbh);

    

  };
};

#endif /*__DenseArray_h__*/

