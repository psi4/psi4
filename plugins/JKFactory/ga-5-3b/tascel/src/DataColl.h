#ifndef __DataColl_h__
#define __DataColl_h__

namespace tascel {
  
  /**
   *  Abstract data collection.
   *
   * The tasking library deals with data collections that inherit this
   * base class (a.k.a. implement this interface). The data structure
   * is treated as a collection of blocks with opaque indices.
   */
  class DataColl {
    public:

    /**
     * Size of a data block.
     *
     * Return the size (in bytes) of data block referred to by the
     * index.
     * @param[in] idx    Pointer to index
     * @param[in] idxlen Length of index (in bytes)
     * @return           Size of data block (in bytes)
     */
    virtual int getSize(const void *idx, int idxlen) = 0;

    /**
     * Process holding the data block.
     *
     * This is a hint for locality optimizations. The return value does
     * not affect correctmess.
     * @param[in] idx    Pointer to index
     * @param[in] idxlen Length of index (in bytes)
     * @return           Process holding the data block
     */
    virtual int getProc(const void *idx, int idxlen) = 0;

    /**
     * Get the data block.
     *
     * Gets the data block pointed to by the index. The buffer is
     * expected to be of length sufficient to fit the block.
     * @param[in] idx    Pointer to index
     * @param[in] idxlen Length of index (in bytes)
     * @param[out] buf   Pointer to place the data
     * @param[in] buflen Length of buffer
     * @return           none
     */
    virtual void get(const void *idx, int idxlen, void *buf, int buflen) = 0;

    /**
     * Put the data block.
     *
     * Puts the data block in buffer at location pointed to by the
     * index. The buffer is expected to be of length sufficient to fit
     * the block. 
     * @param[in] idx    Pointer to index
     * @param[in] idxlen Length of index (in bytes)
     * @param[in] buf    Pointer to place the data
     * @param[in] buflen Length of buffer
     * @return           none
     */
    virtual void put(const void *idx, int idxlen, const void *buf, int buflen) = 0;

    /**
     * Add the data block.
     *
     * Adds the data block in buffer at location pointed to by the
     * index. Adds are atomic w.r.t. other add operations. The buffer
     * is expected to be of length sufficient to fit the block. 
     * @param[in] idx    Pointer to index
     * @param[in] idxlen Length of index (in bytes)
     * @param[in] buf    Pointer to place the data
     * @param[in] buflen Length of buffer
     * @return           none
     */
    virtual void add(const void *idx, int idxlen, const void *buf, int buflen) = 0;
    
    /**
     * Opaque non-blocking wait handle. 
     */
    typedef void *nbh_t;

    /**
     * Get the data block.
     *
     * Gets the data block pointed to by the index. The buffer is
     * expected to be of length sufficient to fit the block.
     * @param[in] idx    Pointer to index
     * @param[in] idxlen Length of index (in bytes)
     * @param[out] buf   Pointer to place the data
     * @param[in] buflen Length of buffer
     * @param[inout] nbh A specific handle to be used (if possible).
     * @return           Handle to initiated communication. 
     */
    virtual nbh_t nbGet(const void *idx, int idxlen, void *buf, int buflen, nbh_t nbh=0) = 0;

    /**
     * Put the data block.
     *
     * Puts the data block in buffer at location pointed to by the
     * index. The buffer is expected to be of length sufficient to fit
     * the block. 
     * @param[in] idx    Pointer to index
     * @param[in] idxlen Length of index (in bytes)
     * @param[in] buf    Pointer to place the data
     * @param[in] buflen Length of buffer
     * @param[inout] nbh A specific handle to be used (if possible).
     * @return           Handle to initiated communication. 
     */
    virtual nbh_t nbPut(const void *idx, int idxlen, const void *buf, int buflen, nbh_t nbh=0) = 0;

    /**
     * Add the data block.
     *
     * Adds the data block in buffer at location pointed to by the
     * index. Adds are atomic w.r.t. other add operations. The buffer
     * is expected to be of length sufficient to fit the block. 
     * @param[in] idx    Pointer to index
     * @param[in] idxlen Length of index (in bytes)
     * @param[in] buf    Pointer to place the data
     * @param[in] buflen Length of buffer
     * @param[inout] nbh A specific handle to be used (if possible).
     * @return           Handle to initiated communication. 
     */
    virtual void *nbAdd(const void *idx, int idxlen, const void *buf, int buflen, nbh_t nbh=0) = 0;

    /**
     * Wait for communication.
     *
     * Blocks till communication corresponding to the handle is
     * complete. The use of the handle after the call's return results
     * in undefined behavior.
     * @param[in] nbh Non-blocking handle
     * @return        None
     */
    virtual void waitHandle(nbh_t nbh)=0;

  }; /* DataColl*/
}; /**tascel*/

#endif /*__DataColl_h__*/

