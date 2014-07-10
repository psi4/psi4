#ifndef __tascel_TslFuncReg_h__
#define __tascel_TslFuncReg_h__

#include <vector>

namespace tascel {
  class UniformTaskCollection;

  /**
   * Handle for a registered tascel function.
   */
  typedef int TslFunc;

  /**
   * Tascel function type.
   *
   * Takes a pointer to a task collection and a user-defined data structure.
   * The function is expected to cast the given data to an appropriate
   * type. Some of the arguments depend on the task collection that is
   * invoking the function and may be unitialized.   
   *
   * @param[in]     coll       Task collection
   * @param[in]     dscr       User-defined task descriptor
   * @param[in]     dscr_len   Length of the task descriptor
   * @param[in]     pldata     Pointer to process-local data
   * @param[in]     pldata_len Length of process-local data
   * @param[in,out] data_bufs  Vector of buffers for use by the task
   */
  typedef void (*TslFunc_t)(UniformTaskCollection *coll,
                            void *dscr, int dscr_len,
                            void *pldata, int pldata_len,
                            std::vector<void *> data_bufs);

  /**
   * Function Registration Table for tascel functions.
   *
   * A function registration table is necessary because the addresses of the
   * functions are specific to the address space of each process.  For this
   * table to work correctly, the table must be constructed collectively and
   * functions added to the table collectively.  Retrieving a function need
   * not be (and will most often not be) collective.
   *
   * The table is required to live as long as it is being used in a task
   * collection.
   */
  class TslFuncRegTbl {
    private:
      std::vector<TslFunc_t> ftbl; /**< the function lookup 'table' */
    public:
      /**
       * Constructs the function lookup table.
       */
      TslFuncRegTbl();

      /**
       * Copy constructor.
       */
      TslFuncRegTbl(const TslFuncRegTbl &that);

      /**
       * Destroys the function lookup table.
       */
      ~TslFuncRegTbl();

      /**
       * Registers a function within the table and returns a handle to the
       * function.
       *
       * This is a collective operation.
       *
       * @return handle to the registered function
       */
      TslFunc add(TslFunc_t f);

      /**
       * Retrieves the actual function within the address space of the caller.
       *
       * @param[in] fn handle to the function to retrieve
       * @return the tascel function
       */
      TslFunc_t get(TslFunc fn) const;

      /**
       * Assignment operator.
       */
      TslFuncRegTbl& operator = (const TslFuncRegTbl &that);

  }; /*TslFuncRegTbl*/

}; /*tascel*/

#endif /*__tascel_TslFuncReg_h__*/
