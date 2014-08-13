/* begin_generated_IBM_copyright_prolog                             */
/*                                                                  */
/* ---------------------------------------------------------------- */
/* (C)Copyright IBM Corp.  2007, 2008                               */
/* IBM BSD License                                                  */
/* ---------------------------------------------------------------- */
/*                                                                  */
/* end_generated_IBM_copyright_prolog                               */
/**
 * \file armci/src/x/armcix.h
 * \brief ARMCI Extension interface.
 */
#ifndef __armci_src_x_armcix_h
#define __armci_src_x_armcix_h

#include "armci.h"
#include "armcip.h"
#include "request.h"
#include "memlock.h"


/**
 * \brief Creates a compile error if the condition is false.
 *
 * This macro must be used within a function for the compiler to process it.
 * It is suggested that C++ classes and C files create an inline function
 * similar to the following example. The inline function is never used at
 * runtime and should be optimized out by the compiler. It exists for the sole
 * purpose of moving runtime \c assert calls to compile-time errors.
 *
 * \code
 * static inline void compile_time_assert ()
 * {
 *   // This compile time assert will succeed.
 *   COMPILE_TIME_ASSERT(sizeof(char) <= sizeof(double));
 *
 *   // This compile time assert will fail.
 *   COMPILE_TIME_ASSERT(sizeof(double) <= sizeof(char));
 * }
 * \endcode
 *
 * Compile time assert errors will look similar to the following:
 *
 * \code
 * foo.h: In function compile_time_assert:
 * foo.h:43: error: duplicate case value
 * foo.h:43: error: previously used here
 * \endcode
 */
#define COMPILE_TIME_ASSERT(expr) switch(0){case 0:case expr:;}

/**
 * \brief Assert during compile if certain conditions are not met.
 */
static inline void armcix_compile_time_assert ()
{
  /*
   * Assert that the size of the internal armci handle data structure is less
   * than or equal to the size of the public armci handle data structure.
   */
  COMPILE_TIME_ASSERT(sizeof(armci_ireq_t)<=sizeof(armci_hdl_t));
}

/**
 * \brief Initialize the ARMCI Extension.
 *
 * \todo Define return values.
 * \return ?
 */
int ARMCIX_Init ();

/**
 * \brief Initialize the ARMCI Extention lock resources.
 *
 * \param[in]  local_memlock_table memlock table
 */
void ARMCIX_init_memlock (memlock_t * local_memlock_table);

/**
 * \brief ARMCI Extension blocking memory lock operation.
 *
 * Send a lock request to the remote node and block until the lock has been
 * acquired on the remote node.
 *
 * \param[in] pstart The start virtual address of the range of memory to lock.
 * \param[in] pend   The end virtual address of the range of memory to lock.
 * \param[in] proc   Remote process(or) ID
 */
void ARMCIX_Lockmem (void * pstart, void * pend, int proc);

/**
 * \brief ARMCI Extension release memory lock operation.
 *
 * Send a lock release message to the remote node. This is a \e fire-and-forget
 * operation because the node does not block for an acknowledgement that the
 * lock release was successful.
 *
 * \param[in] proc   Remote process rank
 */
void ARMCIX_Unlockmem (int proc);

/**
 * \brief ARMCI Extension blocking wait operation for a specifc request
 *
 * The armcix_opaque_t structure is a field in the armci_ireq_t structure
 * and is used to maintain ARMCIX state information for an operation in
 * progress.
 *
 * \param[in] cmpl_info Pointer to the ARMCIX opaque object
 *
 * \todo define return values
 * \return ???
 *
 * \see armci_ireq_t
 */
int ARMCIX_Wait (armcix_opaque_t * cmpl_info);

/**
 * \brief ARMCI Extension blocking wait operation for all requests to a specific process
 *
 * All existing requests to the remote process are compelte after this function returns.
 *
 * \param[in] proc Remote process rank
 *
 * \todo define return values
 * \return ???
 */
int ARMCIX_WaitProc (int proc);

/**
 * \brief ARMCI Extension blocking wait operation for all requests to all processes
 *
 * All existing requests to all processes are completed after this function returns.
 *
 * \todo define return values
 * \return ???
 */
int ARMCIX_WaitAll ();


/**
 * \brief Point-to-point fence operation.
 *
 * Blocks until all active messages between the local node and the remote
 * node have completed and acknowledged by the remote node.
 *
 * \param[in] proc       Rank of the remote node to fence
 *
 * \see ARMCIX_AllFence
 */
void ARMCIX_Fence (int proc);

/**
 * \brief Global fence operation.
 *
 * Blocks until all active messages between the local node and all remote
 * nodes have completed and acknowledged by the remote node.
 *
 * \see ARMCIX_Fence
 */
void ARMCIX_AllFence ();

/**
 * \brief ARMCI Extension blocking read-modify-write operation.
 *
 * \param[in] op
 * \param[in] ploc
 * \param[in] prem
 * \param[in] extra
 * \param[in] proc
 *
 * \todo define return code and input parameters; add detailed doxygen description
 * \retval ???
 */
int ARMCIX_Rmw (int op, int * ploc, int * prem, int extra, int proc);

/**
 * \brief ARMCI Extension blocking get operation.
 *
 * \param[in] src       Source buffer on the remote node
 * \param[in] dst       Destination buffer on the local node
 * \param[in] bytes     Number of bytes to transfer
 * \param[in] proc      Remote node rank
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_Get (void * src, void * dst, int bytes, int proc);

/**
 * \brief ARMCI Extension blocking vector get operation.
 *
 * \param[in] darr      Descriptor array
 * \param[in] len       Length of descriptor array
 * \param[in] proc      Remote process(or) ID
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_GetV (armci_giov_t * darr, int len, int proc);

/**
 * \brief ARMCI Extension blocking strided get operation.
 *
 * \param[in] src_ptr        pointer to 1st segment at source
 * \param[in] src_stride_arr array of strides at source
 * \param[in] dst_ptr        pointer to 1st segment at destination
 * \param[in] dst_stride_arr array of strides at destination
 * \param[in] seg_count      number of segments at each stride levels: count[0]=bytes
 * \param[in] stride_levels  number of stride levels
 * \param[in] proc           remote process(or) ID
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_GetS (void * src_ptr, int * src_stride_arr, 
                 void * dst_ptr, int * dst_stride_arr, 
                 int * seg_count, int stride_levels, int proc);

/**
 * \brief ARMCI Extension non-blocking get operation.
 *
 * \param[in] src       Source buffer on the remote node
 * \param[in] dst       Destination buffer on the local node
 * \param[in] bytes     Number of bytes to transfer
 * \param[in] proc      Remote node rank
 * \param[in] nb_handle ARMCI non-blocking handle
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_NbGet (void * src, void * dst, int bytes, int proc, armci_ihdl_t nb_handle);

/**
 * \brief ARMCI Extension non-blocking vector get operation.
 *
 * \param[in] darr      Descriptor array
 * \param[in] len       Length of descriptor array
 * \param[in] proc      Remote process(or) ID
 * \param[in] nb_handle ARMCI non-blocking handle
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_NbGetV (armci_giov_t * darr, int len, int proc, armci_ihdl_t nb_handle);

/**
 * \brief ARMCI Extension non-blocking strided get operation.
 *
 * \param[in] src_ptr        pointer to 1st segment at source
 * \param[in] src_stride_arr array of strides at source
 * \param[in] dst_ptr        pointer to 1st segment at destination
 * \param[in] dst_stride_arr array of strides at destination
 * \param[in] seg_count      number of segments at each stride levels: count[0]=bytes
 * \param[in] stride_levels  number of stride levels
 * \param[in] proc           remote process(or) ID
 * \param[in] nb_handle      ARMCI non-blocking handle
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_NbGetS (void * src_ptr, int * src_stride_arr, 
                   void * dst_ptr, int * dst_stride_arr, 
                   int * seg_count, int stride_levels, int proc,
                   armci_ihdl_t nb_handle);

/**
 * \brief ARMCI Extension blocking put operation.
 *
 * \param[in] src       Source buffer on the local node
 * \param[in] dst       Destination buffer on the remote node
 * \param[in] bytes     Number of bytes to transfer
 * \param[in] proc      Remote node rank
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_Put (void * src, void * dst, int bytes, int proc);

/**
 * \brief ARMCI Extension blocking vector put operation.
 *
 * \param[in] darr      Descriptor array
 * \param[in] len       Length of descriptor array
 * \param[in] proc      Remote process(or) ID
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_PutV (armci_giov_t * darr, int len, int proc);

/**
 * \brief ARMCI Extension blocking strided put operation.
 *
 * \param[in] src_ptr        pointer to 1st segment at source
 * \param[in] src_stride_arr array of strides at source
 * \param[in] dst_ptr        pointer to 1st segment at destination
 * \param[in] dst_stride_arr array of strides at destination
 * \param[in] seg_count      number of segments at each stride levels: count[0]=bytes
 * \param[in] stride_levels  number of stride levels
 * \param[in] proc           remote process(or) ID
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_PutS (void * src_ptr, int * src_stride_arr, 
                 void * dst_ptr, int * dst_stride_arr, 
                 int * seg_count, int stride_levels, int proc);

/**
 * \brief ARMCI Extension non-blocking put operation.
 *
 * \param[in] src       Source buffer on the local node
 * \param[in] dst       Destination buffer on the remote node
 * \param[in] bytes     Number of bytes to transfer
 * \param[in] proc      Remote node rank
 * \param[in] nb_handle ARMCI non-blocking handle
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_NbPut (void * src, void * dst, int bytes, int proc, armci_ihdl_t nb_handle);

/**
 * \brief ARMCI Extension non-blocking vector put operation.
 *
 * \param[in] darr      Descriptor array
 * \param[in] len       Length of descriptor array
 * \param[in] proc      Remote process(or) ID
 * \param[in] nb_handle ARMCI non-blocking handle
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_NbPutV (armci_giov_t * darr, int len, int proc, armci_ihdl_t nb_handle);

/**
 * \brief ARMCI Extension non-blocking strided put operation.
 *
 * \param[in] src_ptr        pointer to 1st segment at source
 * \param[in] src_stride_arr array of strides at source
 * \param[in] dst_ptr        pointer to 1st segment at destination
 * \param[in] dst_stride_arr array of strides at destination
 * \param[in] seg_count      number of segments at each stride levels: count[0]=bytes
 * \param[in] stride_levels  number of stride levels
 * \param[in] proc           remote process(or) ID
 * \param[in] nb_handle      ARMCI non-blocking handle
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_NbPutS (void * src_ptr, int * src_stride_arr, 
                   void * dst_ptr, int * dst_stride_arr, 
                   int * seg_count, int stride_levels, int proc,
                   armci_ihdl_t nb_handle);

/**
 * \brief ARMCI Extension blocking accumulate operation.
 *
 * \param[in] datatype  accumulate datatype (operation code)
 * \param[in] scale     opaque pointer to the scaling factor for accumulate
 * \param[in] src       Source buffer on the local node
 * \param[in] dst       Destination buffer on the remote node
 * \param[in] bytes     Number of bytes to transfer
 * \param[in] proc      Remote node rank
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_Acc (int datatype, void * scale, void * src, void * dst, int bytes, int proc);

/**
 * \brief ARMCI Extension blocking vector accumulate operation.
 *
 * \param[in] datatype accumulate datatype (operation code)
 * \param[in] scale    opaque pointer to the scaling factor for accumulate
 * \param[in] darr     descriptor array
 * \param[in] len      length of the descriptor array
 * \param[in] proc     process(or) ID
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_AccV (int datatype, void * scale, armci_giov_t * darr, int len, int proc);

/**
 * \brief ARMCI Extension blocking strided accumulate operation.
 *
 * \param[in] datatype       accumulate datatype (operation code)
 * \param[in] scale          opaque pointer to the scaling factor for accumulate
 * \param[in] src_ptr        pointer to 1st segment at source
 * \param[in] src_stride_arr array of strides at source
 * \param[in] dst_ptr        pointer to 1st segment at destination
 * \param[in] dst_stride_arr array of strides at destination
 * \param[in] seg_count      number of segments at each stride levels: count[0]=bytes
 * \param[in] stride_levels  number of stride levels
 * \param[in] proc           remote process(or) ID
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_AccS (int datatype, void * scale,
                 void * src_ptr, int * src_stride_arr, 
                 void * dst_ptr, int * dst_stride_arr, 
                 int * seg_count, int stride_levels, int proc);

/**
 * \brief ARMCI Extension non-blocking accumulate operation.
 *
 * \param[in] datatype  accumulate datatype (operation code)
 * \param[in] scale     opaque pointer to the scaling factor for accumulate
 * \param[in] src       Source buffer on the local node
 * \param[in] dst       Destination buffer on the remote node
 * \param[in] bytes     Number of bytes to transfer
 * \param[in] proc      Remote node rank
 * \param[in] nb_handle ARMCI non-blocking handle
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_NbAcc (int datatype, void * scale, void * src, void * dst, int bytes, int proc, armci_ihdl_t nb_handle);

/**
 * \brief ARMCI Extension non-blocking vector accumulate operation.
 *
 * \param[in] datatype  accumulate datatype (operation code)
 * \param[in] scale     opaque pointer to the scaling factor for accumulate
 * \param[in] darr      Descriptor array
 * \param[in] len       Length of descriptor array
 * \param[in] proc      Remote process(or) ID
 * \param[in] nb_handle ARMCI non-blocking handle
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_NbAccV (int datatype, void * scale, armci_giov_t * darr, int len, int proc, armci_ihdl_t nb_handle);

/**
 * \brief ARMCI Extension non-blocking strided accumulate operation.
 *
 * \param[in] datatype       accumulate datatype (operation code)
 * \param[in] scale          opaque pointer to the scaling factor for accumulate
 * \param[in] src_ptr        pointer to 1st segment at source
 * \param[in] src_stride_arr array of strides at source
 * \param[in] dst_ptr        pointer to 1st segment at destination
 * \param[in] dst_stride_arr array of strides at destination
 * \param[in] seg_count      number of segments at each stride levels: count[0]=bytes
 * \param[in] stride_levels  number of stride levels
 * \param[in] proc           remote process(or) ID
 * \param[in] nb_handle      ARMCI non-blocking handle
 *
 * \todo define return code; add detailed doxygen description
 * \return ???
 */
int ARMCIX_NbAccS (int datatype, void * scale,
                   void * src_ptr, int * src_stride_arr, 
                   void * dst_ptr, int * dst_stride_arr, 
                   int * seg_count, int stride_levels, int proc,
                   armci_ihdl_t nb_handle);


/**
 * \page get_page Get APIs
 *
 * This is a description of the ARMCI Extension Get APIs
 *
 * \section get_blocking Blocking APIs
 * - ARMCIX_Get()
 * - ARMCIX_GetV()
 * - ARMCIX_GetS()
 * \section get_nonblocking Non-blocking APIs
 * - ARMCIX_NbGet()
 * - ARMCIX_NbGetV()
 * - ARMCIX_NbGetS()
 */

/**
 * \page put_page Put APIs
 *
 * This is a description of the ARMCI Extension Put APIs
 *
 * \section put_blocking Blocking APIs
 * - ARMCIX_Put()
 * - ARMCIX_PutV()
 * - ARMCIX_PutS()
 * \section put_nonblocking Non-blocking APIs
 * - ARMCIX_NbPut()
 * - ARMCIX_NbPutV()
 * - ARMCIX_NbPutS()
 */

/**
 * \page acc_page Accumulate APIs
 *
 * This is a description of the ARMCI Extension Accumulate APIs
 *
 * \section acc_blocking Blocking APIs
 * - ARMCIX_Acc()
 * - ARMCIX_AccV()
 * - ARMCIX_AccS()
 * \section acc_nonblocking Non-blocking APIs
 * - ARMCIX_NbAcc()
 * - ARMCIX_NbAccV()
 * - ARMCIX_NbAccS()
 */

/**
 * \page blocking_page Blocking APIs
 *
 * This is a description of the \b blocking ARMCI Extension APIs
 *
 * \section transfer Data transfer APIs
 * - ARMCIX_Get()
 * - ARMCIX_GetV()
 * - ARMCIX_GetS()
 * - ARMCIX_Put()
 * - ARMCIX_PutV()
 * - ARMCIX_PutS()
 * - ARMCIX_Acc()
 * - ARMCIX_AccV()
 * - ARMCIX_AccS()
 * - ARMCIX_Rmw()
 * \section sync Syncronization APIs
 * - ARMCIX_Fence()
 * - ARMCIX_AllFence()
 * - ARMCIX_Wait()
 * - ARMCIX_WaitProc()
 * - ARMCIX_WaitAll()
 * - ARMCIX_Lockmem()
 */

/**
 * \page nonblocking_page Non-blocking APIs
 *
 * This is a description of the \b non-blocking ARMCI Extension APIs
 *
 * \section transfer Data transfer APIs
 * - ARMCIX_NbGet()
 * - ARMCIX_NbGetV()
 * - ARMCIX_NbGetS()
 * - ARMCIX_NbPut()
 * - ARMCIX_NbPutV()
 * - ARMCIX_NbPutS()
 * - ARMCIX_NbAcc()
 * - ARMCIX_NbAccV()
 * - ARMCIX_NbAccS()
 */

/**
 * \page vector_page Vector APIs
 *
 * This is a description of the ARMCI Extension vector APIs
 *
 * \section vector_blocking Blocking APIs
 * - ARMCIX_GetV()
 * - ARMCIX_PutV()
 * - ARMCIX_AccV()
 * \section vector_nonblocking Non-blocking APIs
 * - ARMCIX_NbGetV()
 * - ARMCIX_NbPutV()
 * - ARMCIX_NbAccV()
 */

/**
 * \page strided_page Strided APIs
 *
 * This is a description of the ARMCI Extension strided APIs
 *
 * \section strided_blocking Blocking APIs
 * - ARMCIX_GetS()
 * - ARMCIX_PutS()
 * - ARMCIX_AccS()
 * \section strided_nonblocking Non-blocking APIs
 * - ARMCIX_NbGetS()
 * - ARMCIX_NbPutS()
 * - ARMCIX_NbAccS()
 */



#endif
