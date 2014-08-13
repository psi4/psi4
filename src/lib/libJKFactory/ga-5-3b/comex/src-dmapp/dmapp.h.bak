/*
 * Copyright (c) 2008 Cray Inc. All Rights Reserved.
 *
 * The contents of this file is proprietary information of Cray Inc.
 * and may not be disclosed without prior written consent.
 *
 * $HeadURL: http://svn.us.cray.com/svn/baker/packages/dmapp/branches/RB-3.1/include/dmapp.h $
 * $LastChangedRevision: 2964 $
 *
 *
 * Header file for API for Distributed Memory Applications on Gemini
 *
 */


#ifndef _DMAPP_H
#define _DMAPP_H

#include <stddef.h>
#include <stdint.h>
#include "gni_pub.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* Specifies input _and_ output arguments to DMAPP functions */
#define INOUT

/* include version info */
#include "dmapp_rev.h"

/*
 * ******** Defines/Macros ****************************************
 */

/* Default threshold, in bytes. Transfers smaller than the threshold
   are transfered using cpu-based mechanism, transfers of this size
   or larger are transfered using off-load mechanism. */

#define DMAPP_OFFLOAD_THRESHOLD              (4096)

/* Maximum and minimun number of outstanding non-blocking requests 
   supported; this includes explicit _and_ implicit non-blocking 
   requests. */

#define DMAPP_MAX_OUTSTANDING_NB             (2048)
#define DMAPP_DEF_OUTSTANDING_NB              (512)
#define DMAPP_MIN_OUTSTANDING_NB                (5) /* coordinated w/ GNI */


/*
 * ******** Typedefs, Structs *************************************
 */

/*
 * Function return codes
 */

typedef enum dmapp_return {
        DMAPP_RC_SUCCESS = 0,
        DMAPP_RC_INVALID_PARAM,
        DMAPP_RC_ALIGNMENT_ERROR,
        DMAPP_RC_TRANSACTION_ERROR,
        DMAPP_RC_RESOURCE_ERROR,
        DMAPP_RC_PERMISSION_ERROR,
        DMAPP_RC_NO_SPACE, /* transaction could not be completed due to
                              insufficient resources; user should synchronize 
                              more often or increase max_outstanding_nb */
        DMAPP_RC_NOT_DONE, /* used in _dmapp_nbi_syncid->status */
        DMAPP_RC_NOT_SUPPORTED,
        DMAPP_RC_NOT_FOUND,
        DMAPP_RC_BUSY,
        DMAPP_RC_NOT_USED
} dmapp_return_t;

/*  
 * These are the valid types which can be supplied via the type input 
 * parameter to all data motion funtions. 
 */

typedef enum dmapp_type {
        DMAPP_DQW = 0, /* double quad word (16 byte) */
        DMAPP_QW,      /* quad word (8 byte) */
        DMAPP_DW,      /* double word (4 byte) */
        DMAPP_BYTE     /* byte. Don't use this if you want performance ! */
} dmapp_type_t;

/* Type of routing to be performed for request packet.
   These are valid options for the relaxed_ordering field in
   the dmapp_rma_attrs_t structure below. */

typedef enum dmapp_routing_type {
        DMAPP_ROUTING_IN_ORDER = 0,   /* hash off, adapt off; poor performance */
        DMAPP_ROUTING_DETERMINISTIC,  /* hash on,  adapt off */
        DMAPP_ROUTING_ADAPTIVE        /* hash off, adapt on  */
} dmapp_routing_type_t;

/* Mode of PI access ordering to be used during GNI memory registration
 * For STRICT ordering both posted and non-posted writes arrive in strict order.
 * For DEFAULT and RELAXED these ordering constraints are not honoured and hence
 * extra synchronisation is required when global visibility of data is required
 * (e.g. during a gsync/fence or after a blocking put)
 * These modes do not affect Get operations
 */
typedef enum dmapp_pi_reg_type {
        DMAPP_PI_ORDERING_STRICT = 0, /* Strict PI ordering on HT interface (P_PASS_PW=0, NP_PASS_PW=0) */
        DMAPP_PI_ORDERING_DEFAULT,    /* Default GNI PI ordering (P_PASS_PW=0, NP_PASS_PW=1) */
        DMAPP_PI_ORDERING_RELAXED     /* Relaxed PI ordering (P_PASS_PW=1, NP_PASS_PW=1) */
} dmapp_pi_reg_type_t;


/*  
 * This is a memory segment descriptor, with an address and length.
 * The "len" field contains the currently mapped size of the segment.
 */

typedef struct dmapp_seg_desc {
        void             *addr;
        size_t           len;        /* mapped size of segment in bytes */
        gni_mem_handle_t memhndl;
        uint16_t         flags;      /* DMAPP internal use only */
        void             *reserved;  /* DMAPP internal use only */
} dmapp_seg_desc_t;

/*
 * DMAPP Processing element, aka. PE
 */

typedef int dmapp_pe_t; 

/*  
 * This is the application and memory layout information for a DMAPP job.
 */

typedef struct dmapp_jobinfo {
        int              version;          /* DMAPP version */
        int              hw_version;       /* GNI hardware version */
        int              npes;             /* Number of PEs in entire job */
        dmapp_pe_t       pe;               /* PE number, in [0, npes-1] */
        dmapp_seg_desc_t data_seg;         /* Data segment memory */
        dmapp_seg_desc_t sheap_seg;        /* Symmetric heap memory */
} dmapp_jobinfo_t;

/*  
 * RMA attributes can be used to control the way in which DMAPP handles
 * various RMA requests.
 * Some attributes can be set during initialization only, others can be
 * set multiple times over the course of a job.
 */

typedef struct dmapp_rma_attrs {
        /* max nr of outstanding non-blocking explicit requests supported;
           can be specified during initialization only; 
           [DMAPP_MIN_OUTSTANDING_NB, .., DMAPP_MAX_OUTSTANDING_NB] is legal */
        uint32_t max_outstanding_nb;
        /* treshhold in bytes for switch from cpu-based mechanisms to 
           cpu offload mechanisms; can be specified at any time;
           any value is legal */
        uint32_t offload_threshold;
        /* flag to indicate whether relaxed ordering of requests is allowed
           and if so, which specific routing option to use. Valid options
           are specified by dmapp_routing_type_t above. The routing mode
           can be controlled separately for PUTs (including all AMOs) and
           GETs. The default for PUTs is DMAPP_ROUTING_DETERMINISTIC, the
           default for GETs is DMAPP_ROUTING_ADAPTIVE. DMAPP_ROUTING_IN_ORDER
           guarantees the requests arrive in order and is expected to result
           in poor performance. */
        uint8_t put_relaxed_ordering;
        uint8_t get_relaxed_ordering;
        /* max nr of threads accessing DMAPP; only used with thread-safety 
           enabled! default is 1; specified during init only; must be >= 1 */
        uint8_t max_concurrency;
} dmapp_rma_attrs_t; 

/*  
 * EXTENDED
 * RMA attributes can be used to control the way in which DMAPP handles
 * various RMA requests.
 * Some attributes can be set during initialization only, others can be
 * set multiple times over the course of a job.
 */

typedef struct dmapp_rma_attrs_ext {
        /* max nr of outstanding non-blocking explicit requests supported;
           can be specified during initialization only; 
           [DMAPP_MIN_OUTSTANDING_NB, .., DMAPP_MAX_OUTSTANDING_NB] is legal */
        uint32_t max_outstanding_nb;
        /* treshhold in bytes for switch from cpu-based mechanisms to 
           cpu offload mechanisms; can be specified at any time;
           any value is legal */
        uint32_t offload_threshold;
        /* flag to indicate whether relaxed ordering of requests is allowed
           and if so, which specific routing option to use. Valid options
           are specified by dmapp_routing_type_t above. The routing mode
           can be controlled separately for PUTs (including all AMOs) and
           GETs. The default for PUTs is DMAPP_ROUTING_DETERMINISTIC, the
           default for GETs is DMAPP_ROUTING_ADAPTIVE. DMAPP_ROUTING_IN_ORDER
           guarantees the requests arrive in order and is expected to result
           in poor performance. */
        uint8_t put_relaxed_ordering;
        uint8_t get_relaxed_ordering;
        /* max nr of threads accessing DMAPP; only used with thread-safety 
           enabled! default is 1; specified during init only; must be >= 1 */
        uint8_t max_concurrency;

        /* defines the PI ordering registration flags used by DMAPP when
         * registering all memory regions with GNI. Hence applies to the DATA,
         * SYMMETRIC HEAP and any user or dynamically mapped regions.
         *
         * Should be set to a value defined in dmapp_pi_reg_type_t the default
         * is DMAPP_PI_ORDERING_STRICT
         */
        uint8_t PI_ordering;

        uint8_t unused[32];	/* reserved for future expansion */

} dmapp_rma_attrs_ext_t;

/*
 * Synchronization ID
 */

typedef struct dmapp_syncid *dmapp_syncid_handle_t;

/*
 * DMAPP structures for collective extension support
 */

#define DMAPP_C_PSET_MODE_CONCAT                  0x00000001
/* set for specifying use of collective algorithms which
   are tolerant to network faults, may be slower than
   non-fault tolerant algorithms */
#define DMAPP_C_PSET_MODE_FAULT_TOLERANT          0x00000002
/* try to use any available hardware offload capabilities
   for handling collective operations */
#define DMAPP_C_PSET_MODE_HW_OFFLOAD              0x00000004

typedef enum dmapp_c_pset_delimiter_type {
           DMAPP_C_PSET_DELIMITER_VEC = 1,
           DMAPP_C_PSET_DELIMITER_STRIDED,
           DMAPP_C_PSET_DELIMITER_LAST
} dmapp_c_pset_delimiter_type_t;

typedef struct {
        uint32_t n_pes;
        uint32_t *vec_pes;
} dmapp_c_pset_delimiter_vec_t;

typedef struct {
        uint32_t n_pes;
        uint32_t base_pe;
        uint32_t stride_pe;
} dmapp_c_pset_delimiter_strided_t;

typedef struct {
        void *concat_buf;
        uint64_t concat_buf_size;
        dmapp_c_pset_delimiter_type_t type;
        union {
                dmapp_c_pset_delimiter_vec_t vec_type;
                dmapp_c_pset_delimiter_strided_t stride_type;
       } u;
} dmapp_c_pset_desc_t;

typedef struct dmapp_c_pset *dmapp_c_pset_handle_t;

typedef enum dmapp_c_pset_attrs_type {
           DMAPP_C_PSET_ATTRS_VER1 = 1,
           DMAPP_C_PSET_ATTRS_LAST
} dmapp_c_pset_attrs_type_t;

typedef struct {
        uint32_t radix;
} dmapp_c_pset_attrs_ver1_t; 

typedef struct {
        dmapp_c_pset_attrs_type_t type; 
        union {
                dmapp_c_pset_attrs_ver1_t type1;
       } u;
} dmapp_c_pset_attrs_t;

typedef enum dmapp_c_type {
        DMAPP_C_INT32 = 101,
        DMAPP_C_UINT32,      
        DMAPP_C_INT64,
        DMAPP_C_UINT64,
        DMAPP_C_FLOAT,
        DMAPP_C_DOUBLE, 
        DMAPP_C_COMPLEX8, 
        DMAPP_C_COMPLEX16, 
        DMAPP_C_FLOAT_UINT64,
        DMAPP_C_DOUBLE_UINT64,
        DMAPP_C_INT32_UINT64,
        DMAPP_C_INT64_UINT64,     
        DMAPP_C_TYPE_LAST
} dmapp_c_type_t;

typedef enum dmapp_c_op {
        DMAPP_C_SUM = 201,
        DMAPP_C_MAX,
        DMAPP_C_MIN,      
        DMAPP_C_PROD,
        DMAPP_C_LAND,
        DMAPP_C_BAND,
        DMAPP_C_LOR,      
        DMAPP_C_BOR,
        DMAPP_C_LXOR,
        DMAPP_C_BXOR,
        DMAPP_C_MINLOC,
        DMAPP_C_MAXLOC,
        DMAPP_C_OP_LAST
} dmapp_c_op_t;

/*
 * special data types for minloc/maxloc and complex
 * global reduction operations
 */

typedef struct dmapp_c_float_uint64 {
        float value;
        uint64_t loc;
} dmapp_c_float_uint64_t;

typedef struct dmapp_c_double_uint64 {
        double value;
        uint64_t loc;
} dmapp_c_double_uint64_t;

typedef struct dmapp_c_int_uint64 {
        int value;
        uint64_t loc;
} dmapp_c_int_uint64_t;

typedef struct dmapp_c_int64_uint64 {
        int64_t value;
        uint64_t loc;
} dmapp_c_int64_uint64_t;

typedef struct dmapp_c_complex8 {
        float re;
        float im;
} dmapp_c_complex8_t;

typedef struct dmapp_c_complex16 {
        double re;
        double im;
} dmapp_c_complex16_t;

/*
 * ******** Function Declarations ********************************
 */

/* 
 * Initialization and Query Functions 
 */

/* dmapp_init - Initialize resources for a DMAPP job
 *              Must be called by all DMAPP applications!
 *
 * Parameters:
 * IN  requested_attrs     Desired job attributes
 * OUT actual_attrs        Actual job attributes
 *
 * Returns: 
 * DMAPP_RC_SUCCESS        Operation completed successfully
 * DMAPP_RC_INVALID_PARAM  One or more argument is invalid
 * DMAPP_RC_RESOURCE_ERROR An error occurred during initialization
 *
 */

extern dmapp_return_t
dmapp_init(IN  dmapp_rma_attrs_t *requested_attrs,
           OUT dmapp_rma_attrs_t *actual_attrs);


/* dmapp_init_ext - Initialize resources for a DMAPP job
 *                  Must be called by all DMAPP applications!
 *
 * Version that takes dmapp_rma_attrs_ext_t
 * Either dmapp_init() or dmapp_init_ext() should be called
 *
 */
extern dmapp_return_t
dmapp_init_ext(IN  dmapp_rma_attrs_ext_t *requested_attrs,
               OUT dmapp_rma_attrs_ext_t *actual_attrs);

/* dmapp_finalize - Synchronizes and cleans up DMAPP resources
 *                  Must be called by all DMAPP applications!
 *
 * Returns:
 * DMAPP_RC_SUCCESS
 */

extern dmapp_return_t
dmapp_finalize(void);


/* dmapp_get_jobinfo - Retrieve general job information
 *
 * Parameters:
 * OUT info   Current information about the job
 *
 * Returns: 
 * DMAPP_RC_SUCCESS    Operation completed successfully
 * DMAPP_RC_INVALID_PARAM Input parameter is invalid
 */

extern dmapp_return_t 
dmapp_get_jobinfo(OUT dmapp_jobinfo_t *info);


/* dmapp_get_rma_attrs - Retrieve RMA attributes of a DMAPP job
 *
 * Parameters:
 * OUT attrs   Current RMA attributes of the job
 *
 * Returns: 
 * DMAPP_RC_SUCCESS    Operation completed successfully
 * DMAPP_RC_INVALID_PARAM Input parameter is invalid
 */

extern dmapp_return_t 
dmapp_get_rma_attrs(OUT dmapp_rma_attrs_t *attrs);

/* dmapp_get_rma_attrs_ext - Retrieve RMA attributes of a DMAPP job
 *
 * Version that returns the dmapp_rma_attrs_ext_t extended structure
 *
 */
extern dmapp_return_t
dmapp_get_rma_attrs_ext(OUT dmapp_rma_attrs_ext_t *attrs);


/* dmapp_set_rma_attrs - Set dynamic RMA attributes for a DMAPP job
 *
 * Parameters:
 * IN  requested_attrs Desired job attributes
 * OUT actual_attrs    Actual job attributes
 *
 * Returns: 
 * DMAPP_RC_SUCCESS    Operation completed successfully
 * DMAPP_RC_INVALID_PARAM One or more input parameters is invalid
 *
 * Description:
 * A process can set RMA attributes to control the way that DMAPP handles
 * various RMA requests. Some attributes can be set only during initialization.
 * They will be referred to as static attributes. Others can be set multiple 
 * times over the course of the job, and will be referred to as dynamic attributes. 
 * Settings dynamic attributes does not affect RMA requests previously issued by
 * the PE, only subsequent RMA requests. Dynamic attributes include when to 
 * switch from CPU-based mechanisms for handling RMA requests to using CPU
 * offload mechanisms.
 */

extern dmapp_return_t 
dmapp_set_rma_attrs(IN  dmapp_rma_attrs_t *requested_attrs,
                    OUT dmapp_rma_attrs_t *actual_attrs);

/* dmapp_set_rma_attrs - Set dynamic RMA attributes for a DMAPP job
 *
 * Version that accepts and returns the extended dmapp_rma_attrs_ext_t
 * structure
 *
 */
extern dmapp_return_t
dmapp_set_rma_attrs_ext(IN  dmapp_rma_attrs_ext_t *requested_attrs,
                        OUT dmapp_rma_attrs_ext_t *actual_attrs);


/* 
 * One-sided RMA Functions:
 *
 * All one-sided RMA functions (PUT type, GET type and atomic memory
 * operations) can be subdivided into three categories:
 * - blocking (no suffix):
 *   The process returns from the function only after the side-effects
 *   of the remote memory access are globally visible in the system.
 * - non-blocking explicit (suffix _nb):
 *   A synchronization ID (syncid) is returned to the process,
 *   the side-effects of the remote memory access are only assured
 *   to be globally visible in the system after the application has
 *   determined via an explicit synchronization call
 *   (dmapp_syncid_test/dmapp_syncid_wait) that the syncid has been
 *   retired.
 * - non-blocking implicit (suffix _nbi):
 *   No explicit synchronization ID is returned to the process, the
 *   side-effects of the remote memory access are only assured to be
 *   globally visible in the system following a call to
 *   dmapp_gsync_test or dmapp_gsync_wait.
 *   This mode is recommended for performance reasons for applications
 *   with lots of small messages, where blocking calls or using
 *   individual synids would be expensive.
 */



/* The PUT functions stores a contiguous block of data starting at
 * address source_addr in local memory into a contiguous block at 
 * a remote address where this remote address is specified by the
 * triplet virtual address target_addr, segment descriptor target_seg
 * and target process target_pe. nelems specifies the number of 
 * elements of type type to be transferred. The memory region defined 
 * by target_addr and nelems must be within an exported memory segement
 * of target_pe.
 */


/* dmapp_put_nb - Non-blocking explicit PUT
 *
 * Parameters:
 * IN  target_addr     Address of target buffer
 *     target_seg      Segment descriptor of target buffer
 *     target_pe       Target PE
 *     source_addr     Address of source buffer
 *     nelems          Number of elements to be transferred
 *     type            Type of elements to be transferred
 * OUT syncid          Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS        Operation completed successfully
 * DMAPP_RC_INVALID_PARAM  One or more input parameters is invalid
 * DMAPP_RC_RESOURCE_ERROR A resource error occurred
 * DMAPP_RC_NO_SPACE       The transaction request could not be completed
 *                         due to insufficient resources; user should
 *                         increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_put_nb(IN  void                  *target_addr,
             IN  dmapp_seg_desc_t      *target_seg,
             IN  dmapp_pe_t            target_pe,
             IN  void                  *source_addr,
             IN  uint64_t              nelems,
             IN  dmapp_type_t          type,
             OUT dmapp_syncid_handle_t *syncid);


/* dmapp_put_nbi - Non-blocking implicit PUT
 *
 * Parameters:
 * IN  target_addr     Address of target buffer
 *     target_seg      Segment descriptor of target buffer
 *     target_pe       Target PE
 *     source_addr     Address of source buffer
 *     nelems          Number of elements to be transferred
 *     type            Type of elements to be transferred
 *
 * Returns: 
 * DMAPP_RC_SUCCESS        Operation completed successfully
 * DMAPP_RC_INVALID_PARAM  One or more input parameters is invalid
 * DMAPP_RC_RESOURCE_ERROR A resource error occurred
 * DMAPP_RC_NO_SPACE       The transaction request could not be completed
 *                         due to insufficient resources; user should
 *                         increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_put_nbi(IN  void             *target_addr,
              IN  dmapp_seg_desc_t *target_seg,
              IN  dmapp_pe_t       target_pe,
              IN  void             *source_addr,
              IN  uint64_t         nelems,
              IN  dmapp_type_t     type);


/* dmapp_put - Blocking PUT
 *
 * Parameters:
 * IN  target_addr     Address of target buffer
 *     target_seg      Segment descriptor of target buffer
 *     target_pe       Target PE
 *     source_addr     Address of source buffer
 *     nelems          Number of elements to be transferred
 *     type            Type of elements to be transferred
 *
 * Returns: 
 * DMAPP_RC_SUCCESS        Operation completed successfully
 * DMAPP_RC_INVALID_PARAM  One or more input parameters is invalid
 * DMAPP_RC_RESOURCE_ERROR A resource error occurred
 * DMAPP_RC_NO_SPACE       The transaction request could not be completed
 *                         due to insufficient resources; user should
 *                         increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_put(IN  void             *target_addr,
          IN  dmapp_seg_desc_t *target_seg,
          IN  dmapp_pe_t       target_pe,
          IN  void             *source_addr,
          IN  uint64_t         nelems,
          IN  dmapp_type_t     type);


/* The GET functions load from a contiguous block of data starting 
 * from a remote source address and returning the data into a 
 * contiguous block starting at address target_addr in local memory. 
 * The remote address is specified by the triplet virtual address 
 * source_adddr, segment descriptor source_seg and source process 
 * source_pe. nelems specifies the number of elements of type type 
 * to be transferred. The memory region described by the remote
 * address and nelems must reside in an exported memory of source_pe.
 * Note that zero-length GETs are not supported. 
 */


/* dmapp_get_nb - Non-blocking explicit GET
 * 
 * Parameters:
 * IN  target_addr    Address of target buffer
 *     source_addr    Address of source buffer
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     nelems         Number of elements to be transferred
 *     type           Type of elements to be transferred
 * OUT syncid         Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS        Operation completed successfully
 * DMAPP_RC_INVALID_PARAM  One or more input parameters is invalid
 * DMAPP_RC_RESOURCE_ERROR A resource error occurred
 * DMAPP_RC_NO_SPACE       The transaction request could not be completed
 *                         due to insufficient resources; user should
 *                         increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_get_nb(IN  void                  *target_addr,
             IN  void                  *source_addr,
             IN  dmapp_seg_desc_t      *source_seg,
             IN  dmapp_pe_t            source_pe,
             IN  uint64_t              nelems,
             IN  dmapp_type_t          type,
             OUT dmapp_syncid_handle_t *syncid);


/* dmapp_get_nbi - Non-blocking implicit GET
 *
 * Parameters:
 * IN  target_addr    Address of target buffer
 *     source_addr    Address of source buffer
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE 
 *     nelems         Number of elements to be transferred
 *     type           Type of elements to be transferred
 *
 * Returns: 
 * DMAPP_RC_SUCCESS        Operation completed successfully
 * DMAPP_RC_INVALID_PARAM  One or more input parameters is invalid
 * DMAPP_RC_RESOURCE_ERROR A resource error occurred
 * DMAPP_RC_NO_SPACE       The transaction request could not be completed
 *                         due to insufficient resources; user should
 *                         increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_get_nbi(IN  void             *target_addr,
              IN  void             *source_addr,
              IN  dmapp_seg_desc_t *source_seg,
              IN  dmapp_pe_t       source_pe,
              IN  uint64_t         nelems,
              IN  dmapp_type_t     type);


/* dmapp_get - Blocking GET
 *
 * Parameters:
 * IN  target_addr    Address of target buffer
 *     source_addr    Address of source buffer
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     nelems         Number of elements to be transferred
 *     type           Type of elements to be transferred
 *
 * Returns: 
 * DMAPP_RC_SUCCESS        Operation completed successfully
 * DMAPP_RC_INVALID_PARAM  One or more input parameters is invalid
 * DMAPP_RC_RESOURCE_ERROR A resource error occurred
 * DMAPP_RC_NO_SPACE       The transaction request could not be completed
 *                         due to insufficient resources; user should
 *                         increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_get(IN  void             *target_addr,
          IN  void             *source_addr,
          IN  dmapp_seg_desc_t *source_seg,
          IN  dmapp_pe_t       source_pe,
          IN  uint64_t         nelems,
          IN  dmapp_type_t     type);



/* The Strided PUT functions deliver data starting at address 
 * source_addr in local memory to a remote address using a stride
 * specified by sst at the source and by tst at the target. The
 * remote address is specified by the triplet virtual address 
 * target_addr, segment descriptor target_seg and target process 
 * target_pe. nelems specifies the number of elements of type type 
 * to be transferred. The memory region defined by target_addr, tst
 * and nelems must be within an exported memory segement of target_pe.
 */


/* dmapp_iput_nb - Non-blocking explicit Strided PUT
 *
 * Parameters:
 * IN  target_addr     Address of target buffer
 *     target_seg      Segment descriptor of target buffer
 *     target_pe       Target PE
 *     source_addr     Address of source buffer
 *     tst             Target stride (>= 1)
 *     sst             Source stride (>= 1)
 *     nelems          Number of elements to be transferred
 *     type            Type of elements to be transferred
 * OUT syncid          Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS       Operation completed successfully
 * DMAPP_RC_INVALID_PARAM One or more input parameter is invalid
 * DMAPP_RC_NO_SPACE      The transaction request could not be completed
 *                        due to insufficient resources; user should
 *                        increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_iput_nb(IN  void                  *target_addr,
              IN  dmapp_seg_desc_t      *target_seg,
              IN  dmapp_pe_t            target_pe,
              IN  void                  *source_addr,
              IN  ptrdiff_t             tst,
              IN  ptrdiff_t             sst,
              IN  uint64_t              nelems,
              IN  dmapp_type_t          type,
              OUT dmapp_syncid_handle_t *syncid);


/* dmapp_iput_nbi - Non-blocking implicit Strided PUT
 *
 * Parameters:
 * IN  target_addr     Address of target buffer
 *     target_seg      Segment descriptor of target buffer
 *     target_pe       Target PE
 *     source_addr     Address of source buffer
 *     tst             Target stride
 *     sst             Source stride
 *     nelems          Number of elements to be transferred
 *     type            Type of elements to be transferred
 *
 * Returns: 
 * DMAPP_RC_SUCCESS       Operation completed successfully
 * DMAPP_RC_INVALID_PARAM One or more input parameters is invalid
 * DMAPP_RC_NO_SPACE      The transaction request could not be completed
 *                        due to insufficient resources; user should
 *                        increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_iput_nbi(IN  void             *target_addr,
               IN  dmapp_seg_desc_t *target_seg,
               IN  dmapp_pe_t       target_pe,
               IN  void             *source_addr,
               IN  ptrdiff_t        tst,
               IN  ptrdiff_t        sst,
               IN  uint64_t         nelems,
               IN  dmapp_type_t     type);


/* dmapp_iput - Blocking Strided PUT
 *
 * Parameters:
 * IN  target_addr     Address of target buffer
 *     target_seg      Segment descriptor of target buffer
 *     target_pe       Target PE
 *     source_addr     Address of source buffer
 *     tst             Target stride
 *     sst             Source stride
 *     nelems          Number of elements to be transferred
 *     type            Type of elements to be transferred
 *
 * Returns: 
 * DMAPP_RC_SUCCESS       Operation completed successfully
 * DMAPP_RC_INVALID_PARAM One or more input parameters is invalid
 * DMAPP_RC_NO_SPACE      The transaction request could not be completed
 *                        due to insufficient resources; user should
 *                        increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_iput(IN  void             *target_addr,
           IN  dmapp_seg_desc_t *target_seg,
           IN  dmapp_pe_t       target_pe,
           IN  void             *source_addr,
           IN  ptrdiff_t        tst,
           IN  ptrdiff_t        sst,
           IN  uint64_t         nelems,
           IN  dmapp_type_t     type);


/* The Strided GET functions load data starting from a remote source 
 * address using the source side stride sst and returning the data to
 * address target_addr in local memory using the target side stride tst.
 * The remote address is specified by the triplet virtual address
 * source_adddr, segment descriptor source_seg and source process
 * source_pe. nelems specifies the number of elements of type type
 * to be transferred. The memory region described by the remote
 * address, sst and nelems must reside in an exported memory of source_pe.
 */

/* dmapp_iget_nb - Non-blocking explicit Strided GET
 *
 * Parameters:
 * IN  target_addr    Address of target buffer
 *     source_addr    Address of source buffer
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     tst            Target stride
 *     sst            Source stride
 *     nelems         Number of elements to be transferred
 *     type           Type of elements to be transferred
 * OUT syncid         Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Source or target buffer or length not
 *                          properly (Dword (4 byte)) aligned
 * DMAPP_RC_RESOURCE_ERROR  A resource error occured
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_iget_nb(IN  void                  *target_addr,
              IN  void                  *source_addr,
              IN  dmapp_seg_desc_t      *source_seg,
              IN  dmapp_pe_t            source_pe,
              IN  ptrdiff_t             tst,
              IN  ptrdiff_t             sst,
              IN  uint64_t              nelems,
              IN  dmapp_type_t          type,
              OUT dmapp_syncid_handle_t *syncid);

/* dmapp_iget_nbi - Non-blocking implicit Strided GET
 *
 * Parameters:
 * IN  target_addr    Address of target buffer
 *     source_addr    Address of source buffer
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     tst            Target stride
 *     sst            Source stride
 *     nelems         Number of elements to be transferred
 *     type           Type of elements to be transferred
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Source or target buffer or length not
 *                          properly (Dword (4 byte)) aligned
 * DMAPP_RC_RESOURCE_ERROR  A resource error occured
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_iget_nbi(IN  void             *target_addr,
               IN  void             *source_addr,
               IN  dmapp_seg_desc_t *source_seg,
               IN  dmapp_pe_t       source_pe,
               IN  ptrdiff_t        tst,
               IN  ptrdiff_t        sst,
               IN  uint64_t         nelems,
               IN  dmapp_type_t     type);

/* dmapp_iget - Blocking Strided GET
 *
 * Parameters:
 * IN  target_addr    Address of target buffer
 *     source_addr    Address of source buffer
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     tst            Target stride (>= 1)
 *     sst            Source stride (>= 1)
 *     nelems         Number of elements to be transferred
 *     type           Type of elements to be transferred
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Source or target buffer or length not
 *                          properly (Dword (4 byte)) aligned
 * DMAPP_RC_RESOURCE_ERROR  A resource error occured
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_iget(IN  void             *target_addr,
           IN  void             *source_addr,
           IN  dmapp_seg_desc_t *source_seg,
           IN  dmapp_pe_t       source_pe,
           IN  ptrdiff_t        tst,
           IN  ptrdiff_t        sst,
           IN  uint64_t         nelems,
           IN  dmapp_type_t     type);


/* The Indexed PUT functions scatter a contiguous block 
 * of data starting at address source_addr in local memory to a remote 
 * address using offsets specified by the tidx array. The remote 
 * address is specified by the triplet virtual address target_addr, 
 * segment descriptor target_seg and target process target_pe. nelems 
 * specifies the number of elements of type type to be transferred. 
 * Offsets into the tidx array are in units of type. The memory region 
 * defined by target_addr, tidx and nelems must be within an exported 
 * memory segement of target_pe.
 */

/* dmapp_ixput_nb - Non-blocking explicit Indexed PUT
 *
 * Parameters:
 * IN  target_addr     Address of target buffer
 *     target_seg      Segment descriptor of target buffer
 *     target_pe       Target PE
 *     source_addr     Address of source buffer
 *     tidx            Array of positive offsets into target buffer
 *     nelems          Number of elements to be transferred
 *     type            Type of elements to be transferred
 * OUT syncid          Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS       Operation completed successfully
 * DMAPP_RC_INVALID_PARAM One or more input parameters is invalid
 * DMAPP_RC_NO_SPACE      The transaction request could not be completed
 *                        due to insufficient resources; user should
 *                        increase max_outstanding_nb or sync more often
 *
 */

extern dmapp_return_t 
dmapp_ixput_nb(IN  void                  *target_addr,
               IN  dmapp_seg_desc_t      *target_seg,
               IN  dmapp_pe_t            target_pe,
               IN  void                  *source_addr,
               IN  ptrdiff_t             *tidx,
               IN  uint64_t              nelems,
               IN  dmapp_type_t          type,
               OUT dmapp_syncid_handle_t *syncid);

/* dmapp_ixput_nbi - Non-blocking implicit Indexed PUT
 *
 * Parameters:
 * IN  target_addr     Address of target buffer
 *     target_seg      Segment descriptor of target buffer
 *     target_pe       Target PE
 *     source_addr     Address of source buffer
 *     tidx            Array of positive offsets into target buffer
 *     nelems          Number of elements to be transferred
 *     type            Type of elements to be transferred
 *
 * Returns: 
 * DMAPP_RC_SUCCESS       Operation completed successfully
 * DMAPP_RC_INVALID_PARAM One or more input parameters is invalid
 * DMAPP_RC_NO_SPACE      The transaction request could not be completed
 *                        due to insufficient resources; user should
 *                        increase max_outstanding_nb or sync more often
 *
 */

extern dmapp_return_t 
dmapp_ixput_nbi(IN  void             *target_addr,
                IN  dmapp_seg_desc_t *target_seg,
                IN  dmapp_pe_t       target_pe,
                IN  void             *source_addr,
                IN  ptrdiff_t        *tidx,
                IN  uint64_t         nelems,
                IN  dmapp_type_t     type);

/* dmapp_ixput - Blocking Indexed PUT
 *
 * Parameters:
 * IN  target_addr     Address of target buffer
 *     target_seg      Segment descriptor of target buffer
 *     target_pe       Target PE
 *     source_addr     Address of source buffer
 *     tidx            Array of positive offsets into target buffer
 *     nelems          Number of elements to be transferred
 *     type            Type of elements to be transferred
 *
 * Returns: 
 * DMAPP_RC_SUCCESS       Operation completed successfully
 * DMAPP_RC_INVALID_PARAM One or more input parameters is invalid
 * DMAPP_RC_NO_SPACE      The transaction request could not be completed
 *                        due to insufficient resources; user should
 *                        increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_ixput(IN  void             *target_addr,
            IN  dmapp_seg_desc_t *target_seg,
            IN  dmapp_pe_t       target_pe,
            IN  void             *source_addr,
            IN  ptrdiff_t        *tidx,
            IN  uint64_t         nelems,
            IN  dmapp_type_t     type);

 
/* The Indexed GET functions gather data starting
 * from a remote source address using offsets specified by the
 * sidx array and returning the data into a contiguous block 
 * starting at address target_addr in local memory. The remote 
 * address is specified by the triplet virtual address source_adddr, 
 * segment descriptor source_seg and source process source_pe. nelems 
 * specifies the number of elements of type type to be transferred. 
 * Offsets in the sidx array are in units of type. The memory region 
 * described by the remote address, sidx and nelems must reside in 
 * an exported memory of source_pe.
 * Note that the three Indexed Get functions are not supported for
 * type DMAPP_BYTE.
 */

/* dmapp_ixget_nb - Non-blocking explicit Indexed GET
 *
 * Parameters:
 * IN  target_addr    Address of target buffer
 *     source_addr    Address of source buffer
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     sidx           Array of positive offsets into source buffer
 *     nelems         Number of elements to be transferred
 *     type           Type of elements to be transferred,
                      DMAPP_BYTE not supported
 * OUT syncid         Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_NOT_IMPLEMENTED Type DMAPP_BYTE was used
 * DMAPP_RC_ALIGNMENT_ERROR Source or target buffer or length not
 *                          properly (Dword (4 byte)) aligned
 * DMAPP_RC_RESOURCE_ERROR  A resource error occured
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_ixget_nb(IN  void                  *target_addr,
               IN  void                  *source_addr,
               IN  dmapp_seg_desc_t      *source_seg,
               IN  dmapp_pe_t            source_pe,
               IN  ptrdiff_t             *sidx,
               IN  uint64_t              nelems,
               IN  dmapp_type_t          type,
               OUT dmapp_syncid_handle_t *syncid);

/* dmapp_ixget_nbi - Non-blocking implicit Indexed GET
 *
 * Parameters:
 * IN  target_addr    Address of target buffer
 *     source_addr    Address of source buffer
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     sidx           Array of positive offsets into source buffer
 *     nelems         Number of elements to be transferred
 *     type           Type of elements to be transferred,
                      DMAPP_BYTE not supported
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_NOT_IMPLEMENTED Type DMAPP_BYTE was used
 * DMAPP_RC_ALIGNMENT_ERROR Source or target buffer or length not
 *                          properly (Dword (4 byte)) aligned
 * DMAPP_RC_RESOURCE_ERROR  A resource error occured
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_ixget_nbi(IN  void             *target_addr,
                IN  void             *source_addr,
                IN  dmapp_seg_desc_t *source_seg,
                IN  dmapp_pe_t       source_pe,
                IN  ptrdiff_t        *sidx,
                IN  uint64_t         nelems,
                IN  dmapp_type_t     type);

/* dmapp_ixget - Blocking Indexed GET
 *
 * Parameters:
 * IN  target_addr    Address of target buffer
 *     source_addr    Address of source buffer
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     sidx           Array of positive offsets into source buffer
 *     nelems         Number of elements to be transferred
 *     type           Type of elements to be transferred,
                      DMAPP_BYTE not supported
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_NOT_IMPLEMENTED Type DMAPP_BYTE was used
 * DMAPP_RC_ALIGNMENT_ERROR Source or target buffer or length not
 *                          properly (Dword (4 byte)) aligned
 * DMAPP_RC_RESOURCE_ERROR  A resource error occured
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_ixget(IN  void             *target_addr,
            IN  void             *source_addr,
            IN  dmapp_seg_desc_t *source_seg,
            IN  dmapp_pe_t       source_pe,
            IN  ptrdiff_t        *sidx,
            IN  uint64_t         nelems,
            IN  dmapp_type_t     type);


/* The PUT with Indexed PE Stride functions deliver data starting 
 * at source_addr in local memory to a list of target PEs target_pe_list 
 * starting at a remote address in their memories. The remote address 
 * is specified by the target virtual address target_addr and the
 * segment descriptor target_seg. nelems specifies the number of
 * elements of type type to be transferred to each target PE. 
 * When the transfer is complete, each target PE will have a copy
 * of the contents of the original source buffer. The address
 * range specified by target_addr and nelems must be within an 
 * exported memory segement of each PE in target_pe_list.
 * Note that the remote address must be a symmetric address.
 * Also note that this function is not a collective operation. It 
 * is best used when a small amount of data needs to be transfered
 * to a set of PEs.
 */

/* dmapp_put_ixpe_nb - Non-blocking explicit PUT with Indexed PE Stride
 *
 * Parameters:
 * IN  target_addr     Address of target buffer
 *     target_seg      Segment descriptor of target buffer
 *     target_pe_list  List of target PEs
 *     num_target_pes  Number of target PEs
 *     source_addr     Address of source buffer
 *     nelems          Number of elements to be transferred
 *     type            Type of elements to be transferred
 * OUT syncid          Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS       Operation completed successfully
 * DMAPP_RC_INVALID_PARAM One or more input parameters is invalid
 * DMAPP_RC_NO_SPACE      The transaction request could not be completed
 *                        due to insufficient resources; user should
 *                        increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_put_ixpe_nb(IN  void                  *target_addr,
                  IN  dmapp_seg_desc_t      *target_seg,
                  IN  dmapp_pe_t            *target_pe_list,
                  IN  uint32_t              num_target_pes,
                  IN  void                  *source_addr,
                  IN  uint64_t              nelems,
                  IN  dmapp_type_t          type,
                  OUT dmapp_syncid_handle_t *syncid);

/* dmapp_put_ixpe_nbi - Non-blocking implicit PUT with Indexed PE Stride
 *
 * Parameters:
 * IN  target_addr     Address of target buffer
 *     target_seg      Segment descriptor of target buffer
 *     target_pe_list  List of target PEs
 *     num_target_pes  Number of target PEs
 *     source_addr     Address of source buffer
 *     nelems          Number of elements to be transferred
 *     type            Type of elements to be transferred
 *
 * Returns: 
 * DMAPP_RC_SUCCESS       Operation completed successfully
 * DMAPP_RC_INVALID_PARAM One or more input parameters is invalid
 * DMAPP_RC_NO_SPACE      The transaction request could not be completed
 *                        due to insufficient resources; user should
 *                        increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_put_ixpe_nbi(IN  void             *target_addr,
                   IN  dmapp_seg_desc_t *target_seg,
                   IN  dmapp_pe_t       *target_pe_list,
                   IN  uint32_t         num_target_pes,
                   IN  void             *source_addr,
                   IN  uint64_t         nelems,
                   IN  dmapp_type_t     type);

/* dmapp_put_ixpe - Blocking PUT with Indexed PE Stride
 *
 * Parameters:
 * IN  target_addr     Address of target buffer
 *     target_seg      Segment descriptor of target buffer
 *     target_pe_list  List of target PEs
 *     num_target_pes  Number of target PEs
 *     source_addr     Address of source buffer
 *     nelems          Number of elements to be transferred
 *     type            Type of elements to be transferred
 *
 * Returns: 
 * DMAPP_RC_SUCCESS       Operation completed successfully
 * DMAPP_RC_INVALID_PARAM One or more input parameters is invalid
 * DMAPP_RC_NO_SPACE      The transaction request could not be completed
 *                        due to insufficient resources; user should
 *                        increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_put_ixpe(IN  void             *target_addr,
               IN  dmapp_seg_desc_t *target_seg,
               IN  dmapp_pe_t       *target_pe_list,
               IN  uint32_t         num_target_pes,
               IN  void             *source_addr,
               IN  uint64_t         nelems,
               IN  dmapp_type_t     type);


/* The Scatter with Indexed PE Stride functions deliver data starting 
 * at source_addr in local memory to a list of target PEs target_pe_list
 * starting at a remote address in their memories. The remote address is
 * specified by the target virtual address target_addr and the
 * segment descriptor target_seg. nelems specifies the number of
 * elements of type type to be transferred to each target PE.
 * Unlike the dmapp_put_ixpe* functions, the source array 
 * specifies an array of size num_target_pes * nelems * sizeof(type).
 * The target PE at index i into the target_pe_list will receive
 * elements i * nelems to (i+1) * nelems - 1. The address range 
 * specified by target_addr and nelems must be within an exported 
 * memory segement of each PE in target_pe_list.
 * Note that the remote address must be a symmetric address.
 * Also note that this function is not a collective operation. It
 * is best used when a small amount of data needs to be transfered
 * to a set of PEs.
 */

/* dmapp_scatter_ixpe_nb - Non-blocking explicit Scatter with Indexed PE Stride
 *
 * Parameters:
 * IN  target_addr     Address of target buffer
 *     target_seg      Segment descriptor of target buffer
 *     target_pe_list  List of target PEs
 *     num_target_pes  Number of target PEs
 *     source_addr     Address of source buffer
 *     nelems          Number of elements to be transferred
 *     type            Type of elements to be transferred
 * OUT syncid          Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS    Operation completed successfully
 * DMAPP_RC_INVALID_PARAM One or more input parameters is invalid
 * DMAPP_RC_NO_SPACE      The transaction request could not be completed
 *                        due to insufficient resources; user should
 *                        increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_scatter_ixpe_nb(IN  void                   *target_addr,
                      IN  dmapp_seg_desc_t       *target_seg,
                      IN  dmapp_pe_t             *target_pe_list,
                      IN  uint32_t               num_target_pes,
                      IN  void                   *source_addr,
                      IN  uint64_t               nelems,
                      IN  dmapp_type_t           type,
                      OUT dmapp_syncid_handle_t *syncid);

/* dmapp_scatter_ixpe_nbi - Non-blocking implicit Scatter with Indexed PE Stride
 *
 * Parameters:
 * IN  target_addr     Address of target buffer
 *     target_seg      Segment descriptor of target buffer
 *     target_pe_list  List of target PEs
 *     num_target_pes  Number of target PEs
 *     source_addr     Address of source buffer
 *     nelems          Number of elements to be transferred
 *     type            Type of elements to be transferred
 *
 * Returns: 
 * DMAPP_RC_SUCCESS    Operation completed successfully
 * DMAPP_RC_INVALID_PARAM One or more input parameters is invalid
 * DMAPP_RC_NO_SPACE      The transaction request could not be completed
 *                        due to insufficient resources; user should
 *                        increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_scatter_ixpe_nbi(IN  void             *target_addr,
                       IN  dmapp_seg_desc_t *target_seg,
                       IN  dmapp_pe_t       *target_pe_list,
                       IN  uint32_t         num_target_pes,
                       IN  void             *source_addr,
                       IN  uint64_t         nelems,
                       IN  dmapp_type_t     type);

/* dmapp_scatter_ixpe - Blocking Scatter with Indexed PE Stride
 *
 * Parameters:
 * IN  target_addr     Address of target buffer
 *     target_seg      Segment descriptor of target buffer
 *     target_pe_list  List of target PEs
 *     num_target_pes  Number of target PEs
 *     source_addr     Address of source buffer
 *     nelems          Number of elements to be transferred
 *     type            Type of elements to be transferred
 *
 * Returns: 
 * DMAPP_RC_SUCCESS       Operation completed successfully
 * DMAPP_RC_INVALID_PARAM One or more input parameters is invalid
 * DMAPP_RC_NO_SPACE      The transaction request could not be completed
 *                        due to insufficient resources; user should
 *                        increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_scatter_ixpe(IN  void             *target_addr,
                   IN  dmapp_seg_desc_t *target_seg,
                   IN  dmapp_pe_t       *target_pe_list,
                   IN  uint32_t         num_target_pes,
                   IN  void             *source_addr,
                   IN  uint64_t         nelems,
                   IN  dmapp_type_t     type);


/* The Gather with Indexed PE Stride functions gather data starting 
 * at a remote source address from a list of PEs specified by 
 * source_pe_list and concatenates the data in a buffer specified by 
 * target_addr in local memory.
 * The remote address is specified by the virtual address source_addr
 * and segment descriptor source_seg. nelems specifies the number of 
 * elements of type type to be collected from each PE. The address
 * range described by the remote address and nelems must reside in
 * an exported memory of each PE in source_pe_list.
 * Note that the remote address must be a symmetric address.
 * Also note that this function is not a collective operation. It
 * is best used when a small amount of data needs to be transfered
 * to a set of PEs.
 */

/* dmapp_gather_ixpe_nb - Non-blocking explicit Gather with Indexed PE Stride
 *
 * Parameters:
 * IN  target_addr    Address of target buffer
 *     source_addr    Address of source buffer
 *     source_seg     Segment descriptor of source buffer
 *     source_pe_list List of source PEs
 *     num_source_pes Number of source PEs
 *     nelems         Number of elements to be transferred
 *     type           Type of elements to be transferred
 * OUT syncid         Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS       Operation completed successfully
 * DMAPP_RC_INVALID_PARAM One or more input parameters is invalid
 * DMAPP_RC_NO_SPACE      The transaction request could not be completed
 *                        due to insufficient resources; user should
 *                        increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_gather_ixpe_nb(IN  void                  *target_addr,
                     IN  void                  *source_addr,
                     IN  dmapp_seg_desc_t      *source_seg,
                     IN  dmapp_pe_t            *source_pe_list,
                     IN  uint32_t              num_source_pes,
                     IN  uint64_t              nelems,
                     IN  dmapp_type_t          type,
                     OUT dmapp_syncid_handle_t *syncid);

/* dmapp_gather_ixpe_nbi - Non-blocking implicit Gather with Indexed PE Stride
 *
 * Parameters:
 * IN  target_addr    Address of target buffer
 *     source_addr    Address of source buffer
 *     source_seg     Segment descriptor of source buffer
 *     source_pe_list List of source PEs
 *     num_source_pes Number of source PEs
 *     nelems         Number of elements to be transferred
 *     type           Type of elements to be transferred
 *
 * Returns: 
 * DMAPP_RC_SUCCESS       Operation completed successfully
 * DMAPP_RC_INVALID_PARAM One or more input parameters is invalid
 * DMAPP_RC_NO_SPACE      The transaction request could not be completed
 *                        due to insufficient resources; user should
 *                        increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_gather_ixpe_nbi(IN  void             *target_addr,
                      IN  void             *source_addr,
                      IN  dmapp_seg_desc_t *source_seg,
                      IN  dmapp_pe_t       *source_pe_list,
                      IN  uint32_t         num_source_pes,
                      IN  uint64_t         nelems,
                      IN  dmapp_type_t     type);

/* dmapp_gather_ixpe - Blocking Gather with Indexed PE Stride
 *
 * Parameters:
 * IN  target_addr    Address of target buffer
 *     source_addr    Address of source buffer
 *     source_seg     Segment descriptor of source buffer
 *     source_pe_list List of source PEs
 *     num_source_pes Number of source PEs
 *     nelems         Number of elements to be transferred
 *     type           Type of elements to be transferred
 *
 * Returns: 
 * DMAPP_RC_SUCCESS       Operation completed successfully
 * DMAPP_RC_INVALID_PARAM One or more input parameters is invalid
 * DMAPP_RC_NO_SPACE      The transaction request could not be completed
 *                        due to insufficient resources; user should
 *                        increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_gather_ixpe(IN  void             *target_addr,
                  IN  void             *source_addr,
                  IN  dmapp_seg_desc_t *source_seg,
                  IN  dmapp_pe_t       *source_pe_list,
                  IN  uint32_t         num_source_pes,
                  IN  uint64_t         nelems,
                  IN  dmapp_type_t     type);



/* 
 *  Scalar-Style AMO RMA Functions 
 */

/* A set of scalar-type AMO functions is provided. Note that
 * AMOs only operate on quad-word entities. For atomic functions 
 * with PUT semantics (AADD, AAND, AOR, AXOR) the target memory
 * location must reside in an exported memory segment of remote 
 * PE target_pe. For atomic functions with GET semantics (AFADD,
 * AFAND, AFOR, AFXOR) the source memory location must reside in 
 * an exported memory segment of remote PE source_pe.
 */

/* dmapp_aadd_qw_nb - Non-blocking explicit atomic ADD
 *
 * Parameters:
 * IN  target_addr    Address of target buffer (qw only)
 *     target_seg     Segment descriptor for target buffer
 *     target_pe      Target PE
 *     operand        Value to be added
 * OUT syncid         Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target buffer not properly (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_aadd_qw_nb(IN  void                  *target_addr,
                 IN  dmapp_seg_desc_t      *target_seg,
                 IN  dmapp_pe_t            target_pe,
                 IN  int64_t               operand,
                 OUT dmapp_syncid_handle_t *syncid);

/* dmapp_aadd_qw_nbi - Non-blocking implicit atomic ADD
 *
 * Parameters:
 * IN  target_addr    Address of target buffer (qw only)
 *     target_seg     Segment descriptor for target buffer
 *     target_pe      Target PE
 *     operand        Value to be added
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target buffer not properly (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_aadd_qw_nbi(IN  void             *target_addr,
                  IN  dmapp_seg_desc_t *target_seg,
                  IN  dmapp_pe_t       target_pe,
                  IN  int64_t          operand);

/* dmapp_aadd_qw - Blocking atomic ADD
 *
 * Parameters:
 * IN  target_addr    Address of target buffer (qw only)
 *     target_seg     Segment descriptor for target buffer
 *     target_pe      Target PE
 *     operand        Value to be added
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target buffer not properly (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_aadd_qw(IN  void             *target_addr,
              IN  dmapp_seg_desc_t *target_seg,
              IN  dmapp_pe_t       target_pe,
              IN  int64_t          operand);


/* dmapp_aand_qw_nb - Non-blocking explicit atomic AND
 *
 * Parameters:
 * IN  target_addr    Address of target buffer (qw only)
 *     target_seg     Segment descriptor for target buffer
 *     target_pe      Target PE
 *     operand        Operand for the AND
 * OUT syncid         Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target buffer not properly (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_aand_qw_nb(IN  void                  *target_addr,
                 IN  dmapp_seg_desc_t      *target_seg,
                 IN  dmapp_pe_t            target_pe,
                 IN  int64_t               operand,
                 OUT dmapp_syncid_handle_t *syncid);

/* dmapp_aand_qw_nbi - Non-blocking implicit atomic AND
 *
 * Parameters:
 * IN  target_addr    Address of target buffer (qw only)
 *     target_seg     Segment descriptor for target buffer
 *     target_pe      Target PE
 *     operand        Operand for the AND
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target buffer not properly (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_aand_qw_nbi(IN  void             *target_addr,
                  IN  dmapp_seg_desc_t *target_seg,
                  IN  dmapp_pe_t       target_pe,
                  IN  int64_t          operand);

/* dmapp_aand_qw - Blocking atomic AND
 *
 * Parameters:
 * IN  target_addr    Address of target buffer (qw only)
 *     target_seg     Segment descriptor for target buffer
 *     target_pe      Target PE
 *     operand        Operand for the AND
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target buffer not properly (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_aand_qw(IN  void             *target_addr,
              IN  dmapp_seg_desc_t *target_seg,
              IN  dmapp_pe_t       target_pe,
              IN  int64_t          operand);


/* dmapp_aor_qw_nb - Non-blocking explicit atomic OR
 *
 * Parameters:
 * IN  target_addr    Address of target buffer (qw only)
 *     target_seg     Segment descriptor for target buffer
 *     target_pe      Target PE
 *     operand        Operand for the OR
 * OUT syncid         Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target buffer not properly (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_aor_qw_nb(IN  void                  *target_addr,
                IN  dmapp_seg_desc_t      *target_seg,
                IN  dmapp_pe_t            target_pe,
                IN  int64_t               operand,
                OUT dmapp_syncid_handle_t *syncid);

/* dmapp_aor_qw_nbi - Non-blocking implicit atomic OR
 *
 * Parameters:
 * IN  target_addr    Address of target buffer (qw only)
 *     target_seg     Segment descriptor for target buffer
 *     target_pe      Target PE
 *     operand        Operand for the OR
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target buffer not properly (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_aor_qw_nbi(IN  void             *target_addr,
                 IN  dmapp_seg_desc_t *target_seg,
                 IN  dmapp_pe_t       target_pe,
                 IN  int64_t          operand);

/* dmapp_aor_qw - Blocking atomic OR
 *
 * Parameters:
 * IN  target_addr    Address of target buffer (qw only)
 *     target_seg     Segment descriptor for target buffer
 *     target_pe      Target PE
 *     operand        Operand for the ORD
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target buffer not properly (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_aor_qw(IN  void             *target_addr,
             IN  dmapp_seg_desc_t *target_seg,
             IN  dmapp_pe_t       target_pe,
             IN  int64_t          operand);


/* dmapp_axor_qw_nb - Non-blocking explicit atomic XOR
 *
 * Parameters:
 * IN  target_addr    Address of target buffer (qw only)
 *     target_seg     Segment descriptor for target buffer
 *     target_pe      Target PE
 *     operand        Operand for the XOR
 * OUT syncid         Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target buffer not properly (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_axor_qw_nb(IN  void                  *target_addr,
                 IN  dmapp_seg_desc_t      *target_seg,
                 IN  dmapp_pe_t            target_pe,
                 IN  int64_t               operand,
                 OUT dmapp_syncid_handle_t *syncid);

/* dmapp_axor_qw_nbi - Non-blocking implicit atomic XOR
 *
 * Parameters:
 * IN  target_addr    Address of target buffer (qw only)
 *     target_seg     Segment descriptor for target buffer
 *     target_pe      Target PE
 *     operand        Operand for the XOR
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target buffer not properly (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_axor_qw_nbi(IN  void             *target_addr,
                  IN  dmapp_seg_desc_t *target_seg,
                  IN  dmapp_pe_t       target_pe,
                  IN  int64_t          operand);

/* dmapp_axor_qw - Blocking atomic XOR
 *
 * Parameters:
 * IN  target_addr    Address of target buffer (qw only)
 *     target_seg     Segment descriptor for target buffer
 *     target_pe      Target PE
 *     operand        Operand for the XOR
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target buffer not properly (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_axor_qw(IN  void             *target_addr,
              IN  dmapp_seg_desc_t *target_seg,
              IN  dmapp_pe_t       target_pe,
              IN  int64_t          operand);


/* dmapp_afadd_qw_nb - Non-blocking explicit atomic FADD
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     operand        Operand for the FADD
 * OUT syncid         Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly 
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_afadd_qw_nb(IN  void                  *target_addr,
                  IN  void                  *source_addr,
                  IN  dmapp_seg_desc_t      *source_seg,
                  IN  dmapp_pe_t            source_pe,
                  IN  int64_t               operand,
                  OUT dmapp_syncid_handle_t *syncid);

/* dmapp_afadd_qw_nbi - Non-blocking implicit atomic FADD
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result 
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     operand        Operand for the FADD
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_afadd_qw_nbi(IN  void             *target_addr,
                   IN  void             *source_addr,
                   IN  dmapp_seg_desc_t *source_seg,
                   IN  dmapp_pe_t       source_pe,
                   IN  int64_t          operand);

/* dmapp_afadd_qw - Blocking atomic FADD
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result 
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     operand        Operand for the FADD
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_afadd_qw(IN  void             *target_addr,
               IN  void             *source_addr,
               IN  dmapp_seg_desc_t *source_seg,
               IN  dmapp_pe_t       source_pe,
               IN  int64_t          operand);


/* dmapp_afand_qw_nb - Non-blocking explicit atomic FAND
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result 
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     operand        Operand for the FAND
 * OUT syncid         Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_afand_qw_nb(IN  void                  *target_addr,
                  IN  void                  *source_addr,
                  IN  dmapp_seg_desc_t      *source_seg,
                  IN  dmapp_pe_t            source_pe,
                  IN  int64_t               operand,
                  OUT dmapp_syncid_handle_t *syncid);

/* dmapp_afand_qw_nbi - Non-blocking implicit atomic FAND
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result 
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     operand        Operand for the FAND
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_afand_qw_nbi(IN  void             *target_addr,
                   IN  void             *source_addr,
                   IN  dmapp_seg_desc_t *source_seg,
                   IN  dmapp_pe_t       source_pe,
                   IN  int64_t          operand);

/* dmapp_afand_qw - Blocking atomic FAND
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result 
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     operand        Operand for the FAND
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_afand_qw(IN  void             *target_addr,
               IN  void             *source_addr,
               IN  dmapp_seg_desc_t *source_seg,
               IN  dmapp_pe_t       source_pe,
               IN  int64_t          operand);


/* dmapp_afxor_qw_nb - Non-blocking explicit atomic FXOR
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result 
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     operand        Operand for the FXOR
 * OUT syncid         Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_afxor_qw_nb(IN  void                  *target_addr,
                  IN  void                  *source_addr,
                  IN  dmapp_seg_desc_t      *source_seg,
                  IN  dmapp_pe_t            source_pe,
                  IN  int64_t               operand,
                  OUT dmapp_syncid_handle_t *syncid);

/* dmapp_afxor_qw_nbi - Non-blocking implicit atomic FXOR
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result 
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     operand        Operand for the FXOR
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_afxor_qw_nbi(IN  void             *target_addr,
                   IN  void             *source_addr,
                   IN  dmapp_seg_desc_t *source_seg,
                   IN  dmapp_pe_t       source_pe,
                   IN  int64_t          operand);

/* dmapp_afxor_qw - Blocking atomic FXOR
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result 
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     operand        Operand for the FXOR
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_afxor_qw(IN  void             *target_addr,
               IN  void             *source_addr,
               IN  dmapp_seg_desc_t *source_seg,
               IN  dmapp_pe_t       source_pe,
               IN  int64_t          operand);


/* dmapp_afor_qw_nb - Non-blocking explicit atomic FOR
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result 
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     operand        Operand for the FOR
 * OUT syncid         Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_afor_qw_nb(IN  void                  *target_addr,
                 IN  void                  *source_addr,
                 IN  dmapp_seg_desc_t      *source_seg,
                 IN  dmapp_pe_t            source_pe,
                 IN  int64_t               operand,
                 OUT dmapp_syncid_handle_t *syncid);

/* dmapp_afor_qw_nbi - Non-blocking implicit atomic FOR
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result 
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     operand        Operand for the FOR
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_afor_qw_nbi(IN  void             *target_addr,
                  IN  void             *source_addr,
                  IN  dmapp_seg_desc_t *source_seg,
                  IN  dmapp_pe_t       source_pe,
                  IN  int64_t          operand);

/* dmapp_afor_qw - Blocking atomic FOR
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result 
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     operand        Operand for the FOR
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_afor_qw(IN  void             *target_addr,
              IN  void             *source_addr,
              IN  dmapp_seg_desc_t *source_seg,
              IN  dmapp_pe_t       source_pe,
              IN  int64_t          operand);



/* 
 *  Two-Operand AMO RMA Functions 
 */

/* A set of two-operand AMO functions is provided. Note that
 * AMOs only operate on quad-word entities.
 * For two-operand atomic functions with GET semantics (AFAX,
 * ACSWAP) the source memory location must reside in an exported 
 * memory segment of remote PE source_pe.
 */

/* dmapp_afax_qw_nb - Non-blocking explicit atomic FAX
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result 
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     andMask        mask for AND operation
 *     xorMask        mask for XOR operation
 * OUT syncid         Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly 
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_afax_qw_nb(IN  void                  *target_addr,
                 IN  void                  *source_addr,
                 IN  dmapp_seg_desc_t      *source_seg,
                 IN  dmapp_pe_t            source_pe,
                 IN  int64_t               andMask,
                 IN  int64_t               xorMask,
                 OUT dmapp_syncid_handle_t *syncid);

/* dmapp_afax_qw_nbi - Non-blocking implicit atomic FAX
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result 
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     andMask        mask for AND operation
 *     xorMask        mask for XOR operation
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_afax_qw_nbi(IN  void             *target_addr,
                  IN  void             *source_addr,
                  IN  dmapp_seg_desc_t *source_seg,
                  IN  dmapp_pe_t       source_pe,
                  IN  int64_t          andMask,
                  IN  int64_t          xorMask);

/* dmapp_afax_qw - Blocking atomic FAX
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result 
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     andMask        mask for AND operation
 *     xorMask        mask for XOR operation
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_afax_qw(IN  void             *target_addr,
              IN  void             *source_addr,
              IN  dmapp_seg_desc_t *source_seg,
              IN  dmapp_pe_t       source_pe,
              IN  int64_t          andMask,
              IN  int64_t          xorMask);


/* dmapp_acswap_qw_nb - Non-blocking explicit atomic CSWAP
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result 
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     comperand      Operand against which to compare
 *     swaperand      Operand which to swap in
 * OUT syncid         Synchronization ID
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_acswap_qw_nb(IN  void                  *target_addr,
                   IN  void                  *source_addr,
                   IN  dmapp_seg_desc_t      *source_seg,
                   IN  dmapp_pe_t            source_pe,
                   IN  int64_t               comperand,
                   IN  int64_t               swaperand,
                   OUT dmapp_syncid_handle_t *syncid);

/* dmapp_acswap_qw_nbi - Non-blocking implicit atomic CSWAP
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result 
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     comperand      Operand against which to compare
 *     swaperand      Operand which to swap in
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 */

extern dmapp_return_t 
dmapp_acswap_qw_nbi(IN  void             *target_addr,
                    IN  void             *source_addr,
                    IN  dmapp_seg_desc_t *source_seg,
                    IN  dmapp_pe_t       source_pe,
                    IN  int64_t          comperand,
                    IN  int64_t          swaperand);

/* dmapp_acswap_qw - Blocking atomic CSWAP
 *
 * Parameters:
 * IN  target_addr    Address of target buffer where result 
 *                    is returned (qw only)
 *     source_addr    Address of source buffer (qw only)
 *     source_seg     Segment descriptor of source buffer
 *     source_pe      Source PE
 *     comperand      Operand against which to compare
 *     swaperand      Operand which to swap in
 *
 * Returns: 
 * DMAPP_RC_SUCCESS         Operation completed successfully
 * DMAPP_RC_INVALID_PARAM   One or more input parameters is invalid
 * DMAPP_RC_ALIGNMENT_ERROR Target or source buffer not properly
 *                          (Qword (8 byte)) aligned
 * DMAPP_RC_NO_SPACE        The transaction request could not be completed
 *                          due to insufficient resources; user should
 *                          increase max_outstanding_nb or sync more often
 * DMAPP_RC_TRANSACTION_ERROR A transaction error has occured
 */

extern dmapp_return_t 
dmapp_acswap_qw(IN  void             *target_addr,
                IN  void             *source_addr,
                IN  dmapp_seg_desc_t *source_seg,
                IN  dmapp_pe_t       source_pe,
                IN  int64_t          comperand,
                IN  int64_t          swaperand);



/* 
 *  Synchronization Functions
 */

/* DMAPP applications use synchronization functions to determine
 * when side-effects of locally initiated, non-blocking RMA requests
 * are globally visible. A process can determine when side-effects
 * of a non-blocking explicit RMA function are globally visible in
 * the system using the following functions.
 */

/* dmapp_syncid_test - Test syncid for completion
 *
 * Parameters:
 * IN/OUT syncid Syncid to be tested for completion
 * OUT    flag   Flag indicating global visibility
 *
 * Returns: 
 * DMAPP_RC_SUCCESS           Operation completed successfully
 * DMAPP_RC_INVALID_PARAM     One or more input parameters is invalid
 * DMAPP_RC_TRANSACTION_ERROR The operation encountered a transaction error
 * 
 * Description:
 * Sets *flag to 1 if the remote memory accesses associated with syncid
 * are globally visible in the system. If the RMA request associated 
 * with the syncid has not completed, *flag is set to 0.
 */

extern dmapp_return_t 
dmapp_syncid_test(INOUT dmapp_syncid_handle_t *syncid,
                  OUT   int                 *flag);

/* dmapp_syncid_wait - Wait for completion of request associated with syncid
 *
 * Parameters:
 * IN/OUT syncid Syncid to be tested for completion
 *
 * Returns: 
 * DMAPP_RC_SUCCESS           Operation completed successfully
 * DMAPP_RC_INVALID_PARAM     One or more input parameters is invalid
 * DMAPP_RC_TRANSACTION_ERROR The operation encountered a transaction error
 *
 * Description:
 * This is the blocking version of dmapp_syncid_test. The function
 * only returns when all remote memory accesses associated with syncid
 * are globally visible in the system.
 */

extern dmapp_return_t 
dmapp_syncid_wait(INOUT dmapp_syncid_handle_t *syncid);


/* DMAPP applications use synchronization functions to determine
 * when side-effects of locally initiated, non-blocking RMA requests
 * are globally visible. A process can determine when side-effects
 * of one or more non-blocking implicit RMA requests are globally
 * visible in the system by using one of the following functions.
 * Invoking gsync style functions does not free resources associated
 * with non-blocking explicit RMA requests. dmapp_syncid_test or
 * dmapp_syncid_wait must still be called for each previously issued
 * non-blocking explicit RMA request in order to free up DMAPP internal
 * resources.
 */

/* dmapp_gsync_test - Test for completion of issued nb implicit requests
 *
 * Parameters:
 * OUT flag   Flag indicating global visibility
 *
 * Returns: 
 * DMAPP_RC_SUCCESS           Operation completed successfully
 * DMAPP_RC_TRANSACTION_ERROR An operation encountered a transaction error
 *
 * Description:
 * Sets *flag to 1 if remote memory accesses associated with previously
 * issued non-blocking implicit RMA requests are globally visible in the 
 * system. Otherwise, *flag is set to 0.
 */

extern dmapp_return_t 
dmapp_gsync_test(OUT int *flag);

/* dmapp_gsync_wait - Wait for completion of issued nb implicit requests
 *
 * Returns: 
 * DMAPP_RC_SUCCESS           Operation completed successfully
 * DMAPP_RC_TRANSACTION_ERROR An operation encountered a transaction error
 *
 * Description:
 * This is the blocking version of dmapp_gsync_test. The function
 * only returns when all remote memory accesses associated with 
 * previously isused non-blocking implicit RMA requests are globally 
 * visible in the system.
 */

extern dmapp_return_t 
dmapp_gsync_wait(void);



/*
 * Symmetric Heap Functions
 */

/* Allocate size bytes of memory from the symmetric heap. The space 
 * returned is left uninitialized. It cannot be assume that the memory
 * returned is zeroed out. There are no address equality guarantees 
 * across ranks.
 */

extern void *
dmapp_sheap_malloc(IN size_t size);

/* The dmapp_sheap_malloc function changes the size of the block to which 
 * ptr points to the size (in bytes) specified by size. The contents of the 
 * block are unchanged up to the lesser of the new and old sizes. If the 
 * new size is larger, the value of the newly allocated portion of the block 
 * is indeterminate. If ptr is a null pointer, dmapp_sheap_malloc behaves 
 * like dmapp_sheap_malloc for the specified size. If size is 0 and ptr is 
 * not a null pointer, the block to which it points is freed. Otherwise, if 
 * ptr does not match a pointer earlier returned by a symmetric heap function, 
 * or if the space has already been deallocated, dmapp_sheap_realloc returns 
 * a null pointer. If the space cannot be allocated, the block to which ptr 
 * points is unchanged.
 */

extern void *
dmapp_sheap_realloc(IN void   *ptr, 
                    IN size_t size);

/* Free a block allocated by dmapp_sheap_malloc, dmapp_sheap_realloc. */

extern void 
dmapp_sheap_free(IN void *ptr);


/*
 * Memory Registration Functions.
 */

/* dmapp_mem_register - Register memory other than the statically linked
 * data segment or the symmetric heap on the fly with the NIC.
 *
 * Parameters:
 * IN    addr      Starting address of the memory region to be registered.
 *                 Must be non-NULL
 * IN    length    Length of the memory region in bytes. Must be > 0
 * INOUT seg_desc  Pointer to segment descriptor.
 *                 Must be non-NULL
 *
 * Returns:
 * DMAPP_RC_SUCCESS   Upon success
 *
 * Description:
 * The memory region described by starting address addr and length 
 * length is registered with the NIC on the fly. The content of the
 * segment descriptor is updated to reflect actual starting address
 * and length of the region which was registered. These values can
 * differ from the input values due to, for instance, rounding.
 */

extern dmapp_return_t
dmapp_mem_register(IN    void             *addr, 
                   IN    uint64_t         length,
                   INOUT dmapp_seg_desc_t *seg_desc);


/* dmapp_mem_unregister - Unregister memory other than the statically linked
 * data segment or the symmetric heap on the fly from the NIC. 
 *
 * Parameters:
 * IN seg_desc   Pointer to segment descriptor for the region to be unregistered.
 *
 * Returns:
 * DMAPP_RC_SUCCESS   Upon success
 *
 * Description:
 * The memory region described by the segment descriptor will be deregistered
 * from the NIC. The memory region must previously have been registered via
 * a call to dmapp_mem_register.
 */

extern dmapp_return_t
dmapp_mem_unregister(INOUT dmapp_seg_desc_t *seg_desc);


/* dmapp_segdesc_compare - Compares two segment descriptors.
 *
 * Parameters:
 * IN    seg_desc1 Pointer to one segment descriptor.
 * IN    seg_desc1 Pointer to second segment descriptor.
 * INOUT flag      Pointer to int
 *
 * Returns:
 * DMAPP_RC_SUCCESS if successfuly
 * In this case, content of flag will be set to 1 if segment descriptors
 * are the same. If they are different, content of flag will be set to zero.
 *
 * Description:
 * The function compares two segment descriptors. If they describe
 * the same memory region, 1 is returned in the flag. If they describe
 * different memory regions, 0 is returned in the flag.
 */

extern dmapp_return_t
dmapp_segdesc_compare(IN    dmapp_seg_desc_t *seg_desc1,
                      IN    dmapp_seg_desc_t *seg_desc2,
                      INOUT int              *flag);


/* dmapp_checkpoint - DMAPP checkpoint function
 *
 * Parameters: None
 *
 * Returns:
 * DMAPP_RC_SUCCESS   Upon success
 * DMAPP_RC_NOT_DONE  Upper-level software has not quieced the network
 *
 * Description:
 * This function destroys GNI resources and DMAPP internal resources
 * that cannot be checkpointed. This is a checkpoint-restart specific
 * function and should not be called outside of that feature.
 * In a threaded environment, this function should only be called once
 * per process.
 */

extern dmapp_return_t
dmapp_checkpoint(void);


/* dmapp_restart - DMAPP restart function
 *
 * Parameters:
 * IN    restart_modes
 *
 * Returns:
 * DMAPP_RC_SUCCESS   Upon success
 *
 * Description:
 * This function will reinitialize GNI resources and DMAPP internal
 * resources which could not be checkpointed. This is a checkpoint-restart
 * specific function and should not be called outside of that feature.
 * In a threaded environment, this function should only be called once
 * per process.
 */

extern dmapp_return_t
dmapp_restart(IN uint32_t restart_modes);


/* dmapp_register_process_cb - Register a progress callback function which DMAPP
 *
 * Parameters:
 * IN progress_cb  Callback function to be registered with DMAPP
 * IN data         Pointer to optional data supplied as an argument
 *                 to the callback function.  Can be NULL.
 *
 * Returns:
 * DMAPP_RC_SUCCESS           Operation completed successfully
 * DMAPP_RC_INVALID_PARAM     One or more input parameters is invalid
 * DMAPP_RC_NO_SPACE          Too many callback functions already registered
 *
 * Description:
 * This function provides a means for software with special progress
 * requirements to be able to use DMAPP even when invoking DMAPP calls
 * which may block internally within DMAPP, e.g. dmapp_get.
 *
 * A progress callback function should return 0 upon success, -1 upon
 * failure.
 */

extern dmapp_return_t
dmapp_register_progress_cb(IN int (*progress_cb)(void *),IN void *data);


/* dmapp_deregister_process_cb - Deregisters a progress callback function which DMAPP
 *                               which had previously been registered using
 *                               dmapp_register_progress_cb
 * 
 * Parameters:
 * IN progress_cb  Callback function to be deregistered
 * 
 * Returns:
 * DMAPP_RC_SUCCESS           Operation completed successfully
 * DMAPP_RC_INVALID_PARAM     One or more input parameters is invalid
 * DMAPP_RC_NOT_FOUND         Supplied callback was not found in callback list
 *
 * Description:
 * Deregisters a progress callback function which had previously been
 * registered with DMAPP.
 */

extern dmapp_return_t
dmapp_deregister_progress_cb(IN int (*progress_cb)(void *));

/* dmapp_c_pset_create - Create a pset(processs set) to use with DMAPP collective
 *                       operations
 * 
 * Parameters:
 * IN pdesc        Pointer to a dmapp_c_pset_desc_t structure 
 *                 previously initialized to describe the ranks
 *                 contained in the pset
 * IN identifier   Unique, non-zero 64-bit identifier to be associated with
 *                 this pset.  The same value must be supplied by
 *                 all ranks creating the pset.
 * IN modes        Used to specify additional information about the
 *                 pset.  
 * IN attrs        Pointer to a previously initialized attribute structure.
 *                 Can be NULL.
 * OUT pset        Pointer to an opaque dmapp_c_pset_handle_t structure 
 *                 to be used for subsequent DMAPP collective calls.
 *
 * Returns:
 * DMAPP_RC_SUCCESS           Operation completed successfully
 * DMAPP_RC_INVALID_PARAM     One or more input parameters is invalid
 * DMAPP_RC_NO_SPACE          Too many psets have already been created
 * DMAPP_RC_BUSY              The supplied identifier is already in use
 * 
 * Description:
 * This is a purely local function which is used to delimit ranks to be
 * involved in particular collective operations. This must be done by all ranks
 * in the job which are going to be used in the collective operations involving
 * the pset. All ranks must provide the same unique idenitifier value to
 * be associated with the pset.  It is an error for ranks not enumerated in
 * the pdesc to call this function.  DMAPP_RC_INVALID_PARAM will be returned
 * in this case.
 *
 * Note:
 * For psets that are to be used for concatentate operations (dmapp_c_concat_start)
 * the DMAPP_C_PSET_MODE_CONCAT mode bit must be set in the mode argument.
 *
 * Space pointed to by the pdesc and attrs parameters can be reused upon
 * return from this function.
 * 
 */

extern dmapp_return_t 
dmapp_c_pset_create(IN dmapp_c_pset_desc_t *pdesc, 
                    IN uint64_t identifier,
                    IN uint64_t modes,
                    IN dmapp_c_pset_attrs_t *attrs,
                    OUT dmapp_c_pset_handle_t *pset_hndl);


/* dmapp_c_pset_export - Export a pset in preparation for use with DMAPP collective
 *                       operations
 *
 * Parameters:
 * IN pset_hndl   The dmapp_c_pset_handle_t structure previously
 *                returned by a call to dmapp_c_pset_create
 *
 * Returns:
 * DMAPP_RC_SUCCESS           Operation completed successfully
 * DMAPP_RC_INVALID_PARAM     One or more input parameters is invalid
 * DMAPP_RC_NOT_FOUND         One or more of the input ranks described
 *                            by the pset has not invoked dmapp_c_pset_create
 *                            for the pset associated with the identifier
 *                            supplied in the local dmapp_c_pset_create function
 *                            that returned this pset handle
 *
 * Description:
 * Before being used for DMAPP collective operations, the pset must be exported
 * using this function. The application must use some out-of-band synchronization mechanism
 * after calling dmapp_c_pset_create and before invoking this function.  Without
 * this out-of-band synchronization, there is a high probability that for one or
 * more ranks, this function will return DMAPP_RC_NOT_FOUND.
 *
 * Notes:
 * This may be a heavy-weight operation and should not be invoked frequently
 * by the application.   This is not a synchronizing function.
 *
 */

extern dmapp_return_t dmapp_c_pset_export(IN dmapp_c_pset_handle_t pset_hndl);

/* dmapp_c_pset_destroy -  Destroy a pset in previously created using
 *                         dmapp_c_pset_create
 *
 * Parameters:
 * IN pset_hndl   The dmapp_c_pset_handle_t structure previously
 *                returned by a call to dmapp_c_pset_create
 *
 * Returns:
 * DMAPP_RC_SUCCESS           Operation completed successfully
 * DMAPP_RC_INVALID_PARAM     Input handle is invalid
 * DMAPP_RC_BUSY              The handle is currently in use for
 *                            a collective operation
 *
 * Description:
 * This function frees any dmapp internal resources associated with the pset.
 * The pset handle cannot be used for any susbsequent DMAPP collective calls.
 * DMAPP_RC_BUSY is returned if either the pset is still involved in a
 * collective operation, or there are outstanding network transactions
 * associated with previous collective operations on the pset.
 * In such cases, an application can invoke this function at a later
 * time to insure resources are released.
 *
 * Notes:
 * This can be considered a local operation.  This is not an explicitly
 * synchronizing function.
 *
 */

extern dmapp_return_t dmapp_c_pset_destroy(IN dmapp_c_pset_handle_t pset_hndl);

/* dmapp_c_barrier_join -  Initiate a barrier join operation on a pset
 *
 * Parameters:
 * IN pset_hndl            The dmapp_c_pset_handle_t structure from a previously
 *                         exported pset
 *
 * Returns:
 * DMAPP_RC_SUCCESS           The join operation completed successfully
 * DMAPP_RC_INVALID_PARAM     Input handle is invalid
 * DMAPP_RC_BUSY              The handle is currently in use for
 *                            another collective operation
 *
 * Description:
 * This function initiates a barrier join operation on a previously
 * exported pset.  A successful return from this operation only means
 * the barrier operation is proceeding, not that it is complete.
 * A return value of DMAPP_RC_BUSY indicates the pset is currently
 * in use for another collective operation.
 *
 */

extern dmapp_return_t dmapp_c_barrier_join(IN dmapp_c_pset_handle_t pset_hndl);

/* dmapp_c_greduce_start -  Initiate a global reduction operation over the pset.
 *
 * Parameters:
 * IN pset_hndl            The dmapp_c_pset_handle_t structure from a previously
 *                         exported pset
 * IN in                   Pointer to buffer containing this PE's contribution
 *                         to the global reduction over the pset
 * OUT out                 Pointer to buffer where the result of the global reduction
 *                         operation is to be returned
 * IN nelems               Number of elements in the global operation  
 *
 * IN type                 Data type of the elements in the global operation
 *
 * IN op                   Type of global reduction operation
 *
 * Returns:
 * DMAPP_RC_SUCCESS           The operation was successfully initiated
 * DMAPP_RC_INVALID_PARAM     One or more input arguments is invalid
 * DMAPP_RC_BUSY              The handle is currently in use for
 *                            another collective operation
 *
 * Description:
 * This function initiates a global reduction operation on a previously
 * exported pset.  A successful return from this operation only means
 * the operation is proceeding, not that it is complete.
 * A return value of DMAPP_RC_BUSY indicates the pset is currently
 * in use for another collective operation.
 *
 * Notes:
 * The number of elements specified in the nelems argument cannot exceed
 * the value returned by dmapp_c_greduce_nelems_max for the given
 * data type and operation type.
 */

extern dmapp_return_t dmapp_c_greduce_start(IN dmapp_c_pset_handle_t pset,
                                            IN void *in,
                                            OUT void *out, 
                                            IN uint32_t nelems,
                                            IN dmapp_c_type_t type, 
                                            IN dmapp_c_op_t op);

/* dmapp_c_greduce_nelems_max -  Return the maximum number of elements of a given
 *                               data type which can be used for a single global
 *                               reduction operation
 *
 * Parameters:
 * IN type                 Data type of the element to be used in a global operation
 * OUT nelems_max          Maximum number of elements of the supplied data type
 *                         which can be used for a single global reduction operation
 *
 * Returns:
 * DMAPP_RC_SUCCESS           The operation was successfully initiated
 * DMAPP_RC_INVALID_PARAM     One or more input arguments is invalid
 *
 * Description:
 * This utility function can be used to determine the maximum number of
 * elements of a given data type which can be used for a single global
 * reduction operation
 *
 * Notes:
 * This function can be used before dmapp_init or other dmapp initialization
 * call is made.
 */

extern dmapp_return_t dmapp_c_greduce_nelems_max(IN dmapp_c_type_t type, 
                                                 OUT uint32_t *nelems_max);

/* dmapp_c_pset_test - Test for completion of a collective operation
 *                     for a given pset
 *
 * Parameters:
 * IN pset_hndl        The dmapp_c_pset_handle_t structure for which
 *                     a collective operation has been initiated.
 *
 * Returns:
 * DMAPP_RC_SUCCESS           Collective operation completed successfully
 * DMAPP_RC_INVALID_PARAM     One or more input parameters is invalid
 * DMAPP_RC_NOT_DONE          Collective operation has not completed
 * DMAPP_RC_TRANSACTION_ERROR An unrecoverable network error was encountered
 *
 * Description:
 * This function can be used to test in a non-blocking manner for completion
 * of an outstanding collective operation for the given pset.
 *
 */

extern dmapp_return_t dmapp_c_pset_test(IN dmapp_c_pset_handle_t pset_hndl);

/* dmapp_c_pset_wait - Wait for completion of a collective operation
 *                     for a given pset
 *
 * Parameters:
 * IN pset_hndl        The dmapp_c_pset_handle_t structure for which
 *                     a collective operation has been initiated.
 *
 * Returns:
 * DMAPP_RC_SUCCESS           Collective operation completed successfully
 * DMAPP_RC_INVALID_PARAM     One or more input parameters is invalid
 * DMAPP_RC_NOT_DONE          Collective operation has not completed
 * DMAPP_RC_TRANSACTION_ERROR An unrecoverable network error was encountered
 *
 * Description:
 * This function is the blocking equivalent of dmapp_c_pset_test.
 */

extern dmapp_return_t dmapp_c_pset_wait(IN dmapp_c_pset_handle_t pset_hndl);

/* dmapp_c_pset_cancel_op - Cancel an outstanding collective operation on a pset
 *
 * Parameters:
 * IN pset_hndl   The dmapp_c_pset_handle_t structure previously
 *                returned by a call to dmapp_c_pset_create
 *
 * Returns:
 * DMAPP_RC_SUCCESS           Operation completed successfully
 * DMAPP_RC_INVALID_PARAM     Input handle is invalid
 *
 * Description:
 * This function cancels an outstanding collective operation on a given
 * pset.
 *
 * Notes:
 * The action of this function is purely local.  Successful return
 * from this function only means that any local DMAPP internal resources
 * associated with the outstanding operation have been released.
 * This function is a no-op if there is no outstanding collective
 * operation on the pset. Once an operation on a pset has been cancelled,
 * the pset cannot be used for any further collective operations.  The
 * only allowed operation on the pset after cancellation is to destroy
 * the pset.
 *
 */

dmapp_return_t dmapp_c_pset_cancel_op(IN dmapp_c_pset_handle_t pset_hndl);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* dmapp.h */
