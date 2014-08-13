#ifndef _SERVICES_H
#define _SERVICES_H

namespace GA {

class GlobalArray;

/**  
 * Creates an ndim-dimensional array using the regular distribution model 
 * and returns integer handle representing the array. 
 *
 * The array can be distributed evenly or not. The control over the 
 * distribution is accomplished by specifying chunk (block) size for all or 
 * some of array dimensions.
 *
 * For example, for a 2-dimensional array, setting chunk[0]=dim[0] gives 
 * distribution by vertical strips (chunk[0]*dims[0]); 
 * setting chunk[1]=dim[1] gives distribution by horizontal strips 
 * (chunk[1]*dims[1]). Actual chunks will be modified so that they are at 
 * least the size of the minimum and each process has either zero or one 
 * chunk. Specifying chunk[i] as <1 will cause that dimension to be 
 * distributed evenly. 
 *
 * As a convenience, when chunk is specified as NULL, the entire array is 
 * distributed evenly.
 *
 * This is a collective operation. 
 *
 * @param[in] type        data type(MT_F_DBL,MT_F_INT,MT_F_DCPL)
 * @param[in] ndim        number of array dimensions
 * @param[in] dims[ndim]  array of dimensions
 * @param[in] arrayname   a unique character string
 * @param[in] chunk[ndim] array of chunks, each element specifies 
 *                        minimum size that given dimensions should be
 *                        chunked up into
 *
 * @return pointer to GlobalArray object created; NULL if it fails
 */
GlobalArray* createGA(int type, int ndim, int dims[], char *arrayname, 
        int chunk[]);

/**
 * Creates an array by following the user-specified distribution and 
 * returns integer handle representing the array. 
 *
 * The distribution is specified as a Cartesian product of distributions 
 * for each dimension. The array indices start at 0. For example, the 
 * following figure demonstrates distribution of a 2-dimensional array 8x10 
 * on 6 (or more) processors. nblock[2]={3,2}, the size of map array is s=5 
 * and array map contains the following elements map={0,2,6, 0, 5}. The 
 * distribution is nonuniform because, P1 and P4 get 20 elements each and 
 * processors P0,P2,P3, and P5 only 10 elements each. 
 *        
 * <TABLE>
 * <TR> <TD>5</TD>  <TD>5</TD>  </TR>
 * <TR> <TD>P0</TD> <TD>P3</TD> <TD>2</TD> </TR>
 * <TR> <TD>P1</TD> <TD>P4</TD> <TD>4</TD> </TR>
 * <TR> <TD>P2</TD> <TD>P5</TD> <TD>2</TD> </TR>
 * </TABLE>
 *
 * This is a collective operation. 
 *
 * @param[in] arrayname   a unique character string
 * @param[in] type        MA data type (MT_F_DBL,MT_F_INT,MT_F_DCPL)
 * @param[in] ndim        number of array dimensions
 * @param[in] dims        array of dimension values
 * @param[in] block       [ndim] no. of blocks each dimension is divided into
 * @param[in] maps        [s] starting index for for each block;
 *                        the size s is a sum all elements of nblock array
 *
 * @return pointer to GlobalArray object created; NULL if it fails
 */
GlobalArray * createGA(int type, int ndim, int dims[], char *arrayname, 
        int block[], int maps[]);

/**
 * Creates a new array by applying all the properties of another existing 
 * array.
 *
 * This is a collective operation. 
 *
 * @param[in] arrayname a character string
 * @param[in] g_b       integer handle for reference array
 *
 * @return pointer to GlobalArray object created; NULL if it fails
 */
GlobalArray * createGA(const GlobalArray *g_b, char *arrayname);

/**
 * Creates a new array by applying all the properties of another existing 
 * array.
 *
 * This is a collective operation. 
 *
 * @param[in] g_b integer handle for reference array
 *
 * @return pointer to GlobalArray object created; NULL if it fails
 */
GlobalArray * createGA(const GlobalArray &g_b);

/**
 * Creates a 10x10 global array of type "double"(default).
 *
 * This is a collective operation. 
 *
 * @return pointer to GlobalArray object created; NULL if it fails
 */
GlobalArray * createGA();  

/**
 * Creates an ndim-dimensional array with a layer of ghost cells around 
 * the visible data on each processor using the regular distribution 
 * model and returns an integer handle representing the array. 
 * The array can be distributed evenly or not evenly. The control over 
 * the distribution is accomplished by specifying chunk (block) size for 
 * all or some of the array dimensions. For example, for a 2-dimensional 
 * array, setting chunk(1)=dim(1) gives distribution by vertical strips 
 * (chunk(1)*dims(1)); setting chunk(2)=dim(2) gives distribution by 
 * horizontal strips (chunk(2)*dims(2)). Actual chunks will be modified 
 * so that they are at least the size of the minimum and each process 
 * has either zero or one chunk. Specifying chunk(i) as <1 will cause
 * that dimension (i-th) to be distributed evenly. The  width of the 
 * ghost cell layer in each dimension is specified using the array 
 * width().  The local data of the global array residing on each 
 * processor will have a layer width[n] ghosts cells wide on either 
 * side of the visible data along the dimension n. 
 *
 * This is a collective operation. 
 * 
 * @param[in] array_name a unique character string
 * @param[in] type       data type (MT_DBL,MT_INT,MT_DCPL)
 * @param[in] ndim       number of array dimensions
 * @param[in] dims       [ndim] array of dimensions
 * @param[in] width      [ndim] array of ghost cell widths
 * @param[in] chunk      [ndim] array of chunks, each element specifies
 *                       minimum size that given dimensions should be
 *                       chunked up into
 *
 * @returns pointer to GlobalArray object created; NULL if it fails
 */
GlobalArray * createGA_Ghosts(int type, int ndim, int dims[], 
        int width[], char *array_name, int chunk[]);

/**
 * Creates an array with ghost cells by following the user-specified 
 * distribution and returns integer handle representing the array. 
 * The distribution is specified as a Cartesian product of distributions 
 * for each dimension. For example, the following figure demonstrates 
 * distribution of a 2-dimensional array 8x10 on 6 (or more) processors. 
 * nblock(2)={3,2}, the size of map array is s=5 and array map contains 
 * the following elements map={1,3,7, 1, 6}. The distribution is 
 * nonuniform because, P1 and P4 get 20 elements each and processors 
 * P0,P2,P3, and P5 only 10 elements each. 
 *
 * <TABLE>
 * <TR> <TD>5</TD>  <TD>5</TD>  </TR>
 * <TR> <TD>P0</TD> <TD>P3</TD> <TD>2</TD> </TR>
 * <TR> <TD>P1</TD> <TD>P4</TD> <TD>4</TD> </TR>
 * <TR> <TD>P2</TD> <TD>P5</TD> <TD>2</TD> </TR>
 * </TABLE>
 *
 * The array width[] is used to control the width of the ghost cell 
 * boundary around the visible data on each processor. The local data 
 * of the global array residing on each processor will have a layer 
 * width[n] ghosts cells wide on either side of the visible data along 
 * the dimension n. 
 *
 * This is a collective operation. 
 *
 * @param[in] array_name a unique character string
 * @param[in] type       data type (MT_DBL,MT_INT,MT_DCPL)
 * @param[in] ndim       number of array dimensions
 * @param[in] dims       [ndim] array of dimensions
 * @param[in] width      [ndim] array of ghost cell widths
 * @param[in] nblock     [ndim] no. of blocks each dimension is divided into
 * @param[in] map        [s] starting index for for each block;
 *                       the size s is a sum of all elements of nblock array
 *
 * @return pointer to GlobalArray object created; NULL if it fails
 */
GlobalArray * createGA_Ghosts(int type, int ndim, int dims[], 
        int width[], char *array_name, int map[], 
        int nblock[]);

/**
 * Broadcast from process root to all other processes a message of 
 * length lenbuf. This is operation is provided only for convenience 
 * purposes: it is available regardless of the message-passing library 
 * that GA is running with. 
 *
 * This is a collective operation. 
 *
 * @param[in]     lenbuf length of buffer
 * @param[in,out] buf    [lenbuf] data
 * @param[in]     root   root process
 */
void brdcst(void *buf, int lenbuf, int root);

/**
 * Returns the current value of the internal debug flag.
 *
 * This is a local operation.
 *
 * @return 0 if the debug flag is false, 1 if it is true.
 */
int getDebug();

/**
 * This functions returns the total number of nodes that the program is 
 * running on.
 *
 * On SMP architectures, this will be less than or equal to the total
 * number of processors. 
 *
 * This is a  local operation. 
 *
 * @return the number of nodes the program is running on
 */
int clusterNnodes();

/**  
 * This function returns the node ID of the process.
 *
 * On SMP architectures with more than one processor per node, several
 * processes may return the same node id. 
 *
 * This is a  local operation. 
 *
 * @return the node ID of the process
 */
int clusterNodeid();

/**  
 * This function returns the cluster node ID of the specified process.
 *
 * On SMP architectures with more than one processor per node, several
 * processes may return the same node id. 
 *
 * This is a  local operation. 
 *
 * @return the cluster node ID of the specified process
 */
int clusterProcNodeid(int iproc);

/**
 * This function returns the number of processors available on node inode. 
 *
 * This is a  local operation. 
 *
 * @param[in] inode
 *
 * @return the number of processors available on the given node
 */
int clusterNprocs(int inode);

/**
 * This function returns the processor id associated with node inode and 
 * the local processor id iproc.
 *
 * If node inode has N processors, then the value of iproc lies between
 * 0 and N-1. 
 *
 * This is a  local operation. 
 *
 * @param[in] inode
 * @param[in] iproc
 *
 * @return the processor ID associated with the given node and local processor
 * ID
 */
int clusterProcid(int inode, int iproc);

/**
 * Creates a set containing the number of mutexes.
 *
 * Mutex is a simple synchronization object used to protect Critical
 * Sections. Only one set of mutexes can exist at a time. Array of mutexes
 * can be created and destroyed as many times as needed. 
 * Mutexes are numbered: 0, ..., number -1.
 *
 * This is a collective operation. 
 *
 * @param[in] number of mutexes in mutex array
 *
 * @return 0 if the opereation succeeded or 1 when failed. 
 */
int createMutexes(int number);

/**
 * Remove a user defined data type from GA
 *
 * @param[in] type - user defined data type
 *
 * @return 0 is operation is successful
 *         -2 if type not registered
 *         -1 if type reserved
 */
int deregisterType(int type);

/** 
 * Destroys the set of mutexes created with ga_create_mutexes.
 *
 * This is a collective operation. 
 *
 * @return 0 if the operation succeeded or 1 when failed. 
 */
int destroyMutexes();

/**
 * Double Global OPeration. 
 *
 * X(1:N) is a vector present on each process. DGOP 'sums' elements of 
 * X accross all nodes using the commutative operator OP. The result is 
 * broadcast to all nodes. Supported operations include '+', '*', 'max', 
 * 'min', 'absmax', 'absmin'. The use of lowerecase for operators is 
 * necessary. This is operation is provided only for convenience purposes: 
 * it is available regardless of the message-passing library that GA is 
 * running with.
 *
 * This is a collective operation. 
 *
 * @param[in]     n  number of elements
 * @param[in,out] x  [n] array of elements
 * @param[in]     op operator
 */
void dgop(double x[], int n, char *op);

/**
 * Creates a new array by applying all the properties of another existing 
 * array.
 *
 * This is a collective operation. 
 *
 * @param[in] array_name a character string
 * @param[in] g_a        integer handle for reference array
 *
 * @return array handle; a non-zero array handle means the call was succesful. 
 */
int duplicate(int g_a, char* array_name);

/**
 * To be called in case of an error.
 *
 * Print an error message and an integer value that represents error code.
 * Releases some system resources. 
 * This is the required way of aborting the program execution. 
 *
 * This operation is local. 
 *
 * @param[in] message string to print
 * @param[in] code    code to print
 */
void error(const char *message, int code);

/**
 * Blocks the calling process until all the data transfers corresponding to 
 * GA operations called after ga_init_fence complete.
 *
 * For example, since ga_put might return before the data reaches the final
 * destination, ga_init_fence and ga_fence allow process to wait until the
 * data tranfer is fully completed: 
 *
 * @code
 *   ga_init_fence();
 *   ga_put(g_a, ...);
 *   ga_fence();
 * @endcode
 * 
 * ga_fence must be called after ga_init_fence. A barrier, ga_sync, assures 
 * completion of all data transfers and implicitly cancels all outstanding
 * ga_init_fence calls. ga_init_fence and ga_fence must be used in pairs, 
 * multiple calls to ga_fence require the same number of corresponding
 * ga_init_fence calls. ga_init_fence/ga_fence pairs can be nested. 
 * 
 * ga_fence works for multiple GA operations. For example: 
 * 
 * @code
 *   ga_init_fence();
 *   ga_put(g_a, ...);
 *   ga_scatter(g_a, ...);
 *   ga_put(g_b, ...);
 *   ga_fence();
 * @endcode
 *
 * The calling process will be blocked until data movements initiated by 
 * two calls to ga_put and one ga_scatter complete. 
 */
void fence();

/**
 * Integer Global OPeration.
 *
 * The integer version of ga_dgop described above, also include the bitwise OR
 * operation.  This is operation is provided only for convenience purposes: it
 * is available regardless of the message-passing library that GA is running
 * with.
 *
 * This is a collective operation. 
 *
 * @param[in]     n  number of elements
 * @param[in,out] x  [n] array of elements
 * @param[in]     op operator
 */
void gop(int x[], int n, char *op);

/**
 * Long Global OPeration. 
 *
 * X(1:N) is a vector present on each process. LGOP 'sums' elements of 
 * X accross all nodes using the commutative operator OP. The result is 
 * broadcast to all nodes. Supported operations include '+', '*', 'max', 
 * 'min', 'absmax', 'absmin'. The use of lowerecase for operators is 
 * necessary. This is operation is provided only for convenience purposes: 
 * it is available regardless of the message-passing library that GA is 
 * running with.
 *
 * This is a collective operation. 
 *
 * @param[in]     n  number of elements
 * @param[in,out] x  [n] array of elements
 * @param[in]     op operator
 */
void gop(long x[], int n, char *op);

/**
 * Float Global OPeration. 
 *
 * X(1:N) is a vector present on each process. FGOP 'sums' elements of 
 * X accross all nodes using the commutative operator OP. The result is 
 * broadcast to all nodes. Supported operations include '+', '*', 'max', 
 * 'min', 'absmax', 'absmin'. The use of lowerecase for operators is 
 * necessary. This is operation is provided only for convenience purposes: 
 * it is available regardless of the message-passing library that GA is 
 * running with.
 *
 * This is a collective operation. 
 *
 * @param[in]     n  number of elements
 * @param[in,out] x  [n] array of elements
 * @param[in]     op operator
 */
void gop(float x[], int n, char *op);

/**
 * Double Global OPeration. 
 *
 * X(1:N) is a vector present on each process. DGOP 'sums' elements of 
 * X accross all nodes using the commutative operator OP. The result is 
 * broadcast to all nodes. Supported operations include '+', '*', 'max', 
 * 'min', 'absmax', 'absmin'. The use of lowerecase for operators is 
 * necessary. This is operation is provided only for convenience purposes: 
 * it is available regardless of the message-passing library that GA is 
 * running with.
 *
 * This is a collective operation. 
 *
 * @param[in]     n  number of elements
 * @param[in,out] x  [n] array of elements
 * @param[in]     op operator
 */
void gop(double x[], int n, char *op);

/**
 * Integer Global OPeration.
 *
 * The integer (more precisely long) version of ga_dgop described above,
 * also include the bitwise OR operation. 
 * This is operation is provided only for convenience purposes: it is 
 * available regardless of the message-passing library that GA is running 
 * with.
 *
 * This is a collective operation. 
 *
 * @param[in]     n  number of elements
 * @param[in,out] x  [n] array of elements
 * @param[in]     op operator
 */
void igop(int x[], int n, char *op);

/**
 * Initializes tracing of completion status of data movement operations. 
 *
 * This operation is local. 
 */
void initFence();

/**
 * Returns amount of memory (in bytes) used in the allocated global 
 * arrays on the calling processor.
 *
 * This operation is local. 
 *
 * @return amount of memory (in bytes) used in the allocated global arrays on
 *         the calling processor
 */
size_t inquireMemory();

/**
 * Long Global OPeration. 
 *
 * X(1:N) is a vector present on each process. LGOP 'sums' elements of 
 * X accross all nodes using the commutative operator OP. The result is 
 * broadcast to all nodes. Supported operations include '+', '*', 'max', 
 * 'min', 'absmax', 'absmin'. The use of lowerecase for operators is 
 * necessary. This is operation is provided only for convenience purposes: 
 * it is available regardless of the message-passing library that GA is 
 * running with.
 *
 * This is a collective operation. 
 *
 * @param[in]     n  number of elements
 * @param[in,out] x  [n] array of elements
 * @param[in]     op operator
 */
void lgop(long x[], int n, char *op);

/**
 * Locks a mutex object identified by the mutex number. It is a fatal 
 * error for a process to attempt to lock a mutex which was already 
 * locked by this process. 
 *
 * @param[in] mutex object id
 */
void lock(int mutex);

/**
 * Mask the intrinsic sync operations during collective calls.
 *
 * GA Collective calls has Sync calls at the begining and ending of
 * of the call. Sometimes there may be some redundacy in sync calls, which
 * can be avoided by masking the sync operations. 
 *
 * Setting the parameters as zero will mask (disable) the call. Any non-zero 
 * value will enable the call. Initially these params are set to non-zero 
 * value.
 *
 * @param[in] first masks the sync at the begining of the collective call.
 * @param[in] last  masks the sync at the end of the collective call.
 */
void maskSync(int first, int last);

/**
 * If GA_uses_ma returns true, then GA_Memory_avail returns the 
 * lesser of the amount available under the GA limit and the amount 
 * available from MA (according to ma_inquire_avail operation). 
 * If no GA limit has been set, it returns what MA says is available. 
 * If ( ! GA_Uses_ma() && ! GA_Memory_limited() ) returns < 0, indicating 
 * that the bound on currently available memory cannot be determined. 
 *
 * This operation is local. 
 *
 * @return amount of memory (in bytes) left for allocation of new 
 *         global arrays on the calling processor. 
 *
 */
int memoryAvailable() ;

/**
 * Indicates if limit is set on memory usage in Global Arrays on the 
 * calling processor.
 * 
 * This operation is local. 
 *
 * @return 1 means "yes", "0" means "no".
 */
int memoryLimited();

/**
 * Force completion of a nonblocking operation locally.
 *
 * Waiting on a nonblocking put or an accumulate operation assures that data
 * was injected into the network and the user buffer can be now be reused.
 * Completing a get operation assures data has arrived into the user memory
 * and is ready for use. Wait operation ensures only local completion. Unlike
 * their blocking counterparts, the nonblocking operations are not ordered
 * with respect to the destination. Performance being one reason, the other
 * reason is that by ensuring ordering we incur additional and possibly
 * unnecessary overhead on applications that do not require their operations
 * to be ordered. For cases where ordering is necessary, it can be done by
 * calling a fence operation.  The fence operation is provided to the user to
 * confirm remote completion if needed.
 *
 * This is a local operation.
 *
 * @param[in] nbhandle nonblocking handle
 */
void nbWait(GANbhdl *nbhandle);

/**
 * Returns the GA process id (0, ..., ga_Nnodes()-1) of the requesting 
 * compute process.
 *
 * This operation is local. 
 *
 * @return the GA process ID of the requesting process
 */
int nodeid();

/**
 * Returns the number of the GA compute (user) processes. 
 *
 * This operation is local. 
 *
 * @return the number of GA processes
 */
int nodes();

/** 
 * Print statistical information on GA use.
 *
 * This non-collective (MIMD) operation prints information about:
 *   - number of calls to
 *     - create
 *     - duplicate
 *     - destroy
 *     - get
 *     - put
 *     - scatter
 *     - gather
 *     - read_and_inc operations
 *   - total amount of data moved in the primitive operations
 *   - amount of data moved in the primitive operations to logicaly remote
 *     locations
 *   - maximum memory consumption in global arrays
 *   - number of requests serviced in the interrupt-driven implementations
 *     by the calling process.
 *
 * This operation is local. 
 */
void printStats();

/**
 * Add a user defined data type to GA
 *
 * @param[in] size - size (in bytes) of user defined data type
 *
 * @return  handle for new data type
 */
int registerType(size_t size);

/**
 * This function sets an internal flag in the GA library to either true or
 * false.
 *
 * The value of this flag can be recovered at any time using the
 * getDebug function. The flag is set to false when the the GA library
 * is initialized. This can be useful in a number of debugging situations,
 * especially when examining the behavior of routines that are called in
 * multiple locations in a code. 
 *
 * This is a local operation.
 *
 * @param[in] dbg value to set internal flag
 */
void setDebug(int dbg);

/**
 * Sets the amount of memory to be used (in bytes) per process.
 * 
 * This is a local operation. 
 *
 * @param[in] limit the amount of memory in bytes per process
 */
void setMemoryLimit(size_t limit);

/** 
 * Prints info about allocated arrays. 
 *
 * @param[in] verbose if true print distribution info
 */
void summarize(int verbose);

/**
 * Synchronize processes (a barrier) and ensure that all GA operations 
 * completed. 
 *
 * This is a collective operation. 
 */
void sync();

/**
 * Unlocks a mutex object identified by the mutex number.
 *
 * It is a fatal error for a process to attempt to unlock a mutex which has
 * not been locked by this process. 
 *
 * @param[in] mutex object id
 */
void unlock(int mutex);

/**
 * Returns whether memory comes from internal or external allocator.
 *
 * This operation is local. 
 *
 * @return "1" if memory comes from MA;
 *         "0" if memory comes from another source e.g. System V shared memory
 */
int usesMA();

/**
 * Returns whether GA is using Fortran indexing.
 *
 * @return "1" if uses fortran API, else returns "0"
 */
int usesFAPI();

/**
 * This function return a wall (or elapsed) time on the calling processor.
 *
 * Returns time in seconds representing elapsed wall-clock time
 * since an arbitrary time in the past. Example:
 *
 * @code
 *     double starttime, endtime;
 *     starttime = GA::wtime();
 *     // {{.... code snippet to be timed ....}}
 *     endtime   = GA::wtime();
 *     printf("Time taken = %lf seconds\n", endtime-starttime);
 * @endcode
 *
 * This is a local operation.
 *
 * @note This function is only available in release 4.1 or greater.
 */
double wtime();

}

#endif /* _SERVICES_H */
