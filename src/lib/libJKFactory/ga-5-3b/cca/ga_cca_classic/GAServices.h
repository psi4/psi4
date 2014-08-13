#ifndef _GA_SERVICES_H
#define _GA_SERVICES_H

/**              
 * GAServices : Global Arrays Services class.
 * 
 * Author: Manoj Kumar Krishnan, PNNL. 
 * 
 * Collecting the global information: who am I, and how many processors 
 * are being used. Initialize the communication library (either MPI or 
 * TCSMSG) and Global array. Allocate momory to be used by GA by calling MA
 * and create global arrays.
 */

class GAServices : public virtual ::classic::gov::cca::Component, 
		public virtual GAClassicPort, 
		public virtual ::classic::gov::cca::DistArrayTemplFactoryPort,
		public virtual ::classic::gov::cca::DistArrayDescrFactoryPort {
  
 public:
  /**  
   * Null-constructor. The component won't really be 'alive'
   * much at all until after setServices is called on it.
   */
  GAServices();
  
  /** Destructor. */
  virtual ~GAServices();
  
  /**  
   * Creates an ndim-dimensional array using the regular distribution model 
   * and returns integer handle representing the array. 
   
   * The array can be distributed evenly or not. The control over the 
   * distribution is accomplished by specifying chunk (block) size for all or 
   * some of array dimensions.
   
   * For example, for a 2-dimensional array, setting chunk[0]=dim[0] gives 
   * distribution by vertical strips (chunk[0]*dims[0]); 
   * setting chunk[1]=dim[1] gives distribution by horizontal strips 
   * (chunk[1]*dims[1]). Actual chunks will be modified so that they are at 
   * least the size of the minimum and each process has either zero or one 
   * chunk. Specifying chunk[i] as <1 will cause that dimension to be 
   * distributed evenly. 
   
   * As a convenience, when chunk is specified as NULL, the entire array is 
   * distributed evenly.
   
   * \n This is a collective operation. 
   
   * @param arrayname  - a unique character string               [input]
   * @param type        - data type(MT_F_DBL,MT_F_INT,MT_F_DCPL)  [input]
   * @param ndim        - number of array dimensions              [input]
   * @param dims[ndim]  - array of dimensions                     [input]
   * @param chunk[ndim] - array of chunks, each element specifies 
   * minimum size that given dimensions should be chunked up into [input]
   
   * @return Returns pointer to GlobalArray object created. Returns
   * NULL if it fails to create a GA object.
   */
  GlobalArray * createGA(int type, int ndim, int dims[], char *arrayname, 
			 int chunk[]);

  /**
   * Creates an array by following the user-specified distribution and 
   * returns integer handle representing the array. 
     
   * The distribution is specified as a Cartesian product of distributions 
   * for each dimension. The array indices start at 0. For example, the 
   * following figure demonstrates distribution of a 2-dimensional array 8x10 
   * on 6 (or more) processors. nblock[2]={3,2}, the size of map array is s=5 
   * and array map contains the following elements map={0,2,8, 0, 5}. The 
   * distribution is nonuniform because, P1 and P4 get 20 elements each and 
   * processors P0,P2,P3, and P5 only 10 elements each. 
   *        
   * <TABLE>
   * <TR> <TD>5</TD>  <TD>5</TD>  </TR>
   * <TR> <TD>P0</TD> <TD>P3</TD> <TD>2</TD> </TR>
   * <TR> <TD>P1</TD> <TD>P4</TD> <TD>4</TD> </TR>
   * <TR> <TD>P2</TD> <TD>P5</TD> <TD>2</TD> </TR>
   *  </TABLE>
   *
   * \n This is a collective operation. 
   * @param arrayname    - a unique character string          [input]
   * @param type  - MA data type (MT_F_DBL,MT_F_INT,MT_F_DCPL) [input]
   * @param ndim  - number of array dimensions                 [input]
   * @param  dims - array of dimension values                 [input]
   * @param block[ndim] - no. of blocks each dimension is divided into [input]
   * @param maps[s]  - starting index for for each block; the size s is a sum 
   * all elements of nblock array      [input]
   * @return Returns pointer to GlobalArray object created. Returns
   * NULL if it fails to create a GA object.
   */
  GlobalArray * createGA(int type, int ndim, int dims[], char *arrayname, 
			 int maps[], int block[]);

  /**
   * Creates a new array by applying all the properties of another existing 
   * array.
   * \n This is a collective operation. 
   * @param arrayname    - a character string                 [input]
   * @param g_b           - integer handle for reference array [input]
   * @return Returns pointer to GlobalArray object created. Returns
   * NULL if it fails to create a GA object.
   */
  GlobalArray * createGA(const GlobalArray *g_b, char *arrayname);
  
  /**
   * Creates a new array by applying all the properties of another existing 
   * array.
   * \n This is a collective operation. 
   * @param g_b           - integer handle for reference array [input]
   * @return Returns pointer to GlobalArray object created. Returns
   * NULL if it fails to create a GA object.
   */
  GlobalArray * createGA(const GlobalArray &g_b);
  
  /**
   * Creates a 10x10 global array of type "double"(default).
   * @return Returns pointer to GlobalArray object created. Returns
   * NULL if it fails to create a GA object.
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
   * @param array_name   - a unique character string                [input]
   * @param type         - data type (MT_DBL,MT_INT,MT_DCPL)        [input]
   * @param ndim         - number of array dimensions               [input]
   * @param dims[ndim]   - array of dimensions                      [input]
   * @param width[ndim]  - array of ghost cell widths               [input]
   * @param chunk[ndim]  - array of chunks, each element specifies
   *                       minimum size that given dimensions should be
   *                       chunked up into                          [input]
   *
   * @returns Returns pointer to GlobalArray object created. Returns
   * NULL if it fails to create a GA object.
   */
   GlobalArray * createGA_Ghosts(int type, int ndim, int dims[], 
				 int width[], char *array_name, 
				 int chunk[]);
   
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
   *  </TABLE>
   *
   * The array width[] is used to control the width of the ghost cell 
   * boundary around the visible data on each processor. The local data 
   * of the global array residing on each processor will have a layer 
   * width[n] ghosts cells wide on either side of the visible data along 
   * the dimension n. 
   *
   * @param array_name   - a unique character string                [input]
   * @param type         - data type (MT_DBL,MT_INT,MT_DCPL)        [input]
   * @param ndim         - number of array dimensions               [input]
   * @param dims[ndim]   - array of dimensions                      [input]
   * @param width[ndim]  - array of ghost cell widths               [input]
   * @param nblock[ndim] - no. of blocks each dimension is divided into[input]
   * @param  map[s]      - starting index for for each block; the size     
   *                       s is a sum of all elements of nblock array[input]
   *
   * @return Returns pointer to GlobalArray object created. Returns
   * NULL if it fails to create a GA object.
   * \n This is a collective operation. 
   */
   GlobalArray * createGA_Ghosts(int type, int ndim, int dims[], 
				 int width[], char *array_name, int map[], 
				 int nblock[]);
  
  /**
   * @param lenbuf      - length of buffer      [input]
   * @param buf[lenbuf] - data                  [input/output]
   * @param root        - root process          [input]
   *
   * Broadcast from process root to all other processes a message of 
   * length lenbuf. This is operation is provided only for convenience 
   * purposes: it is available regardless of the message-passing library 
   * that GA is running with. 
   * \n This is a collective operation. 
   */
  void brdcst(void *buf, int lenbuf, int root);

  /**
   * This functions returns the total number of nodes that the program is 
   * running on. On SMP architectures, this will be less than or equal to 
   * the total number of processors. 
   * \n This is a  local operation. 
   */
  int clusterNnodes();
  
  /**  
   * This function returns the node ID of the process.  On SMP architectures 
   * with more than one processor per node, several processes may return the
   * same node id. 
   * \n This is a  local operation. 
   */
  int clusterNodeid();
  
  /**
   * This function returns the number of processors available on node inode. 
   * \n This is a  local operation. 
   * @param inode                [input]
   */
  int clusterNprocs(int inode);
  
  /**
   * This function returns the processor id associated with node inode and 
   * the local processor id iproc. If node inode has N processors, then the 
   * value of iproc lies between 0 and N-1. 
   * @param inode,iproc          [input]
   * \n This is a  local operation. 
   */
  int clusterProcid(int inode, int iproc);

   /**
   * Creates a set containing the number of mutexes. Returns 0 if the 
   * opereation succeeded or 1 when failed. Mutex is a simple 
   * synchronization object used to protect Critical Sections. Only one 
   * set of mutexes can exist at a time. Array of mutexes can be created 
   * and destroyed as many times as needed. 
   * Mutexes are numbered: 0, ..., number -1. 
   * \n This is a collective operation. 
   * @param number  - number of mutexes in mutex array   [input]
   */
  int createMutexes(int number);
    
  /** 
   * Destroys the set of mutexes created with ga_create_mutexes. Returns 0 
   * if the operation succeeded or 1 when failed. 
   * \n This is a collective operation. 
   */
  int destroyMutexes();
  
  /**
   * @param n     - number of elements      [input]
   * @param x[n]  - array of elements       [input/output]
   * @param op    - operator                [input]
   * 
   * Double Global OPeration. 
   *
   * X(1:N) is a vector present on each process. DGOP 'sums' elements of 
   * X accross all nodes using the commutative operator OP. The result is 
   * broadcast to all nodes. Supported operations include '+', '*', 'max', 
   * 'min', 'absmax', 'absmin'. The use of lowerecase for operators is 
   * necessary. This is operation is provided only for convenience purposes: 
   * it is available regardless of the message-passing library that GA is 
   * running with. \n This is a collective operation. 
   */
  void dgop(double x[], int n, char *op);

  /**
   * @param array_name    - a character string                 [input]
   * @param g_a           - integer handle for reference array [input]
   *
   * Creates a new array by applying all the properties of another existing 
   * array. It returns array handle. 
   * Return value: a non-zero array handle means the call was succesful. 
   * \n This is a collective operation. 
   */
  int duplicate(int g_a, char* array_name);
  
  /**
   * To be called in case of an error. Print an error message and an integer 
   * value that represents error code. Releases some system resources. 
   * This is the required way of aborting the program execution. 
   * This operation is local. 
   * @param message  - string to print          [input]
   * @param code     - code to print            [input]
   */
  void error(const char *message, int code);

  /**
  * Blocks the calling process until all the data transfers corresponding to 
  * GA operations called after ga_init_fence complete. For example, since 
  * ga_put might return before the data reaches the final destination, 
  * ga_init_fence and ga_fence allow process to wait until the data tranfer 
  * is fully completed: 
  *
  *          ga_init_fence();
  *          ga_put(g_a, ...);
  *          ga_fence();
  * 
  * ga_fence must be called after ga_init_fence. A barrier, ga_sync, assures 
  * completion of all data transfers and implicitly cancels all outstanding
  * ga_init_fence calls. ga_init_fence and ga_fence must be used in pairs, 
  * multiple calls to ga_fence require the same number of corresponding
  * ga_init_fence calls. ga_init_fence/ga_fence pairs can be nested. 
  * 
  * ga_fence works for multiple GA operations. For example: 
  * 
  *         ga_init_fence();
  *         ga_put(g_a, ...);
  *         ga_scatter(g_a, ...);
  *         ga_put(g_b, ...);
  *         ga_fence();
  *
  * The calling process will be blocked until data movements initiated by 
  *two calls to ga_put and one ga_scatter complete. 
  */
  void fence();

  /**
   * @param n     - number of elements      [input]
   * @param x[n]  - array of elements       [input/output]
   * @param op    - operator                [input]
   * 
   * Integer Global OPeration. The integer (more precisely long) version 
   * of ga_dgop described above, also include the bitwise OR operation. 
   * This is operation is provided only for convenience purposes: it is 
   * available regardless of the message-passing library that GA is running 
   * with. \n This is a collective operation. 
   */
  void igop(Integer x[], int n, char *op);
  
  /**
   * Initializes tracing of completion status of data movement operations. 
   * This operation is local. 
   */
  void initFence();

  /**
   * Returns amount of memory (in bytes) used in the allocated global 
   * arrays on the calling processor. This operation is local. 
   */
  size_t inquireMemory();

   /**
   * 
   * Long Global OPeration. 
   *
   * X(1:N) is a vector present on each process. LGOP 'sums' elements of 
   * X accross all nodes using the commutative operator OP. The result is 
   * broadcast to all nodes. Supported operations include '+', '*', 'max', 
   * 'min', 'absmax', 'absmin'. The use of lowerecase for operators is 
   * necessary. This is operation is provided only for convenience purposes: 
   * it is available regardless of the message-passing library that GA is 
   * running with. \n This is a collective operation. 
   * @param n     - number of elements      [input]
   * @param x[n]  - array of elements       [input/output]
   * @param op    - operator                [input]
   */
  void lgop(long x[], int n, char *op);

  /**
   * @param mutex - mutex object id  [input]
   *
   * Locks a mutex object identified by the mutex number. It is a fatal 
   * error for a process to attempt to lock a mutex which was already 
   * locked by this process. 
   */
  void lock(int mutex);

  /**
   * GA Collective calls has Sync calls at the begining and ending of
   * of the call. Sometimes there may be some redundacy in sync calls, which
   * can be avoided by masking the sync operations. 
   * @ param first - masks the sync at the begining of the collective call.
   * @ param last  - masks the sync at the end of the collective call.
   * setting the parameters as zero will mask (disable) the call. Any non-zero 
   * value will enable the call. Initially these params are set to non-zero 
   * value.
   */
  void maskSync(int first, int last);

  /**
   * @return Returns amount of memory (in bytes) left for allocation of new 
   * global arrays on the calling processor. 
   
   * @note If GA_uses_ma returns true, then GA_Memory_avail returns the 
   * lesser of the amount available under the GA limit and the amount 
   * available from MA (according to ma_inquire_avail operation). 
   * If no GA limit has been set, it returns what MA says is available. 
   * If ( ! GA_Uses_ma() && ! GA_Memory_limited() ) returns < 0, indicating 
   * that the bound on currently available memory cannot be determined. 
   * This operation is local. 
   */
  int memoryAvailable() ;
  
  /**
   * Indicates if limit is set on memory usage in Global Arrays on the 
   * calling processor. "1" means "yes", "0" means "no". This operation 
   * is local. 
   */
  int memoryLimited();

  /**
   * Returns the GA process id (0, ..., ga_Nnodes()-1) of the requesting 
   * compute process. This operation is local. 
   */
  int nodeid();
  
  /**
   * Returns the number of the GA compute (user) processes. 
   * This operation is local. 
   */
  int nodes();
  
  /** 
   * This non-collective (MIMD) operation prints information about: 
   *
   * number of calls to the GA create/duplicate, destroy, get, put, scatter,
   * gather, and read_and_inc operations total amount of data moved in the 
   * GA primitive operations amount of data moved in GA primitive operations 
   * to logicaly remote locations maximum memory consumption in global 
   * arrays, and number of requests serviced in the interrupt-driven 
   * implementations by the calling process. This operation is local. 
   */
  void printStats();
  
  /**
   * @param limit    - the amount of memory in bytes per process    [input]
   * 
   * Sets the amount of memory to be used (in bytes) per process. 
   * \n This is a local operation. 
   */
  void setMemoryLimit(size_t limit);
  
  /** 
   * @param verbose     - If true print distribution info [input]
   * Prints info about allocated arrays. 
   */
  void summarize(int verbose);
  
  /**
   * Synchronize processes (a barrier) and ensure that all GA operations 
   * completed. 
   * \n This is a collective operation. 
   */
  void sync();
  
   /**
   * @param mutex  - mutex object id [input]
   * 
   * Unlocks a mutex object identified by the mutex number. It is a fatal 
   * error for a process to attempt to unlock a mutex which has not been 
   * locked by this process. 
   */
  void unlock(int mutex);

  /**
   * Returns "1" if memory in global arrays comes from the Memory Allocator 
   * (MA). "0"means that memory comes from another source, for example 
   * System V shared memory is used. This operation is local. 
   */
  int usesMA();

  /**
   * Returns "1" if uses fortran API, else returns "0"
   */
  int usesFAPI();
  
  /**
   * *****************************
   *        DADF Interfaces
   * *****************************
   */
  /****************************************************************************
   * DistArrayTemplFactoryPort
   ***************************************************************************/
 
  /** Return an uninitialized template object. */
  virtual DistArrayTemplate * createTemplate(std::string name);
 
  /** Return a template object initialized with the contents of
      another, but not frozen against modification. */
  virtual DistArrayTemplate * cloneTemplate(
					    DistArrayTemplate * original, std::string cloneName);
  
  /** Destroy an existing template. */
  virtual int destroyTemplate(DistArrayTemplate * & victim);
 
  /****************************************************************************
   * DistArrayDescrFactoryPort
   ***************************************************************************/
 
  /** Return an uninitialized descriptor object. */
  virtual DistArrayDescriptor * createDescriptor(std::string name);
 
  /** Return a descriptor object initialized with the contents of
      another, but not frozen against modification. */
  virtual DistArrayDescriptor * cloneDescriptor(
						DistArrayDescriptor * original, std::string cloneName);
 
  /** Destroy an existing descriptor. */
  virtual int destroyDescriptor(DistArrayDescriptor * & victim);
  
  /** Return an uninitialized ga-dadf distributed array object. */
  virtual DistArray * createArray(std::string name);
  
  /** Return an distributed array object initialized with the contents of
      another, but not frozen against modification. */
  virtual DistArray * cloneArray(DistArray* original, 
				 std::string cloneName);
  
  /** Destroy an existing distributed array. */
  virtual int destroyArray(DistArray* & victim);
   
  
  /**
   * The components containing framework provides services through
   * the Services interface. This will be called with a Services
   * when the component is created and with 0/NULL when the component
   * is about to be destroyed.
   * If the framework unrecoverably misbehaves (a port definition
   * call on 'cc' fails, then this call will not return (abort() is invoked)
   * This component is not compatible with such a framework.
   */
  void setServices(::classic::gov::cca::Services *cc);

 private:
  /** The services we receive from the outer framework or container. */
  ::classic::gov::cca::Services *svc;
  
  /** Map to track templates we've handed out and destroyed */
  std::map<DistArrayTemplate *, int> templateRecord;
 
  /** Map to track descriptors we've handed out and destroyed */
  std::map<DistArrayDescriptor *, int> descriptorRecord;
 
  /** Map to track arrays we've handed out and destroyed */
  std::map<DistArray *, int> arrayRecord;
 
  /** List outstanding templates on stderr
      @param label (in) String to label output
  */
  void listTemplates(std::string label);
 
  /** List outstanding descriptors on stderr
      @param label (in) String to label output
  */
  void listDescriptors(std::string label);

  /** List outstanding arrays on stderr
      @param label (in) String to label output
  */
  void listArrays(std::string label);  
};


#endif // _GA_SERVICES_H
