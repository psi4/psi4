#ifndef  _GLOBALARRAY_H
#define  _GLOBALARRAY_H

namespace GA {

class PGroup;

/**
 * This is the GlobalArray class.
 */
class GlobalArray { 

 public:
  
  /**  
   * Creates an ndim-dimensional array using the regular distribution
   * model and returns integer handle representing the array. 
   
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
   
   * This is a collective operation. 
   
   * @param[in] type      data type(MT_F_DBL,MT_F_INT,MT_F_DCPL)
   * @param[in] ndim      number of array dimensions
   * @param[in] dims      [ndim] array of dimensions
   * @param[in] arrayname a unique character string
   * @param[in] chunk     [ndim] array of chunks, each element specifies
   *                      minimum size that given dimensions should be chunked
   *                      up into
   */
  GlobalArray(int type, int ndim, int dims[], char *arrayname, int chunk[]);

  /**  
   * @copydoc GlobalArray::GlobalArray(int,int,int[],char*,int[])
   * @param[in] p_handle processor group handle
   */
  GlobalArray(int type, int ndim, int dims[], char *arrayname, int chunk[],
              PGroup* p_handle);

  /**
   * @copydoc GlobalArray::GlobalArray(int,int,int[],char*,int[])
   */
  GlobalArray(int type, int ndim, int64_t dims[], char *arrayname,
              int64_t chunk[]);

  /**
   * @copydoc GlobalArray::GlobalArray(int,int,int[],char*,int[])
   * @param[in] p_handle processor group handle
   */
  GlobalArray(int type, int ndim, int64_t dims[], char *arrayname,
              int64_t chunk[], PGroup* p_handle);
      
  /**
   * Creates an array by following the user-specified distribution.
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
   * @param[in] type      MA data type (MT_F_DBL,MT_F_INT,MT_F_DCPL)
   * @param[in] ndim      number of array dimensions
   * @param[in] dims      array of dimension values
   * @param[in] arrayname a unique character string
   * @param[in] block     [ndim] no. of blocks each dimension is divided into
   * @param[in] maps      [s] starting index for for each block;
   *                      the size s is a sum all elements of nblock array
   */
  GlobalArray(int type, int ndim, int dims[], char *arrayname, int block[],
	      int maps[]);

  /**
   * @copydoc GlobalArray::GlobalArray(int,int,int[],char*,int[],int[])
   * @param[in] p_handle  processor group handle
   */
  GlobalArray(int type, int ndim, int dims[], char *arrayname, int block[],
	      int maps[], PGroup* p_handle);
  
  /**
   * @copydoc GlobalArray::GlobalArray(int,int,int[],char*,int[],int[])
   */
  GlobalArray(int type, int ndim, int64_t dims[], char *arrayname,
              int64_t block[], int64_t maps[]);

  /**
   * @copydoc GlobalArray::GlobalArray(int,int,int[],char*,int[],int[])
   * @param[in] p_handle  processor group handle
   */
  GlobalArray(int type, int ndim, int64_t dims[], char *arrayname,
              int64_t block[], int64_t maps[], PGroup* p_handle);
      
  /**
   * Creates an ndim-dimensional array with a layer of ghost cells around 
   * the visible data on each processor using the regular distribution model.
   *
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
   * @param[in] type       data type (MT_DBL,MT_INT,MT_DCPL)
   * @param[in] ndim       number of array dimensions
   * @param[in] dims       [ndim] array of dimensions
   * @param[in] width      [ndim] array of ghost cell widths
   * @param[in] arrayname  a unique character string
   * @param[in] chunk      [ndim] array of chunks, each element specifies
   *                       minimum size that given dimensions should be
   *                       chunked up into
   * @param[in] ghosts     this is a dummy parameter: added to increase the 
   *                       number of arguments, inorder to avoid the conflicts
   *                       among constructors. (ghosts = 'g' or 'G')
   */
  GlobalArray(int type, int ndim, int dims[], int width[], char *arrayname, 
	      int chunk[], char ghosts);

  /**
   * @copydoc GlobalArray::GlobalArray(int,int,int[],int[],char*,int[],char)
   * @param[in] p_handle   processor group handle
   */
  GlobalArray(int type, int ndim, int dims[], int width[], char *arrayname, 
	      int chunk[], PGroup* p_handle, char ghosts);
  
  /**
   * @copydoc GlobalArray::GlobalArray(int,int,int[],int[],char*,int[],char)
   */
  GlobalArray(int type, int ndim, int64_t dims[], int64_t width[],
          char *arrayname, int64_t chunk[], char ghosts);

  /**
   * @copydoc GlobalArray::GlobalArray(int,int,int[],int[],char*,int[],char)
   * @param[in] p_handle   processor group handle
   */
  GlobalArray(int type, int ndim, int64_t dims[], int64_t width[],
          char *arrayname, int64_t chunk[], PGroup* p_handle, char ghosts);
  
  /**
   * Creates an array with ghost cells by following the user-specified 
   * distribution.
   *
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
   * the dimension n. This is a collective operation. 
   *
   * @param[in] type      data type (MT_DBL,MT_INT,MT_DCPL)
   * @param[in] ndim      number of array dimensions
   * @param[in] dims      [ndim] array of dimensions
   * @param[in] width     [ndim] array of ghost cell widths
   * @param[in] arrayname a unique character string
   * @param[in] block     [ndim] no. of blocks each dimension is divided into
   * @param[in] maps      [s] starting index for for each block;
   *                      the size s is a sum of all elements of nblock array
   * @param[in] ghosts    this is a dummy parameter: added to increase the 
   *                      number of arguments, inorder to avoid the conflicts
   *                      among constructors. (ghosts = 'g' or 'G')
   */
  GlobalArray(int type, int ndim, int dims[], int width[], char *arrayname, 
	      int block[], int maps[], char ghosts);

  /**
   * @copydoc GlobalArray::GlobalArray(int,int,int[],int[],char*,int[],int[],char)
   * @param[in] p_handle  processor group handle
   */
  GlobalArray(int type, int ndim, int dims[], int width[], char *arrayname, 
	      int block[], int maps[], PGroup* p_handle, char ghosts);

  /**
   * @copydoc GlobalArray::GlobalArray(int,int,int[],int[],char*,int[],int[],char)
   */
  GlobalArray(int type, int ndim, int64_t dims[], int64_t width[],
          char *arrayname, int64_t block[], int64_t maps[], char ghosts);

  /**
   * @copydoc GlobalArray::GlobalArray(int,int,int[],int[],char*,int[],int[],char)
   * @param[in] p_handle  processor group handle
   */
  GlobalArray(int type, int ndim, int64_t dims[], int64_t width[],
          char *arrayname, int64_t block[], int64_t maps[], PGroup* p_handle,
          char ghosts);
  
  /**
   * Creates a new array by applying all the properties of another existing 
   * array.
   *
   * This is a collective operation. 
   *
   * @param[in] arrayname a character string
   * @param[in] g_a       integer handle for reference array
   */
  GlobalArray(const GlobalArray &g_a, char *arrayname); 
  
  /**
   * Creates a new array by applying all the properties of another existing 
   * array.
   *
   * This is a collective operation. 
   *
   * @param[in] g_a integer handle for reference array
   */
  GlobalArray(const GlobalArray &g_a);

  /**
   * Creates a new array with no existing attributes.
   *
   * @note All attributes must subsequently be set using the "set" methods.
   *
   * This is a collective operation. 
   */
  GlobalArray();
  
  /** Destructor */
  ~GlobalArray();

  /* access the data */

  /** @return the array handle */
  int handle() const { return mHandle; }
  
  /* Global Array operations */
  
  /** 
   * Combines data from local array buffer with data in the global array 
   * section.
   *
   * @note The local array is assumed to be have the same number of dimensions
   * as the global array. 
   *
   * global array section (lo[],hi[]) += *alpha * buffer
   *
   * This is a one-sided and atomic operation.  
   *
   * @param[in] lo    [ndim] array of starting indices for array section
   * @param[in] hi    [ndim] array of ending indices for array section
   * @param[in] buf   pointer to the local buffer array
   * @param[in] ld    [ndim-1] array specifying leading
   *                  dimensions/strides/extents for buffer array
   * @param[in] alpha scale factor (double/DoubleComplex/long *)
   */
  void acc(int lo[], int hi[], void *buf, int ld[], void *alpha) const;

  /**
   * @copydoc GlobalArray::acc(int[],int[],void*,int[],void*)const
   */
  void acc(int64_t lo[], int64_t hi[], void *buf, int64_t ld[], void *alpha) const;
  
  /**
   * Provides access to the specified patch of a global array.
   *
   * Returns array of leading dimensions ld and a pointer to the first element 
   * in the patch. This routine allows to access directly, in place 
   * elements in the local section of a global array. It useful for 
   * writing new GA operations. A call to ga_access normally follows a 
   * previous call to ga_distribution that returns coordinates of the 
   * patch associated with a processor. You need to make sure that the 
   * coordinates of the patch are valid (test values returned from 
   * ga_distribution). 
   *
   * Each call to ga_access has to be followed by a call to either 
   * ga_release or ga_release_update. You can access in this fashion only 
   * local data. Since the data is shared with other processes, you need 
   * to consider issues of mutual exclusion. This operation is local. 
   * 
   * @param[in]  lo  [ndim] array of starting indices for array section
   * @param[in]  hi  [ndim] array of ending indices for array section
   * @param[out] ptr points to location of first element in patch
   * @param[out] ld  [ndim-1] leading dimensions for the pacth elements
   */
  void access(int lo[], int hi[], void *ptr, int ld[]) const;

  /**
   * @copydoc GlobalArray::access(int[],int[],void*,int[])const
   */
  void access(int64_t lo[], int64_t hi[], void *ptr, int64_t ld[]) const;

  /**
   * Provides access to the specified block of a global array that is using
   * simple block-cyclic data distribution. Returns array of leading
   * dimensions ld and a pointer to the first element in the patch. This
   * routine allows user to access directly, in-place * elements in the
   * local section of a global array. It useful for writing new GA
   * operations. A call to ga_access normally follows a previous call to
   * ga_distribution that returns coordinates of the patch associated with
   * a processor. You need to make sure that the coordinates of the patch
   * are valid (test values returned from * ga_distribution). 
   *
   * Each call to ga_access_block has to be followed by a call to either 
   * ga_release_block or ga_release_block_update. You can access in this
   * fashion only local data. Since the data is shared with other processes,
   * you need to consider issues of mutual exclusion. This operation is
   * local. 
   * 
   * @param[in]  idx index of block
   * @param[out] ptr points to location of first element in patch
   * @param[out] ld  [ndim-1] leading dimensions for the pacth elements
   */
  void accessBlock(int idx, void *ptr, int ld[]) const;

  /**
   * @copydoc GlobalArray::accessBlock(int,void*,int[])const
   */
  void accessBlock(int64_t idx, void *ptr, int64_t ld[]) const;

  /**
   * Provides access to the specified block of a global array that is using
   * SCALAPACK type block-cyclic data distribution. Returns array of leading
   * dimensions ld and a pointer to the first element in the patch. This
   * routine allows user to access directly, in-place * elements in the
   * local section of a global array. It useful for writing new GA
   * operations. A call to ga_access_block normally follows a previous call to
   * ga_distribution that returns coordinates of the patch associated with
   * a processor. You need to make sure that the coordinates of the patch
   * are valid (test values returned from * ga_distribution). 
   *
   * Each call to ga_access_block_grid has to be followed by a call to either 
   * ga_release_block_grid or ga_release_block_grid_update. You can access in
   * this fashion only local data. Since the data is shared with other
   * processes, you need to consider issues of mutual exclusion. This
   * operation is local. 
   * 
   * @param[in]  index [ndim] indices of block in processor grid
   * @param[out] ptr   points to location of first element in patch
   * @param[out] ld    [ndim-1] leading dimensions for the pacth elements
   */
  void accessBlockGrid(int index[], void *ptr, int ld[]) const;

  /**
   * @copydoc GlobalArray::accessBlockGrid(int[],void*,int[])const
   */
  void accessBlockGrid(int64_t index[], void *ptr, int64_t ld[]) const;

  /**
   * Provides access to the local data of a global array that is using
   * either the simple or SCALAPACK type block-cyclic data distribution.
   * Returns the length of the local data block and a pointer to the first
   * element. This routine allows user to access directly, in-place
   * elements in the local section of a global array. It useful for writing
   * new GA operations.
   *
   * Each call to ga_access_segment has to be followed by a call to either 
   * ga_release_segment or ga_release_segmentupdate. You can access in
   * this fashion only local data. Since the data is shared with other
   * processes, you need to consider issues of mutual exclusion. This
   * operation is local. 
   * 
   * @param[in]  index processor ID
   * @param[out] ptr   points to location of first element
   * @param[out] len   length of locally held data
   */
  void accessBlockSegment(int index, void *ptr, int *len) const;

  /**
   * @copydoc GlobalArray::accessBlockSegment(int,void*,int*)const
   */
  void accessBlockSegment(int index, void *ptr, int64_t *len) const;

  /**
   * Provides access to the local patch of the  global array. Returns 
   * leading dimension ld and and pointer for the data.  This routine 
   * will provide access to the ghost cell data residing on each processor. 
   * Calls to accessGhosts should normally follow a call to  
   * distribution that returns coordinates of the visible data patch 
   * associated with a processor. You need to make sure that the coordinates 
   * of the patch are valid (test values returned from distribution). 
   *    
   * You can only access local data. 
   * This is a local operation. 
   * 
   * @param[out] dims [ndim] array of dimensions of local patch,
   *                  including ghost cells
   * @param[out] ptr  returns an index corresponding to the origin the global
   *                  array patch held locally on the processor
   * @param[out] ld   [ndim-1] physical dimensions of the local array patch,
   *                  including ghost cells
   */
  void accessGhosts(int dims[], void *ptr, int ld[]) const;

  /**
   * @copydoc GlobalArray::accessGhosts(int[],void*,int[])const
   */
  void accessGhosts(int64_t dims[], void *ptr, int64_t ld[]) const;
  
  /**
   * This function can be used to return a pointer to any data element 
   * in the locally held portion of the global array and can be used to 
   * directly access ghost cell data. The array subscript refers to the 
   * local index of the  element relative to the origin of the local 
   * patch (which is assumed to be indexed by (0,0,...)). 
   *
   * This is a  local operation. 
   *
   * @param[out] ptr       index pointing to location of element
   *                       indexed by subscript[]
   * @param[in]  subscript [ndim] array of integers that index desired element
   * @param[out] ld        [ndim-1] array of strides for local data patch.
   *                       These include ghost cell widths.
   */
  void accessGhostElement(void *ptr, int subscript[], int ld[]) const;

  /**
   * @copydoc GlobalArray::accessGhostElement(void*,int[],int[])const
   */
  void accessGhostElement(void *ptr, int64_t subscript[], int64_t ld[]) const;
  
  /**
   * The arrays are aded together elemet-wise:
   * [for example: g_c.add(...,g_a, .., g_b);]
   * c = alpha * a + beta * b
   * The result c may replace one of he input arrays(a/b).
   * This is a collective operation.
   *
   * @param[in] alpha scale factor
   * @param[in] g_a   array
   * @param[in] beta  scale factor
   * @param[in] g_b   array
   */
  void add(void *alpha, const GlobalArray * g_a,
           void *beta,  const GlobalArray * g_b) const;
  
  /**
   * Patches of arrays (which must have the same number of elements) are 
   * added together element-wise. 
   * c[ ][ ] = alpha * a[ ][ ] + beta * b[ ][ ]. 
   *
   * This is a collective operation. 
   *
   * @param[in] alpha scale factor
   * @param[in] g_a   global array
   * @param[in] alo   patch of g_a
   * @param[in] ahi   patch of g_a
   * @param[in] beta  scale factor
   * @param[in] g_b   global array
   * @param[in] blo   patch of g_b
   * @param[in] bhi   patch of g_b
   * @param[in] clo   patch of this GlobalArray
   * @param[in] chi   patch of this GlobalArray
   */
  void addPatch(void *alpha, const GlobalArray * g_a, int alo[], int ahi[],
		 void *beta,  const GlobalArray * g_b, int blo[], int bhi[],
		 int clo[], int chi[]) const;

  /**
   * @copydoc GlobalArray::addPatch(void*,const GlobalArray*,int[],int[],void*,const GlobalArray*,int[],int[],int[],int[])const
   */
  void addPatch(
          void *alpha, const GlobalArray * g_a, int64_t alo[], int64_t ahi[],
		  void *beta,  const GlobalArray * g_b, int64_t blo[], int64_t bhi[],
		  int64_t clo[], int64_t chi[]) const;
 
  /**
   * Allocate internal memory etc. to create a global array
   *
   * @return TODO
   */
  int allocate() const;

  /**
   * This function can be used to preallocate internal buffers that are used by
   * the gather, scatter and scatter accumulate calls. This avoids repeated
   * memory allocations in these calls that can reduce performance. The value of
   * nelems should be set to the maximum number of elements that will be moved
   * in any single call.
   *
   * This is a  local operation. 
   *
   * @param[in]  nelems    The maximum number of elements that will be moved in
   *                       any gather, scatter, scatter-accumulate call
   */
  void allocGatscatBuf(int nelems) const;

  /**
   * Check that the global array handle g_a is valid ... if not call 
   * ga_error with the string provided and some more info. 
   *
   * This operation is local. 
   *
   * @param[in] string message
   */
  void checkHandle(char* string) const;
  
  /**   
   * Compares distributions of two global arrays.
   *
   * This is a collective operation. 
   *
   * @param[in] g_a GlobalArray to compare
   *
   * @return 0 if distributions are identical and 1 when they are not. 
   */
  int compareDistr(const GlobalArray *g_a) const;

  /**   
   * Copies elements in array represented by g_a into the array 
   * represented by g_b [say for example: g_b.copy(g_a);].
   * The arrays must be the same type, shape, and identically aligned. 
   *
   * This is a collective operation.
   *
   * @param[in] g_a GlobalArray to copy
   */
  void copy(const GlobalArray *g_a) const; 
  
  /**
   * Copies elements in a patch of one array (ga) into another one (say for 
   * example:gb.copyPatch(...,ga,....); ).
   *
   * The patches of arrays may be of different shapes but must have the same
   * number of elements. Patches must be nonoverlapping (if gb=ga). 
   *
   * trans = 'N' or 'n' means that the transpose operator should not be 
   * applied. trans = 'T' or 't' means that transpose operator should be 
   * applied. This is a collective operation. 
   *
   * @param[in] trans see above
   * @param[in] ga    global array
   * @param[in] alo   ga patch coordinates
   * @param[in] ahi   ga patch coordinates
   * @param[in] blo   this GlobalArray's patch coordinates
   * @param[in] bhi   this GlobalArray's patch coordinates
   */
  void copyPatch(char trans, const GlobalArray* ga, int alo[], int ahi[], 
		 int blo[], int bhi[]) const;

  /**
   * @copydoc GlobalArray::copyPatch(char,const GlobalArray*,int[],int[],int[],int[])const
   */
  void copyPatch(
          char trans, const GlobalArray* ga, int64_t alo[], int64_t ahi[], 
		  int64_t blo[], int64_t bhi[]) const;
  
  /**
   * Computes element-wise dot product of the two arrays which must be of
   * the same types and same number of elements.
   * return value = SUM_ij a(i,j)*b(i,j)
   *
   * This is a collective operation. 
   *
   * @param[in] g_a GlobalArray operand
   */
  double ddot(const GlobalArray * g_a) const; 
  
  /**
   * Computes the element-wise dot product of the two (possibly transposed) 
   * patches which must be of the same type and have the same number of 
   * elements. 
   *
   * @param[in] ta  transpose flags
   * @param[in] alo g_a patch coordinates
   * @param[in] ahi g_a patch coordinates
   * @param[in] g_a global array
   * @param[in] tb  transpose flags
   * @param[in] blo g_b patch coordinates
   * @param[in] bhi g_b patch coordinates
   */
  double ddotPatch(char ta, int alo[], int ahi[], const GlobalArray * g_a, 
		   char tb, int blo[], int bhi[]) const;

  /**
   * @copydoc GlobalArray::ddotPatch(char,int[],int[],const GlobalArray*,char,int[],int[]const
   */
  double ddotPatch(
          char ta, int64_t alo[], int64_t ahi[], const GlobalArray * g_a, 
		  char tb, int64_t blo[], int64_t bhi[]) const;
  
  /**
   * Deallocates the array and frees any associated resources.
   */
  void destroy();
  
  /**
   * Performs one of the matrix-matrix operations: 
   * [say: g_c.dgemm(..., g_a, g_b,..);]
   *
   *     C := alpha*op( A )*op( B ) + beta*C, \n 
   * where op( X ) is one of \n 
   *     op( X ) = X   or   op( X ) = X', \n 
   * alpha and beta are scalars, and A, B and C are matrices, with op( A ) 
   * an m by k matrix, op( B ) a k by n matrix and C an m by n matrix. 
   * On entry, transa specifies the form of op( A ) to be used in the 
   * matrix multiplication as follows:\n  
   *         ta = 'N' or 'n', op( A ) = A.  \n 
   *         ta = 'T' or 't', op( A ) = A'. \n 
   *
   * This is a collective operation. 
   *
   * @param[in] ta    transpose operators
   * @param[in] tb    transpose operators
   * @param[in] m     number of rows of op(A) and of matrix C
   * @param[in] n     number of columns of op(B) and of matrix C
   * @param[in] k     number of columns of op(A) and rows of matrix op(B)
   * @param[in] alpha scale factors
   * @param[in] g_a   input arrays
   * @param[in] g_b   input arrays
   * @param[in] beta  scale factors
   */
  void dgemm(char ta, char tb, int m, int n, int k, double alpha,  
	     const GlobalArray *g_a, const GlobalArray *g_b,double beta) const;
  /**
   * @copydoc GlobalArray::dgemm(char,char,int,int,int,double,const GlobalArray*,const GlobalArray*,double)const
   */
  void dgemm(char ta, char tb, int64_t m, int64_t n, int64_t k, double alpha,  
	     const GlobalArray *g_a, const GlobalArray *g_b,double beta) const;
  
  /**
   * Solve the generalized eigen-value problem returning all eigen-vectors 
   * and values in ascending order. The input matrices are not overwritten 
   * or destroyed. 
   *
   * This is a collective operation. 
   *
   * @param[in]  g_s  Metric
   * @param[out] g_v  Global matrix to return evecs
   * @param[out] eval Local array to return evals
   * 
   */
  void diag(const GlobalArray *g_s, GlobalArray *g_v, void *eval) const;

  /**
   * Solve the generalized eigen-value problem returning all eigen-vectors 
   * and values in ascending order. Recommended for REPEATED calls if g_s 
   * is unchanged. Values of the control flag: 
   * 
   *          value       action/purpose 
   * 
   *           0          indicates first call to the eigensolver
   * 
   *          >0          consecutive calls (reuses factored g_s) 
   *
   *          <0          only erases factorized g_s; g_v and eval unchanged 
   *                      (should be called after previous use if another 
   *                      eigenproblem, i.e., different g_a and g_s, is to 
   *                      be solved) 
   *
   * The input matrices are not destroyed. 
   *
   * This is a collective operation. 
   *
   * @param[in]  control Control flag
   * @param[in]  g_s     Metric
   * @param[out] g_v     Global matrix to return evecs
   * @param[out] eval    Local array to return evals
   */
  void diagReuse(int control, const GlobalArray *g_s, GlobalArray *g_v, 
		 void *eval) const;
  
  /**
   * Solve the standard (non-generalized) eigenvalue problem returning 
   * all eigenvectors and values in the ascending order. The input matrix 
   * is neither overwritten nor destroyed. 
   *
   * This is a collective operation. 
   *
   * @param[out] g_v  Global matrix to return evecs
   * @param[out] eval Local array to return evals
   */
  void diagStd(GlobalArray *g_v, void *eval) const;

  /**
   * TODO
   */
  void diagSeq(const GlobalArray * g_s, const GlobalArray * g_v, 
	       void *eval) const;
  
  /**
   * TODO
   */
  void diagStdSeq(const GlobalArray * g_v, void *eval) const;
  
  /** 
   * If no array elements are owned by process 'me', the range is returned
   * as lo[]=-1 and hi[]=-2 for all dimensions.
   *
   * The operation is local.
   *
   * @param[in] me process number
   * @param[in] lo [ndim] array of starting indices for array section
   * @param[in] hi [ndim] array of ending indices for array section
   */
  void distribution(int me, int* lo, int* hi) const;
      
  /**
   * @copydoc GlobalArray::distribution(int,int*,int*)const
   */
  void distribution(int me, int64_t* lo, int64_t* hi) const;

  /**
   * TODO
   */
  float fdot(const GlobalArray * g_a) const;

  /**
   * TODO
   */
  float fdotPatch(
          char t_a, int alo[], int ahi[], const GlobalArray * g_b, 
		  char t_b, int blo[], int bhi[]) const;
  /**
   * @copydoc GlobalArray::fdotPatch(char,int[],int[],const GlobalArray*,char,int[],int[])const
   */
  float fdotPatch(
          char t_a, int64_t alo[], int64_t ahi[], const GlobalArray * g_b, 
		  char t_b, int64_t blo[], int64_t bhi[]) const;

  /**
   * Assign a single value to all elements in the array.
   *
   * This is a collective operation. 
   *
   * @param[in] value pointer to the value of appropriate type 
   *            (double/DoubleComplex/long) that matches array type.
   */
  void fill(void *value) const;
  
  /**
   * Fill the patch with  value of 'val' 
   *
   * This is a collective operation. 
   *
   * @param[in] lo  patch of this GlobalArray
   * @param[in] hi  patch of this GlobalArray
   * @param[in] val value to fill
   *
   */
  void fillPatch (int lo[], int hi[], void *val) const;

  /**
   * @copydoc GlobalArray::fillPatch(int[],int[],void*)const
   */
  void fillPatch (int64_t lo[], int64_t hi[], void *val) const;
  
  /**
   * This function can be used to free preallocate internal buffers that were
   * set using the allocGatscatBuf call.
   *
   * This is a  local operation. 
   */
  void freeGatscatBuf();

  /** 
   * Gathers array elements from a global array into a local array. 
   * The contents of the input arrays (v, subscrArray) are preserved, 
   * but their contents might be (consistently) shuffled on return. 
   *
   * @code
   * for(k=0; k<= n; k++){
   *     v[k] = a[subsArray[k][0]][subsArray[k][1]][subsArray[k][2]]...;    
   * }
   * @endcode
   *
   * This is a one-sided operation.  
   *
   * @param[in] n         number of elements
   * @param[in] v         [n] array containing values
   * @param[in] subsarray [n][ndim] array of subscripts for each element
   */
  void gather(void *v, int * subsarray[], int n) const;

  /**
   * @copydoc GlobalArray::gather(void*,int*[],int)const
   */
  void gather(void *v, int64_t * subsarray[], int64_t n) const;

  /**
   * Copies data from global array section to the local array buffer. The 
   * local array is assumed to be have the same number of dimensions as the 
   * global array. Any detected inconsitencies/errors in the input arguments
   * are fatal. 
   * 
   * Example: For ga_get operation transfering data from the [10:14,0:4] 
   * section of 2-dimensional 15x10 global array into local buffer 5x10 
   * array we have: lo={10,0}, hi={14,4}, ld={10}  
   *
   * One-side operation. 
   *
   * @param[in]  lo  [ndim] array of starting indices for global array section
   * @param[in]  hi  [ndim] array of ending indices for global array section
   * @param[out] buf pointer to the local buffer array where the data goes
   * @param[in]  ld  [ndim-1] array specifying leading
   *                 dimensions/strides/extents for buffer array
   */
  void get(int lo[], int hi[], void *buf, int ld[]) const;

  /**
   * @copydoc GlobalArray::get(int[],int[],void*,int[])const
   */
  void get(int64_t lo[], int64_t hi[], void *buf, int64_t ld[]) const;

  /**
   * The function retrieves the number of blocks along each coordinate dimension
   * and the dimensions of the individual blocks for a global array with a
   * block-cyclic data distribution.
   *
   * This is a local operation.
   *
   * @param[out] num_blocks [ndim] array containing number of blocks along each
   *                        coordinate direction
   * @param[out] block_dims [ndim] array containing block dimensions
   */
  void getBlockInfo(int num_blocks[], int block_dims[]);

  /**
   * This function returns 1 if the global array has some dimensions for 
   * which the ghost cell width is greater than zero, it returns 0 otherwise. 
   *
   * This is a local operation. 
   *
   * @return 1 if this GlobalArray has some dimensions for which teh ghost
   *         cell width is greater than zero; 0 otherwise
   */
  int hasGhosts() const;
  
  /**
   * Computes element-wise dot product of the two arrays which must be of
   * the same types and same number of elements.
   *
   * This is a collective operation. 
   *
   * @param[in] g_a GlobalArray
   * @return value = SUM_ij a(i,j)*b(i,j)
   */
  int idot(const GlobalArray * g_a) const; 

  /**
   * Computes the element-wise dot product of the two (possibly transposed) 
   * patches which must be of the same type and have the same number of 
   * elements. 
   *
   * @param[in] ta  transpose flags
   * @param[in] alo g_a patch coordinates
   * @param[in] ahi g_a patch coordinates
   * @param[in] g_a global array
   * @param[in] tb  transpose flags
   * @param[in] blo this GlobalArray's patch coordinates
   * @param[in] bhi this GlobalArray's patch coordinates
   */
  long idotPatch(
          char ta, int alo[], int ahi[], const GlobalArray * g_a, 
		  char tb, int blo[], int bhi[]) const;

  /**
   * @copydoc GlobalArray::idotPatch(char,int[],int[],const GlobalArray*,char,int[],int[])const
   */
  long idotPatch(
          char ta, int64_t alo[], int64_t ahi[], const GlobalArray * g_a, 
		  char tb, int64_t blo[], int64_t bhi[]) const;

  
  /** 
   * Returns data type and dimensions of the array. 
   *
   * This operation is local.   
   *
   * @param[out] type data type
   * @param[out] ndim number of dimensions
   * @param[out] dims array of dimensions
   */
  void inquire(int *type, int *ndim, int dims[]) const;

  /**
   * @copydoc GlobalArray::inquire(int*,int*,int[])const
   */
  void inquire(int *type, int *ndim, int64_t dims[]) const;

  /** 
   * Returns the name of an array represented by the handle g_a. 
   *
   * This operation is local. 
   *
   * @return copy of the name of this GlobalArray
   */
  char* inquireName() const;

  /**
   * Computes element-wise dot product of the two arrays which must be of
   * the same types and same number of elements.
   *          
   *
   * This is a collective operation. 
   *
   * @param[in] g_a array handle
   *
   * @return value = SUM_ij a(i,j)*b(i,j)
   */
   long ldot(const GlobalArray * g_a) const; 

  /**
   * Computes the element-wise dot product of the two (possibly transposed) 
   * patches which must be of the same type and have the same number of 
   * elements. 
   *
   * @param[in] ta  transpose flags
   * @param[in] alo g_a patch coordinates
   * @param[in] ahi g_a patch coordinates
   * @param[in] g_a global array
   * @param[in] tb  transpose flags
   * @param[in] blo this GlobalArray's patch coordinates
   * @param[in] bhi this GlobalArray's patch coordinates
   */
  long ldotPatch(
          char ta, int alo[], int ahi[], const GlobalArray * g_a, 
		  char tb, int blo[], int bhi[]) const;

  /**
   * @copydoc GlobalArray::ldotPatch(char,int[],int[],const GlobalArray*,char,int[],int[])const
   */
  long ldotPatch(
          char ta, int64_t alo[], int64_t ahi[], const GlobalArray * g_a, 
		    char tb, int64_t blo[], int64_t bhi[]) const;

  /**
   * Solves a system of linear equations 
   * 
   *            A * X = B 
   *
   * using the Cholesky factorization of an NxN double precision symmetric 
   * positive definite matrix A (epresented by handle g_a). On successful 
   * exit B will contain the solution X. 
   *
   * This is a collective operation. 
   *
   * @param[in] g_a coefficient matrix
   *
   * @return = 0 : successful exit\n
   *         > 0 : the leading minor of this order is not positive 
   *               definite and the factorization could not be completed 
   */
  int lltSolve(const GlobalArray * g_a) const;

  /**
   * Return in owner the GA compute process id that 'owns' the data. If any 
   * element of subscript[] is out of bounds "-1" is returned.
   *
   * This operation is local. 
   *
   * @param[in] subscript [ndim] element subscript
   *
   * @return ID of compute process which owns the data
   */
  int locate(int subscript[]) const;

  /**
   * @copydoc GlobalArray::locate(int[])const
   */
  int locate(int64_t subscript[]) const;
  
  /**
   * Return the list of the GA processes id that 'own' the data. Parts of the 
   * specified patch might be actually 'owned' by several processes. If lo/hi 
   * are out of bounds "0" is returned, otherwise return value is equal to the 
   * number of processes that hold the data. This operation is local.
   * 
   *   map[i][0:ndim-1]       - lo[i]
   * 
   *   map[i][ndim:2*ndim-1]  - hi[i]
   * 
   *   procs[i]               - processor id that owns data in patch 
   *                            lo[i]:hi[i] 
   * 
   * @param[in]  lo    [ndim] array of starting indices for array section
   * @param[in]  hi    [ndim] array of ending indices for array section
   * @param[out] map   [][2*ndim] array with mapping information
   * @param[out] procs [nproc] list of processes that own a part of selection
   *
   * @return 0 if lo/hi are out of bounds, otherwise the number of processes
   *         holding data
   */
  int locateRegion(int lo[], int hi[], int map[], int procs[]) const;

  /**
   * @copydoc GlobalArray::locateRegion(int[],int[],int[],int[])const
   */
  int locateRegion(int64_t lo[], int64_t hi[], int64_t map[], int procs[]) const;

  /**
   * Solve the system of linear equations op(A)X = B based on the LU 
   * factorization. 
   * 
   * op(A) = A or A' depending on the parameter trans: 
   * 
   * trans = 'N' or 'n' means that the transpose operator should not 
   *          be applied. 
   * 
   * trans = 'T' or 't' means that the transpose operator should be applied. 
   * 
   * Matrix A is a general real matrix. Matrix B contains possibly multiple 
   * rhs vectors. The array associated with the handle g_b is overwritten 
   * by the solution matrix X. 
   * This is a collective operation. 
   *
   * @param[in] trans transpose or not transpose
   * @param[in] g_a   coefficient matrix
   */
  void luSolve(char trans, const GlobalArray * g_a) const;
  
  /**
   * ga_matmul_patch is a patch version of ga_dgemm: 
   * 
   *      C[cilo:cihi,cjlo:cjhi] := alpha* AA[ailo:aihi,ajlo:ajhi] *
   *                                BB[bilo:bihi,bjlo:bjhi] ) + 
   *                                beta*C[cilo:cihi,cjlo:cjhi],
   * 
   * where AA = op(A), BB = op(B), and op( X ) is one of 
   *      op( X ) = X   or   op( X ) = X',
   * 
   * Valid values for transpose arguments: 'n', 'N', 't', 'T'. It works 
   * for both double and DoubleComplex data tape. 
   * This is a collective operation. 
   *
   * @param[in] g_a    global array
   * @param[in] g_b    global array
   * @param[in] ailo   patch of g_a
   * @param[in] aihi   patch of g_a
   * @param[in] ajlo   patch of g_a
   * @param[in] ajhi   patch of g_a
   * @param[in] bilo   patch of g_b
   * @param[in] bihi   patch of g_b
   * @param[in] bjlo   patch of g_b
   * @param[in] bjhi   patch of g_b
   * @param[in] cilo   patch of g_c
   * @param[in] cihi   patch of g_c
   * @param[in] cjlo   patch of g_c
   * @param[in] cjhi   patch of g_c
   * @param[in] alpha  scale factors
   * @param[in] beta   scale factors
   * @param[in] transa transpose operators
   * @param[in] transb transpose operators
   */
  void matmulPatch(char transa, char transb, void* alpha, void *beta,
		   const GlobalArray *g_a, 
		   int ailo, int aihi, int ajlo, int ajhi,
		   const GlobalArray *g_b, 
		   int bilo, int bihi, int bjlo, int bjhi,
		   int cilo, int cihi, int cjlo, int cjhi) const;

  /**
   * @copydoc GlobalArray::matmulPatch(char,char,void*,void*,const GlobalArray*,int,int,int,int,const GlobalArray*,int,int,int,int,int,int,int,int)const
   */
  void matmulPatch(char transa, char transb, void* alpha, void *beta,
		   const GlobalArray *g_a, 
		   int64_t ailo, int64_t aihi, int64_t ajlo, int64_t ajhi,
		   const GlobalArray *g_b, 
		   int64_t bilo, int64_t bihi, int64_t bjlo, int64_t bjhi,
		   int64_t cilo, int64_t cihi, int64_t cjlo, int64_t cjhi) const;

  /**
   * nga_matmul_patch is a n-dimensional patch version of ga_dgemm: 
   * 
   *      C[clo[]:chi[]] := alpha* AA[alo[]:ahi[]] *
   *                               BB[blo[]:bhi[]]) + 
   *                               beta*C[clo[]:chi[]],
   * 
   * where AA = op(A), BB = op(B), and op( X ) is one of 
   *      op( X ) = X   or   op( X ) = X',
   * 
   * Valid values for transpose arguments: 'n', 'N', 't', 'T'. It works 
   * for both double and DoubleComplex data tape. 
   *
   * This is a collective operation. 
   *
   * @param[in] g_a    global array
   * @param[in] g_b    global array
   * @param[in] alo    array of patch of g_a
   * @param[in] ahi    array of patch of g_a
   * @param[in] blo    array of patch of g_b
   * @param[in] bhi    array of patch of g_b
   * @param[in] clo    array of patch of g_c
   * @param[in] chi    array of patch of g_c
   * @param[in] alpha  scale factors
   * @param[in] beta   scale factors
   * @param[in] transa transpose operators
   * @param[in] transb transpose operators
   */
  void matmulPatch(char transa, char transb, void* alpha, void *beta,
		   const GlobalArray *g_a, int *alo, int *ahi,
		   const GlobalArray *g_b, int *blo, int *bhi,
		   int *clo, int *chi) const;
  /**
   * @copydoc GlobalArray::matmulPatch(char,char,void*,void*,const GlobalArray*,int*,int*,const GlobalArray*,int*,int*,int*,int*)const
   */
  void matmulPatch(char transa, char transb, void* alpha, void *beta,
		   const GlobalArray *g_a, int64_t *alo, int64_t *ahi,
		   const GlobalArray *g_b, int64_t *blo, int64_t *bhi,
		   int64_t *clo, int64_t *chi) const;

  /**
   * This function merges all values in a patch of a mirrored array into
   * a patch in another global array g_b.
   *
   * This is a collective operation.
   *
   * @param[in]  alo [ndim] patch indices of mirrored array
   * @param[in]  ahi [ndim] patch indices of mirrored array
   * @param[in]  blo [ndim] patch indices of result array
   * @param[in]  bhi [ndim] patch indices of result array
   * @param[out] g_a global array containing result
   */
  void mergeDistrPatch(int alo[], int ahi[], GlobalArray *g_a,
                       int blo[], int bhi[]);

  /**
   * @copydoc GlobalArray::mergeDistrPatch(int[],int[],GlobalArray*,int[],int[])
   */
  void mergeDistrPatch(int64_t alo[], int64_t ahi[], GlobalArray *g_a,
                       int64_t blo[], int64_t bhi[]);

  /**
   * This function returns 0 if a global array is not mirrored and 1 if it is.
   */
  int isMirrored();

  /**
   * This function adds together all copies of a mirrored array so that all
   * copies are the same.
   *
   * This is a collective operation.
   */
  void mergeMirrored();

  /**
   * Non-blocking accumalate operation. This is function performs an
   * accumulate operation and returns a nblocking handle. Completion of the
   * operation can be forced by calling the nbwait method on the handle.
   *
   * This is a onesided operation.
   *
   * @param[in]  lo       [ndim] patch coordinates of block
   * @param[in]  hi       [ndim] patch coordinates of block
   * @param[in]  buf      local buffer containing data
   * @param[in]  ld       [ndim-1] array of strides for local data
   * @param[in]  alpha    multiplier for data before adding to existing results
   * @param[out] nbhandle nonblocking handle
   */
  void nbAcc(int lo[], int hi[], void *buf, int ld[], void *alpha,
             GANbhdl *nbhandle);

  /**
   * @copydoc GlobalArray::nbAcc(int[],int[],void*,int[],void*,GANbhdl*)
   */
  void nbAcc(int64_t lo[], int64_t hi[], void *buf, int64_t ld[], void *alpha,
             GANbhdl *nbhandle);
  
  /**
   * Non-blocking get operation. This is function gets a data block from a
   * global array, copies it into a local buffer, and returns a nonblocking
   * handle. Completion of the operation can be forced by calling the nbwait
   * method on the handle.
   *
   * This is a onesided operation.
   *
   * @param[in]  lo       [ndim] patch coordinates of block
   * @param[in]  hi       [ndim] patch coordinates of block
   * @param[in]  buf      local buffer to receive data
   * @param[in]  ld       [ndim-1] array of strides for local data
   * @param[out] nbhandle nonblocking handle
   */
  void nbGet(int lo[], int hi[], void *buf, int ld[], GANbhdl *nbhandle);

  /**
   * @copydoc GlobalArray::nbGet(int[],int[],void*,int[],GANbhdl*)
   */
  void nbGet(int64_t lo[], int64_t hi[], void *buf, int64_t ld[], GANbhdl *nbhandle);

  /**
   * Non-blocking update operation for arrays with ghost cells. Ghost cells
   * along the coordinates specified in the mask array are updated with
   * non-blocking get calls. The mask array must contain either 0's or 1's.
   *
   * This is a onesided operation.
   *
   * @param[in]  mask     [ndim] array with flags for directions that are
   *                      to be updated
   * @param[out] nbhandle nonblocking handle
   */
  void nbGetGhostDir(int mask[], GANbhdl *nbhandle);

  /**
   * @copydoc GlobalArray::nbGetGhostDir(int[],GANbhdl*)
   */
  void nbGetGhostDir(int64_t mask[], GANbhdl *nbhandle);
  
  /**
   * Given a distribution of an array represented by the handle g_a, 
   * returns the number of partitions of each array dimension. 
   *
   * This operation is local. 
   *
   * @param[out] nblock [ndim] number of partitions for each dimension
   */
  void nblock(int nblock[]) const;

  /**
   * Non-blocking put operation. This is function puts a data block from a
   * local array, copies it into a global array, and returns a nonblocking
   * handle. Completion of the operation can be forced by calling the nbwait
   * method on the handle.
   *
   * This is a onesided operation.
   *
   * @param[in]  lo       [ndim] patch coordinates of block
   * @param[in]  hi       [ndim] patch coordinates of block
   * @param[in]  buf      local buffer that supplies data
   * @param[in]  ld       [ndim-1] array of strides for local data
   * @param[out] nbhandle nonblocking handle
   */
  void nbPut(int lo[], int hi[], void *buf, int ld[], GANbhdl *nbhandle);

  /**
   * @copydoc GlobalArray::nbPut(int[],int[],void*,int[],GANbhdl*)
   */
  void nbPut(int64_t lo[], int64_t hi[], void *buf, int64_t ld[], GANbhdl *nbhandle);

  /**
   * Returns the number of dimensions in this GlobalArray.
   *
   * This operation is local. 
   *
   * @return number of dimensions aka rank
   */
  int ndim() const;
  
  /**
   * The pack subroutine is designed to compress the values in the source vector
   * g_src into a smaller destination array g_dest based on the values in an
   * integer mask array g_mask. The values lo and hi denote the range of
   * elements that should be compressed and icount is a variable that on output
   * lists the number of values placed in the compressed array. This operation
   * is the complement of the ga_unpack operation. An example is shown below
   *
   *  g_src->pack(g_dest, g_mask, 1, n, icount)
   *  g_mask:   1  0  0  0  0  0  1  0  1  0  0  1  0  0  1  1  0
   *  g_src:    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
   *  g_dest:   1  7  9 12 15 16
   *  icount:   6
   *
   * The calling array is the source array.
   *
   * This is a collective operation.
   *
   * @param[out] g_dest destination array
   * @param[in]  g_mask mask array
   * @param[in]  lo     coordinate interval to pack
   * @param[in]  hi     coordinate interval to pack
   * @param[out] icount number of packed elements
   */
  void pack(const GlobalArray *g_dest, const GlobalArray *g_mask,
            int lo, int hi, int *icount) const;

  /**
   * @copydoc GlobalArray::pack(const GlobalArray*,const GlobalArray*,int,int,int*)const
   */
  void pack(const GlobalArray *g_dest, const GlobalArray *g_mask,
            int64_t lo, int64_t hi, int64_t *icount) const;

  /**
   * This subroutine enumerates the values of an array between elements lo and
   * hi starting with the value istart and incrementing each subsequent value by
   * inc. This operation is only applicable to 1-dimensional arrays. An example
   * of its use is shown below:
   *
   * call g_a->patch_enum(g_a, 1, n, 7, 2)
   * g_a:  7  9 11 13 15 17 19 21 23 ...
   *
   * This is a collective operation.
   *
   * @param[in] lo     coordinate interval to enumerate
   * @param[in] hi     coordinate interval to enumerate
   * @param[in] istart starting value of enumeration
   * @param[in] inc    increment value
   */
   void patchEnum(int lo, int hi, void *istart, void *inc);

  /**
   * @copydoc GlobalArray::patchEnum(int,int,int,int)
   */
   void patchEnum(int64_t lo, int64_t hi, void *start, void *inc);

  /**  
   * Same as nga_acc except the indices can extend beyond the array 
   * boundary/dimensions in which case the library wraps them around. 
   *
   * This is a one-sided and atomic operation.     
   *
   * @param[in] lo    [ndim] array of starting indices for array section
   * @param[in] hi    [ndim] array of ending indices for array section
   * @param[in] buf   pointer to the local buffer array
   * @param[in] ld    [ndim-1] array specifying leading 
   *                  dimensions/strides/extents for buffer array
   * @param[in] alpha double/DoubleComplex/long scale factor 
   */
  void periodicAcc(int lo[], int hi[], void* buf, int ld[], void* alpha) const;

  /**
   * @copydoc GlobalArray::periodicAcc(int[],int[],void*,int[],void*)const
   */
  void periodicAcc(int64_t lo[], int64_t hi[], void* buf, int64_t ld[], void* alpha) const;
  
  /**
   * Same as nga_get except the indices can extend beyond the array 
   * boundary/dimensions in which case the library wraps them around. 
   *
   * This is a one-sided operation.
   *
   * @param[in]  lo  [ndim] array of starting indices for global array section
   * @param[in]  hi  [ndim] array of ending indices for global array section
   * @param[out] buf pointer to the local buffer array where the data goes
   * @param[in]  ld  [ndim-1] array specifying leading
   *                 dimensions/strides/extents for buffer array
   */
  void periodicGet(int lo[], int hi[], void* buf, int ld[]) const;

  /**
   * @copydoc GlobalArray::periodicGet(int[],int[],void*,int[])const
   */
  void periodicGet(int64_t lo[], int64_t hi[], void* buf, int64_t ld[]) const;
  
  /**
   * Same as nga_put except the indices can extend beyond the array 
   * boundary/dimensions in which case the library wraps them around. 
   *
   * This is a one-sided operation.
   *
   * @param[in] lo  [ndim] array of starting indices for global array section
   * @param[in] hi  [ndim] array of ending indices for global array section
   * @param[in] buf pointer to the local buffer array where the data goes
   * @param[in] ld  [ndim-1] array specifying leading
   *                dimensions/strides/extents for buffer array
   */
  void periodicPut(int lo[], int hi[], void* buf, int ld[]) const;

  /**
   * @copydoc GlobalArray::periodicPut(int[],int[],void*,int[])const
   */
  void periodicPut(int64_t lo[], int64_t hi[], void* buf, int64_t ld[]) const;
  
  /** 
   * Prints an entire array to the standard output. 
   *
   * This is a collective operation. 
   */
  void print() const ; 
  
  /** 
   * Prints the array distribution. 
   *
   * This is a collective operation. 
   */
  void printDistribution() const ;

  /** 
   * Prints the array distribution to a file. 
   *
   * This is a collective operation. 
   */
   void printFile(FILE *file) const;

  /**
   * Prints a patch of g_a array to the standard output. If pretty has the 
   * value 0 then output is printed in a dense fashion. If pretty has the 
   * value 1 then output is formatted and rows/columns labeled. 
   
   * This is a collective operation.  
   *
   * @param[in] lo     coordinates of the patch
   * @param[in] hi     coordinates of the patch
   * @param[in] pretty formatting flag
   */
  void printPatch(int* lo, int* hi, int pretty)  const;

  /**
   * @copydoc GlobalArray::printPatch(int*,int*,int)const
   */
  void printPatch(int64_t* lo, int64_t* hi, int pretty)  const;
  
  /**
   * Based on the distribution of an array associated with handle g_a, 
   * determines coordinates of the specified processor in the virtual 
   * processor grid corresponding to the distribution of array g_a. The 
   * numbering starts from 0. The values of -1 means that the processor 
   * doesn't 'own' any section of array represented by g_a. 
   *
   * This operation is local. 
   *
   * @param[in]  proc  process id
   * @param[out] coord [ndim] coordinates in processor grid
   * 
   */
  void procTopology(int proc, int coord[]) const;
  
  /*void procTopology(int proc, int *prow, int *pcol);*/
  
  /**
   * Copies data from local array buffer to the global array section . The 
   * local array is assumed to be have the same number of dimensions as the 
   * global array. Any detected inconsitencies/errors in input arguments are 
   * fatal. This is a one-sided operation. 
   *
   * @param[in] lo  [ndim] array of starting indices for global array section
   * @param[in] hi  [ndim] array of ending indices for global array section
   * @param[in] buf pointer to the local buffer array where the data is
   * @param[in] ld  [ndim-1] array specifying leading
   *                dimensions/strides/extents for buffer array
   * @param[in] buf buffer array
   */
  void put(int lo[], int hi[], void *buf, int ld[]) const;
  
  /**
   * @copydoc GlobalArray::put(int[],int[],void*,int[])const
   */
  void put(int64_t lo[], int64_t hi[], void *buf, int64_t ld[]) const;

  /**
   * Atomically read and increment an element in an integer array. 
   * 
   *      *BEGIN CRITICAL SECTION*
   * 
   *       old_value = a(subscript)
   * 
   *       a(subscript) += inc
   * 
   *      *END CRITICAL SECTION*
   * 
   *       return old_value
   *
   * This is a one-sided and atomic operation. 
   *
   * @param[in] subscript [ndim] subscript array for the referenced element
   * @param[in] inc how much to increment by
   * @return the incremented value
   */
  long readInc(int subscript[], long inc) const;

  /**
   * @copydoc GlobalArray::readInc(int[],long)const
   */
  long readInc(int64_t subscript[], long inc) const;

  /**
   * Releases access to a global array when the data was read only. 
   * Your code should look like: 
   * 
   * @code
   * g_a->distribution(myproc, lo,hi);
   * g_a->access(lo, hi, &ptr, ld);
   * // <operate on the data referenced by ptr> 
   * g_a->release(lo, hi);
   * @endcode
   * 
   * @note see restrictions specified for ga_access. 
   *
   * This operation is local. 
   *
   * @param[in] lo [ndim] array of starting indices for array section
   * @param[in] hi [ndim] array of ending indices for array section
   */
  void release(int lo[], int hi[]) const;

  /**
   * @copydoc GlobalArray::release(int[],int[])const
   */
  void release(int64_t lo[], int64_t hi[]) const;

  /**
   * Releases access to the block of data specified by the integer
   * index when data was accessed as read only. This is only applicable to
   * block-cyclic data distributions created using the simple block-cyclic
   * distribution.
   *
   * This is a local operation.
   *
   * @param[in] index block index
   */
  void releaseBlock(int index) const;

  /**
   * Releases access to the block of data specified by the subscript
   * array when data was accessed as read only. This is only applicable to
   * block-cyclic data distributions created using the SCALAPACK data
   * distribution.
   *
   * This is a local operation.
   *
   * @param[in] index [ndim] indices of block in array
   */      
  void releaseBlockGrid(int index[]) const;
      
  /**
   * Releases access to the block of locally held data for a block-cyclic
   * array, when data was accessed as read-only. This is a local operation.
   *
   * @param[in] proc process ID/rank
   */
  void releaseBlockSegment(int proc) const;
      
  /**
   * Releases access to the data. It must be used if the data was accessed 
   * for writing. NOTE: see restrictions specified for ga_access. 
   *
   * This operation is local. 
   *
   * @param[in] lo [ndim] array of starting indices for array section
   * @param[in] hi [ndim] array of ending indices for array section
   */
  void releaseUpdate(int lo[], int hi[]) const;

  /**
   * @copydoc GlobalArray::releaseUpdate(int[],int[])const
   */
  void releaseUpdate(int64_t lo[], int64_t hi[]) const;

  /**
   * Releases access to the block of data specified by the integer index when
   * data was accessed in read-write mode. This is only applicable to
   * block-cyclic data distributions created using the simple block-cyclic
   * distribution.
   *
   * This is a local operation.
   *
   * @param[in] index block index
   */
  void releaseUpdateBlock(int index) const;
      
  /**
   * Releases access to the block of data specified by the subscript
   * array when data was accessed in read-write mode. This is only applicable
   * to block-cyclic data distributions created using the SCALAPACK data
   * distribution.
   *
   * This is a local operation.
   *
   * @param[in] index [ndim] indices of block in array
   */   
  void releaseUpdateBlockGrid(int index[]) const;

  /**
   * Releases access to the block of locally held data for a block-cyclic
   * array, when data was accessed in read-write mode.
   *
   * This is a local operation.
   *
   * @param[in] proc process ID/rank
   */
  void releaseUpdateBlockSegment(int proc) const;    

  /**
   * Releases access to a global array containing ghost cells when the data was
   * read only. 
   * Your code should look like: 
   * 
   * @code
   * g_a->accessGhosts(dims, &ptr, ld)
   * // <operate on the data referenced by ptr> 
   * g_a->releasGhosts();
   * @endcode
   * 
   * This operation is local. 
   *
   */
  void releaseGhosts() const;

  /**
   * Releases access to a global array containing ghost cells when the data was
   * accessed in read-write mode. 
   * 
   * This operation is local. 
   *
   */
  void releaseUpdateGhosts() const;

  /**
   * Releases access to a global array containing ghost cells when the data was
   * read only. 
   * Your code should look like: 
   * 
   * @code
   * g_a->accessGhostElement(&ptr, subscript, ld)
   * // <operate on the data referenced by ptr> 
   * g_a->releaseGhostElement(subscript);
   * @endcode
   * 
   * This operation is local. 
   * @param[in]  indices of element
   *
   */
  void releaseGhostElement(int subscript[]) const;

  /**
   * @copydoc GlobalArray::releaseGhostElement(int subscript[]) const
   */
  void releaseGhostElement(int64_t subscript[]) const;

  /**
   * Releases access to a global array containing ghost cells when the data was
   * accessed in read-write mode. 
   * 
   * This operation is local. 
   * @param[in]  indices of element
   *
   */
  void releaseUpdateGhostElement(int subscript[]) const;

  /**
   * @copydoc GlobalArray::releaseUpdateGhostElement(int subscript[]) const
   */
  void releaseUpdateGhostElement(int64_t subscript[]) const;

  /** 
   * Scales an array by the constant s. Note that the library is unable 
   * to detect errors when the pointed value is of different type than 
   * the array. 
   *
   * This is a collective operation. 
   *
   * @param[in] value pointer to the value of appropriate type 
   */
  void scale(void *value) const;

  /** 
   * Scale an array by the factor 'val'. 
   *
   * This is a collective operation. 
   *
   * @param[in] lo  patch of g_a
   * @param[in] hi  patch of g_a
   * @param[in] val scale factor
   */
  void scalePatch (int lo[], int hi[], void *val) const;

  /**
   * @copydoc GlobalArray::scalePatch(int[],int[],void*)const
   */
  void scalePatch (int64_t lo[], int64_t hi[], void *val) const;
      
  /**
   * This operation will add successive elements in a source vector g_src
   * and put the results in a destination vector g_dest. The addition will
   * restart based on the values of the integer mask vector g_mask. The scan
   * is performed within the range specified by the integer values lo and
   * hi. Note that this operation can only be applied to 1-dimensional
   * arrays. The excl flag determines whether the sum starts with the value
   * in the source vector corresponding to the location of a 1 in the mask
   * vector (excl=0) or whether the first value is set equal to 0
   * (excl=1). Some examples of this operation are given below.
   *
   * g_src->scanAdd(g_dest, g_mask, 1, n, 0);
   * g_mask:   1  0  0  0  0  0  1  0  1  0  0  1  0  0  1  1  0
   * g_src:    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
   * g_dest:   1  3  6 10 16 21  7 15  9 19 30 12 25 39 15 16 33
   * 
   * g_src->scanAdd(g_dest, g_mask, 1, n, 1);
   * g_mask:   1  0  0  0  0  0  1  0  1  0  0  1  0  0  1  1  0
   * g_src:    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
   * g_dest:   0  1  3  6 10 15  0  7  0  9 19  0 12 25  0  0 16
   * 
   * This is a collective operation.
   *
   * @param[out] g_dest handle for destination array
   * @param[in]  g_mask handle for integer array representing mask
   * @param[in]  lo     low and high values of range on which operation
   *                    is performed
   * @param[in]  hi     low and high values of range on which operation
   *                    is performed
   * @param[in]  excl   value to signify if masked values are included in in add
   */
  void scanAdd(const GlobalArray *g_dest, const GlobalArray *g_mask,
               int lo, int hi, int excl) const;

  /**
   * @copydoc GlobalArray::scanAdd(const GlobalArray*,const GlobalArray*,int,int,int)const
   */
  void scanAdd(const GlobalArray *g_dest, const GlobalArray *g_mask,
               int64_t lo, int64_t hi, int excl) const;

  /**
   * This subroutine does a segmented scan-copy of values in the
   * source array g_src into a destination array g_dest with segments
   * defined by values in the integer mask array g_mask. The scan-copy
   * operation is only applied to the range between the lo and hi
   * indices. This operation is restriced to 1-dimensional arrays. The
   * resulting destination array will consist of segments of consecutive
   * elements with the same value. An example is shown below
   *
   * g_src->scanCopy(g_dest, g_mask, 1, n);
   * g_mask:   1  0  0  0  0  0  1  0  1  0  0  1  0  0  1  1  0
   * g_src:    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
   * g_dest:   1  1  1  1  1  1  7  7  9  9  9 12 12 12 15 16 16
   *
   * This is  a collective operation.
   *
   * @param[out] g_dest handle for destination array
   * @param[in]  g_mask handle for integer array representing mask
   * @param[in]  lo     low and high values of range on which operation
   *                    is performed
   * @param[in]  hi     low and high values of range on which operation
   *                    is performed
   */
  void scanCopy(const GlobalArray *g_dest, const GlobalArray *g_mask,
                int lo, int hi) const;

  /**
   * @copydoc GlobalArray::scanCopy(const GlobalArray*,const GlobalArray*,int,int)const
   */
  void scanCopy(const GlobalArray *g_dest, const GlobalArray *g_mask,
                int64_t lo, int64_t hi) const;
                 
  /** 
   * Scatters array elements into a global array. The contents of the input 
   * arrays (v,subscrArray) are preserved, but their contents might be 
   * (consistently) shuffled on return.    
   *
   * @code
   * for(k=0; k<= n; k++) {
   *    a[subsArray[k][0]][subsArray[k][1]][subsArray[k][2]]... = v[k];    
   * }
   * @endcode
   *
   * This is a one-sided operation.  
   *
   * @param[in] n         number of elements
   * @param[in] v         [n] array containing values
   * @param[in] subsarray [n][ndim] array of subscripts for each element
   */
  void scatter(void *v, int *subsarray[], int n) const;

  /**
   * @copydoc GlobalArray::scatter(void*,int*[],int)const
   */
  void scatter(void *v, int64_t *subsarray[], int64_t n) const;

  /** 
   * Adds element a local array to array elements into a global array after
   * multiplying by alpha.  The contents of the input arrays (v,subscrArray)
   * are preserved, but their contents might be (consistently) shuffled on
   * return.    
   *
   * @code
   * for(k=0; k<= n; k++) {
   *    a[subsArray[k][0]][subsArray[k][1]][subsArray[k][2]]... = v[k];    
   * }
   * @endcode
   *
   * This is a one-sided operation.  
   *
   * @param[in] n         number of elements
   * @param[in] v         [n] array containing values
   * @param[in] subsarray [n][ndim] array of subscripts for each element
   * @param[in] alpha     scale factor
   */
  void scatterAcc(void *v, int *subsarray[], int n, void *alpha) const;

  /**
   * @copydoc GlobalArray::scatterAcc(void*,int*[],int,void*)const
   */
  void scatterAcc(void *v, int64_t *subsarray[], int64_t n, void *alpha) const;

  /**
   * Returns the value and index for an element that is selected by the 
   * specified operator  in a global array corresponding to g_a handle. 
   *
   * This is a collective operation. 
   *
   * @param[in]  op    operator {"min","max"}
   * @param[out] val   address where value should be stored
   * @param[out] index [ndim] array index for the selected element
   */
  void selectElem(char *op, void* val, int index[]) const;

  /**
   * @copydoc GlobalArray::selectElem(char*,void*,int[])const
   */
  void selectElem(char *op, void* val, int64_t index[]) const;

  /**
   * This function can be used to assign a unique character
   * string name to a global array handle that was obtained
   * using the createHandle function.
   *
   * This is a collective operation.
   *
   * @param[in] name array name
   */
  void setArrayName(char *name) const;

  /**
   * This subroutine is used to create a global array with a simple
   * block-cyclic data distribution. The array is broken up into blocks of
   * size dims and each block is numbered sequentially using a column major
   * indexing scheme. The blocks are then assigned in a simple round-robin
   * fashion to processors. This is illustrated in the figure below for an
   * array containing 25 blocks distributed on 4 processors. Blocks at the
   * edge of the array may be smaller than the block size specified in
   * dims. In the example below, blocks 4,9,14,19,20,21,22,23, and 24 might
   * be smaller thatn the remaining blocks. Most global array operations
   * are insensitive to whether or not a block-cyclic data distribution is
   * used, although performance may be slower in some cases if the global
   * array is using a block-cyclic data distribution. Individual data
   * blocks can be accessesed using the block-cyclic access functions.
   *
   * This is a collective operation.
   *
   * @param[in] dims array of block dimensions
   */
  void setBlockCyclic(int dims[]) const;

  /**
   * This subroutine is used to create a global array with a
   * SCALAPACK-type block cyclic data distribution. The user  specifies
   * the dimensions of the processor grid in the array proc_grid. The
   * product of the processor grid dimensions must equal the number of
   * total number of processors  and the number of dimensions in the
   * processor grid must be the same as the number of dimensions in the
   * global array. The data blocks are mapped onto the processor grid
   * in a cyclic manner along each of the processor grid axes. This is
   * illustrated below for an array consisting of 25 data blocks
   * disributed on 6 processors. The 6 processors are configured in a 3
   * by 2 processor grid. Blocks at the edge of the array may be
   * smaller than the block size specified in dims. Most global array
   * operations are insensitive to whether or not a block-cyclic data
   * distribution is used, although performance may be slower in some
   * cases if the global array is using a block-cyclic data
   * distribution. Individual data blocks can be accessesed using the
   * block-cyclic access functions.
   *
   * This is a collective operation.
   *
   * @param[in] dims      array of block dimensions
   * @param[in] proc_grid processor grid dimensions
   */
  void setBlockCyclicProcGrid(int dims[], int proc_grid[]) const;

  /**
   * This function is used to set the chunk array for a global array handle
   * that was obtained using the createHandle function. The chunk array
   * is used to determine the minimum number of array elements assigned to
   * each processor along each coordinate direction.
   *
   * This is a collective operation.
   *
   * @param[in] chunk array of chunk widths
   */
  void setChunk(int chunk[]) const;

  /**
   * @copydoc GlobalArray::setChunk(int[])const
   */
  void setChunk(int64_t chunk[]) const;

  /**
   * This function can be used to set the array dimension, the coordinate
   * dimensions, and the data type assigned to a global array handle obtained
   * using the GA_Create_handle function.
   *
   * This is a collective operation.
   *
   * @param[in] ndim dimension of global array
   * @param[in] dims dimensions of global array
   * @param[in] type data type of global array
   */
  void setData(int ndim, int dims[], int type) const;

  /**
   * @copydoc GlobalArray::setData(int,int[],int)const
   */
  void setData(int ndim, int64_t dims[], int type) const;
      
  /**
   * This function can be used to set the ghost cell widths for a global
   * array handle that was obtained using the createHandle function. The
   * ghosts cells widths indicate how many ghost cells are used to pad the
   * locally held array data along each dimension. The padding can be set
   * independently for each coordinate dimension.
   *
   * This is a collective operation.
   *
   * @param[in] width [ndim] array of ghost cell widths
   */
  void setGhosts(int width[]) const;

  /**
   * @copydoc GlobalArray::setGhosts(int[])const
   */
  void setGhosts(int64_t width[]) const;

  /**
   * This function can be used to partition the array data among the
   * individual processors for a global array handle obtained using the
   * GA_Create_handle function.
   *
   * The distribution is specified as a Cartesian product of distributions
   * for each dimension. For example, the following figure demonstrates
   * distribution of a 2-dimensional array 8x10 on 6 (or more)
   * processors. nblock(2)={3, 2}, the size of mapc array is s=5 and array
   * mapc contains the following elements mapc={1, 3, 7, 1, 6}. The
   * distribution is nonuniform because, P1 and P4 get 20 elements each and
   * processors P0,P2,P3, and P5 only 10 elements each.
   *
   * The array width() is used to control the width of the ghost cell
   * boundary around the visible data on each processor. The local data of
   * the global array residing on each processor will have a layer width(n)
   * ghosts cells wide on either side of the visible data along the dimension
   * n.
   *
   * This is a collective operation.
   *
   * @param[in] mapc   [s] starting index for each block; the size
   *                       s is the sum of all elements of the array nblock
   * @param[in] nblock [ndim] number of blocks that each dimension is
   *                   divided into
   */
  void setIrregDistr(int mapc[], int nblock[]) const;

  /**
   * @copydoc GlobalArray::setIrregDistr(int mapc[], int nblock[]) const
   */
  void setIrregDistr(int64_t mapc[], int64_t nblock[]) const;

  /**
   * This function can be used to set the processor configuration assigned to
   * a global array handle that was obtained using the
   * createHandle function. It can be used to create mirrored arrays by
   * using the mirrored array processor configuration in this function
   * call. It can also be used to create an array on a processor group by
   * using a processor group handle in this call.
   *
   * This is a collective operation.
   *
   * @param[in] pHandle processor group handle
   */
  void setPGroup(PGroup *pHandle) const;

  /**
   * This function is used to restrict the number of processors in a global
   * array that actually contain data. It can also be used to rearrange the
   * layout of data on a processor from the default distribution. Only the
   * processes listed in list[] will actually contain data, the remaining
   * processes will be able to see the data in the global array but they will
   * not contain any of the global array data locally.
   *
   * @param[in] list   list of processors that should contain data
   * @param[in] nprocs number of processors in list
   *
   */
  void setRestricted(int list[], int nprocs) const;

  /**
   * This function is used to restrict the number of processors in a global
   * array that actually contain data. Only the processors in the range
   * [lo_proc:hi_proc] (inclusive) will actually contain data, the remaining
   * processes will be able to see the data in the global array but they will
   * not contain any of the global array data locally.
   *
   * @param[in] lo_proc low end of processor range
   * @param[in] hi_proc high end of processor range
   */
  void setRestrictedRange(int lo_proc, int hi_proc) const;
      
  /**
   * Performs one of the matrix-matrix operations: 
   *
   *     C := alpha*op( A )*op( B ) + beta*C,
   * where op( X ) is one of 
   *     op( X ) = X   or   op( X ) = X',
   * alpha and beta are scalars, and A, B and C are matrices, with op( A ) 
   * an m by k matrix, op( B ) a k by n matrix and C an m by n matrix. 
   * On entry, transa specifies the form of op( A ) to be used in the 
   * matrix multiplication as follows: 
   *
   *         ta = 'N' or 'n', op( A ) = A. 
   *
   *         ta = 'T' or 't', op( A ) = A'. 
   *
   * This is a collective operation. 
   *
   * @param[in] g_a   handles to input arrays
   * @param[in] g_b   handles to input arrays
   * @param[in] ta    transpose operators
   * @param[in] tb    transpose operators
   * @param[in] m     number of rows of op(A) and of matrix C
   * @param[in] n     number of columns of op(B) and of matrix C
   * @param[in] k     number of columns of op(A) and rows of matrix op(B)
   * @param[in] alpha scale factors
   * @param[in] beta  scale factors
   *
   */
  void sgemm(char ta, char tb, int m, int n, int k, float alpha,  
	     const GlobalArray *g_a, const GlobalArray *g_b, float beta) const;

  /**
   * @copydoc GlobalArray::sgemm(char,char,int,int,int,float,const GlobalArray*,const GlobalArray*,float)const
   */
  void sgemm(char ta, char tb, int64_t m, int64_t n, int64_t k, float alpha,  
	     const GlobalArray *g_a, const GlobalArray *g_b, float beta) const;
  
  /**
   * Solves a system of linear equations 
   *            A * X = B 
   * It first will call the Cholesky factorization routine and, if 
   * sucessfully, will solve the system with the Cholesky solver. If 
   * Cholesky will be not be able to factorize A, then it will call the 
   * LU factorization routine and will solve the system with forward/backward 
   * substitution. On exit B will contain the solution X. 
   *
   * This is a collective operation. 
   *
   * @param[in] g_a coefficient matrix
   *
   * @return = 0 : Cholesky factoriztion was succesful\n
   *         > 0 : the leading minor of this order 
   *               is not positive definite, Cholesky factorization 
   *               could not be completed and LU factoriztion was used 
   */
  int solve(const GlobalArray * g_a) const;

  /**
   * It computes the inverse of a double precision using the Cholesky 
   * factorization of a NxN double precision symmetric positive definite 
   * matrix A stored in the global array represented by g_a. On successful 
   * exit, A will contain the inverse. 
   *
   * This is a collective operation. 
   *
   * @return    = 0 : successful exit\n
   *            > 0 : the leading minor of this order is not positive 
   *                  definite and the factorization could not be completed\n
   *            < 0 : it returns the index i of the (i,i) 
   *                  element of the factor L/U that is zero and 
   *                  the inverse could not be computed
   */
  int spdInvert() const;

  /**
   * This operation is the same as "acc", except that the values
   * corresponding to dimension n in buf are accumulated to every skip[n]
   * values of the global array.
   *
   * This is a one-sided operation.
   *
   * @param[in] lo    [ndim] array of starting indices for glob array section
   * @param[in] hi    [ndim] array of ending indices for global array section
   * @param[in] skip  [ndim] array of strides for each dimension
   * @param[in] buf   pointer to local buffer array where data goes
   * @param[in] ld    [ndim-1] rray specifying leading
   *                  dimensions/strides/extents for buffer array
   * @param[in] alpha double/DoublComplex/long scale factor
   */
  void stridedAcc(int lo[], int hi[], int skip[], void*buf, int ld[], void *alpha) const;
  
  /**
   * @copydoc GlobalArray::stridedAcc(int[],int[],int[],void*,int[],void*)const
   */
  void stridedAcc(int64_t lo[], int64_t hi[], int64_t skip[], void*buf, int64_t ld[], void *alpha) const;

  /**
   * This operation is the same as "get", except that the values
   * corresponding to dimension n in buf are accumulated to every skip[n]
   * values of the global array.
   *
   * This is a one-sided operation.
   *
   * @param[in]  lo   [ndim] array of starting indices for glob array section
   * @param[in]  hi   [ndim] array of ending indices for global array section
   * @param[in]  skip [ndim] array of strides for each dimension
   * @param[out] buf  pointer to local buffer array where data goes
   * @param[in]  ld   [ndim-1] array specifying leading
   *                  dimensions/strides/extents for buffer array
   */
  void stridedGet(int lo[], int hi[], int skip[], void*buf, int ld[]) const;

  /**
   * @copydoc GlobalArray::stridedGet(int[],int[],int[],void*,int[])const
   */
  void stridedGet(int64_t lo[], int64_t hi[], int64_t skip[], void*buf, int64_t ld[]) const;
      
  /**
   * This operation is the same as "put", except that the values
   * corresponding to dimension n in buf are accumulated to every skip[n]
   * values of the global array.
   *
   * This is a one-sided operation.
   *
   * @param[in] lo    [ndim] array of starting indices for glob array section
   * @param[in] hi    [ndim] array of ending indices for global array section
   * @param[in] skip  [ndim] array of strides for each dimension
   * @param[in] buf   pointer to local buffer array where data goes
   * @param[in] ld    [ndim-1] array specifying leading
   *                  dimensions/strides/extents for buffer array
   */
  void stridedPut(int lo[], int hi[], int skip[], void*buf, int ld[]) const;
  
  /**
   * "long" interface for stridedPut
   */
  void stridedPut(int64_t lo[], int64_t hi[], int64_t skip[], void*buf, int64_t ld[]) const;

  /**
   * Prints info about allocated arrays.
   *
   * @param[in] verbose If true print distribution info
   */
  void summarize(int verbose) const;
      
  /** 
   * Symmmetrizes matrix A with handle A:=.5 * (A+A').
   *
   * This is a collective operation 
   */
  void symmetrize() const;

  /**
   * This function returns the total number of blocks contained in a global
   * array with a block-cyclic data distribution.
   *
   * This is a local operation.
   *
   * @return number of blocks contained in this block-cyclic distribution
   */
  int totalBlocks() const;
      
  /**
   * Transposes a matrix: B = A', where A and B are represented by 
   * handles g_a and g_b [say, g_b.transpose(g_a);].
   *
   * This is a collective operation.
   *
   * @param[in] g_a GlobalArray to transpose and assign to this GlobalArray
   */
  void transpose(const GlobalArray * g_a) const;
  
  /**
   * The unpack subroutine is designed to expand the values in the source
   * vector g_src into a larger destination array g_dest based on the values
   * in an integer mask array g_mask. The values lo and hi denote the range
   * of elements that should be compressed and icount is a variable that on
   * output lists the number of values placed in the uncompressed array. This
   * operation is the complement of the pack operation. An example is
   * shown below
   *
   * g_src->unpack(g_dest, g_mask, 1, n, &icount);
   * g_src:    1  7  9 12 15 16
   * g_mask:   1  0  0  0  0  0  1  0  1  0  0  1  0  0  1  1  0
   * g_dest:   1  0  0  0  0  0  7  0  9  0  0 12  0  0 15 16  0
   * icount:   6
   *
   * This is a collective operation.
   *
   * @param[out] g_dest handle for destination array
   * @param[in]  g_mask handle for integer array representing mask
   * @param[in]  lo     low value of range on which operation is performed
   * @param[in]  hi     high value of range on which operation is performed
   * @param[out] icount number of values in uncompressed array
   */
  void unpack(GlobalArray *g_dest, GlobalArray *g_mask, int lo, int hi,
              int *icount) const;
  /**
   * @copydoc GlobalArray::unpack(GlobalArray*,GlobalArray*,int,int,int*)const
   */
  void unpack(GlobalArray *g_dest, GlobalArray *g_mask,
            int64_t lo, int64_t hi, int64_t *icount) const;
      
  /**
   * This call updates the ghost cell regions on each processor with the 
   * corresponding neighbor data from other processors. The operation assumes 
   * that all data is wrapped around using periodic boundary data so that 
   * ghost cell data that goes beyound an array boundary is wrapped around to
   * the other end of the array. The updateGhosts call contains two   
   * sync   calls before and after the actual update operation. For some 
   * applications these calls may be unecessary, if so they can be removed 
   * using the maskSync subroutine. 
   *
   * This is a collective operation. 
   */
  void updateGhosts() const;

  /**
   * This operation is similar to the standard updateGhosts operation except
   * that it returns a non-blocking handle after initiating the call. Completion
   * of the operation can be guaranteed by call call the NbWait function on the
   * handle. Data in the local buffers is then ready for use.
   *
   * This is a collective operation. 
   */
  void updateGhostsNb(GANbhdl *nbhandle) const;
 
  /**
   * This function can be used to update the ghost cells along individual 
   * directions. It is designed for algorithms that can overlap updates 
   * with computation. The variable dimension indicates which coordinate 
   * direction is to be updated (e.g. dimension = 1 would correspond to the 
   * y axis in a two or three dimensional system), the variable idir can take
   * the values +/-1 and indicates whether the side that is to be updated lies 
   * in the positive or negative direction, and cflag indicates whether or not
   * the corners on the side being updated are to be included in the update. 
   * The following calls would be equivalent to a call to   updateGhosts 
   * for a 2-dimensional system: 
   *
   * status = g_a->updateGhostDir(0,-1,1);\n 
   * status = g_a->updateGhostDir(0,1,1);\n 
   * status = g_a->updateGhostDir(1,-1,0);\n 
   * status = g_a->updateGhostDir(1,1,0);\n 
   *
   * The variable cflag is set equal to 1 (or non-zero) in the first two
   * calls so that the corner ghost cells are update, it is set equal to 0 in
   * the second two calls to avoid redundant updates of the corners. Note
   * that updating the ghosts cells using several independent calls to the
   * nga_update_ghost_dir functions is generally not as efficient as using
   * updateGhosts  unless the individual calls can be effectively overlapped
   * with computation.
   *
   * This is a  collective operation. 
   *
   * @param[in] dimension array dimension that is to be updated
   * @param[in] idir      direction of update (+/- 1)
   * @param[in] cflag     flag (0/1) to include corners in update
   */
  int updateGhostDir(int dimension, int idir, int cflag) const;

  /**
   * This operation is designed to extract ghost cell data from a global array
   * and copy it to a local array. If the request can be satisfied using
   * completely local data, then a local copy will be used. Otherwise, the
   * method calls periodicGet. The request can be satisfied locally if
   * lo is greater than or equal to the lower bound of data held on the
   * processor minus the ghost cell width and hi is less than or equal to the
   * upper bound of data held on the processor plus the ghost cell width. Cell
   * indices using the global address space should be used for lo and hi. These
   * may exceed the global array dimensions.
   *
   * @param[in]  lo  [ndim] array of starting indices for global array section
   * @param[in]  hi  [ndim] array of ending indices for global array section
   * @param[out] buf pointer to the local buffer array where the data goes
   * @param[in]  ld  [ndim-1] array specifying leading
   *                 dimensions/strides/extents for buffer array
   */
  void getGhostBlock(int lo[], int hi[], void *buf, int ld[]) const;

  /**
   * @copydoc GlobalArray::getGhostBlock(int[],int[],void*,int[])const
   */
  void getGhostBlock(int64_t lo[], int64_t hi[], void *buf, int64_t ld[]) const;

  /**
   * Computes element-wise dot product of the two arrays which must be of
   * the same types and same number of elements.
   *
   * This is a collective operation. 
   *
   * @param[in] g_a array handle
   *
   * @return value = SUM_ij a(i,j)*b(i,j)
   */
  DoubleComplex zdot(const GlobalArray * g_a) const; 

  /**
   * Computes the element-wise dot product of the two (possibly transposed) 
   * patches which must be of the same type and have the same number of 
   * elements. 
   *
   * @param[in] ta  transpose flags
   * @param[in] alo g_a patch coordinates
   * @param[in] ahi g_a patch coordinates
   * @param[in] g_a global array
   * @param[in] tb  transpose flags
   * @param[in] blo g_b patch coordinates
   * @param[in] bhi g_b patch coordinates
   * @return value
   */
  DoubleComplex zdotPatch(char ta, int alo[], int ahi[], 
			  const GlobalArray * g_a, char tb, int blo[], 
			  int bhi[]) const;

  /**
   * @copydoc GlobalArray::zdotPatch(char,int[],int[],const GlobalArray*,char,int[],int[])const
   */
  DoubleComplex zdotPatch(char ta, int64_t alo[], int64_t ahi[], 
			  const GlobalArray * g_a, char tb, int64_t blo[], 
			  int64_t bhi[]) const;
  
  /** 
   * Sets value of all elements in the array to zero. 
   *
   * This is a collective operation. 
   */
  void zero() const;
  
  /**
   * Set all the elements in the patch to zero. 
   * This is a collective operation. 
   *
   * @param[in] lo
   * @param[in] hi
   */
  void zeroPatch (int lo[], int hi[]) const;

  /**
   * @copydoc GlobalArray::zeroPatch(int[],int[])const
   */
  void zeroPatch (int64_t lo[], int64_t hi[]) const;
  
  /**
   * Performs one of the matrix-matrix operations: 
   *     C := alpha*op( A )*op( B ) + beta*C,
   * where op( X ) is one of 
   *     op( X ) = X   or   op( X ) = X',
   * alpha and beta are scalars, and A, B and C are matrices, with op( A ) 
   * an m by k matrix, op( B ) a k by n matrix and C an m by n matrix. 
   * On entry, transa specifies the form of op( A ) to be used in the 
   * matrix multiplication as follows: 
   *
   *         ta = 'N' or 'n', op( A ) = A. 
   *
   *         ta = 'T' or 't', op( A ) = A'. *
   *
   * This is a collective operation. 
   *
   * @param[in] g_a   handles to input arrays
   * @param[in] g_b   handles to input arrays
   * @param[in] ta    transpose operators
   * @param[in] tb    transpose operators
   * @param[in] m     number of rows of op(A) and of matrix C
   * @param[in] n     number of columns of op(B) and of matrix C
   * @param[in] k     number of columns of op(A) and rows of matrix op(B)
   * @param[in] alpha scale factors
   * @param[in] beta  scale factors
   */
  void zgemm(char ta, char tb, int m, int n, int k, DoubleComplex alpha,  
	     const GlobalArray *g_a, const GlobalArray *g_b, 
	     DoubleComplex beta) const;

  /**
   * @copydoc GlobalArray::zgemm(char,char,int,int,int,DoubleComplex,const GlobalArray*,const GlobalArray*,DoubleComplex)const
   */
  void zgemm(char ta, char tb, int64_t m, int64_t n, int64_t k, DoubleComplex alpha,  
	     const GlobalArray *g_a, const GlobalArray *g_b, 
	     DoubleComplex beta) const;
  
  /* New additional functionalities from Limin. */

  /**
   * Take element-wise absolute value of the array. 
   *
   * This is a collective operation. 
   */
   void absValue() const; 

  /**
   * Take element-wise absolute value of the patch. 
   *
   * This is a collective operation.
   * 
   * @param[in] lo patch coordinates
   * @param[in] hi patch coordinates
   */
   void absValuePatch(int *lo, int *hi) const;

  /**
   * @copydoc GlobalArray::absValuePatch(int*,int*)const
   */
   void absValuePatch(int64_t *lo, int64_t *hi) const;

  /**
   * Add the constant pointed by alpha to each element of the array. 
   *
   * This is a collective operation. 
   *
   * @param[in] alpha double/complex/int/long/float
   */
   void addConstant(void* alpha) const;
 
  /**
   * Add the constant pointed by alpha to each element of the patch. 
   *
   * This is a collective operation. 
   *
   * @param[in] lo    g_a patch coordinates
   * @param[in] hi    g_a patch coordinates
   * @param[in] alpha double/complex/int/long/float
   */
   void addConstantPatch(int *lo, int *hi, void *alpha) const;

  /**
   * @copydoc GlobalArray::addConstantPatch(int*,int*,void*)const
   */
   void addConstantPatch(int64_t *lo, int64_t *hi, void *alpha) const;
  
  /**
   * Take element-wise reciprocal of the array. 
   *
   * This is a collective operation. 
   */
   void recip() const;
  
  /**
   * Take element-wise reciprocal of the patch. 
   *
   * This is a collective operation. 
   *
   * @param[in] lo patch coordinates
   * @param[in] hi patch coordinates
   */
   void recipPatch(int *lo, int *hi) const;

  /**
   * @copydoc GlobalArray::recipPatch(int*,int*)const
   */
   void recipPatch(int64_t *lo, int64_t *hi) const;
  
  /**
   * Computes the element-wise product of the two arrays 
   * which must be of the same types and same number of 
   * elements. For two-dimensional arrays, 
   *
   *            c(i, j)  = a(i,j)*b(i,j) 
   *
   * The result (c) may replace one of the input arrays (a/b). 
   * This is a collective operation. 
   *
   * @param[in] g_a GlobalArray
   * @param[in] g_b GlobalArray
   */
   void elemMultiply(const GlobalArray * g_a, const GlobalArray * g_b) const;
  
  /**
   * Computes the element-wise product of the two patches 
   * which must be of the same types and same number of 
   * elements. For two-dimensional arrays, 
   * 
   *             c(i, j)  = a(i,j)*b(i,j) 
   * 
   * The result (c) may replace one of the input arrays (a/b). 
   *
   * This is a collective operation.
   *
   * @param[in] g_a global array
   * @param[in] g_b global array
   * @param[in] alo g_a patch coordinates
   * @param[in] ahi g_a patch coordinates
   * @param[in] blo g_b patch coordinates
   * @param[in] bhi g_b patch coordinates
   * @param[in] clo g_c patch coordinates
   * @param[in] chi g_c patch coordinates
   */ 
   void elemMultiplyPatch(const GlobalArray * g_a,int *alo,int *ahi,
				 const GlobalArray * g_b,int *blo,int *bhi,
				 int *clo,int *chi) const;
  /**
   * @copydoc GlobalArray::elemMultiplyPatch(const GlobalArray*,int*,int*,const GlobalArray*,int*,int*,int*,int*)const
   */
   void elemMultiplyPatch(const GlobalArray * g_a,int64_t *alo,int64_t *ahi,
				 const GlobalArray * g_b,int64_t *blo,int64_t *bhi,
				 int64_t *clo,int64_t *chi) const;

  /**
   * Computes the element-wise quotient of the two arrays 
   * which must be of the same types and same number of 
   * elements. For two-dimensional arrays, 
   * 
   *             c(i, j)  = a(i,j)/b(i,j) 
   * 
   * The result (c) may replace one of the input arrays (a/b). If one of 
   * the elements of array g_b is zero, the quotient for the element of g_c 
   * will be set to GA_NEGATIVE_INFINITY. 
   *
   * This is a collective operation. 
   *
   * @param[in] g_a global array
   * @param[in] g_b global array
   */
   void elemDivide(const GlobalArray * g_a, const GlobalArray * g_b) const;
  
  /**
   * Computes the element-wise quotient of the two patches 
   * which must be of the same types and same number of 
   * elements. For two-dimensional arrays, 
   *
   *            c(i, j)  = a(i,j)/b(i,j) 
   *
   * The result (c) may replace one of the input arrays (a/b). 
   *
   * This is a collective operation. 
   *
   * @param[in] g_a global array
   * @param[in] g_b global array
   * @param[in] alo g_a patch coordinates
   * @param[in] ahi g_a patch coordinates
   * @param[in] blo g_b patch coordinates
   * @param[in] bhi g_b patch coordinates
   * @param[in] clo g_c patch coordinates
   * @param[in] chi g_c patch coordinates
   */ 
   void elemDividePatch(const GlobalArray * g_a,int *alo,int *ahi,
			       const GlobalArray * g_b,int *blo,int *bhi,
			       int *clo,int *chi) const;
  /**
   * @copydoc GlobalArray::elemDividePatch(const GlobalArray*,int*,int*,const GlobalArray*,int*,int*,int*,int*)const
   */
   void elemDividePatch(const GlobalArray * g_a,int64_t *alo,int64_t *ahi,
			       const GlobalArray * g_b,int64_t *blo,int64_t *bhi,
			       int64_t *clo,int64_t *chi) const;

  /**
   * Computes the element-wise maximum of the two arrays 
   * which must be of the same types and same number of 
   * elements. For two dimensional arrays, 
   *
   *                c(i, j)  = max{a(i,j), b(i,j)} 
   *
   * The result (c) may replace one of the input arrays (a/b). 
   *
   * This is a collective operation. 
   *
   * @param[in] g_a global array
   * @param[in] g_b global array
   */
   void elemMaximum(const GlobalArray * g_a, const GlobalArray * g_b) const;
  
  /**
   * Computes the element-wise maximum of the two patches 
   * which must be of the same types and same number of 
   * elements. For two-dimensional of noncomplex arrays, 
   *
   *             c(i, j)  = max{a(i,j), b(i,j)} 
   *
   * If the data type is complex, then 
   *     c(i, j).real = max{ |a(i,j)|, |b(i,j)|} while c(i,j).image = 0. 
   *
   * The result (c) may replace one of the input arrays (a/b). 
   *
   * This is a collective operation. 
   *
   * @param[in] g_a global array
   * @param[in] g_b global array
   * @param[in] alo g_a patch coordinates
   * @param[in] ahi g_a patch coordinates
   * @param[in] blo g_b patch coordinates
   * @param[in] bhi g_b patch coordinates
   * @param[in] clo g_c patch coordinates
   * @param[in] chi g_c patch coordinates
   */
    void elemMaximumPatch(const GlobalArray * g_a,int *alo,int *ahi,
				 const GlobalArray * g_b,int *blo,int *bhi,
				 int *clo,int *chi) const;
   /**
    * @copydoc GlobalArray::elemMaximumPatch(const GlobalArray*,int*,int*,const GlobalArray*,int*,int*,int*,int*)const
    */
    void elemMaximumPatch(const GlobalArray * g_a,int64_t *alo,int64_t *ahi,
				 const GlobalArray * g_b,int64_t *blo,int64_t *bhi,
				 int64_t *clo,int64_t *chi) const;

   /**
    * Computes the element-wise minimum of the two arrays 
    * which must be of the same types and same number of 
    * elements. For two dimensional arrays, 
    * 
    *             c(i, j)  = min{a(i,j), b(i,j)} 
    * 
    * The result (c) may replace one of the input arrays (a/b). 
    *
    * This is a collective operation. 
    *
    * @param[in] g_a global array
    * @param[in] g_b global array
    */
    void elemMinimum(const GlobalArray * g_a, const GlobalArray * g_b) const;
   
   /**
    * Computes the element-wise minimum of the two patches 
    * which must be of the same types and same number of 
    * elements. For two-dimensional of noncomplex arrays, 
    * 
    *             c(i, j)  = min{a(i,j), b(i,j)} 
    * 
    * If the data type is complex, then 
    *             c(i, j).real = min{ |a(i,j)|, |b(i,j)|} while c(i,j).image = 0. 
    * 
    * The result (c) may replace one of the input arrays (a/b). 
    *
    * This is a collective operation. 
    *
    * @param[in] g_a global array
    * @param[in] g_b global array
    * @param[in] alo g_a patch coordinates
    * @param[in] ahi g_a patch coordinates
    * @param[in] blo g_b patch coordinates
    * @param[in] bhi g_b patch coordinates
    * @param[in] clo g_c patch coordinates
    * @param[in] chi g_c patch coordinates
    */
    void elemMinimumPatch(const GlobalArray * g_a,int *alo,int *ahi,
				 const GlobalArray * g_b,int *blo,int *bhi,
				 int *clo,int *chi) const;

    /**
     * @copydoc GlobalArray::elemMinimumPatch(const GlobalArray*,int*,int*,const GlobalArray*,int*,int*,int*,int*)const
     */
    void elemMinimumPatch(const GlobalArray * g_a, int64_t *alo, int64_t *ahi,
				 const GlobalArray * g_b, int64_t *blo, int64_t *bhi,
				 int64_t *clo, int64_t *chi) const;
  
   /** 
    * Calculates the largest multiple of a vector g_b that can be added 
    * to this vector g_a while keeping each element of this vector 
    * nonnegative. 
    *
    * This is a collective operation. 
    *
    * @param[in]  g_b  global array where g_b is the step direction.
    * @param[out] step the maximum step
    */  
    void stepMax(const GlobalArray * g_b, double *step) const;
   
    /**
     * @copydoc GlobalArray::stepMax(const GlobalArray*,double*)const
     * @param[in] alo g_a patch coordinates
     * @param[in] ahi g_a patch coordinates
     * @param[in] blo g_b patch coordinates
     * @param[in] bhi g_b patch coordinates
     */
   void stepMaxPatch(int *alo, int *ahi, 
			    const GlobalArray * g_b, int *blo, int *bhi, 
			    double *step) const;
   /**
    * @copydoc GlobalArray::stepMaxPatch(int*,int*,const GlobalArray*,int*,int*,double*)const
    */
   void stepMaxPatch(int64_t *alo, int64_t *ahi, 
			    const GlobalArray * g_b, int64_t *blo, int64_t *bhi, 
			    double *step) const;
  
  /** Matrix Operations */
  
  /** 
   * Adds this constant to the diagonal elements of the matrix. 
   *
   * This is a collective operation. 
   *
   * @param[in] c double/complex/int/long/float constant to add
   */
   void shiftDiagonal(void *c) const;
  
  /**
   * Sets the diagonal elements of this matrix g_a with the elements of the 
   * vector g_v.
   *
   * This is a collective operation. 
   *
   * @param[in] g_v global array
   */
  void setDiagonal(const GlobalArray * g_v) const;
  
  /**
   * Sets the diagonal elements of this matrix g_a with zeros. 
   *
   * This is a collective operation. 
   */
   void zeroDiagonal() const;
  
  /**
   * Adds the elements of the vector g_v to the diagonal of this matrix g_a. 
   *
   * This is a collective operation. 
   *
   * @param[in] g_v global array
   */
   void addDiagonal(const GlobalArray * g_v) const;
  
  /**
   * Inserts the diagonal elements of this matrix g_a into the vector g_v. 
   *
   * This is a collective operation. 
   *
   * @param[in] g_a global array
   */
   void getDiagonal(const GlobalArray * g_a) const;
  
  /**
   * Scales the rows of this matrix g_a using the vector g_v. 
   *
   * This is a collective operation.  
   *
   * @param[in] g_v global array
   */
   void scaleRows(const GlobalArray * g_v) const;
  
  /** 
   * Scales the columns of this matrix g_a using the vector g_v. 
   *
   * This is a collective operation. 
   *
   * @param[in] g_v global array
   */
   void scaleCols(const GlobalArray * g_v) const;
  
  /**
   * Computes the 1-norm of the matrix or vector g_a. 
   *
   * This is a collective operation. 
   *
   * @param[in] nm matrix/vector 1-norm value 
   */
   void norm1(double *nm) const;
  
  /**
   * Computes the 1-norm of the matrix or vector g_a. 
   *
   * This is a collective operation. 
   *
   * @param[in] nm - matrix/vector 1-norm value 
   */
   void normInfinity(double *nm) const;
  
  /**
   * Computes the componentwise Median of three arrays g_a, g_b, and g_c, and 
   * stores the result in this array g_m.  The result (m) may replace one of
   * the input arrays (a/b/c).
   *
   * This is a collective operation. 
   *
   * @param[in] g_a global array
   * @param[in] g_b global array
   * @param[in] g_c global array
   */
   void median(const GlobalArray * g_a, const GlobalArray * g_b, 
		      const GlobalArray * g_c) const;
  
  /**
   * Computes the componentwise Median of three patches g_a, g_b, and g_c, and 
   * stores the result in this patch g_m.  The result (m) may replace one of
   * the input patches (a/b/c).
   *
   * This is a collective operation. 
   *
   * @param[in] g_a global array
   * @param[in] g_b global array
   * @param[in] g_c global array
   * @param[in] alo g_a patch coordinates
   * @param[in] ahi g_a patch coordinates
   * @param[in] blo g_b patch coordinates
   * @param[in] bhi g_b patch coordinates
   * @param[in] clo g_c patch coordinates
   * @param[in] chi g_c patch coordinates
   * @param[in] mlo g_m patch coordinates
   * @param[in] mhi g_m patch coordinates
   */
   void medianPatch(const GlobalArray * g_a, int *alo, int *ahi, 
			       const GlobalArray * g_b, int *blo, int *bhi, 
			       const GlobalArray * g_c, int *clo, int *chi, 
			       int *mlo, int *mhi) const;
  /**
   * @copydoc GlobalArray::medianPatch(const GlobalArray*,int*,int*,const GlobalArray*,int*,int*,const GlobalArray*,int*,int*,int*,int*)const
   */
   void medianPatch(const GlobalArray * g_a, int64_t *alo, int64_t *ahi, 
			       const GlobalArray * g_b, int64_t *blo, int64_t *bhi, 
			       const GlobalArray * g_c, int64_t *clo, int64_t *chi, 
			       int64_t *mlo, int64_t *mhi) const;

   GlobalArray& operator=(const GlobalArray &g_a);
   int operator==(const GlobalArray &g_a) const;
   int operator!=(const GlobalArray &g_a) const;
  
 private:
  int mHandle; /**<< g_a handle */
};

}

#endif /* _GLOBALARRAY_H */
