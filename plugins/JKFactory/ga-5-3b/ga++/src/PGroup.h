#ifndef _PGROUP_H
#define _PGROUP_H

namespace GA {

/**
 * PGroup class description TODO.
 */
class PGroup {

public:
  /**
   * This constructor creates a processor group.
   *
   * It must be invoked by all processors in the current default processor
   * group. The list of processors use the indexing scheme of the default
   * processor group. If the default processor group is the world group, then
   * these indices are the usual processor indices. This function returns a
   * process group handler that can be used to reference this group by other
   * functions.
   *
   * This is a collective operation on the default processor group.
   *
   * @param[in] plist [size] list of processor IDs in group
   * @param[in] size  number of processors in group
   *
   */
  PGroup(int *plist, int size);

  /**
   * PGroup destructor.
   */
  ~PGroup();

  /* access the data */

  /** @return the array handle */
  int handle() const { return mPHandle; }

  /* Process Group Operations - static */
  
  /**
   * This function will return a handle to the default processor group,
   * which can then be used to create a global array using one of the
   * create_*_config or setPGroup calls.
   *
   * This is a local operation.
   *
   * @return the new PGroup
   */
  static PGroup* getDefault();
  
  /**
   * This function will return a handle to the mirrored processor group,
   * which can then be used to create a global array using one of the
   * GA create_*_config or setPgroup calls.
   *
   * This is a local operation.
   *
   * @return the new PGroup
   */
  static PGroup* getMirror();

  /**
   * This function will return a handle to the world processor group,
   * which can then be used to create a global array using one of the
   * GA create_*_config or GA_Set_pgroup calls.
   *
   * This is a local operation.
   *
   * @return the new PGroup
   */
  static PGroup* getWorld();

  /**
   * This function can be used to reset the default processor group on a
   * collection of processors.
   *
   * All processors in the group referenced by p_handle must make a call to
   * this function. Any standard global array call that is made after
   * resetting the default processor group will be restricted to processors in
   * that group. Global arrays that are created after resetting the default
   * processor group will only be defined on that group and global operations
   * such as sync or igop will be restricted to processors in that group. The
   * pgroupSetDefault call can be used to rapidly convert large applications,
   * written with GA, into routines that run on processor groups.
   *
   * The default processor group can be overridden by using GA calls that
   * require an explicit group handle as one of the arguments.
   *
   * This is a collective operation on the group represented by the handle
   * p_handle.
   *
   * @param[in] p_handle processor group handle
   */
  static void setDefault(PGroup* p_handle);
  
  /* Process Group Operations */
  
  /**
   * Broadcast data from processor specified by root to all other
   * processors in this processor group.
   *
   * The length of the message in bytes is specified by lenbuf. The initial and
   * broadcasted data can be found in the buffer specified by the pointer
   * buf.
   *
   * This is a collective operation on the processor group.
   *
   * @param[in,out] buf    pointer to buffer containing data
   * @param[in]     lenbuf length of data (in bytes)
   * @param[in]     root   processor sending message
   */
  void brdcst(void* buf, int lenbuf, int root);
     
  /**
   * The pgroup gop 'sums' all elements in buf[n] across all processors in the
   * group using the commutative operation specified by the character string op.
   *
   * buf[n] is a double precision array present on each processor in the
   * processor group. The result is broadcast to all processor in this group.
   * Allowed strings are "+", "*", "max", "min", "absmax", "absmin". The use of
   * lowerecase for operators is necessary.
   *
   * This is a collective operation on the processor group.
   *
   * @param[in,out] buf buffer containing data
   * @param[in]     n   number of elements in x
   * @param[in]     op  operation to be performed
   */
  void gop(double *buf, int n, char* op);

  /**
   * The gop 'sums' all elements in buf[n] across all processors in the group
   * using the commutative operation specified by the character string op.
   *
   * buf[n] is an integer(int) array present on each processor in the
   * processor group.   The result is broadcast to all processor in this group.
   * Allowed strings are "+", "*", "max", "min", "absmax", "absmin". The use of
   * lowerecase for operators is necessary.
   *
   * This is a collective operation on the processor group.
   *
   * @param[in,out] buf buffer containing data
   * @param[in]     n   number of elements in x
   * @param[in]     op  operation to be performed
   */
  void gop(int *buf, int n, char* op);

  /**
   * The gop 'sums' all elements in buf[n] across all processors in the group
   * using the commutative operation specified by the character string op.
   *
   * buf[n] is an integer(long) array present on each processor in the
   * processor group.   The result is broadcast to all processor in this group.
   * Allowed strings are "+", "*", "max", "min", "absmax", "absmin". The use of
   * lowerecase for operators is necessary.
   *
   * This is a collective operation on the processor group.
   *
   * @param[in,out] buf buffer containing data
   * @param[in]     n   number of elements in x
   * @param[in]     op  operation to be performed
   */
  void gop(long *buf, int n, char* op);
  
  /**
   * The gop 'sums' all elements in buf[n] across all processors in the group
   * using the commutative operation specified by the character string op.
   *
   * buf[n] is a float array present on each processor in the
   * processor group.   The result is broadcast to all processor in this group.
   * Allowed strings are "+", "*", "max", "min", "absmax", "absmin". The use of
   * lowerecase for operators is necessary.
   *
   * This is a collective operation on the processor group.
   *
   * @param[in,out] buf buffer containing data
   * @param[in] n       number of elements in x
   * @param[in] op      operation to be performed
   */
  void gop(float *buf, int n, char* op);

  /**
   * This function returns the relative index of the processor in the
   * processor group specified by p_handle.
   *
   * This index will generally differ from the absolute processor index
   * returned by GA_Nodeid if the processor group is not the world group.
   *
   * This is a local operation.
   *
   * @return relative index of the processor in the group given by p_handle
   */
  int nodeid();

  /**
   * This function returns the number of processors contained in the
   * group specified by p_handle.
   *
   * This is a local local operation.
   *
   * @return number of processors contained in the group given by p_handle
   */
  int nodes();

  /**
   * This operation executes a synchronization group across the
   * processors in the processor group specified by p_handle.
   *
   * Nodes outside this group are unaffected.
   *
   * This is a collective operation on the processor group.
   */
  void sync();
   
private:
  int mPHandle; /*<< process group handle */
  static PGroup *pgMirror;
  static PGroup *pgDefault;
  static PGroup *pgWorld;
  PGroup(void);
};

}

#endif /* _PGROUP_H */
