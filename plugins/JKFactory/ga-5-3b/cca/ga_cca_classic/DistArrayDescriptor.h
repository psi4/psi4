#ifndef DistArrayDescriptor_h_seen
#define DistArrayDescriptor_h_seen

#include "DistArrayTemplate.h"

/** This is the public interface (abstract base class) for defining
    and querying distributed array descriptors.

    This interface is an experimental first implementation of an
    interface which has been under discussion in the CCA Forum's
    Scientific Data Components Working Group.  This implementation
    does not conform exactly to the interface developed by the Working
    Group, and is intended to be exploratory rather than normative.

    $Id: DistArrayDescriptor.h,v 1.1 2003-08-01 00:41:53 manoj Exp $

    This interface is intended to support the creation of distributed
    data objects structured like dense multi-dimensional distributed
    arrays.  The object is constructed from size and type information,
    pointers to the local data on each process, and a mapping onto a
    distribution template (see DistTemplateCreate).  In this way, it
    should be possible to describe any array-structured distributed
    data object sufficiently to allow the construction of parallel
    communications schedules and other data movement-related
    operations.  This interface is primarily intended to accomodate
    existing applications and perhaps low-level use within new
    components.  We strongly advise developers of new code to use
    "first class" scientific data objects appropriate to their
    problem.

    This and the DistTemplateCreate interfaces have been modeled in
    large measure on the array distribution capability present in High
    Performance Fortran version 2.0, including the "Approved
    Extensions.".  This fact manifests itself most strongly in this
    interface in the generality of mapping (or "aligning") the data
    object to the distribution.  Unfortunately, compilers can do this
    a little more neatly than we can, so please read the documentation
    carefully and if necessary refer to HPF2's alignment
    capabilities.

    Ghost regions are not explicitly supported at this time.  The
    current suite of interfaces can be used with data objects
    including ghost regions, but it may not be as pretty as it could
    be.  It is also possible to build a new set of interfaces on top
    of these which include explicit support for ghosts.  We should
    revisit the question of how to support ghosts most effectively
    once we have a little actual experience using these interfaces and
    derived ghost-supporting ones.

    This interface is intended for use in a parallel environment.
    Functions are identified by which processes are expected to call
    them:
    - "collectively" by the entire parallel cohort
    - by individual processes ("MIMD style") active in this process
    topology.
    The term "collective" is used in quotes because all callers must
    invoke the function with the same arguments, but except where
    noted, no synchronization of the processes is implied or required.
    The commit() function is truly collective, in that it does imply
    synchronization of the calling processes.

    To construct a data object, the functions of this interface must
    be called in an appropriate sequence.  Functions with the same
    position in the calling sequence can be called in any order.  The
    sequence is as follows (note that setName() and setDataType() can
    be called any time before commit()):
    -# setTemplate()
    -# setIdentityAlignmentMap()
    -# setMyProcCoords()
    -# setRegionDataPointer()
    -# commit()

    @returns In general, return codes >= 0 signify success and those <
    0 denote error conditions.  Errors in this interface should be
    user recoverable.  Here is a list of error codes used in this
    interface:

    @retval  -2 Invalid or incompatible distribution type specification
    @retval  -6 Invalid process or process topology
    @retval  -8 Method used out of sequence or internal state invalid
    @retval  -9 Internal failure (i.e. memory allocation)
    @retval -11 Attempt to change commit()ed descriptor
    @retval -12 Template not commit()ed
    @retval -13 Template is not of expected type
    @retval -14 Template not defined
    @retval -15 Invalid region ID
    @retval -16 Invalid strides

    @note Error code definitions should be in harmony with those for
    DistArrayTemplate.

    @note This interface miuses the term <i>stride</i>.  The intent of
    the Working Group was that actual memory strides be used in the
    interface. These are the number of memory units to traverse to get
    to the next location in a given dimension of an array.  Using
    strides allows you to express any storage order (i.e. row-major,
    column-major) unambiguously.  Because of time limitations, we have
    used leading dimensions in this implementation, where you specify
    just the <i>size</i> of the dimension and it gives no information
    as to storage order, thus requiring agreement between the creator
    and the user of the array descriptor as to what storage order is
    being used.  For the purposes of these demonstrations, the only
    user was the CUMULVS-based MxN component, which (because it is
    implemented in C) expects C storage order.  This was done for
    expediency and <i>will</i> change in subsequent versions of this
    interface.

*/

class DistArrayDescriptor {
 public:

  /****************************************************************************
   * Constructors and destructors
   ***************************************************************************/

  virtual ~DistArrayDescriptor(){}

  /****************************************************************************
   * Define the descriptor
   ***************************************************************************/

  /** Support just CUMULVS's data types for now.

      Retain stv prefix to make them easier to locate in user code
      because this _will_ be changed in the future.
  */
  enum DataType {
    stv_Int, stv_Float, stv_Cplx, stv_Double, stv_Dcplx, stv_Long,
    stv_Short, stv_Str, stv_Ushort, stv_Uint, stv_Ulong, stv_Byte
  };

  /** Set data type.

      @param type (In) Data type specification.

      @retval   0 Success
      @retval -11 Attempt to change commit()ed descriptor

      @note Called any time prior to commit(). Called by: cohort.

      @note Need a much more general typing mechanism
  */
  virtual int setDataType(const enum DataType type) = 0;

  /** Associate this data object with a distribution template.

      @param template (In) Distribution template

      @retval   0 Success
      @retval -11 Attempt to change commit()ed descriptor

      @note Calling sequence: 1. Called by: cohort.
  */
  virtual int setTemplate(DistArrayTemplate * & templ) = 0;

  /** Sets this process's location in the process topology.

      Example: Consider a 2-d process topology composed of 6 processes
      in a 3x2 arrangement.  The coordinates might be labeled:
      <pre>
      {0,0} {0,1} {0,2}
      {1,0} {1,1} {1,2}
      </pre>
      and might be assigned to processes as follows:
      - proc 0: {0,0}
      - proc 1: {0,1}
      - proc 2: {0,2}
      - proc 3: {1,0}
      - proc 4: {1,1}
      - proc 5: {1,2}
      Note that the association of processes with coordinates in the
      process topology is entirely up to you, the user.  The entire
      purpose of having this routine in the interface is so that we do
      not have to make assumptions about it.

      @param location (In) array containing the coordinates of this
      process in the process topology.  Coordinates are in the
      range 0..N-1 for N processes.

      @retval   0 Success
      @retval  -6 Invalid process (not within declared topology)
      @retval -11 Attempt to change commit()ed descriptor

      @todo Should probably throw an exception instead of returning
      an error code.

      @note Calling sequence: 2. Called by: active processes, MIMD-style */
  virtual int setMyProcCoords(const int procCoords[] ) = 0;

  /** Align object to template with identity mapping.
      
      This alignment operation specifies how the actual data object
      relates to the referenced distribution template.  The data
      object then inherits all the decomposition characteristics from
      the template according to the specified alignment. This allows
      many data objects to use the same distribution template, even if
      the data objects are not the same size as each other or the
      template, and even if they are positioned somewhere other than
      the upperleft corner of the template.  In fact, some rather
      esoteric mappings are possible through the
      setGeneralAlignmentMap() function.  Note that alignment uses the
      actual index space of the data object and template, based on
      their declared global bounds.  This means that the index space
      of the template must always be at least as large as the index
      space of the data object.

      This function specifies a simple identity mapping of data object
      axes and elements to those of the template, in other words,
      data[i,j,...] maps to template[i,j,...].  

      @retval   0 Success
      @retval  -8 Internal state invalid (probably in template)
      @retval -11 Attempt to change commit()ed descriptor
      @retval -13 Template not defined

      @todo Should probably throw an exception rather than returning
      an error code.

      @note Calling sequence: 3. Called by: cohort.
  */
  virtual int setIdentityAlignmentMap() = 0;

  /** Set pointer for the local region of the data object.

      For everything other than explicit distributions, we assume
      there is one contiguous block of memory per process which
      contains all that process's elements of the array.  This implies
      certain storage arrangements in the case of block-cyclic and
      implicit which may seem a bit strange.  I <i>think</i> this is
      what other packages typically do anyway, but this needs further
      investigation.  As noted below, this routine is largely a
      stopgap measure, so we might be better off rethinking the whole
      thing instead of that little factor.

      @param data (In) Pointer to local memory holding region
      @param strides (In) Array of leading dimensions of the actual
      memory storage.  These are the actual lengths of each axis in
      memory.

      @note The combination of an arbitrary pointer and leading
      dimensions allow support of windowing into other arrays.  For
      example, suppose a process holds a 14x14 array as part of a
      distributed array.  This array consists of a 10x10 patch of
      "real" data and a 2 element wide ring of "ghost" elements, which
      duplicate elements on adjacent processes.  One can create a data
      object representing the entire array, explicitly exposing the
      ghost regions, by specifying (in the distribution template) that
      the process's region is 14x14 and registering the pointer to the
      [0,0] element (in C array indexing) with leading dimensions of
      [14,14].  One can also create a data object which hides the
      ghost regions, giving access to only the "real" data by
      registering the region as 10x10 and registering the pointer to
      the [2,2] element with the leading dimensions [14,14].  In other
      words, the beginning of each row/column are 14 elements apart
      instead of 10.

      @retval   0 Success
      @retval  -2 Incompatible distribution type (only valid for
      non-Explicit distributions)
      @retval -11 Attempt to change commit()ed descriptor
      @retval -16 Invalid strides

      @todo Should probably throw an exception instead of returning an
      error code.

      @note Calling sequence: 4. Called by: active processes once per
      process.

      @todo This is just a stopgap until I get a better understanding
      of SIDL and related interlanguage issues.  We already know that
      using pointers is not general, we just have to figure out the
      right way to handle needs like this.  */
  virtual int setLocalDataPointer(void * data, const int
				    strides[]) = 0;

  /** Set pointer for a local region of an explicitly distributed
      data object.

      Each process must call this function for every region of the
      data object they own.  The region's bounds must be provided in
      the data object's global coordinate system, but of course they
      must map to the regions defined for the underlying distribution
      template, taking into account the alignment of the data object
      with the template.  This is clearly a little ugly and perhaps
      more than a little error-prone.  Alternatives welcome!

      @param lower (In) Lower bounds of region, in the data object's
      global coordinates
      @param upper (In) Upper bounds of region, in the data object's
      global coordinates
      @param data (In) Pointer to local memory holding region
      @param strides (In) Array of leading dimensions of the actual
      memory storage.  These are the actual lengths of each axis in
      memory.

      @note The combination of an arbitrary pointer and leading
      dimensions allow support of windowing into other arrays.  For
      example, suppose a process holds a 14x14 array as part of a
      distributed array.  This array consists of a 10x10 patch of
      "real" data and a 2 element wide ring of "ghost" elements, which
      duplicate elements on adjacent processes.  One can create a data
      object representing the entire array, explicitly exposing the
      ghost regions, by specifying (in the distribution template) that
      the process's region is 14x14 and registering the pointer to the
      [0,0] element (in C array indexing) with leading dimensions of
      [14,14].  One can also create a data object which hides the
      ghost regions, giving access to only the "real" data by
      registering the region as 10x10 and registering the pointer to
      the [2,2] element with the leading dimensions [14,14].  In other
      words, the beginning of each row/column are 14 elements apart
      instead of 10.

      @retval   0 Success
      @retval  -2 Incompatible distribution type (only valid for
      non-Explicit distributions)
      @retval -11 Attempt to change commit()ed descriptor
      @retval -16 Invalid strides

      @todo Should probably throw an exception instead of returning an
      error code.

      @note Calling sequence: 4. Called by: active processes once per
      region owned, MIMD-style.

      @todo This is just a stopgap until I get a better understanding
      of SIDL and related interlanguage issues.  We already know that
      using pointers is not general, we just have to figure out the
      right way to handle needs like this.  */
  virtual int setRegionDataPointer(const int lower[], const int 
				    upper[], void * data, const int
				    strides[]) = 0;

  /** Signal that data object is completely defined.  This asserts
      that pointers for all regions have been properly registered, and
      gives the implementation a chance to verify that the alignment
      mappings, etc. are valid.

      @retval   0 Success
      @retval -11 Attempt to change commit()ed descriptor

      @note Calling sequence: 5. Called by: cohort.

      @todo Should probably throw an exception instead of returning
      an error code.
  */
  virtual int commit() = 0;

  /****************************************************************************
   * Query the descriptor
   ***************************************************************************/
  /** Return the name given to the descriptor. */
  virtual std::string getName() = 0;

  /** Has the descriptor been commit()ed? */
  virtual bool isDefined() = 0;

  /** Get data type.
  */
  virtual DataType getDataType() = 0;

  /** Return pointer to distribution template associated with this
      data object.
  */
  virtual DistArrayTemplate * getTemplate() = 0;

  /** Get this process's location in the process topology.

      @param location (Out) array containing the coordinates of this
      process in the process topology.  Coordinates are in the
      range 0..N-1 for N processes.

      @retval   0 Success
      @retval  -8 Method used out of sequence or internal state invalid
  */
  virtual int getMyProcCoords(int procCoords[] ) = 0;

  /** Part of a kludge for debugging.  Supports only GenBlock and
      Explicit distributions at present.

  */
  virtual int getNumLocalRegions() = 0;

  /** Part of a kludge for debugging.  Supports only Explicit
      distributions at present.

      @retval   0 Success
      @retval -15 Invalid region ID
  */
  virtual int getLocalRegionInfo(int region, int lower[], int upper[],
				 void * & data, int strides[]) = 0;

  /** Mainly for testing and debugging */
  virtual void printDescriptor() = 0;
};
#endif // DistArrayDescriptor_h_seen

