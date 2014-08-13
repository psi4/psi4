#ifndef DistArrayTemplate_h_seen
#define DistArrayTemplate_h_seen

#include <string>

/** This is the public interface (abstract base class) for defining
    and querying array distribution templates for dense
    multi-dimensional rectangular arrays.

    This interface is an experimental first implementation of an
    interface which has been under discussion in the CCA Forum's
    Scientific Data Components Working Group.  This implementation
    does not conform exactly to the interface developed by the Working
    Group, and is intended to be exploratory rather than normative.

   $Id: DistArrayTemplate.h,v 1.1 2003-08-01 00:41:54 manoj Exp $

    This interface is intended to support the creation of distribution
    templates to which actual data objects can later be
    aligned. Following the model of High Performance Fortran,
    distribution templates are conceptually arrays which have no real
    existence, but for which the decomposition on the processor array
    is specified. Data objects can be created referencing the template
    for their decomposition.  Many data objects can be "aligned" to
    the same template.

    This interface is intended to be used in several ways.  Many
    existing distributed array and related packages do not have the
    distribution as a first-class, standlone object.  In this case,
    the distribution template can be created separately, allowing
    other components to access the information through the
    complementary DataDistQuery interface.  In newer packages, a
    distribution template created here might be handed to a factory
    to create an actual data object.   This interface can also serve
    as a guide for development of data distributions specification
    interfaces for existing or new distibuted array packages which may
    not support the full range of flexibility embodied here.

    This interface is desgned to support most of the array
    distribution capability present in High Performance Fortran
    version 2.0, including the "Approved Extensions.".  It also
    supports completely general "explicit" or table-based
    distributions. Based on an informal survey of distributed array
    and related packages, this seemed to be the best target.

    Data distributions are independent of the data type and storage
    order (i.e. row-major vs column-major).

    Data distributions are essentially assumed to be static.  It is
    important that the user be able to make a set of queries about a
    distribution and be confident of getting a consistent set of
    results.  Therefore, if an implementation wants to support dynamic
    distributions, it should provide mechanisms/rules so that the user
    can be confident of getting consistent results from properly
    formulated inquiries.  One way to do this might be by making any
    redistribution operations explicit user calls, so that while not
    in such a call, the data is consistent.

    Distribution templates do not allow the direct specification of
    replication, but in the manner of HPF2, it is possible to align
    real data to a template in such a way that the real data is
    replicated across many elements of the template.

    Ghost regions are not explicitly supported at this time.  The
    current suite of interfaces can be used with data objects
    including ghost regions, but it may not be as "pretty" as it could
    be.  It is also possible to build a new set of interfaces on top
    of these which include explicit support for ghosts.  We should
    revisit the question of how to support ghosts most effectively
    once we have a little actual experience using these interfaces and
    derived ghost-supporting ones.

    This interface is intended for use in a parallel environment.
    Functions are identified by which processes are expected to call
    them:
    - "collectively" by the entire parallel cohort
    - by individual participating processes ("MIMD style").
    The term "collective" is used in quotes because all callers must
    invoke the function with the same arguments, but except where
    noted, no synchronization of the processes is implied or required.

    To construct a distribution template, the functions of this
    interface must be called in an appropriate sequence.  Functions
    with the same position in the calling sequence can be called in
    any order.  The sequence is as follows:
    -# setRank()
    -# setGlobalBounds(), setProcTopology(), setDistType()
    -# setDistParameters(), setGenBlock(), setImplicitMap(),
    addExplicitRegion()
    -# commit()

    @returns In general, return codes >= 0 signify success and those <
    0 denote error conditions.  Errors in this interface should be
    user recoverable.  Here is a list of error codes used in this
    interface:

    @retval  -1 Invalid bounds
    @retval  -2 Invalid or incompatible distribution type specification
    @retval  -3 Invalid rank
    @retval  -4 Invalid axis
    @retval  -5 Invalid blockSize
    @retval  -6 Invalid process or process topology
    @retval  -7 Region overlaps existing region
    @retval  -8 Method used out of sequence or internal state invalid
    @retval  -9 Internal failure (i.e. memory allocation)
    @retval -10 Attempt to change commit()ed template */

class DistArrayTemplate {
 public:

  /****************************************************************************
   * Constructors and destructors
   ***************************************************************************/

  /** The usual destructor */
  virtual ~DistArrayTemplate(){}

  /****************************************************************************
   * Define the template
   ***************************************************************************/

  /** Name associated with this distribution. (default value is
      "_UNNAMED")

      @param name (In) Name of distribution template, for convenience
             of debugging.

      @retval   0 Success
      @retval -10 Attempt to modify commit()ed template

      @note May be used at any time before commit() is called.  Called
      by: cohort.
 */
  virtual int setName(const std::string name) = 0;

  /** Set rank (number of dimensions) of distribution template. 

      @param rank (In) Rank (number of dimensions) in array template.

      @retval   0 Success
      @retval  -3 Invalid rank, must be a positive integer.
      @retval -10 Attempt to modify commit()ed template

      @note Calling sequence: 1. Called by: cohort.

      @todo Should probably throw an exception instead of returning an
      error code.

      @note Implementation note: the current implementation tries to
      do sensible things of the rank is changed before the template is
      commit()ed. (It is open to debate whether/under what
      circumstances it is desirable to change the rank.)  If the rank
      is reduced, all rank-based data structures will be truncated to
      the new rank.  If the rank is increased, the new elements are
      filled with something mildly invalid (i.e. array bounds of 0:-1)
      or something innocuous (i.e collapsed distributions) with the
      expectation that they will be properly set later on.  For
      explicit distributions, a change in rank will cause the region
      list to be cleared.
  */
  virtual int setRank(const int rank) = 0;

  /** Set the global upper and lower bounds of the array index space.

      @param lower (In) array of global lower bounds of array indices
      @param upper (In) array of global upper bounds of array indices

      @retval   0 Success
      @retval  -1 Bounds invalid (i.e. upper bound smaller than lower)
      @retval  -8 Method used out of sequence or internal state invalid
      @retval -10 Attempt to modify commit()ed template

      @note Calling sequence: 2. Called by: cohort.

      @todo Should probably throw an exception instead of returning an
      error code.
  */
  virtual int setGlobalBounds(const int lower[], const int upper[]) = 0;

  /** Sets process topology.  This is the number of processes in each
      in axis of the distribution.  The topology array is of length
      getRank().  Axes not distributed have an entry of 1 in the
      topology array.  The product of the elements of the topology
      array should be the total number of processes participating in
      the distribution template.

      For example, a 3-d array on a 6-process system might have a
      topology with three processes in the first dimension, one in
      the second (this dimension not distributed) and two in the
      third, which would be a topology of [3, 1, 2]. 3x1x2 = 6.

      @param topology (In) array containing the number of processes
      involved in each axis of the data distribution.

      @retval   0 Success
      @retval  -6 Invalid process topology (must be at least one
      process in every dimension)
      @retval  -8 Method used out of sequence or internal state invalid
      @retval -10 Attempt to modify commit()ed template

      @note Explicit distributions are by convention described by
      effectively one-dimensional process topologies.  In other
      words, the topology array will be [getNumProcs(), 1, 1, ...].

      @note Calling sequence: 2. Called by: cohort.  
  */

    virtual int setProcTopology(const int topology[] ) = 0;

  /** The types of distributions supported. */

  enum DistType {
    Collapsed, /**< Designates a dimension held entirely on a single
		  process */

    Block, /**< HPF-style block distribution.  There may be zero or
	      more blocks on each process, depending on the block
	      size. Smaller block sizes lead to multiple blocks per
	      process laid out in cyclic fashion -- the block-cyclic
	      distribution widely used in linear algebra
	      packages. This type also describes "cyclic"
	      distributions. */

    GenBlock, /**< HPF-style generalized block distribution.

	A generalization of the simple block distribution which allows
	blocks to be of arbitrary size on each process, one block per
	process.  This DistType can be combined with any other
	DistType in the other axes except for Explicit.  However since
	this DistType is irregular, it cannot be represented in the
	simple parameters of getDistParameters -- the process-based
	inquiry functions or the getMap function must be used instead.

	Here are some examples of 5x5 arrays distributed over 4
	processes using various combinations of GenBlock and other
	DistTypes with with block sizes of 4+1 for the GenBlock.
	Values indicate process number.
	<pre>
	Row DistType: GenBlock      GenBlock      Cyclic
	Col DistType: GenBlock      Block         GenBlock
	              0 0 0 0 2     0 0 0 2 2     0 0 0 0 2
		      0 0 0 0 2     0 0 0 2 2     1 1 1 1 3
		      0 0 0 0 2     0 0 0 2 2     0 0 0 0 2
		      0 0 0 0 2     0 0 0 2 2     1 1 1 1 3
		      1 1 1 1 3     1 1 1 3 3     0 0 0 0 2
	</pre> */

    Implicit, /**< HPF-style implicitly mapped distribution.

	Allows a completely general mapping of elements to processes
	on a per-axis basis (overall distribution must be decomposible
	as a cartesian product of the per-axis mapping).  This
	DistType can be combined with any other DistType in the other
	axes except for Explicit.  However since this DistType is
	irregular, it cannot be represented in the simple parameters
	of getDistParameters -- the process-based inquiry functions or
	the getMaping function must be used instead.

	Each position in an implicit mapping array holds the process
	number (on that axis) on which elements at that position of
	the array live.  By way of example, consider a 9x9 matrix laid
	out on a process grid that has three processes across the top
	and two down the side.  We can label these processes A-F. The
	matrix, shown on the left, have been filled in with the
	label of the process on which that element lives, and the rows
	and columns of the matrix have been labeled with process
	number on that axis.  The corresponding implicit maps for the
	rows and columns are shown on the right.
        <pre>
           0 1 1 2 2 2 1 1 0    Row map array: [0, 1, 1, 2, 2, 2, 1, 1, 0]
           - - - - - - - - -    Column map array: [0, 1, 0, 1, 1, 0, 0]
        0: A B B C C C B B A
        1: D E E F F F E E D    Process grid:
        0: A B B C C C B B A           0 1 2
        1: D E E F F F E E D        0: A B C
        1: D E E F F F E E D        1: D E F
        0: A B B C C C B B A
        0: A B B C C C B B A
	</pre> */

    Explicit /**< Designates a distribution that is completely user
        specified and in general the full distribution cannot be
        represented as the cartesian product of the axes.  The
        Explicit DistType must be used on all dimensions of a
        distribution.

	This DistType is a step beyond the GenBlock because this one
	cannot in general be expressed one axis at a time.  Here
	are some examples of Explicit distributions of a 5x5 array
	over 4 processes (values indicate process number):
	<pre>
	0 0 0 0 2     0 0 0 0 0     0 0 0 0 0
	0 0 0 0 2     0 0 0 0 0     0 0 2 2 0 
	1 1 3 3 3     1 1 2 2 2     1 1 2 2 3
	1 1 3 3 3     1 1 2 2 2     1 1 3 3 3 
	1 1 3 3 3     1 1 2 2 2     1 1 3 3 3
	</pre>
	Note that in the third case, process 0 has 3 distinct patches
	and process 3 has 2 (there are of course several choices for
	the exact layout of these patches).  */
  };
  
  /** Sets distribution type on each axis.

      @param dist (In) array of distribution types of size getRank().
      Note that all axes or none must be specified as Explicit, or it
      will return an error.

      @retval   0 Success
      @retval  -2 Invalid distribution type specification (either all
      axes must be Explicit or none)
      @retval  -8 Method used out of sequence or internal state invalid
      @retval -10 Attempt to modify commit()ed template

      @note Calling sequence: 2. Called by: cohort.

      @note An alternative (more object-oriented) approach would
      involve the construction of objects defining each axis and then
      using an array of them to define the distribution.  Note that
      this still wouldn't accomodate Explicit distributions, at least
      not without further work.  Comments welcome!
  */
  virtual int setDistType(const enum DistType dist[]) = 0;

  /** Set distribution parameters for an axis with a regular distributions.

      Uses a compact set of parameters which, together with
      information about the process topology, fully define the
      distribution (only for regular distributions, of course).  If
      blockSize is smaller than that required to give one block per
      process, blocks are laid out in cyclic fashion, allowing the
      specification of block-cyclic distributions.

      @param axis (In) Axis for which descriptor is being defined
      (0..getRank()-1)
      @param blockSize (In) Block size for this axis
      @param first (In) Process coordinate (0..N-1) to which first
      block is assigned in this axis

      @retval   0 Success
      @retval  -2 Incompatible distribution type (not Block)
      @retval  -4 Invalid axis
      @retval  -5 Invalid blockSize (block size larger than declared bounds)
      @retval  -6 Invalid process (first does not fit declared topology)
      @retval -10 Attempt to modify commit()ed template

      @note Calling sequence: 3. Called by: cohort

      @todo Should probably throw an exception instead of returning
      an error code.  
  */
  virtual int setDistParameters(int axis, int blockSize,
				int first) = 0;


  /** Set distribution parameters for a GenBlock axis.

      @param axis (In) axis for which descriptor is being defined
      (0..getRank()-1)
      @param blockSizes (In) Array of block size on each process (in
      this axis).  Length of array is the number of processes on this
      axis of the process topology.

      @retval   0 Success
      @retval  -2 Incompatible distribution type (not GenBlock)
      @retval  -4 Invalid axis
      @retval  -5 Invalid blockSize (Size of individual block or total
      size inconsistent with declared bounds)
      @retval -10 Attempt to modify commit()ed template

      @note Calling sequence: 3. Called by: cohort.

      @todo Should probably throw an exception instead of returning an
      error code.
  */
  virtual int setGenBlock(int axis, int blockSizes[]) = 0;

  /** Set distribution parameters for an Implicit axis.

      @param axis (In) axis for which descriptor is being defined
      (0..getRank()-1)
      @param map (In) Array mapping elements to processes.  Length of
      array is the number of elements in this axis of the template,
      values refer processes on this axis of the process topology
      (0..N-1).

      @retval   0 Success
      @retval  -2 Incompatible distribution type (not Implicit)
      @retval  -4 Invalid axis
      @retval  -6 Invalid process (in map array)
      @retval -10 Attempt to modify commit()ed template

      @note Calling sequence: 3. Called by: cohort.

      @todo Should probably throw an exception instead of returning an
      error code.
  */
  virtual int setImplicitMap(int axis, int map[]) = 0;


  /** Add a region to an Explicit distribution.

      Specifies the ownership of a given sub-array of the template.
      Applicable only to templates with Explicit distribution on all
      axes.

      Such regions must tile the template through repeated calls to
      this function.  This is an MIMD-style call -- each region must
      be specified exactly once across the entire cohort.  Typical
      usage would involve each process computing its own portion of
      the distribution, however it is allowed for one process to add
      regions that live on another (in case, for example, the entire
      distribution is computed on one node).

      @param lower (In) lower bounds of sub-array region
      @param upper (In) upper bounds of sub-array region

      @retval   0 Success
      @retval  -1 Invalid bounds
      @retval  -7 Region overlaps with existing region
      @retval -10 Attempt to modify commit()ed template

      @note Calling sequence: 3. Called by: individual (MIMD
      style).

      @note At one point in its design, this method had an additional
      argument: "procCoords (In) Coordinates of process of interest.
      Recall that Explicit distributions are by convention described
      by effectively one-dimensional process topologies."  There is
      some question in my mind as to whether this argument is
      necessary or desirable.  Without it, the owning process must
      register the region, while with it any process can register any
      region.  This may be useful, because it allows an arbitrary
      process to compute and register the distribution.

      @note There is no way to remove a region.  Perhaps there should
      be? This particularly effects cloned templates.

      @todo Should probably throw an exception instead of returning
      an error code.  

  */
  virtual int addExplicitRegion(int lower[], int upper[]) = 0;

  /** Signal that template is completely defined. Up to this point,
      things can be added or changed (as long as there is no other
      conflict. After calling commit(), the distribution template is
      set in stone and cannot be changed. 

      @retval   0 Success
      @retval -10 Attempt to modify commit()ed template

      @note Calling sequence: 4. Called by: cohort.

      @todo Should probably throw an exception instead of returning
      an error code.

      @note Implementation note: There are a lot of consistency checks
      that could/should be done at this point but aren't.  If you
      build up a template in the normal fashion, many of those checks
      will have been applied piecemeal.  But if you clone and then
      modify a template, or just abuse the calling sequence too much,
      it is possible you could have set some things that are not
      globally consistent.

*/
  virtual int commit() = 0;

  /****************************************************************************
   * Query the template
   *
   * Should these be available only after the template is commit()ed?
   ***************************************************************************/

  /** Name associated with this distribution. (default value is "_UNNAMED") */
  virtual std::string getName() = 0;

  /** Get rank (number of dimensions) of distributed object. */
  virtual int getRank() = 0;

  /** The global upper and lower bounds of the array.

      @param lower (Out) array of global lower bounds of array
      @param upper (Out) array of global upper bounds of array

      @retval   0 Success
      @retval  -8 Internal state invalid
  */
  virtual int getGlobalBounds(int lower[], int upper[]) = 0;

  /** Returns process topology.  This is the number of processes in
      each in axis of the distribution.  The topology array is of
      length getRank().  Dimensions not distributed have an entry of 1
      in the topology array.  The product of the elements of the
      topology array should be getNumProcs().

      For example, a 3-d array on a 6-process system might have a
      topology with three processes in the first axis, one in
      the second (this axis not distributed) and two in the
      third, which would be a topology of [3, 1, 2]. 3x1x2 = 6.

      @param topology (Out) array containing the number of
      processes involved in each axis of the data distribution.

      @retval   0 Success
      @retval  -8 Internal state invalid

      @note Explicit distributions are by convention described by
      effectively one-dimensional process topologies.  In other
      words, the topology array will be [getNumProcs(), 1, 1, ...].

  */
  virtual int getProcTopology(int topology[] ) = 0;

  /** Returns distribution type on each axis.

      @param dist (Out) array of distribution types of size getRank()

      @retval   0 Success
      @retval  -8 Internal state invalid
  */
  virtual int getDistType(enum DistType dist[] ) = 0;

  /** Get distribution parameters for an axis with a regular distributions.

      Support block/block-cyclic/cyclic distributions

      @param axis (In) Axis for which descriptor is being defined
      (0..getRank()-1)
      @param blockSize (Out) Block size for this axis
      @param first (Out) Process coordinate (0..N-1) to which first
      block is assigned in this axis

      @retval   0 Success
      @retval  -2 Incompatible distribution type (not Block)
      @retval  -4 Invalid axis
      @retval  -8 Method used out of sequence or internal state invalid
  */
  virtual int getDistParameters(int axis, int blockSize,
				int first) = 0;


  /** Get distribution parameters for a GenBlock axis.

      @param axis (In) axis for which descriptor is being defined
      (0..getRank()-1)
      @param blockSizes (Out) Array of block size on each process (in
      this axis).  Length of array is the number of processes on this
      axis of the process topology.

      @retval   0 Success
      @retval  -2 Incompatible distribution type (not GenBlock)
      @retval  -4 Invalid axis
      @retval  -8 Method used out of sequence or internal state invalid
  */
  virtual int getGenBlock(int axis, int blockSizes[]) = 0;

  /** Get distribution parameters for an Implicit axis.

      @param axis (In) axis for which descriptor is being defined
      (0..getRank()-1)
      @param map (Out) Array mapping elements to processes.  Length of
      array is the number of elements in this axis of the template,
      values refer processes on this axis of the process topology
      (0..N-1).

      @retval   0 Success
      @retval  -2 Incompatible distribution type (not Implicit)
      @retval  -4 Invalid axis
      @retval  -8 Method used out of sequence or internal state invalid
  */
  virtual int getImplicitMap(int axis, int map[]) = 0;

  /** Check if the tempate is fully defined (i.e. has been
      commit()ed).

      @retval true if commit() has been successfully called on this template
      @retval false otherwise
  */
  virtual bool isDefined() = 0;

  /** Mainly for testing and debugging */
  virtual void printTemplate() = 0;

};
#endif // DistArrayTemplate_h_seen
