#ifndef DistArray_h_seen
#define DistArray_h_seen

#include "DistArrayTemplate.h"
#include "DistArrayDescriptor.h"

class DistArray {
 public:

  /****************************************************************************
   * Constructors and destructors
   ***************************************************************************/

  virtual ~DistArray(){}

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
  virtual void printArray() = 0;
  virtual void printArrayDistribution() = 0;

  

};
#endif // DistArray_h_seen

