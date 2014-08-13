#ifndef DADFDescriptor_h_seen
#define DADFDescriptor_h_seen

/** This is our implementation of DistArrayDescriptor.  It is kept
    private within the DistArrayDescriptorFactory.  Only
    DistArrayDescriptor is exposed to the outside. 

    $Id: DADFDescriptor.h,v 1.1 2003-08-01 00:41:53 manoj Exp $
 */

#include <string>
#include <vector>

#include "DistArrayDescriptor.h"
#include "DADFTemplate.h"

class DADFDescriptor : public DistArrayDescriptor {
 public:

  /****************************************************************************
   * Constructors and destructors
   ***************************************************************************/

  /** Simple constructor for internal use. */
  DADFDescriptor() ;

  /** Normal constructor -- forces setting of name */
  DADFDescriptor( const std::string name ) ;

  /** Copy constructor */
  DADFDescriptor( const std::string name, DADFDescriptor & original );

  /** The usual destructor */
  virtual ~DADFDescriptor();

  /** Set data type. */
  virtual int setDataType(const enum DataType type);

  /** Associate this data object with a distribution template. */
  virtual int setTemplate(DistArrayTemplate * & templ);

  /** Sets this process's location in the process topology. */
  virtual int setMyProcCoords(const int procCoords[] );

  /** Align object to template with identity mapping. */
  virtual int setIdentityAlignmentMap();

  /** Set pointer for the local region of the data object. */
  virtual int setLocalDataPointer(void* data, const int strides[]);

  /** Set pointer for a local region of an explicitly distributed
      data object.
  */
  virtual int setRegionDataPointer(const int lower[], const int
				    upper[], void* data, const int
				    strides[]);

  /** Signal that data object is completely defined. */
  virtual int commit();

  /****************************************************************************
   * Query the descriptor
   ***************************************************************************/

  /** Return name given to descriptor */
  virtual std::string getName();

  /** Has commit() been called on this descriptor? */
  virtual bool isDefined();

  /** Get data type. */
  virtual DistArrayDescriptor::DataType getDataType();

  /** Return pointer to distribution template associated with this
      data object. 
  */
  virtual DistArrayTemplate * getTemplate();

  /** Get this process's location in the process topology. */
  virtual int getMyProcCoords(int procCoords[] );

  /** Mainly for testing and debugging */
  virtual void printDescriptor();

  /** Part of a kludge for debugging. */
  virtual int getNumLocalRegions();

  /** Part of a kludge for debugging. */
  virtual int getLocalRegionInfo(int region, int lower[], int upper[],
				 void * & data, int strides[]);

  /****************************************************************************
   * Internals
   ***************************************************************************/
 private:

  // Human-readable name for this descriptor
  std::string _name;           
  
  // Whether or not commit() has been called
  int _frozen;            

  // Rank of array
  int _rank;

  // Type of data
  DataType _type;

  // Shorthand check if this is an explicit distribution
  bool _isExplicitDist;

  // Array distribution template for this descriptor
  DistArrayTemplate * _templ;

  // Global lower bounds of array
  std::vector<int>                 _lowerBounds;

  // Global upper bounds of array
  std::vector<int>                 _upperBounds;

  // Process topology
  std::vector<int>                 _topology;

  // Coordinates of this process in topology
  std::vector<int>                 _procCoords;

  // List of regions associated with this process of array
  std::list<DADFRegionInfo *> _regionList;
};
#endif // DADFDescriptor_h_seen
