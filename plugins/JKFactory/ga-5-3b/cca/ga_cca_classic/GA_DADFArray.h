#ifndef GA_DADFArray_h_seen
#define GA_DADFArray_h_seen

/** This is our implementation of DistArray.  It is kept
    private within the GAServices.  Only  DistArray is exposed 
    to the outside.
*/

#include <string>
#include <vector>
 
#include "DistArray.h"
#include "DADFTemplate.h"
#include "DADFDescriptor.h"

class DADFArray : public DistArray {
 public:
 
  /****************************************************************************
   * Constructors and destructors
   ***************************************************************************/
 
  /** Simple constructor for internal use. */
  DADFArray() ;
 
  /** Normal constructor -- forces setting of name */
  DADFArray( const std::string name ) ;
 
  /** Copy constructor */
  DADFArray( const std::string name, DADFArray & original );
 
  /** The usual destructor */
  virtual ~DADFArray();
 
  /** Set data type. */
  virtual int setDataType(const enum DataType type);
 
  /** Associate this data object with a distribution template. */
  virtual int setTemplate(DistArrayTemplate * & templ);
 
  /** Sets this process's location in the process topology. */
  virtual int setMyProcCoords(const int procCoords[] );
 
  /** Align object to template with identity mapping. */
  virtual int setIdentityAlignmentMap();
 
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
  virtual DistArray::DataType getDataType();
 
  /** Return pointer to distribution template associated with this
      data object.
  */
  virtual DistArrayTemplate * getTemplate();
 
  /** Get this process's location in the process topology. */
  virtual int getMyProcCoords(int procCoords[] );
 
  /** Mainly for testing and debugging */
  virtual void printArrayDistribution();

  virtual void printArray();
 
  /** Part of a kludge for debugging. */
  virtual int getNumLocalRegions();
 
  /** Part of a kludge for debugging. */
  virtual int getLocalRegionInfo(int region, int lower[], int upper[],
                                 void * & data, int strides[]);

  /****************************************************************************
   * Internals
   ***************************************************************************/
 private:
 
  // GA specific array Handle
  int _handle;
  
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

#endif // GA_DADFArray_h_seen
