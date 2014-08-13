#ifndef DADFTemplate_h_seen
#define DADFTemplate_h_seen

/** This is our implementation of DistArrayTemplate.  It is kept
    private within the DistArrayDescriptorFactory.  Only
    DistArrayTemplate is exposed to the outside. 

    $Id: DADFTemplate.h,v 1.1 2003-08-01 00:41:53 manoj Exp $
 */

#include <string>
#include <vector>
#include <list>
#include "DistArrayTemplate.h"
#include "DADFAxisInfo.h"
#include "DADFRegionInfo.h"

class DADFTemplate : public DistArrayTemplate {
 public:

  /****************************************************************************
   * Constructors and destructors
   ***************************************************************************/

  /** Simple constructor for internal use */
  DADFTemplate() ;

  /** Normal constructor.  Forces setting of name. */
  DADFTemplate( const std::string name ) ;

  /** Simple copy constructor for internal use */
  DADFTemplate( DADFTemplate & original );

  /** Normal copy constructor.  Forces setting of name. */
  DADFTemplate( const std::string name, DADFTemplate & original );

  /** The usual destructor */
  virtual ~DADFTemplate();

  /****************************************************************************
   * Define the template
   ***************************************************************************/

  /** Name associated with this distribution. */
  virtual int setName(const std::string name);

  /** Set rank (number of dimensions) of distribution template. */
  virtual int setRank(const int rank);

  /** Set the global upper and lower bounds of the array. */
  virtual int setGlobalBounds(const int lower[], const int upper[]);

  /** Sets process topology. */
  virtual int setProcTopology(const int topology[] );

  /** Sets distribution type on each axis. */
  virtual int setDistType(const enum DistType dist[] );

  /** Set distribution parameters for an axis with a regular distributions. */
  virtual int setDistParameters(int axis, int blockSize,
				int first);


  /** Set distribution parameters for a GenBlock axis. */
  virtual int setGenBlock(int axis, int blockSizes[]);

  /** Set distribution parameters for an Implicit axis. */
  virtual int setImplicitMap(int axis, int map[]);


  /** Add a region to an Explicit distribution. */
  virtual int addExplicitRegion(int lower[], int upper[]);

  /** Signal that template is completely defined. */
  virtual int commit();

  /****************************************************************************
   * Query the template
   ***************************************************************************/

  /** Name associated with this distribution. (default value is "_UNNAMED") */
  virtual std::string getName();

  /** Get rank (number of dimensions) of distributed object. */
  virtual int getRank();

  /** The global upper and lower bounds of the array. */
  virtual int getGlobalBounds(int lower[], int upper[]);

  /** Returns process topology.  */
  virtual int getProcTopology(int topology[] );

  /** Returns distribution type on each axis. */
  virtual int getDistType(enum DistType dist[] );

  /** Get distribution parameters for an axis with a regular distributions. */
  virtual int getDistParameters(int axis, int blockSize,
				int first);

  /** Get distribution parameters for a GenBlock axis. */
  virtual int getGenBlock(int axis, int blockSizes[]);

  /** Get distribution parameters for an Implicit axis. */
  virtual int getImplicitMap(int axis, int map[]);

  /** Has commit() been called on this template? */
  virtual bool isDefined();
  
  /** Mainly for testing and debugging */
  virtual void printTemplate();

  /****************************************************************************
   * Internals
   ***************************************************************************/
 private:
  // Human-readable name for this template
  std::string _name; 
  
  // Whether or not commit() has been called
  int _frozen;

  // Rank of template
  int _rank;

  /** Caution: it is possible the following (and how they're used)
      could give rise to overflows for large problems. I'm inclined to
      think that if the volume is that large, there will be other
      overflow problems too so this is not important, but I might be
      wrong. 
  */

  // Number of elements
  int _volume;
  
  // Running total of volume of explicit regions
  int _volDefined;

  // Global lower bounds of array
  std::vector<int>                 _lowerBounds;

  // Global upper bounds of array
  std::vector<int>                 _upperBounds;

  // Process topology
  std::vector<int>                 _topology;

  // Distribution types
  std::vector<enum DistType>       _dist;

  // Per-axis distribution parameters
  std::vector<DADFAxisInfo *> _axisInfo;

  // List of regions associated with this process of array
  std::list<DADFRegionInfo *> _regionList;
};
#endif // DADFTemplate_h_seen
