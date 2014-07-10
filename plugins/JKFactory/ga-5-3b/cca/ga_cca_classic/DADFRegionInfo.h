#ifndef DADFRegionInfo_h_seen
#define DADFRegionInfo_h_seen

#include <vector>

/** A concrete class to hold a region specification for a
    multidimensional array.

   $Id: DADFRegionInfo.h,v 1.1 2003-08-01 00:41:53 manoj Exp $

   There's nothing fancy here, just lower bounds and upper bounds.
   Maybe eventually we'll need something fancier, with IDs/handles for
   the regions and other stuff.  But for now this suffices.
 */

class DADFRegionInfo {
 public:

  /** Normal constructor. 

      @param size (in) length of the upper/lower bounds arrays
  */
  DADFRegionInfo( int size );

  /** Copy constructor.

      @param original (in) Instance of DADFRegionInfo from which to
      initialize new object.
  */
  DADFRegionInfo( const DADFRegionInfo & original);

  ~DADFRegionInfo();

  /** Set region bounds.

      @param lower (in) lower bounds of region
      @param upper (in) upper bounds of region
  */
  void setBounds(const int lower[], const int upper[]);

  /** Get region bounds.

      @param lower (in) lower bounds of region
      @param upper (in) upper bounds of region
  */
  void getBounds(int lower[], int upper[]);

  /** Compare region bounds.  Determine if object's region and
      argument region are identical, disjoint, or overlapping.

      @param lower (in) lower bounds of region
      @param upper (in) upper bounds of region

      @retval  1 Object region and argument region are identical
      @retval  0 Object region and argument region are disjoint
      @retval -1 Object region and argument region overlap but aren't
      identical
  */
  int compareBounds(const int lower[], const int upper[]);

  /** Set data location.

      @param dataPtr (in) pointer to data
      @param strides (in) stride in each dimension to access data locations

  */
  void setDataLocation(void * dataPtr, const int strides[] );

  /** Get data location.

      @param dataPtr (out) pointer to data
      @param strides (out) stride in each dimension to access data locations

  */
  void getDataLocation(void * & dataPtr, int strides[] );

  void printRegionInfo();

 private:
  std::vector<int> _lowerBounds;
  std::vector<int> _upperBounds;
  std::vector<int> _strides;
  void * _data;
};

#endif // DADFRegionInfo_h_seen
