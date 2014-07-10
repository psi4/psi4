#ifndef DADFAxisInfo_h_seen
#define DADFAxisInfo_h_seen

#include <vector>
#include "DistArrayTemplate.h"

/** A set of classes to hold distribution parameters on a per-axis basis.

    $Id: DADFAxisInfo.h,v 1.1 2003-08-01 00:41:53 manoj Exp $

    DADFAxisInfo is the abstract base class, but includes only two
    functions -- one to report the distribution type, and one to print
    the distribution parameters on the stdout stream (as an aid in
    debugging).

    Each type of distribution has a concrete class that inherits from
    DADFAxisInfo and adds specialized functions to set and get the
    parameters appropriate to the type of distribution. (Exception:
    since Collapsed distributions have no such data, they have no such
    functions.)

    At present, no error checking is done here, on the assumption that
    whoever creates us is better equipped to do it -- our purpose is
    mainly storage.
*/

class DADFAxisInfo {
 public:

  /** The usual destructor. */
  virtual ~DADFAxisInfo(){}

  /** Get the type of distribution this object represents. 
   */
  virtual DistArrayTemplate::DistType getDistType() = 0;

  /** Print distribution parameters to stdout stream in a
      human-readable form. 

      This is mainly intended to facilitate debugging.
  */
  virtual void printAxisInfo() = 0;
};

/** Concrete DADFAxisInfo for collapsed distribution. */

class CollapsedAxisInfo : public DADFAxisInfo {
 public:

  /** The usual destructor. */
  virtual ~CollapsedAxisInfo();

  /** Get the type of distribution this object represents. 
   */
  virtual DistArrayTemplate::DistType getDistType();

  /** Print distribution parameters to stdout stream in a
      human-readable form. 

      This is mainly intended to facilitate debugging.
  */
  virtual void printAxisInfo();
};

/** Concrete DADFAxisInfo for block distribution. */

class BlockAxisInfo : public DADFAxisInfo {
 public:

  /** Normal constructor.
   */
  BlockAxisInfo();

  /** Copy constructor.

      @param original (in) Instance of BlockAxisInfo from which to
      initialize new object.
  */
  BlockAxisInfo(const BlockAxisInfo & original);

  /** The usual destructor */
  virtual ~BlockAxisInfo();

  /** Get the type of distribution this object represents. 
   */
  virtual DistArrayTemplate::DistType getDistType();

  /** Set block distribution parameters.

      @param blockSize (in) Size of block.
      @param first (in) process owning first block of distribution.
  */
  void setDistParameters(const int blockSize, const int first);

  /** Get block distribution parameters.

      @param blockSize (out) Size of block.
      @param first (out) process owning first block of distribution.
  */
  void getDistParameters(int & blockSize, int & first);

  /** Print distribution parameters to stdout stream in a
      human-readable form. 

      This is mainly intended to facilitate debugging.
  */
  virtual void printAxisInfo();

 private:
  /** Block size for block/cyclic distribution. */
  int _blockSize;

  /** Process on which first block resides */
  int _first;
};

/** Concrete DADFAxisInfo for generalized block distribution. */

class GenBlockAxisInfo : public DADFAxisInfo {
 public:

  /** Normal constructor. 

      @param size (in) Size of block size vector (should be same as
      number of processes in this axis of process topology).
  */
  GenBlockAxisInfo(const int size);

  /** Copy constructor.

      @param original (in) Instance of GenBlockAxisInfo from which to
      initialize new object.
  */
  GenBlockAxisInfo(const GenBlockAxisInfo & original);

  /** The usual destructor */
  virtual ~GenBlockAxisInfo();

  /** Get the type of distribution this object represents. 
   */
  virtual DistArrayTemplate::DistType getDistType();

  /** Set generalized block distribution parameters.

      @param blockSizes (in) Array of block sizes
  */
  void setDistParameters(const int blockSizes[]);

  /** Get generalized block distribution parameters.

      @param blockSizes (out) Array of block sizes
  */
  void getDistParameters(int blockSizes[]);

  /** Print distribution parameters to stdout stream in a
      human-readable form. 

      This is mainly intended to facilitate debugging.
  */
  virtual void printAxisInfo();

 private:

  /** Array of block sizes.  Should size as the number of processes. */
  std::vector<int> _blockSizes;
};

/** Concrete DADFAxisInfo for an implicit distribution. */

class ImplicitAxisInfo : public DADFAxisInfo {
 public:

  /** Normal constructor.  

      @param size (in) Size of element to process map (should be same
      as number of elements in this axis of array template).  
  */
  ImplicitAxisInfo(const int size);

  /** Copy constructor.

      @param original (in) Instance of ImplicitAxisInfo from which to
      initialize new object.  
  */
  ImplicitAxisInfo(const ImplicitAxisInfo & original);

  /** The usual destructor */
  virtual ~ImplicitAxisInfo();

  /** Get the type of distribution this object represents. 
   */
  virtual DistArrayTemplate::DistType getDistType();

  /** Set implicit map distribution parameters.

      @param map (in) Array mapping elements to processes
  */
  void setDistParameters(const int map[]);

  /** Get implicit map distribution parameters.

      @param map (out) Array mapping elements to processes
  */
  void getDistParameters(int map[]);

  /** Print distribution parameters to stdout stream in a
      human-readable form. 

      This is mainly intended to facilitate debugging.
  */
  virtual void printAxisInfo();

 private:
  /** Mapping of elements to processes. */
  std::vector<int> _map;
};

#endif // DADFAxisInfo_h_seen
