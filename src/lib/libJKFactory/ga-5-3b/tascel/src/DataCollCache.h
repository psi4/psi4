#ifndef __DataCollCache_h__
#define __DataCollCache_h__

#include "AccessMode.h"

namespae tascel {

  class DataColl;

  /**
   * Cache of data blocks in an abstract data collection.
   */
  class DataCollCache {
  private:
    const int idxlen;
    const int size;
    DataColl *coll;
    AccessMode mode;

  public:
    DataCollCache(int idxlen, int size, DataColl *coll, AccessMode mode);
    
    ~DataCollCache();    

    void *get(void *idx, int idxlen);

    void *allocate(void *idx, int idxlen);
    
    void put(void *idx, int idxlen, void *buf);
    
    void add(void *idx, int idxlen, void *buf);

    void flush();

  }; /*DataCollCache*/
}; /*tascel*/


#endif /*__DataCollCache_h__*/




