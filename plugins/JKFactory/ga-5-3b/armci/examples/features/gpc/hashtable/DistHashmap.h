/* $Id:  */
#ifndef DISTHASHMAP_H_
#define DISTHASHMAP_H_

#include <string>
using std::string;
#define ARMCI_ENABLE_GPC_CALLS
#include "gpc.h"
#include "Hash_common.h"

class DistHashmap {
      
 public:
  /**
   * Constructor
   */
  DistHashmap();
  
  /**
   * Default Destructor
   */
  ~DistHashmap();
  
  /**
   * creates a new distributed hashmap
   */
  void create();
  
  /**
   * destroys this distributed hashmap
   */
  void destroy();
  
  /**
   * str - string to be inserted into the distributed hashmap
   * size - size of strlen array
   * Insert() is just to ensure that there is a new entry to be inserted into
   * the distributed hash map, and this entry is marked in the local buffer
   * as "to be inserted". The reason is, it will be expensive to do the insert
   * for each and every new element. In order to avoid those latency costs,
   * insert aggregates the entries into a single buffer of size BUF_LIMIT and
   * only completes the insertion into the distributed hashh map, if and only
   * if this aggregate local buffer is full. In order to complete this
   * aggregation call (even though the buffer is not full), commit() can be
   * used to flush the local buffer into the distributed hash map. So it is a
   * good practice to call commit after all insert or whenever consistency is
   * required.
   */
  void insert(string str);
   
  /**
   * Call commit() to complete the data insertion into the distributed
   * hashmap.
   */
  void commit();

  /**
   * Call commit() to complete the data insertion into the distributed
   * hashmap of a specified process.
   */
  void commit(int proc);
    
  /**
   * prints this distributed hashmap
   */
  void print();

  void print2();

  void rehash();

  const VocabIntMap * getLocalMapPtr();
  
  /**
   * returns the Global HashMap Size
   */
  int getGlobalHashMapSize();
  
  /**
   * returns true if a hashmap already exists
   */
  static bool isCreated();
  
 private:
  int mGPCHandle;
  int m_me;
  static short int sm_initialized;
  char **mBuf;
  char *mTmpBuf; /* temporary buffer */
  int *mIndex;
  int m_globalHashMapSize;
  VocabIntMap *mGlobalIdMap;
};

#endif /* DISTHASHMAP_H_ */
