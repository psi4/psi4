
#include <cstring>
#include "Util.h"
#include "armci.h"

int armci_hashmap_pack(char *buf, string str) {

    int len = str.length();
    
    // first copy the length of the string
    *((int*)buf) = len;
    
    // then copy the actual string
    int index = sizeof(int);
    char *dst = &buf[index];
    strcpy(dst, str.c_str());
    index += len;
    
    // offset "index" to integer-byte boundary
    int adjust = sizeof(int) - (index % sizeof(int));
    index += adjust;
    
    return (int)index;
}

int armci_hashmap_unpack(const char *buf, string& str) {

    // first get the string length and the corresponding string
    int len = *((int*)buf);
    
    // then get the actual term string
    int index = sizeof(int);
    string termStr(&buf[index], len);
    str = termStr;
    index += len;
    
    // offset "index" to integer-byte boundary
    int adjust = sizeof(int) - (index % sizeof(int));
    index += adjust;

    return (int)index;
}

void armci_hashmap_insert(VocabIntMap *vocabMap,
                          const char *buf, size_t bufsize) {

    string termStr;
    size_t index=0;

    // get the number of strings to be inserted (i.e.first field in the
    // buffer) and increment the index accordingly.
    int numstrings = *((int*)buf);
    index += sizeof(int);

#if DEBUG
    int me; MP_MYID(&me);
    printf("%d: armci_hashmap_insert(): numstrings=%d bufsize=%ld\n",
           me, numstrings, bufsize); fflush(stdout);
#endif
    
    // unpack the buffer
    for(int i=0; i<numstrings; i++) {

       // unpack: get string length and the corresponding string (termStr)
       // armci_hashmap_unpack() returns the string length. 
       index += armci_hashmap_unpack(&buf[index], termStr);
       
       if(index > bufsize)
         ARMCI_Error("GPCHashmap::insert() failed. Buffer overflow.", 0);

       // add the term to the hashmap
       VocabIntMap::const_iterator iter = vocabMap->find(termStr);
       int termID = -1;       
       if (iter != vocabMap->end()) {
          // term already in map
          termID = (*iter).second;
       }
       else {
          // new term. Add to vocab hashmap
          termID = vocabMap->size(); // starts with zero
          (*vocabMap)[termStr] = termID;
       }
    }
}

void armci_hashmap_insert2(VocabIntMap *vocabMap, const char *buf,
                           size_t bufsize, int *globalTermIds, int op) {

    string termStr;
    size_t index=0;

    // get the number of strings to be inserted (i.e.first field in the
    // buffer) and increment the index accordingly.
    int numstrings = *((int*)buf);
    index += sizeof(int);

    // unpack the buffer
    for(int i=0; i<numstrings; i++) {

       // unpack: get string length and the corresponding string (termStr)
       // armci_hashmap_unpack() returns the string length. 
       index += armci_hashmap_unpack(&buf[index], termStr);
       
       if(index > bufsize)
         ARMCI_Error("armci_hashmap_insert2() failed. Buffer overflow.", 0);

       if(op==HASHMAP_INSERT) { // insert the term to the hashmap
         (*vocabMap)[termStr] = globalTermIds[i];
       } else if(op==HASHMAP_GET){ // retrieve a term's global id from hashmap
         globalTermIds[i] = (*vocabMap)[termStr];
       } else {
         ARMCI_Error("armci_hashmap_insert2() Invalid operation", 0);
       }
    }
}
