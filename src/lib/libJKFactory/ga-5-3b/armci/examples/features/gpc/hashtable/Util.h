#ifndef UTIL_H
#define UTIL_H

#include <string>
using std::string;
#include "Hash_common.h"

extern int armci_hashmap_pack(char *buf, string str);
extern int armci_hashmap_unpack(const char *buf, string& str);
extern void armci_hashmap_insert(VocabIntMap *vocabMap,
                                 const char *buf, size_t bufsize);
extern void armci_hashmap_insert2(VocabIntMap *vocabMap, const char *buf,
                                  size_t bufsize, int *globalTermIds, int op);

#endif /* UTIL_H */
