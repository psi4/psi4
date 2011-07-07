#ifndef yeti_utils_h
#define yeti_utils_h

#include "class.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

static void
quicksort(int* item, int* index, int n);

static void
quicksort(usi* item, usi* index, usi n);

}

#ifdef redefine_size_t
#undef size_t
#endif

#endif

