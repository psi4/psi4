#ifndef yeti_OPIMPL_H
#define yeti_OPIMPL_H

#include "dataimpl.h"
#include "class.h"
#include "elementop.h"

#include <libsmartptr/printstream.h>

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {
}

#ifdef redefine_size_t
#undef size_t
#endif

#endif // OPIMPL_H
