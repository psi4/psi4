#include "plugin.h"

namespace psi {

#if HAVE_DLFCN_H == 1

#include <dlfcn.h>

void plugin_close(const plugin_info& info)
{
    dlclose(info.plugin_handle);
}

#else

void plugin_close(const plugin_info& info)
{
}

#endif

}
