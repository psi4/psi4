#include <psiconfig.h>

#include "plugin.h"

namespace psi {

#ifdef HAVE_DLFCN_H

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
