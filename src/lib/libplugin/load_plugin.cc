#include "plugin.h"
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

namespace psi {
#if HAVE_DLFCN_H == 1

#include <dlfcn.h>

plugin_info plugin_load(std::string& plugin_pathname)
{
    plugin_info info;

    info.plugin_handle = dlopen(plugin_pathname.c_str(), RTLD_LAZY);
    if (info.plugin_handle == NULL) {
        std::string msg = "load_plugin: Cannot open library: ";
        msg += dlerror();
        throw PSIEXCEPTION(msg.c_str());
    }

    info.init_plugin = (init_plugin_t) dlsym(info.plugin_handle, "init_plugin");
    const char *dlsym_error1 = dlerror();
    if (dlsym_error1) {
        dlclose(info.plugin_handle);

        std::string msg = "load_plugin: Cannot symbol: ";
        msg += dlsym_error1;
        throw PSIEXCEPTION(msg);
    }

    info.read_options = (read_options_t) dlsym(info.plugin_handle, "read_options");
    const char *dlsym_error2 = dlerror();
    if (dlsym_error2) {
        dlclose(info.plugin_handle);
        std::string msg = "load_plugin: Cannot symbol: ";
        msg += dlsym_error2;
        throw PSIEXCEPTION(msg);
    }

    info.name = boost::filesystem::basename(plugin_pathname);
    info.plugin = (plugin_t) dlsym(info.plugin_handle, info.name.c_str());
    const char *dlsym_error3 = dlerror();
    if (dlsym_error3) {
        dlclose(info.plugin_handle);
        std::string msg = "load_plugin: Cannot symbol: ";
        msg += dlsym_error3;
        throw PSIEXCEPTION(msg);
    }

    // Tell the plugin to initialize itself (found in libplugin)
    info.init_plugin(Communicator::world, Process::environment);

    // Store the name of the plugin for read_options
    boost::to_upper(info.name);

    return info;
}

#else

plugin_info plugin_load(std::string& plugin_path)
{
    fprintf(outfile, "Plugins are not supported on your platform.\n");
    return plugin_info();
}

#endif

}

