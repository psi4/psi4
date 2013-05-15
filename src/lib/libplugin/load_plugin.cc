/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <boost/regex.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/xpressive/regex_actions.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem.hpp>

#include "plugin.h"

#include <libparallel/parallel.h>
#include <libchkpt/chkpt.hpp>
#include <libpsio/psio.hpp>

using namespace boost;

namespace psi {
#ifdef HAVE_DLFCN_H

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

        std::string msg = "load_plugin: Cannot find symbol: init_plugin; ";
        msg += dlsym_error1;
        throw PSIEXCEPTION(msg);
    }

    info.read_options = (read_options_t) dlsym(info.plugin_handle, "read_options");
    const char *dlsym_error2 = dlerror();
    if (dlsym_error2) {
        dlclose(info.plugin_handle);
        std::string msg = "load_plugin: Cannot find symbol: read_options; ";
        msg += dlsym_error2;
        throw PSIEXCEPTION(msg);
    }

    boost::filesystem::path pluginPath(plugin_pathname);
    boost::filesystem::path pluginStem = pluginPath.stem();
    info.name = pluginStem.string();

    // Modify info.name converting things that are allowed
    // filename characters to allowed C++ function names.
    std::string format_underscore("_");
    // Replace all '-' with '_'
    xpressive::sregex match_format = xpressive::as_xpr("-");
    info.name = regex_replace(info.name, match_format, format_underscore);

    info.plugin = (plugin_t) dlsym(info.plugin_handle, info.name.c_str());
    const char *dlsym_error3 = dlerror();
    if (dlsym_error3) {
        dlclose(info.plugin_handle);
        std::string msg = "load_plugin: Cannot find symbol: ";
        msg += info.name;
        msg += dlsym_error3;
        throw PSIEXCEPTION(msg);
    }

    // Store the name of the plugin for read_options
    boost::to_upper(info.name);

    // Get the plugin's options into the global space
    Process::environment.options.set_read_globals(true);
    info.read_options(info.name, Process::environment.options);
    Process::environment.options.set_read_globals(false);

    // Tell the plugin to initialize itself (found in libplugin)
    info.init_plugin();

    return info;
}

#else

plugin_info plugin_load(std::string& plugin_path)
{
    throw PSIEXCEPTION("Plugins are not supported on your platform.\n");
    return plugin_info();
}

#endif

}

