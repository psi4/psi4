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

#include <boost/python.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <libplugin/plugin.h>
#include <libparallel/parallel.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>

#include <string>
#include <vector>

using namespace boost::python;
using namespace psi;

extern void py_psi_prepare_options_for_module(const std::string &name);

std::map<std::string, plugin_info> plugins;

/**************************************************************************
 * Plug-In functions                                                      *
 **************************************************************************/

/**
    Python interface for loading plugins.

        Python:
            plugin_load("integrals.so")

    @param fullpathname Full path and filename of the plugin to load.
    @returns 0 if not loaded, 1 if loaded, 2 if already loaded.
*/

int py_psi_plugin_load(std::string fullpathname)
{
    int ret = 0;

    boost::filesystem::path pluginPath(fullpathname);
    boost::filesystem::path pluginStem = pluginPath.stem();
    std::string uc = boost::algorithm::to_upper_copy(pluginStem.string());

    // Make sure the plugin isn't already loaded.
    if (plugins.count(uc) == 0) {
        plugins[uc] = plugin_load(fullpathname);
        psi::fprintf(outfile, "%s loaded.\n", fullpathname.c_str());
        ret = 1;
    }
    else
        ret = 2;

    return ret;
}

/**
    Python interface for calling plugins.

        Python:
            plugin("integrals.so")

    @param fullpathname Used to identity loaded plugin.
    @returns The result from the plugin.
*/
int py_psi_plugin(std::string fullpathname)
{
    boost::filesystem::path pluginPath(fullpathname);
    boost::filesystem::path pluginStem = pluginPath.stem();
    std::string uc = boost::algorithm::to_upper_copy(pluginStem.string());
    if (plugins.count(uc) == 0) {
        plugins[uc] = plugin_load(fullpathname);
    }
    plugin_info& tmpinfo = plugins[uc];
//    Process::environment.options.set_current_module(name);
    psi::fprintf(outfile, "Reading options from the %s block\n", tmpinfo.name.c_str());
    py_psi_prepare_options_for_module(tmpinfo.name);
    fflush(outfile);
    tmpinfo.read_options(tmpinfo.name, Process::environment.options);

    plugin_info& info = plugins[uc];

    // Call the plugin
    // Should be wrapped in a try/catch block.
    psi::fprintf(outfile, "Calling plugin %s.\n", fullpathname.c_str());
    fflush(outfile);

    // Have the plugin copy the environment to get current options.
    info.init_plugin();

    // Call the plugin
    int ret = info.plugin(Process::environment.options);

    return ret;
}

/**
    Python interface for closing plugin.

        Python:
            plugin_close("integrals.so")

    @param fullpathname Used to identity loaded plugin.
*/
void py_psi_plugin_close(std::string fullpathname)
{
    boost::filesystem::path pluginPath(fullpathname);
    boost::filesystem::path pluginStem = pluginPath.stem();
    std::string uc = boost::algorithm::to_upper_copy(pluginStem.string());
    if (plugins.count(uc) > 0) {
        plugin_info& info = plugins[uc];
        plugin_close(info);
        plugins.erase(uc);
    }
}

/**
    Python interface for closing all plugins.

        Python:
            plugin_close_all()
*/
void py_psi_plugin_close_all()
{
    std::map<std::string, plugin_info>::const_iterator iter = plugins.begin();

    for (; iter != plugins.end(); ++iter)
        plugin_close(plugins[iter->first]);

    plugins.clear();
}

/**************************************************************************
 * End of Plug-In functions                                               *
 **************************************************************************/

void export_plugins()
{
    // plugins
    def("plugin_load",      py_psi_plugin_load, "docstring");
    def("plugin",           py_psi_plugin, "docstring");
    def("plugin_close",     py_psi_plugin_close, "docstring");
    def("plugin_close_all", py_psi_plugin_close_all, "docstring");
}
