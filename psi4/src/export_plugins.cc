/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include <string>
#include <vector>

#include "psi4/pybind11.h"

#include "psi4/libfilesystem/path.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libplugin/plugin.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/process.h"

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

int py_psi_plugin_load(std::string fullpathname) {
    int ret = 0;

    filesystem::path pluginPath(fullpathname);
    std::string uc = to_upper_copy(pluginPath.stem());

    // Make sure the plugin isn't already loaded.
    if (plugins.count(uc) == 0) {
        plugins[uc] = plugin_load(fullpathname);
        // outfile->Printf("%s loaded.\n", fullpathname.c_str());
        ret = 1;
    } else {
        ret = 2;
    }
    return ret;
}

/**
    Python interface for calling plugins.

        Python:
            plugin("integrals.so")

    @param fullpathname Used to identity loaded plugin.
    @returns The result from the plugin.
*/
SharedWavefunction py_psi_plugin(std::string fullpathname, SharedWavefunction ref_wfn) {
    filesystem::path pluginPath(fullpathname);
    std::string uc = to_upper_copy(pluginPath.stem());
    if (plugins.count(uc) == 0) {
        plugins[uc] = plugin_load(fullpathname);
    }
    plugin_info &tmpinfo = plugins[uc];

    outfile->Printf("\nReading options from the %s block\n", tmpinfo.name.c_str());
    py_psi_prepare_options_for_module(tmpinfo.name);

    tmpinfo.read_options(tmpinfo.name, Process::environment.options);

    plugin_info &info = plugins[uc];

    // Call the plugin
    // Should be wrapped in a try/catch block.
    outfile->Printf("Calling plugin %s.\n\n\n", fullpathname.c_str());

    outfile->Printf(
        "Plugins that use gradients: set Da, Db, and Lagrangian for gradient theory on the wavefunction. The old way of passing these will stop working "
        "as soon as 1.8.");
    // Call the plugin
    if (ref_wfn) {
        return info.plugin(ref_wfn, Process::environment.options);
    } else {
        throw PSIEXCEPTION("Psi4::plugin: No wavefunction passed into the plugin, aborting");
    }
}

/**
    Python interface for closing plugin.

        Python:
            plugin_close("integrals.so")

    @param fullpathname Used to identity loaded plugin.
*/
void py_psi_plugin_close(std::string fullpathname) {
    filesystem::path pluginPath(fullpathname);
    std::string uc = to_upper_copy(pluginPath.stem());
    if (plugins.count(uc) > 0) {
        plugin_info &info = plugins[uc];
        plugin_close(info);
        plugins.erase(uc);
    }
}

/**
    Python interface for closing all plugins.

        Python:
            plugin_close_all()
*/
void py_psi_plugin_close_all() {
    std::map<std::string, plugin_info>::const_iterator iter = plugins.begin();

    for (; iter != plugins.end(); ++iter) plugin_close(plugins[iter->first]);

    plugins.clear();
}

/**************************************************************************
 * End of Plug-In functions                                               *
 **************************************************************************/

void export_plugins(py::module &m) {
    // plugins
    m.def("plugin_load", py_psi_plugin_load,
          "Load the plugin of name arg0. Returns 0 if not loaded, 1 if loaded, 2 if already loaded");
    m.def("plugin", py_psi_plugin, "Call the plugin of name arg0. Returns the plugin code result.");
    m.def("plugin_close", py_psi_plugin_close, "Close the plugin of name arg0.");
    m.def("plugin_close_all", py_psi_plugin_close_all, "Close all open plugins.");
}
