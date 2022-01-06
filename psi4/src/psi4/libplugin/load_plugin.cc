/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include "plugin.h"

#include <algorithm>
#include <cctype>

#include "psi4/libfilesystem/path.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsi4util/libpsi4util.h"

namespace psi {
#ifdef HAVE_DLFCN_H

#include <dlfcn.h>

plugin_info plugin_load(std::string &plugin_pathname) {
    plugin_info info;

    info.plugin_handle = dlopen(plugin_pathname.c_str(), RTLD_LAZY);
    if (info.plugin_handle == nullptr) {
        std::string msg = "load_plugin: Cannot open library: ";
        msg += dlerror();
        throw PSIEXCEPTION(msg);
    }

    info.read_options = (read_options_t)dlsym(info.plugin_handle, "read_options");
    const char *dlsym_error2 = dlerror();
    if (dlsym_error2) {
        dlclose(info.plugin_handle);
        std::string msg = "load_plugin: Cannot find symbol: read_options; ";
        msg += dlsym_error2;
        throw PSIEXCEPTION(msg);
    }

    info.name = filesystem::path(plugin_pathname).stem();

    // Modify info.name converting things that are allowed
    // filename characters to allowed C++ function names.
    // Replace all '-' with '_'
    std::transform(info.name.begin(), info.name.end(), info.name.begin(), [](char c) { return c == '-' ? '_' : c; });

    info.plugin = (plugin_t)dlsym(info.plugin_handle, info.name.c_str());
    const char *dlsym_error3 = dlerror();
    if (dlsym_error3) {
        dlclose(info.plugin_handle);
        std::string msg = "load_plugin: Cannot find symbol: ";
        msg += info.name;
        msg += dlsym_error3;
        throw PSIEXCEPTION(msg);
    }

    // Uppercase and store the name of the plugin for read_options
    to_upper(info.name);

    // Get the plugin's options into the global space
    Process::environment.options.set_read_globals(true);
    info.read_options(info.name, Process::environment.options);
    Process::environment.options.set_read_globals(false);

    return info;
}

#else

plugin_info plugin_load(std::string& plugin_path) {
    throw PSIEXCEPTION("Plugins are not supported on your platform.\n");
    return plugin_info();
}

#endif
}
