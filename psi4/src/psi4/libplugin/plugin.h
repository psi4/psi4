/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#ifndef __psi4_src_lib_libplugin_plugin_h
#define __psi4_src_lib_libplugin_plugin_h


#include "psi4/libparallel/parallel.h"
#include "psi4/psi4-dec.h"
#include <string>

namespace psi {
class Chkpt;

class Communicator;

class PSIO;

// Useful typedef's
typedef SharedWavefunction (*plugin_t)(SharedWavefunction, Options&);
typedef int (*read_options_t)(std::string, Options&);

// DEPRECATED but left so that existing plugins will work.
#define INIT_PLUGIN

// The following lines are used by the PSI4 driver.
typedef struct plugin_info_struct
{
    // Uppercased plugin name.
    std::string name;

    // File handle to plugin
    void *plugin_handle;

    // Used to call the user provided plugin code
    plugin_t plugin;

    // Function pointer to the plugin's read_options function
    read_options_t read_options;
} plugin_info;

/**
  Loads a plugin from disk and calls its initialization routine.

  @param plugin_path Absolute or relative path and filename name to the plugin.
  @return A plugin_info structure containing function pointers to entry points in the plugin.
*/
plugin_info plugin_load(std::string& plugin_path);

/**
  Frees the plugin from memory.

  @param plugin A plugin_info struct that contains the plugin to free.
*/
void plugin_close(const plugin_info& plugin);
}

#endif
