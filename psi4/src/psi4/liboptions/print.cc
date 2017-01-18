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

#include <iostream>
#include <vector>
#include <map>
#include <cstddef>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <assert.h>
 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/libpsi4util.h" // Needed for Ref counting, string splitting, and conversions
#include "psi4/libpsi4util/ref.h" // Needed for Ref counting, string splitting, and conversions


#include "liboptions.h"
#include "psi4/psi4-dec.h"

namespace psi {

  std::string Options::to_string() const {
    std::stringstream str;
    int linewidth = 0;
    int largest_key = 0, largest_value = 7;  // 7 for '(empty)'

    std::map<std::string, std::map<std::string, Data> >::const_iterator localmap = locals_.find(current_module_);
    if (localmap == locals_.end()) return str.str(); // Nothing to print
    const std::map<std::string, Data> &keyvals = localmap->second;
    for (const_iterator pos = keyvals.begin(); pos != keyvals.end(); ++pos) {
        pos->first.size()  > largest_key   ? largest_key = pos->first.size() : 0;
        pos->second.to_string().size() > largest_value ? largest_value = pos->second.to_string().size() : 0;
    }

    for (const_iterator local_iter = keyvals.begin(); local_iter != keyvals.end(); ++local_iter) {
      std::stringstream line;
      std::string value;
      bool option_specified;
      const std::string &name = local_iter->first;
      const_iterator global_iter = globals_.find(name);
      if(local_iter->second.has_changed()){
          // Local option was set, use it
          value = local_iter->second.to_string();
          option_specified = true;
      }else if(global_iter->second.has_changed()){
          // Global option was set, get that
          value = global_iter->second.to_string();
          option_specified = true;
      }else{
          // Just use the default local value
          value = local_iter->second.to_string();
          option_specified = false;
      }
      if (value.length() == 0) {
        value = "(empty)";
      }
      line << "  " << std::left << std::setw(largest_key) << local_iter->first << " => " << std::setw(largest_value) << value;
      if (option_specified)
        line << " !";
      else
        line << "  ";

      str << line.str();

      linewidth += line.str().size();
      if (linewidth + largest_key + largest_value + 8 > 80) {
          str << std::endl;
          linewidth = 0;
      }
    }


    return str.str();
  }

  void Options::print() {
    std::string list = to_string();
    outfile->Printf( "\n\n  Options:");
    outfile->Printf( "\n  ----------------------------------------------------------------------------\n");
    outfile->Printf( "%s\n", list.c_str());

  }

  std::string Options::globals_to_string() const {
    std::stringstream str;
    int linewidth = 0;
    int largest_key = 0, largest_value = 7;  // 7 for '(empty)'

    for (const_iterator pos = globals_.begin(); pos != globals_.end(); ++pos) {
        pos->first.size()  > largest_key   ? largest_key = pos->first.size() : 0;
        pos->second.to_string().size() > largest_value ? largest_value = pos->second.to_string().size() : 0;
    }

    for (const_iterator pos = globals_.begin(); pos != globals_.end(); ++pos) {
      std::stringstream line;
      std::string second_tmp = pos->second.to_string();
      if (second_tmp.length() == 0) {
        second_tmp = "(empty)";
      }
      line << "  " << std::left << std::setw(largest_key) << pos->first << " => " << std::setw(largest_value) << second_tmp;
      if (pos->second.has_changed())
        line << " !";
      else
        line << "  ";

      str << line.str();

      linewidth += line.str().size();
      if (linewidth + largest_key + largest_value + 8 > 80) {
          str << std::endl;
          linewidth = 0;
      }
    }


    return str.str();
  }

  void Options::print_globals() {
    std::string list = globals_to_string();
    outfile->Printf( "\n\n  Global Options:");
    outfile->Printf( "\n  ----------------------------------------------------------------------------\n");
    outfile->Printf( "%s\n", list.c_str());
  }

  std::vector<std::string> Options::list_globals() {
    std::vector<std::string> glist(globals_.size());
    int ii = 0;

    for (const_iterator pos = globals_.begin(); pos != globals_.end(); ++pos) {
      glist[ii] = pos->first;
      ii++;
    }
    return glist;
  }

}
