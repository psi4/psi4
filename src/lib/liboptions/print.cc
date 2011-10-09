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

#include <exception.h>
#include <libutil/libutil.h> // Needed for Ref counting, string splitting, and conversions
#include <libutil/ref.h> // Needed for Ref counting, string splitting, and conversions
#include <boost/shared_ptr.hpp>

#include "liboptions.h"
#include <psi4-dec.h>

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
    fprintf(outfile, "\n\n  Options:");
    fprintf(outfile, "\n  ----------------------------------------------------------------------------\n");
    fprintf(outfile, "%s\n", list.c_str());
    fflush(outfile);
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
    fprintf(outfile, "\n\n  Global Options:");
    fprintf(outfile, "\n  ----------------------------------------------------------------------------\n");
    fprintf(outfile, "%s\n", list.c_str());
  }

}
