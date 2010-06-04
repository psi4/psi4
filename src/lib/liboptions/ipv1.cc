#include <vector>
#include <cstdlib>
#include <cstring>
#include <algorithm>

#include <libipv1/ip_lib.h>
#include <libutil/libutil.h>

#include "liboptions.h"

#include <psi4-dec.h>

namespace psi {

  // Read all options from IPV1
  void Options::read_ipv1()
  {
    // Walk through all the options and attempt to read in all the information for it
    // from IPV1.
    for (Options::iterator pos = keyvals_.begin(); pos != keyvals_.end(); ++pos) {
      try {
            // Determine the type of the value
            if (pos->second.type() == "double")
              read_double(pos->second, pos->first);
            else if (pos->second.type() == "int")
              read_int(pos->second, pos->first);
            else if (pos->second.type() == "boolean")
              read_boolean(pos->second, pos->first);
            else if (pos->second.type() == "array")
              read_array(pos->second, pos->first);
            else if (pos->second.type() == "string")
              read_string(pos->second, pos->first);
            else
              throw OptionsException("Unknown data type. [type() == " + pos->second.type() + "]");
      }
      catch (PsiException e) {
        fprintf(stderr, "Key: %s\n%s\n", pos->first.c_str(), e.what());
        abort();
      }
    }
  }

  void Options::read_boolean(Data& data, const std::string& key, int m, int n)
  {
    int int_value = 0;
    int err_code;
    if (m == 0)
      err_code = ip_boolean(const_cast<char*>(key.c_str()), &int_value, m);
    else
      err_code = ip_boolean(const_cast<char*>(key.c_str()), &int_value, m, n);
    if (err_code == IPE_OK)
      data.assign(int_value);
  }

  void Options::read_double(Data& data, const std::string& key, int m, int n)
  {
    double double_value = 0.0;
    int err_code;
    if (m == 0)
      err_code = ip_data(const_cast<char*>(key.c_str()), const_cast<char*>("%lf"), &double_value, m);
    else
      err_code = ip_data(const_cast<char*>(key.c_str()), const_cast<char*>("%lf"), &double_value, m, n);
    if (err_code == IPE_OK)
      data.assign(double_value);
  }

  void Options::read_int(Data& data, const std::string& key, int m, int n)
  {
    int int_value = 0;
    int err_code;
    if (m == 0)
      err_code = ip_data(const_cast<char*>(key.c_str()), const_cast<char*>("%d"), &int_value, m);
    else
      err_code = ip_data(const_cast<char*>(key.c_str()), const_cast<char*>("%d"), &int_value, m, n);
    if (err_code == IPE_OK)
      data.assign(int_value);
  }

  void Options::read_string(Data& data, const std::string& key, int m, int n)
  {
    char *value = NULL;
    int err_code;
    if (m == 0)
      err_code = ip_string(const_cast<char*>(key.c_str()), &value, m);
    else
      err_code = ip_string(const_cast<char*>(key.c_str()), &value, m, n);
    if (err_code == IPE_OK && value != NULL)
      data.assign(std::string(value));
    if (value != NULL)
      free(value);
    value = NULL;
  }

  void Options::read_array(Data& data, const std::string& key)
  {
    // Make sure ipv1 and data agree on the size of the array
    int ipv1_size;
    int err_code = ip_count(const_cast<char*>(key.c_str()), &ipv1_size, 0);
    if (err_code != IPE_OK)
      return;

    data.changed();

    if (data.size() < static_cast<unsigned int>(ipv1_size)) {
      // Add enough StringDataTypes until it is the same size.
      unsigned int diff = static_cast<unsigned int>(ipv1_size) - data.size();
      for (unsigned int i = 0; i < diff; ++i)
        data.add(new StringDataType("", ""));
    }

    // Walk through the array and convert elements
    for (unsigned int i = 0; i < data.size(); ++i) {
      // Check the type of the array element. Only convert int, double, and boolean. Throw on others.
      if (data[i].type() == "double")
        read_double(data[i], key, 1, i);
      else if (data[i].type() == "int")
        read_int(data[i], key, 1, i);
      else if (data[i].type() == "boolean")
        read_boolean(data[i], key, 1, i);
      else if (data[i].type() == "string")
        read_string(data[i], key, 1, i);
      else if (data[i].type() == "array")
        throw OptionsException("Don't nest an array if you expect to use Options with IPV1. (" + key + ")");
      else if (data[i].type() == "map")
        throw OptionsException("Don't nest a map if you expect to use Options with IPV1. (" + key + ")");
      else
        throw OptionsException("Unknown data type detect in Options. (" + data[i].type() + ")");
    }
  }

  std::string Options::to_string() const {
    std::stringstream str;
    int linewidth = 0;
    int largest_key = 0, largest_value = 7;  // 7 for '(empty)'

    for (const_iterator pos = keyvals_.begin(); pos != keyvals_.end(); ++pos) {
        pos->first.size()  > largest_key   ? largest_key = pos->first.size() : 0;
        pos->second.to_string().size() > largest_value ? largest_value = pos->second.to_string().size() : 0;
    }

    for (const_iterator pos = keyvals_.begin(); pos != keyvals_.end(); ++pos) {
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

  void Options::print() {
    std::string list = to_string();
    fprintf(outfile, "\n\n  Options:");
    fprintf(outfile, "\n  ----------------------------------------------------------------------------\n");
    fprintf(outfile, "%s\n", list.c_str());
  }

}
