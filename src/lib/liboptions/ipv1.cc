#include <vector>
#include <cstdlib>
#include <cstring>
#include <algorithm>

#include <libipv1/ip_lib.h>
#include <libutil/libutil.h>

#include "liboptions.hpp"

#include <psi4-dec.h>

namespace psi {
  
  // Read all options from IPV1
  void Options::read_ipv1()
  {
    // Walk through all the options and attempt to read in all the information for it
    // from IPV1.
    for (Options::iterator pos = keyvals_.begin(); pos != keyvals_.end(); ++pos) {
      // Determine the type of the value
      if (pos->second.type() == "double")
        read_double(pos->second, pos->first);
      else if (pos->second.type() == "int")
        read_int(pos->second, pos->first);
      else if (pos->second.type() == "boolean")
        read_boolean(pos->second, pos->first);
      else if (pos->second.type() == "array")
        read_array(pos->second, pos->first);
      else
        throw OptionsException("Unknown data type.");
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
      err_code = ip_data(const_cast<char*>(key.c_str()), const_cast<char*>("%f"), &double_value, m);
    else
      err_code = ip_data(const_cast<char*>(key.c_str()), const_cast<char*>("%f"), &double_value, m, n);
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
  
  void Options::read_array(Data& data, const std::string& key)
  {
    // Make sure ipv1 and data agree on the size of the array
    int ipv1_size;
    int err_code = ip_count(const_cast<char*>(key.c_str()), &ipv1_size, 0);
    if (err_code != IPE_OK)
      return;
      
    if (data.size() != static_cast<unsigned int>(ipv1_size)) 
      throw OptionsException("Options and ipv1 array size for " + key + " do not match!");
    
    // Walk through the array and convert elements
    for (unsigned int i = 0; i < data.size(); ++i) {
      // Check the type of the array element. Only convert int, double, and boolean. Throw on others.
      if (data[i].type() == "double")
        read_double(data[i], key, 1, i);
      else if (data[i].type() == "int")
        read_int(data[i], key, 1, i);
      else if (data[i].type() == "boolean")
        read_boolean(data[i], key, 1, i);
      else if (data[i].type() == "array")
        throw OptionsException("Don't nest an array if you expect to use Options. (" + key + ")");
      else if (data[i].type() == "map")
        throw OptionsException("Don't nest a map if you expect to use Options. (" + key + ")");
      else
        throw OptionsException("Unknown data type detect in Options. (" + data[i].type() + ")");
    }
  }
}