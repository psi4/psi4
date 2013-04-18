#ifndef smartptr_printstream_h_
#define smartptr_printstream_h_

#include <sstream>
#include <stdio.h>
#include <iostream>
#include <map>

namespace std {

/**
    @param fmt A printf format string
    @param va_args The set of arguments
    @return A string with formatted values
*/
std::string 
stream_printf(const char* fmt, ...);

template <
  class Type
> void
print_keys(
    std::map<std::string, Type>& keymap,
    std::stringstream& sstr
)
{
    typename std::map<std::string, Type>::iterator it;
    for (it = keymap.begin(); it != keymap.end(); it++)
        sstr << it->first << std::endl;
}

/**
 @param keymap A string map
 @param os The stream to print the map to
*/
void
print_map(
    std::map<std::string, std::string>& keymap,
    std::ostream& os
);

}

#ifdef redefine_size_t
#undef size_t
#endif

#endif
