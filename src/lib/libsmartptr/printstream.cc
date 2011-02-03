#include <stdarg.h>

#include "printstream.h"

using namespace std;

string
std::stream_printf(const char* fmt, ...)
{
    char str[1024];

    va_list args;
    va_start(args, fmt);
    vsprintf(str, fmt, args);
    va_end(args);

    string strobj = str;
    
    return strobj;
}

void
std::print_map(
    map<string, string>& keymap,
    ostream& os
)
{
    map<string, string>::iterator it;
    for (it = keymap.begin(); it != keymap.end(); ++it)
        os << stream_printf("%20s : %s", it->first.data(), it->second.data()) << endl;
}
