#ifndef smartptr_strop_h
#define smartptr_strop_h

#include <vector>
#include <string>

namespace smartptr {

void
split(std::vector<std::string>& values, const std::string& line, const char* tokens);

}

#endif

