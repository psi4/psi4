#include <boost/algorithm/string.hpp>
#include "strop.h"

void
smartptr::split(std::vector<std::string>& values, const std::string& line, const char* tokens)
{
    boost::split(values, line, boost::is_any_of(tokens), boost::algorithm::token_compress_on);
}
