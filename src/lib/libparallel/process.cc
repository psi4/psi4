#include <psi4-dec.h>

#include <boost/algorithm/string.hpp>

using namespace std;
using namespace psi;
using namespace boost;

Process::Environment Process::environment;
Process::Arguments Process::arguments;

// Need to split each entry by the first '=', left side is key, right the value
void Process::Environment::init(char **envp)
{
    // First set some defaults:
    environment_["PSIDATADIR"] = INSTALLEDPSIDATADIR;

    // Go through user provided environment overwriting defaults if necessary
    int i=0;
    while (envp[i] != NULL) {
        std::vector<std::string> strs;
        boost::split(strs, envp[i], boost::is_any_of("="));
        environment_[strs[0]] = strs[1];
        ++i;
    }

    // For testing, print out the enviroment:
//    for (map<string, string>::const_iterator iter = environment_.begin();
//         iter != environment_.end(); ++iter)
//        printf("Key: %s Value: %s\n", iter->first.c_str(), iter->second.c_str());
}

const string& Process::Environment::operator()(const string& key) const
{
    // Search for the key:
    map<string, string>::const_iterator it = environment_.find(key);

    if (it == environment_.end())
        return string(); // Not found return empty string.
    else
        return it->second; // Found, return the value
}

string Process::Environment::operator()(const string& key)
{
    // Search for the key:
    map<string, string>::const_iterator it = environment_.find(key);

    if (it == environment_.end())
        return string(); // Not found return empty string.
    else
        return it->second; // Found, return the value
}

const string& Process::Environment::set(const std::string &key, const std::string &value)
{
    const string& old = operator()(key);
    environment_[key] = value;
    return old;
}

void Process::Arguments::init(int argc, char **argv)
{
    for (int i=0; i<argc; ++i)
        arguments_.push_back(argv[i]);
}

const string& Process::Arguments::operator ()(int argc) const
{
    if (argc < 0 || argc >= arguments_.size())
        throw SanityCheckError("Process:Arguments::operator ()(int argc): argc is either < 0 or >= the number of arguments.\n", __FILE__, __LINE__);

    return arguments_[argc];
}

string Process::Arguments::operator ()(int argc)
{
    if (argc < 0 || argc >= arguments_.size())
        throw SanityCheckError("Process:Arguments::operator ()(int argc): argc is either < 0 or >= the number of arguments.\n", __FILE__, __LINE__);

    return arguments_[argc];
}
