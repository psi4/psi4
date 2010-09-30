#include <psi4-dec.h>
#include <psiconfig.h>
#include <libmints/molecule.h>
#include <boost/algorithm/string.hpp>

//MKL Header
//#define _MKL
#ifdef _MKL
#include <mkl.h>
#endif

//OpenMP Header
//_OPENMP is defined by the compiler if it exists
#ifdef _OPENMP
#include <omp.h>
#endif


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

    nthread_ = 1;

    #ifdef _OPENMP
        nthread_ = omp_get_max_threads(); 
    #endif

    // For testing, print out the enviroment:
//    for (map<string, string>::const_iterator iter = environment_.begin();
//         iter != environment_.end(); ++iter)
//        printf("Key: %s Value: %s\n", iter->first.c_str(), iter->second.c_str());


    // If madness and or MPI is not set up, COMMUNICATOR is changed to a value
    // that makes sense. (i.e. MPI or LOCAL)
#if HAVE_MADNESS == 0
    #if HAVE_MPI == 1
        if (Process::environment("COMMUNICATOR") != "MPI") {
            environment_["COMMUNICATOR"] = "MPI";
            std::cout << "COMMUNICATOR was changed to MPI" << std::endl;
        }
    #else
        if (Process::environment("COMMUNICATOR") != "LOCAL") {
            environment_["COMMUNICATOR"] = "LOCAL";
            std::cout << "COMMUNICATOR was changed to LOCAL" << std::endl;
        }
    #endif
#else
    if(Process::environment("COMMUNICATOR") != "MPI" && Process::environment("COMMUNICATOR") != "MADNESS") {
        environment_["COMMUNICATOR"] = "MPI";
        std::cout << "COMMUNICATOR was changed to MPI" << std::endl;
    }
#endif


}

void Process::Environment::set_n_threads(int nthread)
{
    nthread_ = nthread;
    #ifdef _OPENMP
        omp_set_num_threads(nthread_);        
    #endif
    #ifdef _MKL
        mkl_set_num_threads(nthread_);
    #endif
}

const string& Process::Environment::operator()(const string& key) const
{
    // Search for the key:
    map<string, string>::const_iterator it = environment_.find(key);

    if (it == environment_.end()) {
        return empty_;      // Not found return empty string.
    }
    else
        return it->second;  // Found, return the value
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

void Process::Environment::set_molecule(boost::shared_ptr<Molecule> molecule)
{
  molecule_ = molecule;
}

boost::shared_ptr<Molecule> Process::Environment::molecule() const
{
  return molecule_;
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
