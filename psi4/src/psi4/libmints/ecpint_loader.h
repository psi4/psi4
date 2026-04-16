#ifndef LIBMINTS_ECPINT_LOADER_H
#define LIBMINTS_ECPINT_LOADER_H

#ifdef USING_ecpint_RUNTIME

#include <string>

namespace psi {
namespace ecpint_runtime {

// Check if libecpint is available at runtime
// Note: First call attempts to load the library (lazy initialization)
bool is_available();

// Get informative error message with install instructions
const char* get_unavailable_message();

// Get path to loaded library (empty if not loaded)
std::string get_library_path();

// Explicitly attempt to load (called by is_available on first use)
bool try_load();

// Unload library (called at psi4 shutdown)
void unload();

} // namespace ecpint_runtime
} // namespace psi

#endif // USING_ecpint_RUNTIME
#endif // header guard
