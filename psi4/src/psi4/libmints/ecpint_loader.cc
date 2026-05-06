#ifdef USING_ecpint_RUNTIME

#include "ecpint_loader.h"
#include "psi4/libpsi4util/exception.h"

#ifdef HAVE_DLFCN_H
#include <dlfcn.h>
#endif

#include <string>

namespace psi {
namespace ecpint_runtime {

namespace {
    void* lib_handle_ = nullptr;
    bool load_attempted_ = false;
    std::string lib_path_;
    
    const char* unavailable_msg_ = 
        "ECP integrals requested but libecpint not available.\n"
        "  Install via: conda install libecpint\n"
        "  Or recompile Psi4 with -DENABLE_ecpint=ON";
}

bool try_load() {
    if (load_attempted_) {
        return lib_handle_ != nullptr;
    }
    
    load_attempted_ = true;
    
    // In runtime mode, the library is preloaded by Python in __init__.py
    // If we got this far (module imported successfully), the library is available
    // Just mark it as loaded
    lib_path_ = "preloaded by Python";
    lib_handle_ = (void*)1;  // Non-null marker to indicate success
    return true;
}

bool is_available() {
    if (!load_attempted_) {
        try_load();
    }
    return lib_handle_ != nullptr;
}

const char* get_unavailable_message() {
    return unavailable_msg_;
}

std::string get_library_path() {
    return lib_path_;
}

void unload() {
#ifdef HAVE_DLFCN_H
    if (lib_handle_) {
        dlclose(lib_handle_);
        lib_handle_ = nullptr;
        lib_path_.clear();
    }
#endif
}

} // namespace ecpint_runtime
} // namespace psi

#endif // USING_ecpint_RUNTIME
