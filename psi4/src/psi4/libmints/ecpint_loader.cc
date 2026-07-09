#ifdef USING_ecpint_RUNTIME

#include "ecpint_loader.h"
#include "psi4/libpsi4util/exception.h"

// dlopen/dlsym/dlclose are POSIX; ENABLE_ecpint_RUNTIME is only supported on
// Linux/macOS (the libecpint.so/.dylib names below are not meaningful on
// Windows), so gate on that rather than depending on a CMake-probed
// HAVE_DLFCN_H that isn't plumbed into this target.
#if !defined(_WIN32)
#define PSI4_ECPINT_HAVE_DLFCN 1
#include <dlfcn.h>
#endif

#include <array>
#include <mutex>
#include <string>

namespace psi {
namespace ecpint_runtime {

namespace {
    void* lib_handle_ = nullptr;
    bool load_attempted_ = false;
    std::string lib_path_;
    std::mutex load_mutex_;

    const char* unavailable_msg_ = 
        "ECP integrals requested but libecpint not available.\n"
        "  Install via: conda install libecpint\n"
        "  Or recompile Psi4 with -DENABLE_ecpint=ON";
}

bool try_load() {
    std::lock_guard<std::mutex> lock(load_mutex_);

    if (load_attempted_) {
        return lib_handle_ != nullptr;
    }
    load_attempted_ = true;

#ifdef PSI4_ECPINT_HAVE_DLFCN
    // libecpint is NOT linked into core.so (that's the whole point of runtime-
    // optionality), so its symbols are undefined in this shared object. We
    // actually dlopen the real library here, with RTLD_GLOBAL so that the
    // dynamic linker can find its symbols when later lazily binding the
    // PLT/JUMP_SLOT relocations left by the ECPInt/ECPSOInt code compiled
    // into core.so. This must be called (and must succeed) before any code
    // that touches libecpint:: symbols runs.
    static const std::array<const char*, 4> candidates{{
        "libecpint.so.1",
        "libecpint.so",
        "libecpint.1.dylib",
        "libecpint.dylib",
    }};

    for (const char* name : candidates) {
        void* handle = dlopen(name, RTLD_NOW | RTLD_GLOBAL);
        if (handle) {
            lib_handle_ = handle;
            lib_path_ = name;
            return true;
        }
    }
#endif

    return false;
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
#ifdef PSI4_ECPINT_HAVE_DLFCN
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
