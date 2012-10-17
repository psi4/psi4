#if !defined(psi_src_lib_libparallel_threaded_storage)
#define psi_src_lib_libparallel_threaded_storage

#include <vector>

namespace psi {

template<typename T>
class threaded_storage
{
    std::vector<T> storage_;

public:

    threaded_storage(const T& value = T())
        : storage_(WorldComm->nthread(), value)
    { }

    void initialize(const T& value) {
        storage_.clear();
        for (int i=0; i<WorldComm->nthread(); ++i) {
            storage_.push_back(value);
        }
    }

    T& operator*() {
        return storage_[WorldComm->thread_id(pthread_self())];
    }

    const T& operator[](size_t ind) const {
        return storage_[ind];
    }
    T& operator[](size_t ind) {
        return storage_[ind];
    }
};

}

#endif

