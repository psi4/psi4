#ifndef _psi_src_lib_libmints_dimension_h_
#define _psi_src_lib_libmints_dimension_h_

#include <string>
#include <cstdio>

namespace psi {

extern FILE *outfile;

class Dimension
{
    std::string name_;
    int n_;
    int *blocks_;

    Dimension();
public:
    Dimension(int n, const std::string& name = "");
    ~Dimension();

    /// Return the dimension
    int n() const { return n_; }

    /// Return the name of the dimension.
    const std::string& name() const { return name_; }

    int& operator[](int i) { return blocks_[i]; }
    const int& operator[](int i) const { return blocks_[i]; }

    void print(FILE* out=outfile) const;
};

}

#endif // _psi_src_lib_libmints_dimension_h_
