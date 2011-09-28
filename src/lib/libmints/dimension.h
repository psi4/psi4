#ifndef _psi_src_lib_libmints_dimension_h_
#define _psi_src_lib_libmints_dimension_h_

#include <string>
#include <cstdio>
#include <vector>

namespace psi {

extern FILE *outfile;

class Dimension
{
    std::string name_;
    int n_;
    int *blocks_;

public:
    Dimension();
    Dimension(const Dimension& other);
    Dimension(int n, const std::string& name = "");
    Dimension(const std::vector<int>& other);
    ~Dimension();

    /// Assignment operator
    Dimension& operator=(const Dimension& other);

    /// Assignment operator, this one can be very dangerous
    Dimension& operator=(const int* other);

    Dimension& operator+=(const Dimension& b);
    Dimension& operator-=(const Dimension& b);

    /// Return the dimension
    int n() const { return n_; }

    /// Return the name of the dimension.
    const std::string& name() const { return name_; }

    /// Blocks access
    int& operator[](int i) { return blocks_[i]; }
    const int& operator[](int i) const { return blocks_[i]; }

    /// Casting operator to int*
    operator int*() const { return blocks_; }
    /// Casting operator to const int*
    operator const int*() const { return blocks_; }

    /// Return the sum of constituent dimensions
    int sum() const;

    int* pointer() const { return blocks_; }

    void print(FILE* out=outfile) const;
};

bool operator==(const Dimension& a, const Dimension& b) {
    if (a.n() != b.n())
        return false;
    for (int i=0; i<a.n(); ++i)
        if (a[i] != b[i])
            return false;
    return true;
}

bool operator!=(const Dimension& a, const Dimension& b) {
    return !operator==(a, b);
}

Dimension operator+(const Dimension& a, const Dimension& b) {
    Dimension result = a;
    for (int i=0; i<a.n(); ++i)
        result[i] += b[i];
    return result;
}

Dimension operator-(const Dimension& a, const Dimension& b) {
    Dimension result = a;
    for (int i=0; i<a.n(); ++i)
        result[i] -= b[i];
    return result;
}

}

#endif // _psi_src_lib_libmints_dimension_h_
