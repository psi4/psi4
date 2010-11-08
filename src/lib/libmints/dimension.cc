#include <string.h>
#include "dimension.h"

using namespace psi;

Dimension::Dimension(int n, const std::string &name)
    : name_(name), n_(n)
{
    blocks_ = new int[n_];
    memset(blocks_, 0, sizeof(int)*n_);
}

Dimension::Dimension(const Dimension &other)
{
    if (n_ < other.n_) {
        delete blocks_;
        blocks_ = new int[other.n_];
    }
    n_ = other.n_;
    for (int i=0; i<n_; ++i)
        blocks_[i] = other.blocks_[i];
}

Dimension::~Dimension()
{
    delete[] blocks_;
}

void Dimension::print(FILE *out) const
{
    fprintf(outfile, "  %s (n = %d): ", name_.c_str(), n_);
    for (int i=0; i<n(); ++i) {
        fprintf(outfile, "%d  ", blocks_[i]);
    }
    fprintf(outfile, "\n");
}

Dimension& Dimension::operator =(const Dimension& other)
{
    if (n_ < other.n_) {
        delete blocks_;
        blocks_ = new int[other.n_];
    }
    n_ = other.n_;
    for (int i=0; i<n_; ++i)
        blocks_[i] = other.blocks_[i];

    return *this;
}
