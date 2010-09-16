#include "dimension.h"

using namespace psi;

Dimension::Dimension(int n, const std::string &name)
    : name_(name), n_(n)
{
    blocks_ = new int[n_];
    memset(blocks_, 0, sizeof(int)*n_);
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
