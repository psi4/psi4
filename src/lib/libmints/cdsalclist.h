#ifndef CDSALCLIST_H
#define CDSALCLIST_H

#include <cstdio>

namespace boost {
template<class T>
class shared_ptr;
}

namespace psi {

extern FILE *outfile;

class Molecule;

class CdSalcList
{
    const boost::shared_ptr<Molecule>& molecule_;

public:
    CdSalcList(const boost::shared_ptr<Molecule>& mol);
};

} // namespace psi

#endif // CDSALCLIST_H
