#ifndef WRITER_H
#define WRITER_H

#include <boost/shared_ptr.hpp>
#include <string>

namespace psi {

class Molecule;
class SimpleMatrix;

class GradientWriter
{
    boost::shared_ptr<Molecule> molecule_;
    const SimpleMatrix& gradient_;

public:
    GradientWriter(boost::shared_ptr<Molecule> mol, const SimpleMatrix& grad);

    void write(const std::string& filename);
};

}

#endif // WRITER_H
