#ifndef WRITER_H
#define WRITER_H

#include <boost/shared_ptr.hpp>
#include <string>

namespace psi {

class Molecule;
class Matrix;
class BasisSet;
class Wavefunction;
class Options;

class GradientWriter
{
    boost::shared_ptr<Molecule> molecule_;
    const Matrix& gradient_;

public:
    GradientWriter(boost::shared_ptr<Molecule> mol, const Matrix& grad);

    void write(const std::string& filename);
};

class MoldenWriter
{
    boost::shared_ptr<Wavefunction> wavefunction_;

public:
    MoldenWriter(boost::shared_ptr<Wavefunction> wavefunction);

    void write(const std::string& filename);
};

class MOWriter
{
    boost::shared_ptr<Wavefunction> wavefunction_;
    Options & options_;
private:
    double * Ca_pointer, * eps;
    int * map, * sym, * occ, nmo, nso;
    void write_mos(Molecule & mol);
public:
    MOWriter(boost::shared_ptr<Wavefunction> wavefunction,Options&options);
    void write();
};

class NBOWriter
{
    boost::shared_ptr<Wavefunction> wavefunction_;

public:
    NBOWriter(boost::shared_ptr<Wavefunction> wavefunction);

    void write(const std::string &filename);
};

}

#endif // WRITER_H
