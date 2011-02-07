#include "writer.h"
#include <libmints/mints.h>

#include <cstdio>

using namespace psi;
using namespace boost;

GradientWriter::GradientWriter(shared_ptr<Molecule> mol, const SimpleMatrix& grad)
    : molecule_(mol), gradient_(grad)
{
}

void GradientWriter::write(const std::string &filename)
{
    FILE *file11 = fopen(filename.c_str(), "a");
    int i;

    if (file11 == NULL)
        throw PSIEXCEPTION("GradientWriter::write: Unable to open file11.dat");

    fprintf(file11, "%-59.59s %-10.10s%-8.8s\n",
            molecule_->name().c_str(),
            "(wfn)",
            "(dertype)");

    fprintf(file11, "%5d%20.10lf\n", molecule_->natom(), Process::environment.globals["CURRENT ENERGY"]);

    for (i=0; i<molecule_->natom(); ++i) {
        fprintf(file11, "%20.10lf%20.10lf%20.10lf%20.10lf\n",
                double(molecule_->Z(i)), molecule_->x(i), molecule_->y(i), molecule_->z(i));
    }

    for (i=0; i<molecule_->natom(); ++i) {
        fprintf(file11, "                    %20.10lf%20.10lf%20.10lf\n",
                gradient_(i, 0), gradient_(i, 1), gradient_(i, 2));
    }

    fclose(file11);
}
