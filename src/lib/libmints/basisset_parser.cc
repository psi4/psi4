#include <cstdio>
#include "basisset_parser.h"

using namespace psi;
using namespace boost;

BasisSetParser::BasisSetParser(const std::string& filename)
   : filename_(filename), file_(NULL)
{
    if ((file_ = fopen(filename.c_str(), "r")) == NULL)
        throw PSIEXCEPTION("BasisSetParser::BasisSetParser: Unable to open basis set file.");
}

BasisSetParser::~BasisSetParser()
{
    fclose(file_);
}

void Gaussian94BasisSetParser::parse(shared_ptr<BasisSet>& basisSet)
{

}

