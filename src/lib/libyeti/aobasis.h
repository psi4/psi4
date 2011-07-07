#ifndef yeti_aobasis_h
#define yeti_aobasis_h

#include "class.h"

#include "index.hpp"
#include "aobasis.hpp"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

class AtomShell :
    public smartptr::Countable
{
    private:
        AtomBasisPtr atom_basis_;

        ShellBasisPtr shell_basis_;

    public:
        typedef std::vector<uli>::const_iterator iterator;

        AtomShell(
            const AtomBasisPtr& atom_basis,
            const ShellBasisPtr& shell_basis
        );

        uli nshell() const;

        iterator begin() const;

        iterator end() const;

};

class ShellBasis :
    public smartptr::Countable
{
    private:
        std::vector<uli> shell_sizes_;

        usi ang_mom_;

    public:
        typedef std::vector<uli>::const_iterator iterator;

        ShellBasis(usi ang_mom);

        void add_shell(uli size);

        usi angular_momentum() const;

        uli nshell() const;

        iterator begin() const;

        iterator end() const;

};

class AtomBasis :
    public smartptr::Countable
{

    private:
        std::vector<ShellBasisPtr> shells_;

        std::string symbol_;

        std::string name_;

        uli nbasis_;

        uli nshell_;

        double x_;

        double y_;

        double z_;

    public:
        typedef std::vector<ShellBasisPtr>::const_iterator iterator;

        AtomBasis(
            const std::string& symbol,
            const std::string& name,
            double x,
            double y,
            double z
        );

        void add_shell(uli size, usi ang_mom);

        iterator begin() const;

        iterator end() const;

        uli nbasis() const;

        uli nshell() const;

        double x() const;

        double y() const;

        double z() const;

        const std::string& symbol() const;

        const std::string& name() const;

};

class AOBasis :
    public smartptr::Countable
{

    private:
        std::vector<AtomBasisPtr> atoms_;

        std::vector<AtomShellPtr> shells_;

        uli nbasis_;

        std::string id_;

        std::string name_;

    public:
        AOBasis(
            const std::string& id,
            const std::string& name
        );

        void add_atom(const AtomBasisPtr& atom);

        IndexRange* get_index_range(uli nranges) const;

        uli nbasis() const;
};

}

#ifdef redefine_size_t
#undef size_t
#endif

#endif


