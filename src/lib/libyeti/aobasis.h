#ifndef yeti_aobasis_h
#define yeti_aobasis_h

#include "class.h"

#include "index.hpp"
#include "aobasis.hpp"


namespace yeti {

class AtomBasis : public smartptr::Countable {

    private:
        std::vector<size_t> shell_sizes_;

        std::string symbol_;

        std::string name_;

        size_t nbasis_;

        size_t nshell_;

    public:
        typedef std::vector<size_t>::const_iterator iterator;

        AtomBasis(
            const std::string& symbol,
            const std::string& name
        );

        void add_shell(size_t size);

        iterator begin() const;

        iterator end() const;

        uli nbasis() const;

        uli nshell() const;

        const std::string& symbol() const;

        const std::string& name() const;

};

class AOBasis : public smartptr::Countable {

    private:
        std::vector<AtomBasisPtr> atoms_;

        size_t nbasis_;

        std::string id_;

        std::string name_;

    public:
        AOBasis(
            const std::string& id,
            const std::string& name
        );

        void add_atom(const AtomBasisPtr& atom);

        IndexRange* get_index_range() const;

        size_t nbasis() const;
};

}


#endif


