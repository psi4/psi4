#ifndef yeti_aobasis_h
#define yeti_aobasis_h

#include "class.h"

#include "index.hpp"
#include "aobasis.hpp"

#include "gigmatrix.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

class Shell :
    public smartptr::Countable
{

    public:
        typedef enum { valence, diffuse } diffuse_type_t;

    private:
        AtomPtr atom_;

        usi ang_mom_;

        uli nfxn_;

        uli nsets_;

        uli program_shell_number_;

        uli program_fxn_start_;

        double min_exponent_;

        diffuse_type_t diffuse_type_;

        static bool statics_done_;

        void init_statics();

        static std::map<std::string, std::map<usi, double> > diffuse_cutoffs_;

    public:
        /**
          @param atom
          @param ang_mom
          @param nsets  The number of sets of functions.  For example, you could have a DZ basis with
                        with 2 sets of p-functions combined in a single shell
          @param nfxn   The total number of basis functions.  For example, you could have a DZ basis
                        with 2 sets of p-functions or 6 total functions
          @param min_exponent   The minimum Gaussian exponent (or minimum ``average'' exponent if
                                a contracted function) for all sets of functions in the shell
          @param program_shell_number   The shell number as stored by the interfacing integral program
        */
        Shell(
            const AtomPtr& atom,
            usi ang_mom,
            uli nsets,
            uli nfxn,
            double min_exponent,
            uli program_shell_number,
            uli program_fxn_start
        );


        diffuse_type_t diffuse_type() const;

        double min_exponent() const;

        usi angular_momentum() const;

        uli nfxn() const;

        const AtomPtr& get_atom() const;

        uli program_shell_number() const;

        uli program_fxn_start() const;

        void print(std::ostream& os = std::cout) const;

};

class Atom :
    public smartptr::Countable
{

    private:
        std::string symbol_;

        std::string name_;

        double x_;

        double y_;

        double z_;

        uli number_;

        uli group_number_;

    public:

        Atom(
            const std::string& symbol,
            const std::string& name,
            double x,
            double y,
            double z,
            uli atom_number
        );

        double x() const;

        double y() const;

        double z() const;

        const std::string& symbol() const;

        const std::string& name() const;

        uli number() const;

        uli group_number() const;

        void set_group_number(uli num);

};

class AOBasis :
    public smartptr::Countable
{

    private:
        std::vector<ShellPtr> shells_;

        std::vector<AtomPtr> atoms_;

        uli nbasis_;

        uli nfxn_min_per_tile_;

        uli* yeti_to_program_shell_map_;

        uli* yeti_to_program_fxn_map_;

        uli* program_to_yeti_shell_map_;

        uli* program_to_yeti_fxn_map_;

        std::string name_;

        std::string descr_;

        IndexRangePtr range_;

        IndexDescrPtr index_descr_;

        MultiShellMapPtr multi_shell_map_;

        void configure_atom_groups_as_atoms();

        void configure_atom_groups_by_number(uli natoms_per_group);

        void configure_many_atoms_large_basis(uli natoms_per_group);

    public:
        AOBasis(const std::string& name, const std::string& descr);

        void add_shell(const ShellPtr& shell);

        void add_atom(const AtomPtr& atom);

        void configure_many_atoms_small_basis();

        void configure_few_atoms_large_basis();

        void configure_many_atoms_large_basis();

        void configure_few_atoms_small_basis();

        void configure_min_nfxn_per_tile(uli nmin);

        void convert_to_yeti_numbering(RectMatrixPtr& coefs);

        IndexRange* get_index_range() const;

        IndexDescr* get_index_descr() const;

        const MultiShellMapPtr& get_multishell_map() const;

        const std::string& name() const;

        const std::string& descr() const;

        uli program_shell_number(uli yeti_shell_number);

        uli program_fxn_number(uli yeti_fxn_number);

        uli yeti_shell_number(uli program_shell_number);

        uli yeti_fxn_number(uli program_fxn_number);

        uli nbasis() const;
};

class MultiShellMap :
    public smartptr::Countable
{

    private:
        uli* shell_sizes_;

        uli* shell_start_on_index_;

        uli* nshells_per_index_;

        uli* nfxn_per_index_;

    public:
        MultiShellMap(
            const std::vector<ShellPtr>& shells,
            const IndexDescrPtr& ao_descr
        );

        uli shell_size(uli shell);

        uli shell_start(uli idx);

        uli nshell(uli idx);

        uli nfxn(uli idx);

};


}

#ifdef redefine_size_t
#undef size_t
#endif

#endif


