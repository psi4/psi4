#ifndef _psi_src_lib_libmints_symmetry_h_
#define _psi_src_lib_libmints_symmetry_h_

#include <cstdio>
#include <libchkpt/chkpt.hpp>

#include <libmints/molecule.h>

namespace psi {

    class SymmetryOperation {
    private:
        double d_[3][3];

    public:
        SymmetryOperation();
        SymmetryOperation(const SymmetryOperation &);
        virtual ~SymmetryOperation();

        double trace() const { return d_[0][0] + d_[1][1] + d_[2][2]; }

        double* operator[](int i) { return d_[i]; }

        const double* operator[](int i) const { return d_[i]; }

        double& operator()(int i, int j) { return d_[i][j]; }

        double operator()(int i, int j) const { return d_[i][j]; }

        void zero() { memset(d_, 0, sizeof(double) * 9); }

        SymmetryOperation operate(const SymmetryOperation& r) const;

        SymmetryOperation transform(const SymmetryOperation& r) const;

        void unit() { zero(); d_[0][0] = d_[1][1] = d_[2][2] = 1.0; }

        void E() { unit(); }

        void i() { zero(); d_[0][0] = d_[1][1] = d_[2][2] = -1.0; }

        void sigma_h() { unit(); d_[2][2] = -1.0; }
        void sigma_xz() { unit(); d_[1][1] = -1.0; }
        void sigma_yz() { unit(); d_[0][0] = -1.0; }

        void rotation(int n);
        void rotation(double theta);

        void c2_x() { i(); d_[0][0] = 1.0; }
        void c2_y() { i(); d_[1][1] = 1.0; }

        void transpose();
    };

    class SymRep {
    private:
        int n_;
        double d_[5][5];

    public:
        SymRep(int =0);
        SymRep(const SymmetryOperation&);
        virtual ~SymRep();

        operator SymmetryOperation() const;

        inline double trace() const;

        void set_dim(int i) { n_ = i; }

        double* operator[](int i) { return d_[i]; }
        const double* operator[](int i) const { return d_[i]; }

        double& operator()(int i, int j) { return d_[i][j]; }
        double operator()(int i, int j) const { return d_[i][j]; }

        void zero() { memset(d_, 0, sizeof(double) * 25); }

        SymRep operate(const SymRep& r) const;
        SymRep transform(const SymRep& r) const;

        void unit() {
            zero();
            d_[0][0] = d_[1][1] = d_[2][2] = d_[3][3] = d_[4][4] = 1.0;
        }

        void E() { unit(); }
        void i() {
            zero();
            d_[0][0] = d_[1][1] = d_[2][2] = d_[3][3] = d_[4][4] = -1.0;
        }

        void sigma_h();
        void sigma_xz();
        void sigma_yz();

        void rotation(int n);
        void rotation(double theta);

        void c2_x();
        void c2_y(); 
    };

    inline double SymRep::trace() const {
        double r=0;
        for (int i=0; i<n_; ++i)
            r += d_[i][i];
        return r;
    }

    class CharacterTable;

    class IrreducibleRepresentation {
        friend class CharacterTable;

    private:
        int g_;
        int degen_;
        int nrot_;
        int ntrans_;
        int complex_;
        char *symb_;
        char *csymb_;
        SymRep *rep_;

    public:
        IrreducibleRepresentation();
        IrreducibleRepresentation(const IrreducibleRepresentation&);
        IrreducibleRepresentation(int,int,const char*,const char*=0);

        virtual ~IrreducibleRepresentation();

        IrreducibleRepresentation& operator=(const IrreducibleRepresentation&);

        void init(int=0, int=0, const char*=0, const char*=0);

        int order() const { return g_; }

        int degeneracy() const { return degen_; }

        int complex() const { return complex_; }

        int nproj() const { return degen_*degen_; }

        int nrot() const { return nrot_; }

        int ntrans() const { return ntrans_; }

        const char* symbol() const { return symb_; }

        const char* symbol_ns() const { return (csymb_?csymb_:symb_); }

        double character(int i) const {
            return complex_ ? 0.5*rep_[i].trace() : rep_[i].trace();
        }

        double p(int x1, int x2, int i) const { return rep_[i](x1,x2); }

        double p(int d, int i) const {
            int dc=d/degen_; int dr=d%degen_;
            return rep_[i](dr,dc);
        }
    };

    class CharacterTable {
    public:
        enum pgroups { C1, CS, CI, CN, CNV, CNH, DN, DND, DNH, SN, T, TH, TD, O, OH, I, IH };

    private:
        int g_;
        int nt_;
        pgroups pg_;
        int nirreps_;
        IrreducibleRepresentation *gamma_;
        SymmetryOperation *symop_;
        int *inv_;
        char *symb_;

        int parse_symbol();
        int make_table();

        void t();
        void th();
        void td();
        void o();
        void oh();
        void i();
        void ih();

    public:
        CharacterTable();
        CharacterTable(const char*);
        CharacterTable(const char*, const SymmetryOperation&);

        CharacterTable(const CharacterTable&);
        virtual ~CharacterTable();

        CharacterTable& operator=(const CharacterTable&);

        int nirrep() const { return nirrep_; }
        int order() const { return g_; }
        const char* symbol() const { return symb_; }

        IrreducibleRepresentation& gamma(int i) { return gamma_[i]; }
        SymmetryOperation& symm_operation(int i) { return symop_[i]; }

        int complex() const {
            if (pg_ == CN || pg_ == SN || pg_ == CNH || 
                pg_ == T || pg_ == TH)
                return 1;
            return 0;
        }

        int inverse(int i) const { return inv_[i]; }

        int ncomp() const {
            int ret=0;
            for (int i=0; i<nirrep_; ++i) {
                int nc = (gamma_[i].complex()) ? 1 : gamma_[i].degen_;
                ret += nc;
            }
            return ret;
        }

        int which_irrep(int i) {
            for (int ir=0, cn=0; ir<nirrep_; ++ir) {
                int nc=(gamma_[ir].complex()) ? 1 : gamma_[ir].degen_;
                for (int c=0; c<nc_; c++,cn++)
                    if (cn==i)
                        return ir;
            }
            return -1;
        }

        int which_compp(int i) {
            for (int ir=0, cn=0; ir<nirrep_; ++ir) {
                int nc=(gamma_[ir].complex()) ? 1 : gamma_[ir].degen_;
                for (int c=0; c<nc_; c++,cn++)
                    if (cn==i)
                        return c;
            }
            return -1;
        }
    };

}

#endif // _basis_symmetry_h_
