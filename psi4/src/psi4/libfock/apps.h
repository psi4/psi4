/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef APPS_H
#define APPS_H
#include <set>
#include <tuple>
#include "psi4/libmints/wavefunction.h"

namespace psi {

class JK;
class VBase;

// => BASE CLASSES <= //

class RBase : public Wavefunction {
   protected:
    int print_;
    int bench_;

    SharedMatrix C_;

    SharedMatrix Cocc_;
    SharedMatrix Cfocc_;
    SharedMatrix Cfvir_;
    SharedMatrix Caocc_;
    SharedMatrix Cavir_;

    std::shared_ptr<Vector> eps_focc_;
    std::shared_ptr<Vector> eps_fvir_;
    std::shared_ptr<Vector> eps_aocc_;
    std::shared_ptr<Vector> eps_avir_;

    SharedMatrix AO2USO_;

    /// How far to converge the two-norm of the residual
    double convergence_;
    /// Global JK object, built in preiterations, destroyed in postiterations
    std::shared_ptr<JK> jk_;
    std::shared_ptr<VBase> v_;

    bool use_symmetry_;
    double Eref_;

   public:
    RBase(SharedWavefunction ref_wfn, Options& options, bool use_symmetry = true);
    // TODO: Remove AS SOON AS POSSIBLE, such a dirty hack
    RBase(bool flag);
    ~RBase() override;

    virtual bool same_a_b_orbs() const { return true; }
    virtual bool same_a_b_dens() const { return true; }

    // TODO: Remove AS SOON AS POSSIBLE, such a dirty hack
    double compute_energy() override { return 0.0; }

    void set_print(int print) { print_ = print; }

    /// Gets a handle to the JK object, if built by preiterations
    std::shared_ptr<JK> jk() const { return jk_; }
    /// Set the JK object, say from SCF
    void set_jk(std::shared_ptr<JK> jk) { jk_ = jk; }
    /// Gets a handle to the VBase object, if built by preiterations
    std::shared_ptr<VBase> v() const { return v_; }
    /// Set the VBase object, say from SCF (except that wouldn't work, right?)
    void set_jk(std::shared_ptr<VBase> v) { v_ = v; }
    /// Builds JK object, if needed
    virtual void preiterations();
    /// Destroys JK object, if needed
    virtual void postiterations();

    /// => Setters <= ///

    void set_use_symmetry(bool usesym) { use_symmetry_ = usesym; }
    /// Set convergence behavior
    void set_convergence(double convergence) { convergence_ = convergence; }

    /// Set reference info
    void set_C(SharedMatrix C) { C_ = C; }
    void set_Cocc(SharedMatrix Cocc) { Cocc_ = Cocc; }
    void set_Cfocc(SharedMatrix Cfocc) { Cfocc_ = Cfocc; }
    void set_Caocc(SharedMatrix Caocc) { Caocc_ = Caocc; }
    void set_Cavir(SharedMatrix Cavir) { Cavir_ = Cavir; }
    void set_Cfvir(SharedMatrix Cfvir) { Cfvir_ = Cfvir; }
    void set_eps_focc(SharedVector eps) { eps_focc_ = eps; }
    void set_eps_aocc(SharedVector eps) { eps_aocc_ = eps; }
    void set_eps_avir(SharedVector eps) { eps_avir_ = eps; }
    void set_eps_fvir(SharedVector eps) { eps_fvir_ = eps; }
    void set_Eref(double Eref) { Eref_ = Eref; }

    /// Update reference info
    void set_reference(std::shared_ptr<Wavefunction> reference);
};

// => APPLIED CLASSES <= //

class RCIS : public RBase {
   protected:
    std::vector<std::tuple<double, int, int, int> > states_;
    std::vector<SharedMatrix> singlets_;
    std::vector<SharedMatrix> triplets_;
    std::vector<double> E_singlets_;
    std::vector<double> E_triplets_;

    void sort_states();

    virtual void print_header();
    virtual void print_wavefunctions();
    virtual void print_amplitudes();
    virtual void print_transitions();
    virtual void print_densities();

    virtual SharedMatrix TDmo(SharedMatrix T1, bool singlet = true);
    virtual SharedMatrix TDso(SharedMatrix T1, bool singlet = true);
    virtual SharedMatrix TDao(SharedMatrix T1, bool singlet = true);

    virtual SharedMatrix Dmo(SharedMatrix T1, bool diff = false);
    virtual SharedMatrix Dso(SharedMatrix T1, bool diff = false);
    virtual SharedMatrix Dao(SharedMatrix T1, bool diff = false);

    virtual std::pair<SharedMatrix, std::shared_ptr<Vector> > Nmo(SharedMatrix T1, bool diff = false);
    virtual std::pair<SharedMatrix, std::shared_ptr<Vector> > Nso(SharedMatrix T1, bool diff = false);
    virtual std::pair<SharedMatrix, std::shared_ptr<Vector> > Nao(SharedMatrix T1, bool diff = false);

    virtual std::pair<SharedMatrix, SharedMatrix> ADmo(SharedMatrix T1);
    virtual std::pair<SharedMatrix, SharedMatrix> ADso(SharedMatrix T1);
    virtual std::pair<SharedMatrix, SharedMatrix> ADao(SharedMatrix T1);

   public:
    RCIS(SharedWavefunction ref_wfn, Options& options);
    ~RCIS() override;

    double compute_energy() override;
};

class RTDHF : public RBase {
   protected:
    std::vector<SharedMatrix> singlets_X_;
    std::vector<SharedMatrix> triplets_X_;
    std::vector<SharedMatrix> singlets_Y_;
    std::vector<SharedMatrix> triplets_Y_;
    std::vector<double> E_singlets_;
    std::vector<double> E_triplets_;

    virtual void print_header();

   public:
    RTDHF(SharedWavefunction ref_wfn, Options& options);
    ~RTDHF() override;

    double compute_energy() override;
};

class RCPHF : public RBase {
   protected:
    // OV-Rotations
    std::map<std::string, SharedMatrix> x_;
    // OV-Perturbations
    std::map<std::string, SharedMatrix> b_;

    virtual void print_header();

    void add_named_tasks();
    void analyze_named_tasks();

    void add_polarizability();
    void analyze_polarizability();

    std::set<std::string> tasks_;

   public:
    RCPHF(SharedWavefunction ref_wfn, Options& options, bool use_symmetry = true);
    ~RCPHF() override;

    /// Solve for all perturbations currently in b
    double compute_energy() override;

    /// Perturbation vector queue, shove tasks onto this guy before compute_energy
    std::map<std::string, SharedMatrix>& b() { return b_; }
    /// Resultant solution vectors, available after compute_energy is called
    std::map<std::string, SharedMatrix>& x() { return x_; }

    /// Add a named task
    void add_task(const std::string& task);
};

class RCPKS : public RCPHF {
   protected:
    void print_header() override;

   public:
    RCPKS(SharedWavefunction ref_wfn, Options& options);
    ~RCPKS() override;

    double compute_energy() override;
};

class RTDA : public RCIS {
   protected:
    void print_header() override;

   public:
    RTDA(SharedWavefunction ref_wfn, Options& options);
    ~RTDA() override;

    double compute_energy() override;
};

class RTDDFT : public RTDHF {
   protected:
    void print_header() override;

   public:
    RTDDFT(SharedWavefunction ref_wfn, Options& options);
    ~RTDDFT() override;

    double compute_energy() override;
};
}
#endif
