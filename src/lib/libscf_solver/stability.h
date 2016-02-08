#ifndef STABILITY_H
#define STABILITY_H

#endif // STABILITY_H

#include <libmints/wavefunction.h>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class BasisSet;
class Matrix;
class TwoBodyAOInt;
class JK;
class VBase;

namespace scf {

class UStab {

protected:

    std::vector<std::pair<SharedMatrix,SharedMatrix> > vecs_;
    std::vector<double> vals_;

    bool unstable = false;
    double unstable_val = 0.0;
    std::pair<SharedMatrix,SharedMatrix> unstable_vec;

    int print_;
    int bench_;
    int debug_;
    long int memory_;

    //SharedMatrix C_;

    SharedMatrix Cocca_;
    SharedMatrix Coccb_;
    SharedMatrix Cvira_;
    SharedMatrix Cvirb_;
    SharedMatrix Ca_;
    SharedMatrix Cb_;

    boost::shared_ptr<Vector> eps_occa_;
    boost::shared_ptr<Vector> eps_vira_;
    boost::shared_ptr<Vector> eps_occb_;
    boost::shared_ptr<Vector> eps_virb_;

    boost::shared_ptr<Wavefunction> reference_wavefunction_;
    boost::shared_ptr<Molecule> molecule_;
    boost::shared_ptr<BasisSet> basis_;

    SharedMatrix AO2USO_;

    Options& options_;

    /// How far to converge the two-norm of the residual
    double convergence_;
    /// Global JK object, built in preiterations, destroyed in postiterations
    boost::shared_ptr<JK> jk_;
    boost::shared_ptr<VBase> v_;

    double Eref_;

    void common_init();
    void print_header();
    void preiterations();

public:

    UStab(SharedWavefunction ref_wfn, Options& options);
    virtual ~UStab();

    /// Gets a handle to the JK object, if built by preiterations
    boost::shared_ptr<JK> jk() const { return jk_;}
    /// Set the JK object, say from SCF
    void set_jk(boost::shared_ptr<JK> jk) { jk_ = jk; }
    /// Gets a handle to the VBase object, if built by preiterations
    boost::shared_ptr<VBase> v() const { return v_;}
    /// Set the VBase object, say from SCF (except that wouldn't work, right?)
    void set_jk(boost::shared_ptr<VBase> v) { v_ = v; }
    /// Is the wavefunction stable ?
    bool is_unstable() const { return unstable;}

    /// Get the eigenvalue for storage and comparison.
    double get_eigval() const {return unstable_val;}

    /// => Setters <= ///

    /// Set convergence behavior
    void set_convergence(double convergence) { convergence_ = convergence; }

    /// Update reference info
    void set_reference(boost::shared_ptr<Wavefunction> reference);

    virtual double compute_energy();
    SharedMatrix analyze();
    void rotate_orbs(double scale);
};

} // namespace scf


} // namespace psi
