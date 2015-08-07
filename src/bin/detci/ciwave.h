/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#ifndef CIWAVE_H
#define CIWAVE_H

// Forward declarations
namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class Wavefunction;
class Options;
class JK;
class DFERI;
class IntegralTransform;
class MOSpace;
struct stringwr;
typedef boost::shared_ptr<Matrix> SharedMatrix;
}

namespace psi { namespace detci {

class CIWavefunction : public Wavefunction
{

public:
    CIWavefunction(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options);
    virtual ~CIWavefunction();

    double compute_energy();
    PsiReturnType cas_update();

    /**!
     * Similar to wavefunction.Ca_subset(); however, this version knows about all of the CI
     * subspaces in the SO basis. We stick to the MCSCF definitions for now.
     * @param  orbital_name FZC, DOCC, ACT, VIR, FZV
     * @return C            Returns the appropriate orbitals in the SO basis.
     */
    SharedMatrix get_orbitals(const std::string& orbital_name);

    /**!
     * Similar to wavefunction.Ca_subset(); however, this version knows about all of the CI
     * subspaces in the SO basis. We stick to the MCSCF definitions for now.
     * @param  orbital_name FZC, DOCC, ACT, VIR, FZV
     * @param  orbitals     SharedMatrix to set
     * @return C            Returns the appropriate orbitals in the SO basis.
     */
    void set_orbitals(const std::string& orbital_name, SharedMatrix orbitals);

    /**
     * Transform the one and two electron integrals.
     * The
     */
    void transform_ci_integrals(void);

    /**!
     * Obtains the OPDM <root| Epq |root> from disk
     * @param root       Root to obtain
     * @param spin       0 returns alpha, 1 returns beta, 2 sums alpha and beta
     * @param transden   Obtain <0 | Epq | root> transition density matrix
     * @return OPDM or TDM shared matrix
     **/
    SharedMatrix get_opdm(int root=0, int spin=2, bool transden=false);

    /**!
     Returns the symmetry block active OPDM
     **/
    SharedMatrix get_active_opdm();

    /**!
     * Loads the OPDM to the wavefunction Da_ and Db_
     * @param use_old_d - If no new OPDM is calculated use the reference density matrices
    **/
    void set_opdm(bool use_old_d=false);

    /**!
     * Obtains the TPDM from disk
     * @param symmetrized Symmetrize the TPDM or not
     * @return TPDM SharedVector
     **/
    SharedVector get_tpdm(bool symmetrized=true, bool act_only=true);

    /**!
     Returns the a full 4D active TPDM
     **/
    SharedMatrix get_active_tpdm();

    /**!
     * Sets the dense TPDM to the wavefunction TPDM_
    **/
    void set_tpdm();

    // Should be private
    void compute_mcscf(struct stringwr **alplist, struct stringwr **betlist);

private:

    // Grabs mo info
    void get_mo_info();

    // Sets the ciwavefunction object
    void common_init();

    // Find out which orbitals belong hwere
    void orbital_locations(const std::string& orbital_name, int* start, int* end);

    // Symmetry block a matrix
    SharedMatrix symm_block(SharedMatrix x, Dimension dim1, Dimension dim2);

    // Integrals
    bool ints_init_;
    boost::shared_ptr<IntegralTransform> ints_; // Non-DF
    boost::shared_ptr<DFERI> dferi_; // For DF
    boost::shared_ptr<JK> jk_;       // For DF

    void tf_onel_ints();
    void form_gmat();

    void setup_mcscf_ints();
    void transform_mcscf_ints();
    void setup_dfmcscf_ints();
    void transform_dfmcscf_ints();
    // void transform_dfmcscf_ints();


};

}}

#endif // CIWAVE_H

