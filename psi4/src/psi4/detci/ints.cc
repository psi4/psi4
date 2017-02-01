/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here
*/


/*
** INTS.C
**
** Sets the one and two-electron integrals
** For DPD X and R will denonte the active and rotatable space, respectively
** For DFERI a and N will denonte the active and rotatable space, respectively
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia
**
** Updated 3/18/95 to exclude frozen virtual orbitals.
** Updated 3/28/95 to exclude frozen core orbitals.
** DGAS: Rewrote 7/15/15 Now uses libtrans and DFERI objects
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "psi4/libiwl/iwl.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/vector.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libthce/thce.h"
#include "psi4/libthce/lreri.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/soscf.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/detci/structs.h"
#include "psi4/detci/ciwave.h"
#include "psi4/detci/globaldefs.h"

namespace psi { namespace detci {

void CIWavefunction::transform_ci_integrals() {
    outfile->Printf("\n   ==> Transforming CI integrals <==\n\n");
    // Grab orbitals
    SharedMatrix Cdrc = get_orbitals("DRC");
    SharedMatrix Cact = get_orbitals("ACT");
    SharedMatrix Cvir = get_orbitals("VIR");
    SharedMatrix Cfzv = get_orbitals("FZV");

    // Build up active space
    std::vector<std::shared_ptr<MOSpace> > spaces;

    // Indices should be empty
    std::vector<int> indices(CalcInfo_->num_ci_orbs, 0);
    std::vector<int> orbitals(CalcInfo_->num_ci_orbs, 0);

    for (int h = 0, cinum = 0, orbnum = 0; h < CalcInfo_->nirreps; h++) {
        orbnum += CalcInfo_->dropped_docc[h];
        for (int i = 0; i < CalcInfo_->ci_orbs[h]; i++) {
            orbitals[cinum++] = orbnum++;
        }
        orbnum += CalcInfo_->dropped_uocc[h];
    }

    std::shared_ptr<MOSpace> act_space(new MOSpace('X', orbitals, indices));
    spaces.push_back(act_space);

    IntegralTransform* ints = new IntegralTransform(
        Cdrc, Cact, Cvir, Cfzv, spaces, IntegralTransform::Restricted,
        IntegralTransform::DPDOnly, IntegralTransform::PitzerOrder,
        IntegralTransform::OccAndVir, true);
    ints_ = std::shared_ptr<IntegralTransform>(ints);
    ints_->set_memory(Process::environment.get_memory() * 0.8);

    // Incase we do two ci runs
    dpd_set_default(ints_->get_dpd_id());
    ints_->set_keep_iwl_so_ints(true);
    ints_->set_keep_dpd_so_ints(true);
    ints_->transform_tei(act_space, act_space, act_space, act_space);

    // Read drc energy
    CalcInfo_->edrc = ints_->get_frozen_core_energy();

    // Read integrals
    read_dpd_ci_ints();

    // Form auxiliary matrices
    tf_onel_ints(CalcInfo_->onel_ints, CalcInfo_->twoel_ints, CalcInfo_->tf_onel_ints);
    form_gmat(CalcInfo_->onel_ints, CalcInfo_->twoel_ints, CalcInfo_->gmat);

    // This is a build and burn call
    ints_.reset();
}
void CIWavefunction::setup_dfmcscf_ints() {
    outfile->Printf("\n   ==> Setting up DF-MCSCF integrals <==\n\n");

    /// Build JK object
    if (options_.get_str("SCF_TYPE") == "DF"){
        jk_ = JK::build_JK(basisset_, get_basisset("DF_BASIS_SCF"), options_);
    } else {
        jk_ = JK::build_JK(basisset_, BasisSet::zero_ao_basis_set(), options_);
    }
    jk_->set_do_J(true);
    jk_->set_do_K(true);
    jk_->initialize();
    jk_->set_memory(Process::environment.get_memory() * 0.8);
    jk_->print_header();

    /// Build DF object
    dferi_ = DFERI::build(get_basisset("ORBITAL"), get_basisset("DF_BASIS_SCF"), options_);
    dferi_->print_header();

    df_ints_init_ = true;
}
void CIWavefunction::transform_mcscf_integrals(bool approx_only) {
    if (Parameters_->mcscf_type == "DF") {
        transform_dfmcscf_ints(approx_only);
    }
    else if (Parameters_->mcscf_type == "AO")
        transform_mcscf_ints_ao(approx_only);
    else {
        transform_mcscf_ints(approx_only);
    }
}
void CIWavefunction::rotate_mcscf_integrals(SharedMatrix k,
                                            SharedVector onel_out,
                                            SharedVector twoel_out) {

    // => Sizing <= //
    Dimension zero_dim(nirrep_, "Zero_dim");
    Dimension av_dim = CalcInfo_->ci_orbs + CalcInfo_->rstr_uocc;
    Dimension oa_dim = CalcInfo_->ci_orbs + CalcInfo_->rstr_docc;
    Dimension rot_dim = get_dimension("ROT");

    if ((k->rowspi() != oa_dim) || (k->colspi() != av_dim)){
        throw PSIEXCEPTION("CIWavefunction::rotate_dfmcscf_ints: Rotation matrix k is not of the correct shape");
    }

    // Build up Uact (act x nmo)
    SharedMatrix Uact(new Matrix("Active U", CalcInfo_->ci_orbs, rot_dim));
    for (int h=0; h<nirrep_; h++){

        for (int i=0; i<CalcInfo_->ci_orbs[h]; i++){
            for (int j=0; j<CalcInfo_->rstr_docc[h]; j++){
                Uact->set(h, i, j, -1 * k->get(h, j, i));
            }

            for (int j=0; j<av_dim[h]; j++){
                Uact->set(h, i, j + CalcInfo_->rstr_docc[h], k->get(h, i + CalcInfo_->rstr_docc[h], j));
            }
        }

    }
    // Uact->print();


    // => Setup <= //
    SharedMatrix Cact = get_orbitals("ACT");
    SharedMatrix Crot = get_orbitals("ROT");

    SharedMatrix H_rot_a = Matrix::triplet(Crot, CalcInfo_->so_onel_ints, Cact, true, false, false);

    // => Rotate onel ints <= //
    SharedMatrix rot_onel = Matrix::doublet(Uact, H_rot_a);
    rot_onel->gemm(true, true, 1.0, H_rot_a, Uact, 1.0);

    SharedVector tmponel(new Vector("Temporary onel storage", CalcInfo_->num_ci_tri));

    pitzer_to_ci_order_onel(rot_onel, tmponel);

    if (Parameters_->mcscf_type == "DF") {
        rotate_dfmcscf_twoel_ints(Uact, twoel_out);
    } else {
        rotate_mcscf_twoel_ints(Uact, twoel_out);
    }

    // Transform the onel ints
    if (Parameters_->fci){
        tf_onel_ints(tmponel, twoel_out, onel_out);
    }
    else{
        form_gmat(tmponel, twoel_out, onel_out);
    }

}
void CIWavefunction::transform_dfmcscf_ints(bool approx_only) {
    if (!df_ints_init_) setup_dfmcscf_ints();
    timer_on("CIWave: DFMCSCF integral transform");

    // => AO C matrices <= //
    // We want a pitzer order C matrix with appended Cact
    SharedMatrix Cocc = get_orbitals("DOCC");
    SharedMatrix Cact = get_orbitals("ACT");
    SharedMatrix Cvir = get_orbitals("VIR");

    int nao = AO2SO_->rowspi()[0];
    int nact = CalcInfo_->num_ci_orbs;
    int nrot = Cocc->ncol() + Cact->ncol() + Cvir->ncol();
    int aoc_rowdim = nrot + Cact->ncol();
    SharedMatrix AO_C = SharedMatrix(new Matrix("AO_C", nao, aoc_rowdim));

    double** AO_Cp = AO_C->pointer();
    for (int h = 0, offset = 0, offset_act = 0; h < nirrep_; h++) {
        int hnso = AO2SO_->colspi()[h];
        if (hnso == 0) continue;
        double** Up = AO2SO_->pointer(h);

        int noccpih = Cocc->colspi()[h];
        int nactpih = Cact->colspi()[h];
        int nvirpih = Cvir->colspi()[h];
        // occupied
        if (noccpih) {
            double** CSOp = Cocc->pointer(h);
            C_DGEMM('N', 'N', nao, noccpih, hnso, 1.0, Up[0], hnso, CSOp[0],
                    noccpih, 0.0, &AO_Cp[0][offset], aoc_rowdim);
            offset += noccpih;
        }
        // active
        if (nactpih) {
            double** CSOp = Cact->pointer(h);
            C_DGEMM('N', 'N', nao, nactpih, hnso, 1.0, Up[0], hnso, CSOp[0],
                    nactpih, 0.0, &AO_Cp[0][offset], aoc_rowdim);
            offset += nactpih;

            C_DGEMM('N', 'N', nao, nactpih, hnso, 1.0, Up[0], hnso, CSOp[0],
                    nactpih, 0.0, &AO_Cp[0][offset_act + nrot], aoc_rowdim);
            offset_act += nactpih;
        }
        // virtual
        if (nvirpih) {
            double** CSOp = Cvir->pointer(h);
            C_DGEMM('N', 'N', nao, nvirpih, hnso, 1.0, Up[0], hnso, CSOp[0],
                    nvirpih, 0.0, &AO_Cp[0][offset], aoc_rowdim);
            offset += nvirpih;
        }
    }

    // => Compute DF ints <= //
    dferi_->clear();
    dferi_->set_C(AO_C);
    dferi_->add_space("R", 0, nrot);
    dferi_->add_space("a", nrot, aoc_rowdim);
    dferi_->add_space("F", 0, aoc_rowdim);

    if (approx_only) {
        dferi_->add_pair_space("aaQ", "a", "a");
        dferi_->add_pair_space("RaQ", "a", "R", -1.0 / 2.0, true);
    } else {
        dferi_->add_pair_space("aaQ", "a", "a");
        dferi_->add_pair_space("RaQ", "a", "R", -1.0 / 2.0, true);
        dferi_->add_pair_space("RRQ", "R", "R");
    }

    dferi_->compute();
    std::map<std::string, std::shared_ptr<Tensor> >& dfints = dferi_->ints();

    // => Compute onel ints <= //
    onel_ints_from_jk();

    // => Compute twoel ints <= //
    int nQ = dferi_->size_Q();

    std::shared_ptr<Tensor> aaQT = dfints["aaQ"];
    SharedMatrix aaQ(new Matrix("aaQ", nact * nact, nQ));

    double* aaQp = aaQ->pointer()[0];
    FILE* aaQF = aaQT->file_pointer();
    fseek(aaQF, 0L, SEEK_SET);
    fread(aaQp, sizeof(double), nact * nact * nQ, aaQF);
    SharedMatrix actMO = Matrix::doublet(aaQ, aaQ, false, true);
    aaQ.reset();

    pitzer_to_ci_order_twoel(actMO, CalcInfo_->twoel_ints);
    actMO.reset();

    tf_onel_ints(CalcInfo_->onel_ints, CalcInfo_->twoel_ints, CalcInfo_->tf_onel_ints);
    form_gmat(CalcInfo_->onel_ints, CalcInfo_->twoel_ints, CalcInfo_->gmat);
    timer_off("CIWave: DFMCSCF integral transform");
}
void CIWavefunction::setup_mcscf_ints() {
    // We need to do a few weird things to make IntegralTransform work for us
    outfile->Printf("\n   ==> Setting up MCSCF integrals <==\n\n");

    // Grab orbitals
    SharedMatrix Cdrc = get_orbitals("DRC");
    SharedMatrix Cact = get_orbitals("ACT");
    SharedMatrix Cvir = get_orbitals("VIR");
    SharedMatrix Cfzv = get_orbitals("FZV");

    // Need active and rot spaces
    std::vector<std::shared_ptr<MOSpace> > spaces;

    std::vector<int> rot_orbitals(CalcInfo_->num_rot_orbs, 0);
    std::vector<int> act_orbitals(CalcInfo_->num_ci_orbs, 0);

    // Indices *should* be zero, DPD does not use it
    std::vector<int> indices(CalcInfo_->num_ci_orbs, 0);

    int act_orbnum = 0;
    int rot_orbnum = 0;
    for (int h = 0, rn = 0, an = 0; h < CalcInfo_->nirreps; h++) {
        act_orbnum += CalcInfo_->dropped_docc[h];
        rot_orbnum += CalcInfo_->frozen_docc[h];

        // Act space
        for (int i = 0; i < CalcInfo_->ci_orbs[h]; i++) {
            act_orbitals[an++] = act_orbnum++;
        }
        act_orbnum += CalcInfo_->dropped_uocc[h];

        int nrotorbs = CalcInfo_->rstr_docc[h] + CalcInfo_->ci_orbs[h] +
                       CalcInfo_->rstr_uocc[h];
        for (int i = 0; i < nrotorbs; i++) {
            rot_orbitals[rn++] = rot_orbnum++;
        }
        rot_orbnum += CalcInfo_->frozen_uocc[h];
    }

    rot_space_ = std::shared_ptr<MOSpace>(new MOSpace('R', rot_orbitals, indices));
    act_space_ = std::shared_ptr<MOSpace>(new MOSpace('X', act_orbitals, indices));
    spaces.push_back(rot_space_);
    spaces.push_back(act_space_);

    // Now the occ space is active, the vir space is our rot space (FZC to FZV)
    IntegralTransform* ints = new IntegralTransform(
        Cdrc, Cact, Cvir, Cfzv, spaces, IntegralTransform::Restricted,
        IntegralTransform::DPDOnly, IntegralTransform::PitzerOrder,
        IntegralTransform::OccAndVir, true);
    ints_ = std::shared_ptr<IntegralTransform>(ints);
    ints_->set_memory(Process::environment.get_memory() * 0.8);

    // Incase we do two ci runs
    dpd_set_default(ints_->get_dpd_id());
    ints_->set_keep_iwl_so_ints(true);
    ints_->set_keep_dpd_so_ints(true);
    ints_->set_print(0);

    // Conventional JK build
    if (options_.get_str("SCF_TYPE") == "DF"){
        jk_ = JK::build_JK(basisset_, get_basisset("DF_BASIS_SCF"), options_);
    } else {
        jk_ = JK::build_JK(basisset_, BasisSet::zero_ao_basis_set(), options_);
    }
    jk_->set_do_J(true);
    jk_->set_do_K(true);
    jk_->set_memory(Process::environment.get_memory() * 0.8);
    jk_->initialize();
    jk_->print_header();

    ints_init_ = true;
}
void CIWavefunction::setup_mcscf_ints_ao()
{
    outfile->Printf("\n   ==> Setting up MCSCF integrals <==\n\n");

    timer_on("CIWave: Setup MCSCF INTS AO");
    std::string scf_type = options_.get_str("SCF_TYPE");
    if(scf_type == "GTFOCK")
    {
 #ifdef HAVE_JK_FACTORY
        Process::environment.set_legacy_molecule(molecule_);
        jk_ = boost::shared_ptr<JK>(new GTFockJK(basisset_));
#else
        throw PSIEXCEPTION("GTFock was not compiled in this version");
#endif
    } else if (scf_type == "DF") {
        jk_ = JK::build_JK(this->basisset(), get_basisset("DF_BASIS_SCF"), options_);
    }
    else if (scf_type == "CD" or scf_type == "PK" or scf_type == "DIRECT" or scf_type == "OUT_OF_CORE")
    {
        jk_ = JK::build_JK(this->basisset(), BasisSet::zero_ao_basis_set(), options_);
    }
    else {
        outfile->Printf("\n Please select GTFock, DF, CD or PK for use with MCSCF_TYPE AO");
        throw PSIEXCEPTION("AO_CASSCF does not work with your SCF_TYPE");
    }
    jk_->set_do_J(true);
    jk_->set_allow_desymmetrization(true);
    jk_->set_do_K(true);
    jk_->set_memory(Process::environment.get_memory() * 0.8);
    jk_->initialize();
    jk_->print_header();
    ints_init_ = true;
    timer_off("CIWave: Setup MCSCF INTS AO");
}
void CIWavefunction::transform_mcscf_ints_ao(bool approx_only)
{
    if(!approx_only)
        throw PSIEXCEPTION("Have not implemented rotated integrals for mcscf_ao");

    ///Perform a integral transformation using jk object
    ///KPH got this idea from Hohenstein's AO-CASSCF paper.
    ///The goal is to use direct scf algorithm for computing the transformed integrals
    ///This takes advantage of sparsity and allows for larger scale casscf calculations.

    /// J(D)_{mu nu} = (mu nu | rho sigma ) D_{rho sigma}
    ///Step 1:  Form a density for every active orbital ie
    /// for every u and v form a density using C_DGER
    ///     D_{mu nu}^{uv} = C_{mu u}C_{nu}^{v}
    ///Step 2:  Fill the JK object with these densities
    ///Step 3:  Build the J matrix
    ///Step 4:  (p u | x y) = C_{mu p}^{T} J(D_{mu nu}^{xy}A)_{mu nu} C_{nu u}
    ///Step 5:  Transfer these integrals to CI and SOMCSCF for CASSCF optimization
    //outfile->Printf("\n   ==> Transforming CI integrals aodirect <==\n");
    timer_on("CIWave: AO MCSCF integral transform");
    if(!ints_init_)
    {
        setup_mcscf_ints_ao();
    }

    int nact = CalcInfo_->num_ci_orbs;
    int nrot = nmo_ - CalcInfo_->num_fzc_orbs - CalcInfo_->num_fzv_orbs;

    // Transform from the SO to the AO basis for the C matrix.
    // just transfroms the C_{mu_ao i} -> C_{mu_so i}
    /// If C1, do not do.  C_{mu_ao i} = C_{mu_so i} if C1
    SharedMatrix Crot;
    if (nirrep_ > 1) {
        Crot = SharedMatrix(new Matrix(nso_, nrot));
        timer_on("CIWave: Transform C matrix from SO to AO");

        SharedMatrix Cso_rot = get_orbitals("ROT");
        int nao = AO2SO_->rowspi()[0];

        int offset = 0;
        double** Crotp = Crot->pointer();
        for (int h = 0, offset = 0, offset_act = 0; h < nirrep_; h++) {
            int hnso = AO2SO_->colspi()[h];
            if (hnso == 0) continue;
            double** Up = AO2SO_->pointer(h);
            int nrotpih = Cso_rot->colspi()[h];

            if (nrotpih) {
                double** CSOp = Cso_rot->pointer(h);
                C_DGEMM('N', 'N', nao, nrotpih, hnso, 1.0, Up[0], hnso, CSOp[0],
                        nrotpih, 0.0, &Crotp[0][offset], nrot);
                offset += nrotpih;
            }
        }

        //std::vector<int> nmo_offset(nirrep_, 0);
        //nmo_offset[0] = 0;

        ///// Have an orbital offset array (nmo_offset[2] = nmopi_[0] + nmopi[1];
        //for (int h = 1; h < nirrep_; h++) {
        //    nmo_offset[h] = nmo_offset[h - 1] + nmopi_[h - 1];
        //}

        //for (size_t h = 0, index = 0; h < nirrep_; ++h) {
        //    #pragma omp parallel for schedule(static)
        //    for (int i = CalcInfo_->frozen_docc[h]; i < nmopi_[h] - CalcInfo_->frozen_uocc[h]; ++i) {
        //        size_t nao = nso_;
        //        size_t nso = nsopi_[h];

        //        if (!nso) continue;
        //        int index = nmo_offset[h] + i;

        //        C_DGEMV('N', nao, nso, 1.0, AO2SO_->pointer(h)[0], nso, &Ca_->pointer(h)[0][i],
        //                nmopi_[h], 0.0, &Crot->pointer()[0][index], nrot);
        //    }
        //}
        timer_off("CIWave: Transform C matrix from SO to AO");
    } else {
        Crot = get_orbitals("ROT");
    }

    std::vector<int> active_abs(nact, 0);
    for (int h = 0, cinum = 0, orbnum = 0; h < CalcInfo_->nirreps; h++) {
        orbnum += CalcInfo_->rstr_docc[h];
        for (int i = 0; i < CalcInfo_->ci_orbs[h]; i++) {
            active_abs[cinum++] = orbnum++;
        }
        orbnum += CalcInfo_->rstr_uocc[h];
    }

    SharedMatrix Cact(new Matrix(nso_, nact));
    for (size_t v = 0; v < nact; v++) {
        SharedVector Crot_vec = Crot->get_column(0, active_abs[v]);
        Cact->set_column(0, v, Crot_vec);
    }

    timer_on("CIWave: Forming Active Psuedo Density");
    /// Step 1:  D_{mu nu} ^{tu} = C_{mu t} C_{nu u} forall t, u in active
    std::vector<std::tuple<int, int, SharedMatrix, SharedMatrix>> D_vec;
    for (size_t i = 0; i < nact; i++) {
        SharedVector C_i = Cact->get_column(0, i);
        SharedMatrix Cmat_i(new Matrix("C_i", nso_, 1));
        C_DCOPY(nso_, C_i->pointer(), 1, Cmat_i->pointer()[0], 1);

        for (size_t j = i; j < nact; j++) {
            SharedVector C_j = Cact->get_column(0, j);
            SharedMatrix Cmat_j(new Matrix("C_j", nso_, 1));
            C_DCOPY(nso_, C_j->pointer(), 1, Cmat_j->pointer()[0], 1);

            D_vec.push_back(std::make_tuple(i, j, Cmat_i, Cmat_j));
        }
    }
    timer_off("CIWave: Forming Active Psuedo Density");

    std::vector<SharedMatrix>  &Cl = jk_->C_left();
    std::vector<SharedMatrix>  &Cr = jk_->C_right();
    Cl.clear();
    Cr.clear();
    for(size_t d  = 0; d < D_vec.size(); d++)
    {
        Cl.push_back(std::get<2>(D_vec[d]));
        Cr.push_back(std::get<3>(D_vec[d]));
    }

    jk_->set_allow_desymmetrization(false);
    jk_->set_do_K(false);
    ///Step 2:  Compute the Coulomb build using these density
    timer_on("CIWave: AO MCSCF Integral Transformation Fock build");
    jk_->compute();
    timer_off("CIWave: AO MCSCF Integral Transformation Fock build");
    SharedMatrix casscf_ints(new Matrix("ALL Active", nrot * nact, nact * nact));

    ///Step 3:  Fill the integrals for SOSCF in_core
    ///TODO: Figure out WTF is going on with ordering.
    /// For LARGE (>40) active spaces, it is possible that N*A^3 might become too big
    timer_on("CIWave: Filling the (pu|xy) integrals");
    for(int D_tasks = 0; D_tasks < D_vec.size(); D_tasks++)
    {
        int i = std::get<0>(D_vec[D_tasks]);
        int j = std::get<1>(D_vec[D_tasks]);
        SharedMatrix J = jk_->J()[D_tasks];
        SharedMatrix half_trans = Matrix::triplet(Crot, J, Cact, true, false, false);
        #pragma omp parallel for schedule(static)
        for(size_t p = 0; p < nrot; p++){
            for(size_t q = 0; q < nact; q++){
                casscf_ints->set(p * nact + q, i * nact + j, half_trans->get(p, q));
                casscf_ints->set(p * nact + q, j * nact + i, half_trans->get(p, q));
            }
        }

    }
    timer_off("CIWave: Filling the (pu|xy) integrals");
    Cl.clear();
    Cr.clear();
    jk_->set_do_K(true);
    jk_->set_allow_desymmetrization(true);

    SharedMatrix actMO(new Matrix("ALL Active", nact * nact, nact * nact));
    timer_on("CIWave: Filling the (zu|xy) integrals");
    for(int u = 0; u < nact; u++){
        for(int v = 0; v < nact; v++){
            for(int x = 0; x < nact; x++){
                for(int y = 0; y < nact; y++){
                    actMO->set(u * nact + v, x * nact + y, casscf_ints->get(active_abs[u]* nact + v, x * nact + y));
                }
            }
        }
    }
    timer_off("CIWave: Filling the (zu|xy) integrals");

    timer_on("CIWave: Computing ^{I}F_{pq}");
    onel_ints_from_jk();
    timer_off("CIWave: Computing ^{I}F_{pq}");

    pitzer_to_ci_order_twoel(actMO, CalcInfo_->twoel_ints);
    tei_raaa_ = casscf_ints;
    tei_aaaa_ = actMO;
    somcscf_->set_eri_tensors(tei_aaaa_, tei_raaa_);

    tf_onel_ints(CalcInfo_->onel_ints, CalcInfo_->twoel_ints, CalcInfo_->tf_onel_ints);
    form_gmat(CalcInfo_->onel_ints, CalcInfo_->twoel_ints, CalcInfo_->gmat);
    timer_off("CIWave: AO MCSCF integral transform");
}
void CIWavefunction::transform_mcscf_ints(bool approx_only) {
    if (!ints_init_) setup_mcscf_ints();
    timer_on("CIWave: MCSCF integral transform");

    // The orbital matrix need to be identical to the previous one
    ints_->set_orbitals(get_orbitals("ALL"));

    if (approx_only) {
        // We only need (aa|aa), (aa|aN) for appoximate update
        ints_->set_keep_ht_ints(true);  // Save the aa half ints here
        ints_->transform_tei_first_half(act_space_, act_space_);
        ints_->transform_tei_second_half(act_space_, act_space_, act_space_, rot_space_);
        ints_->set_keep_ht_ints(false);  // Need to nuke the half ints after building act/act
        ints_->transform_tei_second_half(act_space_, act_space_, act_space_, act_space_);
    } else {
        // We need (aa|aa), (aa|aN), (aa|NN), (aN|aN)
        ints_->set_keep_ht_ints(false);  // Need to nuke the half ints after building
        ints_->transform_tei(act_space_, rot_space_, act_space_, rot_space_);

        // Half trans then work from there
        ints_->set_keep_ht_ints(true);  // Save the aa half ints here
        ints_->transform_tei_first_half(act_space_, act_space_);
        ints_->transform_tei_second_half(act_space_, act_space_, rot_space_, rot_space_);
        ints_->transform_tei_second_half(act_space_, act_space_, act_space_, rot_space_);
        ints_->set_keep_ht_ints(false);  // Need to nuke the half ints after building act/act
        ints_->transform_tei_second_half(act_space_, act_space_, act_space_, act_space_);
    }

    // Read DPD ints into memory
    read_dpd_ci_ints();

    // => Compute onel ints <= //
    // Libtrans does NOT change efzc or MO_FZC unless presort is called.
    // Presort is only called on initialize. As sort is very expensive, this
    // is a good thing.
    onel_ints_from_jk();

    // Form auxiliary matrices
    tf_onel_ints(CalcInfo_->onel_ints, CalcInfo_->twoel_ints, CalcInfo_->tf_onel_ints);
    form_gmat(CalcInfo_->onel_ints, CalcInfo_->twoel_ints, CalcInfo_->gmat);
    timer_off("CIWave: MCSCF integral transform");
}

void CIWavefunction::read_dpd_ci_ints() {
    // => Read one electron integrals <= //
    // Build temporary desired arrays
    int nmotri_full = (CalcInfo_->nmo * (CalcInfo_->nmo + 1)) / 2;
    double* tmp_onel_ints = new double[nmotri_full];

    // Read one electron integrals
    iwl_rdone(PSIF_OEI, PSIF_MO_FZC, tmp_onel_ints, nmotri_full, 0,
              (print_ > 4), "outfile");

    // IntegralTransform does not properly order one electron integrals for
    // whatever reason
    double* onel_ints = CalcInfo_->onel_ints->pointer();
    for (size_t i = 0, ij = 0; i < CalcInfo_->num_ci_orbs; i++) {
        for (size_t j = 0; j <= i; j++) {
            size_t si = i + CalcInfo_->num_drc_orbs;
            size_t sj = j + CalcInfo_->num_drc_orbs;
            size_t order_idx = INDEX(CalcInfo_->order[si], CalcInfo_->order[sj]);
            onel_ints[ij++] = tmp_onel_ints[order_idx];
        }
    }

    delete[] tmp_onel_ints;

    // => Read two electron integral <= //
    // This is not a good algorithm, stores two e's in memory twice!
    // However, this CI is still limited to 256 basis functions

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    dpdbuf4 I;
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ints_->DPD_ID("[X>=X]+"),
                           ints_->DPD_ID("[X>=X]+"), ints_->DPD_ID("[X>=X]+"),
                           ints_->DPD_ID("[X>=X]+"), 0, "MO Ints (XX|XX)");

    // Read everything size_to memory
    for (size_t h = 0; h < CalcInfo_->nirreps; h++) {
        global_dpd_->buf4_mat_irrep_init(&I, h);
        global_dpd_->buf4_mat_irrep_rd(&I, h);
    }

    double* twoel_intsp = CalcInfo_->twoel_ints->pointer();
    for (size_t p = 0; p < CalcInfo_->num_ci_orbs; p++) {
        size_t p_sym = I.params->psym[p];

        for (size_t q = 0; q <= p; q++) {
            size_t q_sym = I.params->qsym[q];
            size_t pq = I.params->rowidx[p][q];
            size_t pq_sym = p_sym ^ q_sym;
            size_t r_pq = INDEX(CalcInfo_->act_reorder[p], CalcInfo_->act_reorder[q]);

            for (size_t r = 0; r <= p; r++) {
                size_t r_sym = I.params->rsym[r];
                size_t smax = (p == r) ? q + 1 : r + 1;

                for (size_t s = 0; s < smax; s++) {
                    size_t s_sym = I.params->ssym[s];
                    size_t rs_sym = r_sym ^ s_sym;
                    if (pq_sym != rs_sym) continue;
                    size_t rs = I.params->colidx[r][s];

                    double value = I.matrix[pq_sym][pq][rs];
                    size_t r_rs = INDEX(CalcInfo_->act_reorder[r], CalcInfo_->act_reorder[s]);
                    size_t r_pqrs = INDEX(r_pq, r_rs);

                    twoel_intsp[r_pqrs] = value;
                }
            }
        }
    }

    // Close everything out
    for (size_t h = 0; h < CalcInfo_->nirreps; h++) {
        global_dpd_->buf4_mat_irrep_close(&I, h);
    }

    global_dpd_->buf4_close(&I);
    psio_->close(PSIF_LIBTRANS_DPD, 1);
}
void CIWavefunction::rotate_dfmcscf_twoel_ints(SharedMatrix Uact,
                                               SharedVector twoel_out) {
    // => Rotate twoel ints <= //
    int nQ = dferi_->size_Q();
    int nrot = CalcInfo_->num_rot_orbs;
    int nact = CalcInfo_->num_ci_orbs;
    int nav = nact + CalcInfo_->num_rsv_orbs;

    // Read RaQ
    std::shared_ptr<Tensor> RaQT = dferi_->ints()["RaQ"];
    SharedMatrix RaQ(new Matrix("RaQ", nrot, nact * nQ));

    double* RaQp = RaQ->pointer()[0];
    FILE* RaQF = RaQT->file_pointer();
    fseek(RaQF, 0L, SEEK_SET);
    fread(RaQp, sizeof(double), nrot * nact * nQ, RaQF);

    // We could slice it or... I like my raw GEMM
    // Uact_av DFERI_R_a_Q - > DFERI_a_aQ

    SharedMatrix dense_Uact = Uact->to_block_sharedmatrix();
    SharedMatrix tmp_rot_aaQ = Matrix::doublet(dense_Uact, RaQ);
    RaQ.reset();


    // Quv += Qvu
    SharedMatrix rot_aaQ(new Matrix("Rotated aaQ Matrix", nact * nact, nQ));
    double* tmp_rot_aaQp = tmp_rot_aaQ->pointer()[0];
    double* rot_aaQp = rot_aaQ->pointer()[0];
    for (size_t u = 0; u < nact; ++u) {
        for (size_t v = 0; v < nact; ++v) {
            for (size_t q = 0; q < nQ; ++q) {
                size_t uvq = u * nact * nQ + v * nQ + q;
                size_t vuq = v * nact * nQ + u * nQ + q;

                rot_aaQp[uvq] += tmp_rot_aaQp[uvq];
                rot_aaQp[uvq] += tmp_rot_aaQp[vuq];

            }
        }
    }

    // Read aaQ
    std::shared_ptr<Tensor> aaQT = dferi_->ints()["aaQ"];
    SharedMatrix aaQ(new Matrix("aaQ", nact * nact, nQ));

    double* aaQp = aaQ->pointer()[0];
    FILE* aaQF = aaQT->file_pointer();
    fseek(aaQF, 0L, SEEK_SET);
    fread(aaQp, sizeof(double), nact * nact * nQ, aaQF);


    // Form ERI's
    SharedMatrix rot_twoel = Matrix::doublet(rot_aaQ, aaQ, false, true);

    rot_aaQ.reset();
    aaQ.reset();

    // Add final symmetry
    rot_twoel->add(rot_twoel->transpose());

    pitzer_to_ci_order_twoel(rot_twoel, twoel_out);

}
void CIWavefunction::rotate_mcscf_twoel_ints(SharedMatrix Uact,
                                             SharedVector twoel_out) {


    // Eh, just read it into memory
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    dpdbuf4 I;
    global_dpd_->buf4_init(&I, PSIF_LIBTRANS_DPD, 0,
                           ints_->DPD_ID("[X>=X]+"), ints_->DPD_ID("[X,R]"),
                           ints_->DPD_ID("[X>=X]+"), ints_->DPD_ID("[X,R]"),
                           0, "MO Ints (XX|XR)");

    // Read everything size_to memory
    for (size_t h = 0; h < CalcInfo_->nirreps; h++) {
        global_dpd_->buf4_mat_irrep_init(&I, h);
        global_dpd_->buf4_mat_irrep_rd(&I, h);
    }

    int nrot = CalcInfo_->num_rot_orbs;
    int nact = CalcInfo_->num_ci_orbs;
    int nav = nact + CalcInfo_->num_rsv_orbs;

    SharedMatrix aaar(new Matrix("Tmp (aa|ar) Matrix", nact*nact*nact, nrot));

    double* aaarp = aaar->pointer()[0];
    for (size_t p = 0; p < nact; p++) {
        size_t p_sym = I.params->psym[p];

        for (size_t q = 0; q < nact; q++) {
            size_t q_sym = I.params->qsym[q];
            size_t pq = I.params->rowidx[p][q];
            size_t pq_sym = p_sym ^ q_sym;
            size_t left_idx = (p * nact + q) * nact * nrot;

            for (size_t r = 0; r < nact; r++) {
                size_t r_sym = I.params->rsym[r];
                size_t pqr_sym = p_sym ^ q_sym ^ r_sym;

                for (size_t s = 0; s < nrot; s++) {
                    size_t s_sym = I.params->ssym[s];
                    size_t rs_sym = r_sym ^ s_sym;
                    if (pq_sym != rs_sym) continue;
                    size_t rs = I.params->colidx[r][s];

                    double value = I.matrix[pq_sym][pq][rs];
                    size_t right_idx = r * nrot + s;

                    // outfile->Printf("%zd %zd %zd %zd | %zd %zd %zd %zd | %lf\n",
                    //                 p, q, r, rel_s, pq_sym, pq, rs, left_idx + right_idx, value);
                    aaarp[left_idx + right_idx] = value;
                }
            }
        }
    }

    // Close everything out
    for (size_t h = 0; h < CalcInfo_->nirreps; h++) {
        global_dpd_->buf4_mat_irrep_close(&I, h);
    }

    global_dpd_->buf4_close(&I);
    psio_->close(PSIF_LIBTRANS_DPD, 1);

    // Rotate the integrals
    SharedMatrix dense_Uact = Uact->to_block_sharedmatrix();
    SharedMatrix half_rot_aaaa = Matrix::doublet(aaar, dense_Uact, false, true);

    aaar.reset();

    // Use symmetry to add the rest
    SharedMatrix rot_twoel(new Matrix("Rotated aaaa Matrix", nact * nact, nact * nact));

    double* half_aaaap = half_rot_aaaa->pointer()[0];
    double* rot_twoelp = rot_twoel->pointer()[0];
    for (size_t u = 0; u < nact; ++u) {
        for (size_t v = 0; v < nact; ++v) {
            size_t left_idx = (u * nact + v) * nact * nact;
            for (size_t x = 0; x < nact; ++x) {
                for (size_t y = 0; y < nact; ++y) {
                    size_t uvxy = left_idx + x * nact + y;
                    size_t uvyx = left_idx + y * nact + x;

                    rot_twoelp[uvxy] += half_aaaap[uvxy];
                    rot_twoelp[uvxy] += half_aaaap[uvyx];
                }
            }
        }
    }
    // rot_twoel->print();

    rot_twoel->add(rot_twoel->transpose());
    // rot_twoel->print();


    pitzer_to_ci_order_twoel(rot_twoel, twoel_out);

}

double CIWavefunction::get_onel(int i, int j) {
    double value;
    if (i > j) {
        size_t ij = ioff[i] + j;
        value = CalcInfo_->onel_ints->get(ij);
    } else {
        size_t ij = ioff[j] + i;
        value = CalcInfo_->onel_ints->get(ij);
    }
    return value;
}

double CIWavefunction::get_twoel(int i, int j, int k, int l) {

    size_t ij = ioff[MAX0(i, j)] + MIN0(i, j);
    size_t kl = ioff[MAX0(k, l)] + MIN0(k, l);
    size_t ijkl = ioff[MAX0(ij, kl)] + MIN0(ij, kl);

    return CalcInfo_->twoel_ints->get(ijkl);
}

/*
** tf_onel_ints(): Function lumps together one-electron contributions
**    so that h'_{ij} = h_{ij} - 1/2 SUM_k (ik|kj)
**    The term h' arises in the calculation of sigma1 and sigma2 via
**    equation (20) of Olsen, Roos, et. al. JCP 1988
**
*/
void CIWavefunction::tf_onel_ints(SharedVector onel, SharedVector twoel, SharedVector output) {

    /* set up some shorthand notation */
    size_t nbf = CalcInfo_->num_ci_orbs;

    if ((output->dim(0) != CalcInfo_->num_ci_tri) || (output->nirrep() != 1)){
        throw PSIEXCEPTION("CIWavefunction::tf_onel_ints: output is not of the correct shape.");
    }

    /* ok, new special thing for CASSCF...if there are *no* excitations
       size_to restricted orbitals, and if Parameters_->fci=TRUE, then we
       do *not* want to sum over the restricted virts in h' or else
       we would need to account for RAS-out-of-space contributions
       (requiring fci=false).
     */
    if (Parameters_->fci && (nbf > Parameters_->ras3_lvl) && (Parameters_->ras34_max == 0)){
        nbf = Parameters_->ras3_lvl;
    }

    /* fill up the new array */
    double tval;
    for (size_t i = 0, ij = 0; i < nbf; i++) {
        for (size_t j = 0; j <= i; j++) {
            tval = onel->get(ij);

            for (size_t k = 0; k < nbf; k++) {
                size_t ik = ioff[MAX0(i, k)] + MIN0(i, k);
                size_t kj = ioff[MAX0(k, j)] + MIN0(k, j);
                size_t ikkj = ioff[ik] + kj;
                tval -= 0.5 * twoel->get(ikkj);
            }

            output->set(ij++, tval);
        }
    }
}

/*
** form_gmat(): Form the g matrix necessary for restriction to the RAS
**    subspaces (i.e. to eliminate contributions of out-of-space terms).
**    See equations (28-29) in Olsen, Roos, et. al. JCP 1988
**
*/
void CIWavefunction::form_gmat(SharedVector onel, SharedVector twoel, SharedVector output) {

    /* set up some shorthand notation */
    double tval;
    size_t i, j, k, ij, ii, ik, kj, ikkj, iiij;
    size_t nbf = CalcInfo_->num_ci_orbs;

    if ((output->dim(0) != (nbf * nbf)) || (output->nirrep() != 1)){
        throw PSIEXCEPTION("CIWavefunction::form_gmat: output is not of the correct shape.");
    }

    /* fill up the new array */
    for (i = 0; i < nbf; i++) {
        for (j = i + 1; j < nbf; j++) {
            ij = ioff[j] + i;
            tval = onel->get(ij);
            for (k = 0; k < i; k++) {
                ik = ioff[i] + k;
                kj = ioff[j] + k;
                ikkj = ioff[kj] + ik;
                tval -= twoel->get(ikkj);
            }
            output->set(i * nbf + j, tval);
        }
    }

    for (i = 0, ij = 0; i < nbf; i++) {
        for (j = 0; j <= i; j++, ij++) {
            tval = onel->get(ij);
            for (k = 0; k < i; k++) {
                ik = ioff[i] + k;
                kj = ioff[MAX0(k, j)] + MIN0(k, j);
                ikkj = ioff[ik] + kj;
                tval -= twoel->get(ikkj);
            }
            ii = ioff[i] + i;
            iiij = ioff[ii] + ij;
            if (i == j)
                tval -= 0.5 * twoel->get(iiij);
            else
                tval -= twoel->get(iiij);
            output->set(i * nbf + j, tval);
        }
    }
}

void CIWavefunction::onel_ints_from_jk() {

    // Compute the frozen fock if needed
    // This is the *only* place where the frozen fock enters
    if ((!fzc_fock_computed_) && (CalcInfo_->frozen_docc.sum() > 0)){
        std::vector<SharedMatrix>& Cl = jk_->C_left();
        std::vector<SharedMatrix>& Cr = jk_->C_right();
        Cl.clear();
        Cr.clear();

        Cl.push_back(get_orbitals("FZC"));
        jk_->compute();
        Cl.clear();

        const std::vector<SharedMatrix>& J = jk_->J();
        const std::vector<SharedMatrix>& K = jk_->K();

        J[0]->scale(2.0);
        J[0]->subtract(K[0]);

        CalcInfo_->fzc_so_onel_ints = J[0]->clone();
        fzc_fock_computed_ = true;

    }

    SharedMatrix Cact = get_orbitals("ACT");
    SharedMatrix Cdocc = get_orbitals("DOCC");

    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    Cl.clear();
    Cr.clear();

    Cl.push_back(Cdocc);
    jk_->compute();
    Cl.clear();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    J[0]->scale(2.0);
    J[0]->subtract(K[0]);
    if (CalcInfo_->frozen_docc.sum() > 0){
        J[0]->add(CalcInfo_->fzc_so_onel_ints);
    }

    J[0]->add(H_);


    CalcInfo_->so_onel_ints = J[0]->clone();
    if (mcscf_object_init_){
        somcscf_->set_AO_IFock(J[0]);
    }
    SharedMatrix onel_ints = Matrix::triplet(Cact, J[0], Cact, true, false, false);

    // Set 1D onel ints
    pitzer_to_ci_order_onel(onel_ints, CalcInfo_->onel_ints);

    // => Compute dropped core energy <= //
    J[0]->add(H_);

    SharedMatrix Cdrc = get_orbitals("DRC");
    SharedMatrix D = Matrix::doublet(Cdrc, Cdrc, false, true);
    CalcInfo_->edrc = J[0]->vector_dot(D);
}

void CIWavefunction::pitzer_to_ci_order_onel(SharedMatrix src, SharedVector dest){
    if ((src->nirrep() != nirrep_) || dest->nirrep() != 1){
        throw PSIEXCEPTION("CIWavefunciton::pitzer_to_ci_order_onel irreps are not of the correct size.");
    }

    if (dest->dim(0) != CalcInfo_->num_ci_tri){
        throw PSIEXCEPTION("CIWavefunciton::pitzer_to_ci_order_onel: Destination vector must be of size ncitri.");
    }


    double* destp = dest->pointer();
    for (size_t h = 0, target = 0, offset = 0; h < nirrep_; h++) {
        size_t nactpih = CalcInfo_->ci_orbs[h];
        if (!nactpih) continue;

        double** srcp = src->pointer(h);
        for (size_t i = 0; i < nactpih; i++) {
            target += offset;
            for (size_t j = 0; j <= i; j++) {
                size_t r_ij = INDEX(CalcInfo_->act_reorder[i + offset],
                                    CalcInfo_->act_reorder[j + offset]);
                destp[r_ij] = srcp[i][j];
            }
        }
        offset += nactpih;
    }
}
void CIWavefunction::pitzer_to_ci_order_twoel(SharedMatrix src, SharedVector dest){
    if ((src->nirrep() != 1) || dest->nirrep() != 1){
        throw PSIEXCEPTION("CIWavefunciton::pitzer_to_ci_order_twoel irreped matrices are not supported.");
    }
    if (dest->dim(0) != CalcInfo_->num_ci_tri2){
        throw PSIEXCEPTION("CIWavefunciton::pitzer_to_ci_order_onel: Destination vector must be of size ncitri2.");
    }

    double** srcp = src->pointer();
    double* destp = dest->pointer();
    size_t nact = CalcInfo_->num_ci_orbs;
    for (size_t i = 0; i < nact; i++) {
        size_t irel = CalcInfo_->act_reorder[i];

        for (size_t j = 0; j <= i; j++) {
            size_t jrel = CalcInfo_->act_reorder[j];
            size_t ijrel = INDEX(irel, jrel);

            for (size_t k = 0; k <= i; k++) {
                size_t krel = CalcInfo_->act_reorder[k];
                size_t lmax = (i == k) ? j + 1 : k + 1;

                for (size_t l = 0; l < lmax; l++) {
                    size_t lrel = CalcInfo_->act_reorder[l];
                    size_t klrel = INDEX(krel, lrel);
                    destp[INDEX(ijrel, klrel)] = srcp[i * nact + j][k * nact + l];
                }
            }
        }
    }
}

}} // namespace psi::detci
