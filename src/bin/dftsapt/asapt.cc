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

#include "dftsapt.h"
#include "vis.h"
#include "atomic.h"
#include <libmints/mints.h>
#include <libmints/sieve.h>
#include <libmints/local.h>
#include <libfock/jk.h>
#include <libfock/points.h>
#include <libfock/cubature.h>
#include <libthce/thce.h>
#include <libthce/lreri.h>
#include <libqt/qt.h>
#include <libpsio/psio.hpp>
#include <libpsio/psio.h>
#include <psi4-dec.h>
#include <psifiles.h>
#include <physconst.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace boost;
using namespace std;

namespace psi {
namespace dftsapt {

ASAPT::ASAPT()
{
    common_init();
}
ASAPT::~ASAPT()
{
}
void ASAPT::common_init()
{
}
boost::shared_ptr<ASAPT> ASAPT::build(boost::shared_ptr<Wavefunction> d,
                                      boost::shared_ptr<Wavefunction> mA,
                                      boost::shared_ptr<Wavefunction> mB)
{
    Options& options = Process::environment.options;

    ASAPT* sapt = new ASAPT();

    sapt->exch_scale_      = options.get_bool("ASAPT_EXCH_SCALE");
    sapt->ind_scale_       = options.get_bool("ASAPT_IND_SCALE");
    sapt->ind_resp_        = options.get_bool("ASAPT_IND_RESPONSE");
    sapt->sep_core_        = options.get_bool("ASAPT_SEPARATE_CORE");
    
    sapt->print_ = options.get_int("PRINT");
    sapt->debug_ = options.get_int("DEBUG");
    sapt->bench_ = options.get_int("BENCH");

    sapt->memory_ = (unsigned long int)(Process::environment.get_memory() * options.get_double("SAPT_MEM_FACTOR") * 0.125);

    sapt->cpks_maxiter_ = options.get_int("MAXITER");
    sapt->cpks_delta_ = options.get_double("D_CONVERGENCE");

    sapt->dimer_     = d->molecule();
    sapt->monomer_A_ = mA->molecule();
    sapt->monomer_B_ = mB->molecule();

    sapt->E_dimer_     = d->reference_energy();
    sapt->E_monomer_A_ = mA->reference_energy();
    sapt->E_monomer_B_ = mB->reference_energy();

    sapt->primary_   = d->basisset();
    sapt->primary_A_ = mA->basisset();
    sapt->primary_B_ = mB->basisset();

    if (sapt->primary_A_->nbf() != sapt->primary_B_->nbf() || sapt->primary_->nbf() != sapt->primary_A_->nbf()) {
        throw PSIEXCEPTION("Monomer-centered bases not allowed in DFT-SAPT");
    }

    sapt->Cocc_A_     = mA->Ca_subset("AO","OCC");
    sapt->Cvir_A_     = mA->Ca_subset("AO","VIR");
    sapt->eps_occ_A_  = mA->epsilon_a_subset("AO","OCC");
    sapt->eps_vir_A_  = mA->epsilon_a_subset("AO","VIR");

    sapt->Cfocc_A_    = mA->Ca_subset("AO","FROZEN_OCC");
    sapt->Caocc_A_    = mA->Ca_subset("AO","ACTIVE_OCC");
    sapt->Cavir_A_    = mA->Ca_subset("AO","ACTIVE_VIR");
    sapt->Cfvir_A_    = mA->Ca_subset("AO","FROZEN_VIR");

    sapt->eps_focc_A_ = mA->epsilon_a_subset("AO","FROZEN_OCC");
    sapt->eps_aocc_A_ = mA->epsilon_a_subset("AO","ACTIVE_OCC");
    sapt->eps_avir_A_ = mA->epsilon_a_subset("AO","ACTIVE_VIR");
    sapt->eps_fvir_A_ = mA->epsilon_a_subset("AO","FROZEN_VIR");

    sapt->Cocc_B_     = mB->Ca_subset("AO","OCC");
    sapt->Cvir_B_     = mB->Ca_subset("AO","VIR");
    sapt->eps_occ_B_  = mB->epsilon_a_subset("AO","OCC");
    sapt->eps_vir_B_  = mB->epsilon_a_subset("AO","VIR");

    sapt->Cfocc_B_    = mB->Ca_subset("AO","FROZEN_OCC");
    sapt->Caocc_B_    = mB->Ca_subset("AO","ACTIVE_OCC");
    sapt->Cavir_B_    = mB->Ca_subset("AO","ACTIVE_VIR");
    sapt->Cfvir_B_    = mB->Ca_subset("AO","FROZEN_VIR");

    sapt->eps_focc_B_ = mB->epsilon_a_subset("AO","FROZEN_OCC");
    sapt->eps_aocc_B_ = mB->epsilon_a_subset("AO","ACTIVE_OCC");
    sapt->eps_avir_B_ = mB->epsilon_a_subset("AO","ACTIVE_VIR");
    sapt->eps_fvir_B_ = mB->epsilon_a_subset("AO","FROZEN_VIR");

    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    sapt->mp2fit_ = BasisSet::construct(parser, sapt->dimer_, "DF_BASIS_SAPT");

    return boost::shared_ptr<ASAPT>(sapt);
}
double ASAPT::compute_energy()
{
    energies_["HF"] = E_dimer_ - E_monomer_A_ - E_monomer_B_; // TODO: get dHF loaded correctly

    print_header();

    fock_terms();

    atomize();

    localize();

    populate();

    elst();
    
    exch();

    ind();

    disp();

    analyze();

    print_trailer();

    Process::environment.globals["SAPT ENERGY"] = 0.0;

    return 0.0;
}
void ASAPT::print_header() const
{
    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\t                      A-SAPT Analysis                    \n");
    fprintf(outfile, "\t                        Rob Parrish                      \n");
    fprintf(outfile, "\t --------------------------------------------------------\n");
    fprintf(outfile, "\n");

    fprintf(outfile, "  ==> Sizes <==\n");
    fprintf(outfile, "\n");

    fprintf(outfile, "   => Resources <=\n\n");

    fprintf(outfile, "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
    fprintf(outfile, "\n");

    fprintf(outfile, "   => Orbital Ranges <=\n\n");

    int nmo_A = eps_focc_A_->dim() + eps_aocc_A_->dim() + eps_avir_A_->dim() + eps_fvir_A_->dim();
    int nmo_B = eps_focc_B_->dim() + eps_aocc_B_->dim() + eps_avir_B_->dim() + eps_fvir_B_->dim();

    int nA = 0;
    for (int A = 0; A < monomer_A_->natom(); A++) {
        if (monomer_A_->Z(A) != 0.0) nA++;
    }

    int nB = 0;
    for (int B = 0; B < monomer_B_->natom(); B++) {
        if (monomer_B_->Z(B) != 0.0) nB++;
    }

    fprintf(outfile, "    ------------------\n");
    fprintf(outfile, "    %-6s %5s %5s\n", "Range", "M_A", "M_B");
    fprintf(outfile, "    ------------------\n");
    fprintf(outfile, "    %-6s %5d %5d\n", "natom", nA, nB);
    fprintf(outfile, "    %-6s %5d %5d\n", "nso", primary_A_->nbf(), primary_B_->nbf());
    fprintf(outfile, "    %-6s %5d %5d\n", "nmo", nmo_A, nmo_B);
    fprintf(outfile, "    %-6s %5d %5d\n", "nocc", eps_aocc_A_->dim() + eps_focc_A_->dim(), eps_aocc_B_->dim() + eps_focc_B_->dim());
    fprintf(outfile, "    %-6s %5d %5d\n", "nvir", eps_avir_A_->dim() + eps_fvir_A_->dim(), eps_avir_B_->dim() + eps_fvir_B_->dim());
    fprintf(outfile, "    %-6s %5d %5d\n", "nfocc", eps_focc_A_->dim(), eps_focc_B_->dim());
    fprintf(outfile, "    %-6s %5d %5d\n", "naocc", eps_aocc_A_->dim(), eps_aocc_B_->dim());
    fprintf(outfile, "    %-6s %5d %5d\n", "navir", eps_avir_A_->dim(), eps_avir_B_->dim());
    fprintf(outfile, "    %-6s %5d %5d\n", "nfvir", eps_fvir_A_->dim(), eps_fvir_B_->dim());
    fprintf(outfile, "    ------------------\n");
    fprintf(outfile, "\n");

    fprintf(outfile, "   => Primary Basis Set <=\n\n");
    primary_->print_by_level(outfile, print_);

    fflush(outfile);
}
void ASAPT::print_trailer()
{
    energies_["delta HF,r (2)"] = 0.0;
    if (energies_["HF"] != 0.0) {
        energies_["delta HF,r (2)"] = energies_["HF"] - energies_["Elst10,r"] - energies_["Exch10"] - energies_["Ind20,r"] - energies_["Exch-Ind20,r"];
    }

    energies_["Electrostatics"] = energies_["Elst10,r"];
    energies_["Exchange"]       = energies_["Exch10"];
    energies_["Induction"]      = energies_["Ind20,r"] + energies_["Exch-Ind20,r"] + energies_["delta HF,r (2)"];
    energies_["Dispersion"]     = energies_["Disp20"] + energies_["Exch-Disp20"];
    energies_["SAPT"]           = energies_["Electrostatics"] + energies_["Exchange"] + energies_["Induction"] + energies_["Dispersion"];

    fprintf(outfile,"\n    SAPT Results  \n");
    fprintf(outfile,"  -----------------------------------------------------------------------\n");
    fprintf(outfile,"    Electrostatics     %16.8lf mH %16.8lf kcal mol^-1\n",
      energies_["Electrostatics"]*1000.0,energies_["Electrostatics"]*pc_hartree2kcalmol);
    fprintf(outfile,"      Elst10,r         %16.8lf mH %16.8lf kcal mol^-1\n\n",
      energies_["Elst10,r"]*1000.0,energies_["Elst10,r"]*pc_hartree2kcalmol);
    fprintf(outfile,"    Exchange           %16.8lf mH %16.8lf kcal mol^-1\n",
      energies_["Exchange"]*1000.0,energies_["Exchange"]*pc_hartree2kcalmol);
    fprintf(outfile,"      Exch10           %16.8lf mH %16.8lf kcal mol^-1\n",
      energies_["Exch10"]*1000.0,energies_["Exch10"]*pc_hartree2kcalmol);
    fprintf(outfile,"      Exch10(S^2)      %16.8lf mH %16.8lf kcal mol^-1\n\n",
      energies_["Exch10(S^2)"]*1000.0,energies_["Exch10(S^2)"]*pc_hartree2kcalmol);
    fprintf(outfile,"    Induction          %16.8lf mH %16.8lf kcal mol^-1\n",
      energies_["Induction"]*1000.0,energies_["Induction"]*pc_hartree2kcalmol);
    fprintf(outfile,"      Ind20,r          %16.8lf mH %16.8lf kcal mol^-1\n",
      energies_["Ind20,r"]*1000.0,energies_["Ind20,r"]*pc_hartree2kcalmol);
    fprintf(outfile,"      Exch-Ind20,r     %16.8lf mH %16.8lf kcal mol^-1\n",
      energies_["Exch-Ind20,r"]*1000.0,energies_["Exch-Ind20,r"]*pc_hartree2kcalmol);
    fprintf(outfile,"      delta HF,r (2)   %16.8lf mH %16.8lf kcal mol^-1\n\n",
      energies_["delta HF,r (2)"]*1000.0,energies_["delta HF,r (2)"]*pc_hartree2kcalmol);
    fprintf(outfile,"    Dispersion         %16.8lf mH %16.8lf kcal mol^-1\n",
      energies_["Dispersion"]*1000.0,energies_["Dispersion"]*pc_hartree2kcalmol);
    fprintf(outfile,"      Disp20           %16.8lf mH %16.8lf kcal mol^-1\n",
      energies_["Disp20"]*1000.0,energies_["Disp20"]*pc_hartree2kcalmol);
    fprintf(outfile,"      Exch-Disp20      %16.8lf mH %16.8lf kcal mol^-1\n\n",
      energies_["Exch-Disp20"]*1000.0,energies_["Exch-Disp20"]*pc_hartree2kcalmol);

    fprintf(outfile,"    Total HF           %16.8lf mH %16.8lf kcal mol^-1\n",
      energies_["HF"]*1000.0,energies_["HF"]*pc_hartree2kcalmol);
    fprintf(outfile,"    Total SAPT0        %16.8lf mH %16.8lf kcal mol^-1\n",
      energies_["SAPT"]*1000.0,energies_["SAPT"]*pc_hartree2kcalmol);
    fprintf(outfile,"\n");

    fprintf(outfile, "  To the pessimist, the glass is half empty.\n");
    fprintf(outfile, "  To the optimist, the glass is half full.\n");
    fprintf(outfile, "  To the engineer, the glass is twice as big as it needs to be.\n");
    fprintf(outfile, "\n");
}
void ASAPT::atomize()
{
    fprintf(outfile," ATOMIZATION:\n\n");
    
    atomic_A_ = AtomicDensity::build("STOCKHOLDER", primary_A_, Process::environment.options);
    atomic_A_->compute(Matrix::doublet(Cocc_A_,Cocc_A_,false,true));
    atomic_A_->compute_charges();

    atomic_B_ = AtomicDensity::build("STOCKHOLDER", primary_B_, Process::environment.options);
    atomic_B_->compute(Matrix::doublet(Cocc_B_,Cocc_B_,false,true));
    atomic_B_->compute_charges();
}
void ASAPT::localize()
{
    fprintf(outfile," LOCALIZATION:\n\n");

    if (sep_core_) {

        fprintf(outfile,"  Local Core Orbitals for Monomer A:\n\n");
        boost::shared_ptr<Localizer> localfA = Localizer::build(primary_, Cfocc_A_, Process::environment.options);
        localfA->localize();
        boost::shared_ptr<Matrix> Lfocc_A = localfA->L();
        boost::shared_ptr<Matrix> Ufocc_A = localfA->U();

        int nfA = eps_focc_A_->dimpi()[0];
        boost::shared_ptr<Matrix> FcfA(new Matrix("FcfA", nfA, nfA));
        FcfA->set_diagonal(eps_focc_A_); 
        boost::shared_ptr<Matrix> FlfA = localfA->fock_update(FcfA);
        FlfA->set_name("FlfA");    
        //FcfA->print();
        //FlfA->print();

        fprintf(outfile,"  Local Valence Orbitals for Monomer A:\n\n");
        boost::shared_ptr<Localizer> localaA = Localizer::build(primary_, Caocc_A_, Process::environment.options);
        localaA->localize();
        boost::shared_ptr<Matrix> Laocc_A = localaA->L();
        boost::shared_ptr<Matrix> Uaocc_A = localaA->U();

        int naA = eps_aocc_A_->dimpi()[0];
        boost::shared_ptr<Matrix> FcaA(new Matrix("FcaA", naA, naA));
        FcaA->set_diagonal(eps_aocc_A_); 
        boost::shared_ptr<Matrix> FlaA = localaA->fock_update(FcaA);
        FlaA->set_name("FlaA");    
        //FcaA->print();
        //FlaA->print();

        std::vector<boost::shared_ptr<Matrix> > LAlist;
        LAlist.push_back(Lfocc_A);
        LAlist.push_back(Laocc_A);
        Locc_A_ = Matrix::horzcat(LAlist); 

        Uocc_A_ = boost::shared_ptr<Matrix>(new Matrix("Uocc A", nfA + naA, nfA + naA));
        double** Uocc_Ap = Uocc_A_->pointer();
        double** Ufocc_Ap = Ufocc_A->pointer();
        double** Uaocc_Ap = Uaocc_A->pointer();
        for (int i = 0; i < nfA; i++) {
            ::memcpy(&Uocc_Ap[i][0],&Ufocc_Ap[i][0],sizeof(double)*nfA);
        }
        for (int i = 0; i < naA; i++) {
            ::memcpy(&Uocc_Ap[i+nfA][nfA],&Uaocc_Ap[i][0],sizeof(double)*naA);
        }

        fprintf(outfile,"  Local Core Orbitals for Monomer B:\n\n");
        boost::shared_ptr<Localizer> localfB = Localizer::build(primary_, Cfocc_B_, Process::environment.options);
        localfB->localize();
        boost::shared_ptr<Matrix> Lfocc_B = localfB->L();
        boost::shared_ptr<Matrix> Ufocc_B = localfB->U();

        int nfB = eps_focc_B_->dimpi()[0];
        boost::shared_ptr<Matrix> FcfB(new Matrix("FcfB", nfB, nfB));
        FcfB->set_diagonal(eps_focc_B_); 
        boost::shared_ptr<Matrix> FlfB = localfB->fock_update(FcfB);
        FlfB->set_name("FlfB");    
        //FcfB->print();
        //FlfB->print();

        fprintf(outfile,"  Local Valence Orbitals for Monomer B:\n\n");
        boost::shared_ptr<Localizer> localaB = Localizer::build(primary_, Caocc_B_, Process::environment.options);
        localaB->localize();
        boost::shared_ptr<Matrix> Laocc_B = localaB->L();
        boost::shared_ptr<Matrix> Uaocc_B = localaB->U();

        int naB = eps_aocc_B_->dimpi()[0];
        boost::shared_ptr<Matrix> FcaB(new Matrix("FcaA", naB, naB));
        FcaB->set_diagonal(eps_aocc_B_); 
        boost::shared_ptr<Matrix> FlaB = localaB->fock_update(FcaB);
        FlaB->set_name("FlaB");    
        //FcaB->print();
        //FlaB->print();

        std::vector<boost::shared_ptr<Matrix> > LBlist;
        LBlist.push_back(Lfocc_B);
        LBlist.push_back(Laocc_B);
        Locc_B_ = Matrix::horzcat(LBlist); 

        Uocc_B_ = boost::shared_ptr<Matrix>(new Matrix("Uocc B", nfB + naB, nfB + naB));
        double** Uocc_Bp = Uocc_B_->pointer();
        double** Ufocc_Bp = Ufocc_B->pointer();
        double** Uaocc_Bp = Uaocc_B->pointer();
        for (int i = 0; i < nfB; i++) {
            ::memcpy(&Uocc_Bp[i][0],&Ufocc_Bp[i][0],sizeof(double)*nfB);
        }
        for (int i = 0; i < naB; i++) {
            ::memcpy(&Uocc_Bp[i+nfB][nfB],&Uaocc_Bp[i][0],sizeof(double)*naB);
        }

    } else {

        fprintf(outfile,"  Local Orbitals for Monomer A:\n\n");
        boost::shared_ptr<Localizer> localA = Localizer::build(primary_, Cocc_A_, Process::environment.options);
        localA->localize();
        Locc_A_ = localA->L();
        Uocc_A_ = localA->U();

        int na = eps_occ_A_->dimpi()[0];
        boost::shared_ptr<Matrix> FcA(new Matrix("FcA", na, na));
        FcA->set_diagonal(eps_occ_A_); 
        boost::shared_ptr<Matrix> FlA = localA->fock_update(FcA);
        FlA->set_name("FlA");    

        //FcA->print();
        //FlA->print();

        fprintf(outfile,"  Local Orbitals for Monomer B:\n\n");
        boost::shared_ptr<Localizer> localB = Localizer::build(primary_, Cocc_B_, Process::environment.options);
        localB->localize();
        Locc_B_ = localB->L();
        Uocc_B_ = localB->U();

        int nb = eps_occ_B_->dimpi()[0];
        boost::shared_ptr<Matrix> FcB(new Matrix("FcB", nb, nb));
        FcB->set_diagonal(eps_occ_B_); 
        boost::shared_ptr<Matrix> FlB = localB->fock_update(FcB);
        FlB->set_name("FlB");    

        //FcB->print();
        //FlB->print();

    }

    fflush(outfile);
}
void ASAPT::populate()
{
    fprintf(outfile,"  POPULATION:\n\n");

    // => Sizing <= //

    int na = Cocc_A_->colspi()[0];
    int nb = Cocc_B_->colspi()[0];

    int nA = 0;
    std::vector<int> cA;
    for (int A = 0; A < monomer_A_->natom(); A++) {
        if (monomer_A_->Z(A) != 0.0) {
            nA++;
            cA.push_back(A);
        }
    }

    int nB = 0;
    std::vector<int> cB;
    for (int B = 0; B < monomer_B_->natom(); B++) {
        if (monomer_B_->Z(B) != 0.0) {
            nB++;
            cB.push_back(B);
        }
    }

    // => Targets <= //

    boost::shared_ptr<Matrix> Q2A(new Matrix("QA (A x a)", nA, na));
    boost::shared_ptr<Matrix> Q2B(new Matrix("QB (B x b)", nB, nb));
    double** Q2Ap = Q2A->pointer();
    double** Q2Bp = Q2B->pointer();

    // => Molecular Grid <= //
    
    boost::shared_ptr<MolecularGrid> grid = atomic_A_->grid();
    int max_points = grid->max_points();
    double* xp = grid->x();
    double* yp = grid->y();
    double* zp = grid->z();
    double* wp = grid->w();
        
    // => Temps <= //
    
    boost::shared_ptr<Matrix> WA(new Matrix("WA", nA, max_points));
    boost::shared_ptr<Matrix> WB(new Matrix("WB", nB, max_points));
    double** WAp = WA->pointer();
    double** WBp = WB->pointer();
    
    // => Local orbital collocation computer <= //

    boost::shared_ptr<UKSFunctions> points = boost::shared_ptr<UKSFunctions>(new UKSFunctions(primary_,grid->max_points(),grid->max_functions()));
    points->set_ansatz(0);
    points->set_Cs(Locc_A_,Locc_B_);
    double** psiap = points->orbital_value("PSI_A")->pointer();
    double** psibp = points->orbital_value("PSI_B")->pointer();

    // => Master Loop <= //

    const std::vector<boost::shared_ptr<BlockOPoints> >& blocks = grid->blocks();
    size_t offset = 0L;
    for (int ind = 0; ind < blocks.size(); ind++) {
        int npoints = blocks[ind]->npoints();
        points->compute_orbitals(blocks[ind]);
        for (int a = 0; a < na; a++) {
            for (int P = 0; P < npoints; P++) {
                psiap[a][P] *= psiap[a][P];
            }
        }
        for (int b = 0; b < nb; b++) {
            for (int P = 0; P < npoints; P++) {
                psibp[b][P] *= psibp[b][P];
            }
        }
        atomic_A_->compute_weights(npoints,&xp[offset],&yp[offset],&zp[offset],WAp,&wp[offset]);
        atomic_B_->compute_weights(npoints,&xp[offset],&yp[offset],&zp[offset],WBp,&wp[offset]);
        C_DGEMM('N','T',nA,na,npoints,1.0,WAp[0],max_points,psiap[0],max_points,1.0,Q2Ap[0],na); 
        C_DGEMM('N','T',nB,nb,npoints,1.0,WBp[0],max_points,psibp[0],max_points,1.0,Q2Bp[0],nb); 
        offset += npoints;
    }

    fprintf(outfile,"    Grid-based orbital charges.\n\n");

    fprintf(outfile,"    Monomer A Grid Errors (Orbital Normalizations):\n");
    for (int a = 0; a < na; a++) {
        double val = 0.0;
        for (int A = 0; A < nA; A++) {
            val += Q2Ap[A][a]; 
        } 
        C_DSCAL(nA,1.0/val,&Q2Ap[0][a],na);
        fprintf(outfile,"    %4d: %11.3E\n", a+1, fabs(1.0 - val));
    }
    fprintf(outfile,"\n");

    fprintf(outfile,"    Monomer B Grid Errors (Orbital Normalizations):\n");
    for (int b = 0; b < nb; b++) {
        double val = 0.0;
        for (int B = 0; B < nB; B++) {
            val += Q2Bp[B][b]; 
        } 
        C_DSCAL(nB,1.0/val,&Q2Bp[0][b],nb);
        fprintf(outfile,"    %4d: %11.3E\n", b+1, fabs(1.0 - val));
    }
    fprintf(outfile,"\n");

    fprintf(outfile,"    Orbital charges renormalized.\n\n");

    // => Globals <= //

    Q_A_ = Q2A; 
    Q_B_ = Q2B; 
    
    R_A_ = Matrix::doublet(Q_A_,Uocc_A_,false,true);
    R_B_ = Matrix::doublet(Q_B_,Uocc_B_,false,true);

    // The ASAPT visualization and analysis container
    vis_ = boost::shared_ptr<ASAPTVis>(new ASAPTVis(primary_,monomer_A_,monomer_B_,Locc_A_,Locc_B_,Q_A_,Q_B_,atomic_A_,atomic_B_));
}
void ASAPT::elst()
{
    fprintf(outfile,"  ELECTROSTATICS:\n\n");

    // ==> Sizing <== //

    int nn  = primary_->nbf();

    int na = Cocc_A_->colspi()[0];
    int nb = Cocc_B_->colspi()[0];

    int nA = 0;
    std::vector<int> cA;
    for (int A = 0; A < monomer_A_->natom(); A++) {
        if (monomer_A_->Z(A) != 0.0) {
            nA++;
            cA.push_back(A);
        }
    }

    int nB = 0;
    std::vector<int> cB;
    for (int B = 0; B < monomer_B_->natom(); B++) {
        if (monomer_B_->Z(B) != 0.0) {
            nB++;
            cB.push_back(B);
        }
    }

    int nr = Cvir_A_->colspi()[0];
    int ns = Cvir_B_->colspi()[0];

    // ==> DF ERI Setup (JKFIT Type, in Full Basis) <== //

    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    boost::shared_ptr<BasisSet> jkfit = BasisSet::construct(parser, primary_->molecule(), "DF_BASIS_SCF");
    int nQ = jkfit->nbf();

    boost::shared_ptr<DFERI> df = DFERI::build(primary_,jkfit,Process::environment.options);
    df->clear();

    std::vector<boost::shared_ptr<Matrix> > Cs;
    Cs.push_back(Cocc_A_);
    Cs.push_back(Cvir_A_);
    Cs.push_back(Cocc_B_);
    Cs.push_back(Cvir_B_);
    boost::shared_ptr<Matrix> Call = Matrix::horzcat(Cs);
    Cs.clear();

    df->set_C(Call);
    df->set_memory(memory_);

    int offset = 0;
    df->add_space("a",offset,offset+Cocc_A_->colspi()[0]); offset += Cocc_A_->colspi()[0];
    df->add_space("r",offset,offset+Cvir_A_->colspi()[0]); offset += Cvir_A_->colspi()[0];
    df->add_space("b",offset,offset+Cocc_B_->colspi()[0]); offset += Cocc_B_->colspi()[0];
    df->add_space("s",offset,offset+Cvir_B_->colspi()[0]); offset += Cvir_B_->colspi()[0];

    df->add_pair_space("Aar", "a", "r");
    df->add_pair_space("Abs", "b", "s");

    df->print_header();
    df->compute();

    std::map<std::string, boost::shared_ptr<Tensor> >& ints = df->ints();
    boost::shared_ptr<Tensor> AarT = ints["Aar"];
    boost::shared_ptr<Tensor> AbsT = ints["Abs"];

    boost::shared_ptr<Matrix> Jm12 = df->Jpow(-1.0/2.0);

    df.reset();

    // ==> Auxiliary basis representation of atomic ESP <== //
    
    // Grid Definition
    boost::shared_ptr<Vector> x = atomic_A_->x();
    boost::shared_ptr<Vector> y = atomic_A_->y();
    boost::shared_ptr<Vector> z = atomic_A_->z();
    boost::shared_ptr<Vector> w = atomic_A_->w();
    boost::shared_ptr<Vector> rhoa = atomic_A_->rho();
    boost::shared_ptr<Vector> rhob = atomic_B_->rho();
    int nP = x->dimpi()[0];
    int max_points = atomic_A_->grid()->max_points();
    double* xp = x->pointer();
    double* yp = y->pointer();
    double* zp = z->pointer();
    double* wp = w->pointer();
    double* rhoap = rhoa->pointer();
    double* rhobp = rhob->pointer();

    // Grid charges (blocked. Nuke 'em Rico!)
    boost::shared_ptr<Matrix> QAP(new Matrix("Q_A_P", nA, max_points));
    boost::shared_ptr<Matrix> QBQ(new Matrix("Q_B_Q", nB, max_points));
    double** QAPp = QAP->pointer();
    double** QBQp = QBQ->pointer();

    // Targets
    boost::shared_ptr<Matrix> QAC(new Matrix("Q_A^C", nA, nQ));
    boost::shared_ptr<Matrix> QBD(new Matrix("Q_B^D", nB, nQ));
    double** QACp = QAC->pointer();
    double** QBDp = QBD->pointer();

    boost::shared_ptr<Matrix> VAB(new Matrix("V_A^B", nA, nB));
    boost::shared_ptr<Matrix> VBA(new Matrix("V_B^A", nB, nA));
    double** VABp = VAB->pointer();
    double** VBAp = VBA->pointer();

    // Threading 
    int nthreads = 0;
    #ifdef _OPENMP
        nthreads = omp_get_max_threads();
    #endif

    // Integral computers and thread-safe targets
    boost::shared_ptr<IntegralFactory> Vfact(new IntegralFactory(jkfit,BasisSet::zero_ao_basis_set()));
    std::vector<boost::shared_ptr<Matrix> > ZxyzT;
    std::vector<boost::shared_ptr<Matrix> > VtempT;
    std::vector<boost::shared_ptr<Matrix> > QACT;
    std::vector<boost::shared_ptr<Matrix> > QBDT;
    std::vector<boost::shared_ptr<PotentialInt> > VintT;
    for (int thread = 0; thread < nthreads; thread++) {
        ZxyzT.push_back(boost::shared_ptr<Matrix>(new Matrix("Zxyz",1,4)));
        VtempT.push_back(boost::shared_ptr<Matrix>(new Matrix("Vtemp",nQ,1)));
        QACT.push_back(boost::shared_ptr<Matrix>(new Matrix("QACT",nA,nQ)));
        QBDT.push_back(boost::shared_ptr<Matrix>(new Matrix("QBDT",nB,nQ)));
        VintT.push_back(boost::shared_ptr<PotentialInt>(static_cast<PotentialInt*>(Vfact->ao_potential())));
        VintT[thread]->set_charge_field(ZxyzT[thread]);
    }

    boost::shared_ptr<Matrix> Zxyz(new Matrix("Zxyz",1,4));
    double** Zxyzp = Zxyz->pointer();
    boost::shared_ptr<PotentialInt> Vint(static_cast<PotentialInt*>(Vfact->ao_potential()));
    Vint->set_charge_field(Zxyz);
    boost::shared_ptr<Matrix> Vtemp(new Matrix("Vtemp",nQ,1));
    double** Vtempp = Vtemp->pointer();

    // Master loop 
    for (int offset = 0; offset < nP; offset += max_points) {
        int npoints = (offset + max_points >= nP ? nP - offset : max_points);
        atomic_A_->compute_weights(npoints,&xp[offset],&yp[offset],&zp[offset],QAPp,&rhoap[offset]);
        atomic_B_->compute_weights(npoints,&xp[offset],&yp[offset],&zp[offset],QBQp,&rhobp[offset]);

        #pragma omp parallel for schedule(dynamic)
        for (int P = 0; P < npoints; P++) {

            // Thread info
            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            // Pointers
            double** ZxyzTp = ZxyzT[thread]->pointer();
            double** VtempTp = VtempT[thread]->pointer();
            double** QACTp = QACT[thread]->pointer();       
            double** QBDTp = QBDT[thread]->pointer();       
 
            // Absolute point indexing
            int Pabs = P + offset;
            
            // => Q_A^P and Q_B^Q <= //

            VtempT[thread]->zero();
            ZxyzTp[0][0] = 1.0;
            ZxyzTp[0][1] = xp[Pabs];
            ZxyzTp[0][2] = yp[Pabs];
            ZxyzTp[0][3] = zp[Pabs]; 
            VintT[thread]->compute(VtempT[thread]);

            // Potential integrals add a spurios minus sign
            C_DGER(nA,nQ,-wp[Pabs],&QAPp[0][P],max_points,VtempTp[0],1,QACTp[0],nQ);
            C_DGER(nB,nQ,-wp[Pabs],&QBQp[0][P],max_points,VtempTp[0],1,QBDTp[0],nQ);

        }

        for (int P = 0; P < npoints; P++) {
        
            // Absolute point indexing
            int Pabs = P + offset;
            
            // => V_A^B <= //

            for (int B = 0; B < nB; B++) {
                double Z  = monomer_B_->Z(cB[B]);
                double xc = monomer_B_->x(cB[B]);
                double yc = monomer_B_->y(cB[B]);
                double zc = monomer_B_->z(cB[B]);
                double R = sqrt((xp[Pabs] - xc) * (xp[Pabs] - xc) + 
                                (yp[Pabs] - yc) * (yp[Pabs] - yc) + 
                                (zp[Pabs] - zc) * (zp[Pabs] - zc));
                if (R < 1.0E-12) continue;
                double Q = Z * wp[Pabs] / R; 
                for (int A = 0; A < nA; A++) {
                    VABp[A][B] += -2.0 * Q * QAPp[A][P]; 
                }
            }
            
            // => V_B^A <= //

            for (int A = 0; A < nA; A++) {
                double Z  = monomer_A_->Z(cA[A]);
                double xc = monomer_A_->x(cA[A]);
                double yc = monomer_A_->y(cA[A]);
                double zc = monomer_A_->z(cA[A]);
                double R = sqrt((xp[Pabs] - xc) * (xp[Pabs] - xc) + 
                                (yp[Pabs] - yc) * (yp[Pabs] - yc) + 
                                (zp[Pabs] - zc) * (zp[Pabs] - zc));
                if (R < 1.0E-12) continue;
                double Q = Z * wp[Pabs] / R; 
                for (int B = 0; B < nB; B++) {
                    VBAp[B][A] += -2.0 * Q * QBQp[B][P]; 
                }
            }

        }
    }

    for (int thread = 0; thread < nthreads; thread++) {
        QAC->add(QACT[thread]);
        QBD->add(QBDT[thread]);
    }

    boost::shared_ptr<Matrix> RAC = Matrix::doublet(QAC,Jm12);
    boost::shared_ptr<Matrix> RBD = Matrix::doublet(QBD,Jm12);
    double** RACp = RAC->pointer();
    double** RBDp = RBD->pointer();

    // ==> Elst <== //

    double Elst10 = 0.0;
    std::vector<double> Elst10_terms;
    Elst10_terms.resize(4);

    boost::shared_ptr<Matrix> Elst_atoms(new Matrix("Elst (A x B)", nA, nB));
    double** Elst_atomsp = Elst_atoms->pointer();

    // => A <-> B <= //
     
    for (int A = 0; A < nA; A++) {
        for (int B = 0; B < nB; B++) {
            double val = monomer_A_->Z(cA[A]) * monomer_B_->Z(cB[B]) / (monomer_A_->xyz(cA[A]).distance(monomer_B_->xyz(cB[B])));
            Elst10_terms[3] += val;
            Elst_atomsp[A][B] += val;
        }
    }

    // => a <-> b <= //

    boost::shared_ptr<Matrix> Elst10_3 = Matrix::doublet(RAC,RBD,false,true);
    double** Elst10_3p = Elst10_3->pointer();
    for (int A = 0; A < nA; A++) {
        for (int B = 0; B < nB; B++) {
            double val = 4.0 * Elst10_3p[A][B];
            Elst10_terms[2] += val;
            Elst_atomsp[A][B] += val;
        }
    }

    // => A <-> b <= //

    for (int A = 0; A < nA; A++) {
        for (int B = 0; B < nB; B++) {
            double val = VBAp[B][A];
            Elst10_terms[1] += val;
            Elst_atomsp[A][B] += val;
        }
    }

    // => a <-> B <= //

    for (int B = 0; B < nB; B++) {
        for (int A = 0; A < nA; A++) {
            double val = VABp[A][B];
            Elst10_terms[0] += val;
            Elst_atomsp[A][B] += val;
        }
    }

    for (int k = 0; k < Elst10_terms.size(); k++) {
        Elst10 += Elst10_terms[k];
    }
    if (debug_) {
        for (int k = 0; k < Elst10_terms.size(); k++) {
            fprintf(outfile,"    Elst10,r (%1d)        = %18.12lf H\n",k+1,Elst10_terms[k]);
        }
    }
    //energies_["Elst10,r"] = Elst10;
    fprintf(outfile,"    Elst10,r            = %18.12lf H\n",Elst10);
    fprintf(outfile,"\n");
    fflush(outfile);

    vis_->vars()["Elst_AB"] = Elst_atoms;

    // ==> Setup Atomic Electrostatic Fields for Induction <== //

    boost::shared_ptr<Tensor> WBarT = DiskTensor::build("WBar", "nB", nB, "na", na, "nr", nr, false, false);
    FILE* WBarf = WBarT->file_pointer();
    boost::shared_ptr<Tensor> WAbsT = DiskTensor::build("WAbs", "nA", nA, "nb", nb, "ns", ns, false, false);
    FILE* WAbsf = WAbsT->file_pointer();

    // => Nuclear Part (PITA) <= //
    
    boost::shared_ptr<Matrix> Zxyz2(new Matrix("Zxyz",1,4));
    double** Zxyz2p = Zxyz2->pointer();
    boost::shared_ptr<IntegralFactory> Vfact2(new IntegralFactory(primary_));
    boost::shared_ptr<PotentialInt> Vint2(static_cast<PotentialInt*>(Vfact2->ao_potential()));
    Vint2->set_charge_field(Zxyz2);
    boost::shared_ptr<Matrix> Vtemp2(new Matrix("Vtemp2",nn,nn));

    for (int A = 0; A < nA; A++) {
        Vtemp2->zero();
        Zxyz2p[0][0] = monomer_A_->Z(cA[A]);
        Zxyz2p[0][1] = monomer_A_->x(cA[A]);
        Zxyz2p[0][2] = monomer_A_->y(cA[A]);
        Zxyz2p[0][3] = monomer_A_->z(cA[A]); 
        Vint2->compute(Vtemp2);
        boost::shared_ptr<Matrix> Vbs = Matrix::triplet(Cocc_B_,Vtemp2,Cvir_B_,true,false,false);
        double** Vbsp = Vbs->pointer();
        fwrite(Vbsp[0],sizeof(double),nb*ns,WAbsf);
    }

    for (int B = 0; B < nB; B++) {
        Vtemp2->zero();
        Zxyz2p[0][0] = monomer_B_->Z(cB[B]);
        Zxyz2p[0][1] = monomer_B_->x(cB[B]);
        Zxyz2p[0][2] = monomer_B_->y(cB[B]);
        Zxyz2p[0][3] = monomer_B_->z(cB[B]); 
        Vint2->compute(Vtemp2);
        boost::shared_ptr<Matrix> Var = Matrix::triplet(Cocc_A_,Vtemp2,Cvir_A_,true,false,false);
        double** Varp = Var->pointer();
        fwrite(Varp[0],sizeof(double),na*nr,WBarf);
    }

    // => Electronic Part (Massive PITA) <= //

    FILE* Absf = AbsT->file_pointer();
    fseek(Absf,0L,SEEK_SET);
    boost::shared_ptr<Matrix> TsQ(new Matrix("TsQ",ns,nQ));
    boost::shared_ptr<Matrix> T1As(new Matrix("T1As",nA,ns));
    boost::shared_ptr<Matrix> T2As(new Matrix("T2As",1,ns));
    double** TsQp = TsQ->pointer(); 
    double** T1Asp = T1As->pointer(); 
    double** T2Asp = T2As->pointer(); 
    for (int b = 0; b < nb; b++) {
        fread(TsQp[0],sizeof(double),ns*nQ,Absf);
        C_DGEMM('N','T',nA,ns,nQ,2.0,RACp[0],nQ,TsQp[0],nQ,0.0,T1Asp[0],ns);
        for (int A = 0; A < nA; A++) {
            fseek(WAbsf,A*nb*ns*sizeof(double) + b*ns*sizeof(double),SEEK_SET);
            fread(T2Asp[0],sizeof(double),ns,WAbsf);
            C_DAXPY(ns,1.0,T1Asp[A],1,T2Asp[0],1);
            fseek(WAbsf,A*nb*ns*sizeof(double) + b*ns*sizeof(double),SEEK_SET);
            fwrite(T2Asp[0],sizeof(double),ns,WAbsf);
        }
    }

    FILE* Aarf = AarT->file_pointer();
    fseek(Aarf,0L,SEEK_SET);
    boost::shared_ptr<Matrix> TrQ(new Matrix("TrQ",nr,nQ));
    boost::shared_ptr<Matrix> T1Br(new Matrix("T1Br",nB,nr));
    boost::shared_ptr<Matrix> T2Br(new Matrix("T2Br",1,nr));
    double** TrQp = TrQ->pointer(); 
    double** T1Brp = T1Br->pointer(); 
    double** T2Brp = T2Br->pointer(); 
    for (int a = 0; a < na; a++) {
        fread(TrQp[0],sizeof(double),nr*nQ,Aarf);
        C_DGEMM('N','T',nB,nr,nQ,2.0,RBDp[0],nQ,TrQp[0],nQ,0.0,T1Brp[0],nr);
        for (int B = 0; B < nB; B++) {
            fseek(WBarf,B*na*nr*sizeof(double) + a*nr*sizeof(double),SEEK_SET);
            fread(T2Brp[0],sizeof(double),nr,WBarf);
            C_DAXPY(nr,1.0,T1Brp[B],1,T2Brp[0],1);
            fseek(WBarf,B*na*nr*sizeof(double) + a*nr*sizeof(double),SEEK_SET);
            fwrite(T2Brp[0],sizeof(double),nr,WBarf);
        }
    }
     
    tensors_["WAbs"] = WAbsT;
    tensors_["WBar"] = WBarT;
}
void ASAPT::exch()
{
    fprintf(outfile, "  EXCHANGE:\n\n");

    // ==> Sizing <== //

    int nn = primary_->nbf();
    int na = Cocc_A_->colspi()[0];
    int nb = Cocc_B_->colspi()[0];

    int nA = 0;
    std::vector<int> cA;
    for (int A = 0; A < monomer_A_->natom(); A++) {
        if (monomer_A_->Z(A) != 0.0) {
            nA++;
            cA.push_back(A);
        }
    }

    int nB = 0;
    std::vector<int> cB;
    for (int B = 0; B < monomer_B_->natom(); B++) {
        if (monomer_B_->Z(B) != 0.0) {
            nB++;
            cB.push_back(B);
        }
    }

    int nr = Cvir_A_->colspi()[0];
    int ns = Cvir_B_->colspi()[0];

    // ==> Stack Variables <== //

    boost::shared_ptr<Matrix> S   = vars_["S"];
    boost::shared_ptr<Matrix> V_A = vars_["V_A"];
    boost::shared_ptr<Matrix> J_A = vars_["J_A"];
    boost::shared_ptr<Matrix> V_B = vars_["V_B"];
    boost::shared_ptr<Matrix> J_B = vars_["J_B"];
    boost::shared_ptr<Matrix> L_A = Locc_A_;
    boost::shared_ptr<Matrix> L_B = Locc_B_;

    // ==> DF ERI Setup (JKFIT Type, in Full Basis) <== //

    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    boost::shared_ptr<BasisSet> jkfit = BasisSet::construct(parser, primary_->molecule(), "DF_BASIS_SCF");
    int nQ = jkfit->nbf();

    boost::shared_ptr<DFERI> df = DFERI::build(primary_,jkfit,Process::environment.options);
    df->clear();

    std::vector<boost::shared_ptr<Matrix> > Cs;
    Cs.push_back(L_A);
    Cs.push_back(Cvir_A_);
    Cs.push_back(L_B);
    Cs.push_back(Cvir_B_);
    boost::shared_ptr<Matrix> Call = Matrix::horzcat(Cs);
    Cs.clear();

    df->set_C(Call);
    df->set_memory(memory_);

    int offset = 0;
    df->add_space("a",offset,offset+Cocc_A_->colspi()[0]); offset += Cocc_A_->colspi()[0];
    df->add_space("r",offset,offset+Cvir_A_->colspi()[0]); offset += Cvir_A_->colspi()[0];
    df->add_space("b",offset,offset+Cocc_B_->colspi()[0]); offset += Cocc_B_->colspi()[0];
    df->add_space("s",offset,offset+Cvir_B_->colspi()[0]); offset += Cvir_B_->colspi()[0];

    df->add_pair_space("Aaa", "a", "a");
    df->add_pair_space("Abb", "b", "b");
    df->add_pair_space("Aar", "a", "r");
    df->add_pair_space("Abs", "b", "s");

    df->print_header();
    df->compute();

    std::map<std::string, boost::shared_ptr<Tensor> >& ints = df->ints();
    boost::shared_ptr<Tensor> AaaT = ints["Aaa"];
    boost::shared_ptr<Tensor> AbbT = ints["Abb"];
    boost::shared_ptr<Tensor> AarT = ints["Aar"];
    boost::shared_ptr<Tensor> AbsT = ints["Abs"];

    df.reset();

    boost::shared_ptr<Matrix> AaQ(new Matrix("AaQ",na,nQ));
    boost::shared_ptr<Matrix> AbQ(new Matrix("AbQ",nb,nQ));
    double** AaQp = AaQ->pointer();
    double** AbQp = AbQ->pointer();
    FILE* Aaaf = AaaT->file_pointer();
    FILE* Abbf = AbbT->file_pointer();

    for (int a = 0; a < na; a++) {
        fseek(Aaaf,(a*na+a)*(size_t)nQ*sizeof(double),SEEK_SET);
        fread(AaQp[a],sizeof(double),nQ,Aaaf);
    }

    for (int b = 0; b < nb; b++) {
        fseek(Abbf,(b*nb+b)*(size_t)nQ*sizeof(double),SEEK_SET);
        fread(AbQp[b],sizeof(double),nQ,Abbf);
    }

    // ==> Electrostatic Potentials <== //

    boost::shared_ptr<Matrix> W_A(J_A->clone());
    W_A->copy(J_A);
    W_A->scale(2.0);
    W_A->add(V_A);

    boost::shared_ptr<Matrix> W_B(J_B->clone());
    W_B->copy(J_B);
    W_B->scale(2.0);
    W_B->add(V_B);

    boost::shared_ptr<Matrix> WAbs = Matrix::triplet(Locc_B_,W_A,Cvir_B_,true,false,false);
    boost::shared_ptr<Matrix> WBar = Matrix::triplet(Locc_A_,W_B,Cvir_A_,true,false,false);
    double** WBarp = WBar->pointer();
    double** WAbsp = WAbs->pointer();
    
    W_A.reset();
    W_B.reset();

    // ==> Exchange S^2 Computation <== //

    double Exch10_2 = 0.0;
    std::vector<double> Exch10_2_terms;
    Exch10_2_terms.resize(3);

    boost::shared_ptr<Matrix> Sab = Matrix::triplet(L_A,S,L_B,true,false,false);
    boost::shared_ptr<Matrix> Sba = Matrix::triplet(L_B,S,L_A,true,false,false);
    boost::shared_ptr<Matrix> Sas = Matrix::triplet(L_A,S,Cvir_B_,true,false,false);
    boost::shared_ptr<Matrix> Sbr = Matrix::triplet(L_B,S,Cvir_A_,true,false,false);
    double** Sabp = Sab->pointer();
    double** Sbap = Sba->pointer();
    double** Sasp = Sas->pointer();
    double** Sbrp = Sbr->pointer();

    boost::shared_ptr<Matrix> WBab(new Matrix("WBab",na,nb));
    double** WBabp = WBab->pointer();
    boost::shared_ptr<Matrix> WAba(new Matrix("WAba",nb,na));
    double** WAbap = WAba->pointer();

    C_DGEMM('N','T',na,nb,nr,1.0,WBarp[0],nr,Sbrp[0],nr,0.0,WBabp[0],nb);
    C_DGEMM('N','T',nb,na,ns,1.0,WAbsp[0],ns,Sasp[0],ns,0.0,WAbap[0],na);

    boost::shared_ptr<Matrix> E_exch1(new Matrix("E_exch [a <x- b]", na, nb));
    double** E_exch1p = E_exch1->pointer();
    boost::shared_ptr<Matrix> E_exch2(new Matrix("E_exch [a -x> b]", na, nb));
    double** E_exch2p = E_exch2->pointer();

    for (int a = 0; a < na; a++) {
        for (int b = 0; b < nb; b++) {
            E_exch1p[a][b] -= 2.0 * Sabp[a][b] * WBabp[a][b];
            E_exch2p[a][b] -= 2.0 * Sbap[b][a] * WAbap[b][a];
        }
    }

    //E_exch1->print();
    //E_exch2->print();

    boost::shared_ptr<Matrix> TrQ(new Matrix("TrQ",nr,nQ));
    double** TrQp = TrQ->pointer();
    boost::shared_ptr<Matrix> TsQ(new Matrix("TsQ",ns,nQ));
    double** TsQp = TsQ->pointer();
    boost::shared_ptr<Matrix> TbQ(new Matrix("TbQ",nb,nQ));
    double** TbQp = TbQ->pointer();
    boost::shared_ptr<Matrix> TaQ(new Matrix("TaQ",na,nQ));
    double** TaQp = TaQ->pointer();

    boost::shared_ptr<Tensor> BabT = DiskTensor::build("BabT","na",na,"nb",nb,"nQ",nQ,false,false);
    FILE* Aarf = AarT->file_pointer();
    FILE* Babf = BabT->file_pointer();
    fseek(Babf,0L,SEEK_SET);
    fseek(Aarf,0L,SEEK_SET);
    for (int a = 0; a < na; a++) {
        fread(TrQp[0],sizeof(double),nr*nQ,Aarf);
        C_DGEMM('N','N',nb,nQ,nr,1.0,Sbrp[0],nr,TrQp[0],nQ,0.0,TbQp[0],nQ);
        fwrite(TbQp[0],sizeof(double),nb*nQ,Babf);
    }

    boost::shared_ptr<Tensor> BbaT = DiskTensor::build("BbaT","nb",nb,"na",na,"nQ",nQ,false,false);
    FILE* Absf = AbsT->file_pointer();
    FILE* Bbaf = BbaT->file_pointer();
    fseek(Bbaf,0L,SEEK_SET);
    fseek(Absf,0L,SEEK_SET);
    for (int b = 0; b < nb; b++) {
        fread(TsQp[0],sizeof(double),ns*nQ,Absf);
        C_DGEMM('N','N',na,nQ,ns,1.0,Sasp[0],ns,TsQp[0],nQ,0.0,TaQp[0],nQ);
        fwrite(TaQp[0],sizeof(double),na*nQ,Bbaf);
    }

    boost::shared_ptr<Matrix> E_exch3(new Matrix("E_exch [a <x-x> b]", na, nb));
    double** E_exch3p = E_exch3->pointer();

    fseek(Babf,0L,SEEK_SET);
    for (int a = 0; a < na; a++) {
        fread(TbQp[0],sizeof(double),nb*nQ,Babf);
        for (int b = 0; b < nb; b++) {
            fseek(Bbaf,(b*na+a)*(size_t)nQ*sizeof(double),SEEK_SET);
            fread(TaQp[0],sizeof(double),nQ,Bbaf);
            E_exch3p[a][b] -= 2.0 * C_DDOT(nQ,TbQp[b],1,TaQp[0],1);
        }
    }

    //E_exch3->print();

    boost::shared_ptr<Matrix> E_exch(new Matrix("E_exch (a x b)", na, nb));
    double** E_exchp = E_exch->pointer();

    for (int a = 0; a < na; a++) {
        for (int b = 0; b < nb; b++) {
            E_exchp[a][b] = E_exch1p[a][b] +
                            E_exch2p[a][b] +
                            E_exch3p[a][b];
            Exch10_2_terms[0] += E_exch1p[a][b];
            Exch10_2_terms[1] += E_exch2p[a][b];
            Exch10_2_terms[2] += E_exch3p[a][b];
        }
    }

    for (int k = 0; k < Exch10_2_terms.size(); k++) {
        Exch10_2 += Exch10_2_terms[k];
    }
    if (debug_) {
        for (int k = 0; k < Exch10_2_terms.size(); k++) {
            fprintf(outfile,"    Exch10(S^2) (%1d)     = %18.12lf H\n",k+1,Exch10_2_terms[k]);
        }
    }
    //energies_["Exch10(S^2)"] = Exch10_2;
    fprintf(outfile,"    Exch10(S^2)         = %18.12lf H\n",Exch10_2);
    fprintf(outfile, "\n");
    fflush(outfile);

    // => Exchange scaling <= //

    if (exch_scale_) {
        double scale = energies_["Exch10"] / energies_["Exch10(S^2)"];
        E_exch->scale(scale);
        fprintf(outfile,"    Scaling ASAPT Exchange by %11.3E to match S^\\infty\n\n", scale);
    }
    
    vis_->vars()["Exch_ab"] = E_exch;
}
void ASAPT::ind()
{
    fprintf(outfile, "  INDUCTION:\n\n");

    // => Sizing <= //

    int nn = primary_->nbf();

    int na = Cocc_A_->colspi()[0];
    int nb = Cocc_B_->colspi()[0];
    int nA = 0;
    std::vector<int> cA;
    for (int A = 0; A < monomer_A_->natom(); A++) {
        if (monomer_A_->Z(A) != 0.0) {
            nA++;
            cA.push_back(A);
        }
    }

    int nB = 0;
    std::vector<int> cB;
    for (int B = 0; B < monomer_B_->natom(); B++) {
        if (monomer_B_->Z(B) != 0.0) {
            nB++;
            cB.push_back(B);
        }
    }

    int nr = Cvir_A_->colspi()[0];
    int ns = Cvir_B_->colspi()[0];
    int nQ = mp2fit_->nbf();
    size_t naQ = na * (size_t) nQ;
    size_t nbQ = nb * (size_t) nQ;

    int nT = 1;
    #ifdef _OPENMP
        nT = omp_get_max_threads();
    #endif

    // ==> Stack Variables <== //

    double*  eap = eps_occ_A_->pointer();
    double*  ebp = eps_occ_B_->pointer();
    double*  erp = eps_vir_A_->pointer();
    double*  esp = eps_vir_B_->pointer();

    FILE* WAbsf = tensors_["WAbs"]->file_pointer();
    FILE* WBarf = tensors_["WBar"]->file_pointer();
    
    boost::shared_ptr<Matrix> S   = vars_["S"];
    boost::shared_ptr<Matrix> D_A = vars_["D_A"];
    boost::shared_ptr<Matrix> V_A = vars_["V_A"];
    boost::shared_ptr<Matrix> J_A = vars_["J_A"];
    boost::shared_ptr<Matrix> K_A = vars_["K_A"];
    boost::shared_ptr<Matrix> D_B = vars_["D_B"];
    boost::shared_ptr<Matrix> V_B = vars_["V_B"];
    boost::shared_ptr<Matrix> J_B = vars_["J_B"];
    boost::shared_ptr<Matrix> K_B = vars_["K_B"];
    boost::shared_ptr<Matrix> J_O = vars_["J_O"];
    boost::shared_ptr<Matrix> K_O = vars_["K_O"];
    boost::shared_ptr<Matrix> J_P_A = vars_["J_P_A"];
    boost::shared_ptr<Matrix> J_P_B = vars_["J_P_B"];

    // ==> MO Amplitudes/Sources (by source atom) <== //

    boost::shared_ptr<Matrix> xA(new Matrix("xA",na,nr));
    boost::shared_ptr<Matrix> xB(new Matrix("xB",nb,ns));
    double** xAp = xA->pointer(); 
    double** xBp = xB->pointer(); 

    boost::shared_ptr<Matrix> wB(new Matrix("wB",na,nr));
    boost::shared_ptr<Matrix> wA(new Matrix("wA",nb,ns));
    double** wBp = wB->pointer(); 
    double** wAp = wA->pointer(); 
    
    // ==> Generalized ESP (Flat and Exchange) <== //

    std::map<std::string, boost::shared_ptr<Matrix> > mapA;
    mapA["Cocc_A"] = Locc_A_;
    mapA["Cvir_A"] = Cvir_A_;
    mapA["Cocc_B"] = Locc_B_;
    mapA["Cvir_B"] = Cvir_B_;
    mapA["S"] = S;
    mapA["D_A"] = D_A;
    mapA["V_A"] = V_A;
    mapA["J_A"] = J_A;
    mapA["K_A"] = K_A;
    mapA["D_B"] = D_B;
    mapA["V_B"] = V_B;
    mapA["J_B"] = J_B;
    mapA["K_B"] = K_B;
    mapA["J_O"] = J_O;
    mapA["K_O"] = K_O;
    mapA["J_P"] = J_P_A; 

    boost::shared_ptr<Matrix> wBT = build_ind_pot(mapA);
    boost::shared_ptr<Matrix> uBT = build_exch_ind_pot(mapA);
    double** wBTp = wBT->pointer();
    double** uBTp = uBT->pointer();

    K_O->transpose_this();

    std::map<std::string, boost::shared_ptr<Matrix> > mapB;
    mapB["Cocc_A"] = Locc_B_;
    mapB["Cvir_A"] = Cvir_B_;
    mapB["Cocc_B"] = Locc_A_;
    mapB["Cvir_B"] = Cvir_A_;
    mapB["S"] = S;
    mapB["D_A"] = D_B;
    mapB["V_A"] = V_B;
    mapB["J_A"] = J_B;
    mapB["K_A"] = K_B;
    mapB["D_B"] = D_A;
    mapB["V_B"] = V_A;
    mapB["J_B"] = J_A;
    mapB["K_B"] = K_A;
    mapB["J_O"] = J_O;
    mapB["K_O"] = K_O;
    mapB["J_P"] = J_P_B; 

    boost::shared_ptr<Matrix> wAT = build_ind_pot(mapB);
    boost::shared_ptr<Matrix> uAT = build_exch_ind_pot(mapB);
    double** wATp = wAT->pointer();
    double** uATp = uAT->pointer();
    
    K_O->transpose_this();

    // ==> Uncoupled Targets <== //

    boost::shared_ptr<Matrix> Ind20u_AB_terms(new Matrix("Ind20 [A<-B] (a x B)", na, nB));
    boost::shared_ptr<Matrix> Ind20u_BA_terms(new Matrix("Ind20 [B<-A] (A x b)", nA, nb));
    double** Ind20u_AB_termsp = Ind20u_AB_terms->pointer();
    double** Ind20u_BA_termsp = Ind20u_BA_terms->pointer();

    double Ind20u_AB = 0.0; 
    double Ind20u_BA = 0.0;

    boost::shared_ptr<Matrix> ExchInd20u_AB_terms(new Matrix("ExchInd20 [A<-B] (a x B)", na, nB));
    boost::shared_ptr<Matrix> ExchInd20u_BA_terms(new Matrix("ExchInd20 [B<-A] (A x b)", nA, nb));
    double** ExchInd20u_AB_termsp = ExchInd20u_AB_terms->pointer();
    double** ExchInd20u_BA_termsp = ExchInd20u_BA_terms->pointer();

    double ExchInd20u_AB = 0.0; 
    double ExchInd20u_BA = 0.0;

    boost::shared_ptr<Matrix> Indu_AB_terms(new Matrix("Ind [A<-B] (a x B)", na, nB));
    boost::shared_ptr<Matrix> Indu_BA_terms(new Matrix("Ind [B<-A] (A x b)", nA, nb));
    double** Indu_AB_termsp = Indu_AB_terms->pointer();
    double** Indu_BA_termsp = Indu_BA_terms->pointer();

    double Indu_AB = 0.0; 
    double Indu_BA = 0.0;
    
    // ==> A <- B Uncoupled <== //

    fseek(WBarf,0L,SEEK_SET);
    for (int B = 0; B < nB; B++) {

        // ESP
        fread(wBp[0],sizeof(double),na*nr,WBarf); 
        
        // Uncoupled amplitude
        for (int a = 0; a < na; a++) {
            for (int r = 0; r < nr; r++) {
                xAp[a][r] = wBp[a][r] / (eap[a] - erp[r]);
            }
        }

        // Backtransform the amplitude to LO
        boost::shared_ptr<Matrix> x2A = Matrix::doublet(Uocc_A_,xA,true,false);
        double** x2Ap = x2A->pointer();

        // Zip up the Ind20 contributions
        for (int a = 0; a < na; a++) {
            double Jval = 2.0 * C_DDOT(nr,x2Ap[a],1,wBTp[a],1);
            double Kval = 2.0 * C_DDOT(nr,x2Ap[a],1,uBTp[a],1);
            Ind20u_AB_termsp[a][B] = Jval;
            Ind20u_AB += Jval;
            ExchInd20u_AB_termsp[a][B] = Kval;
            ExchInd20u_AB += Kval;
            Indu_AB_termsp[a][B] = Jval + Kval;
            Indu_AB += Jval + Kval;
        }

    } 

    // ==> B <- A Uncoupled <== //

    fseek(WAbsf,0L,SEEK_SET);
    for (int A = 0; A < nA; A++) {

        // ESP
        fread(wAp[0],sizeof(double),nb*ns,WAbsf); 
        
        // Uncoupled amplitude
        for (int b = 0; b < nb; b++) {
            for (int s = 0; s < ns; s++) {
                xBp[b][s] = wAp[b][s] / (ebp[b] - esp[s]);
            }
        }

        // Backtransform the amplitude to LO
        boost::shared_ptr<Matrix> x2B = Matrix::doublet(Uocc_B_,xB,true,false);
        double** x2Bp = x2B->pointer();

        // Zip up the Ind20 contributions
        for (int b = 0; b < nb; b++) {
            double Jval = 2.0 * C_DDOT(ns,x2Bp[b],1,wATp[b],1);
            double Kval = 2.0 * C_DDOT(ns,x2Bp[b],1,uATp[b],1);
            Ind20u_BA_termsp[A][b] = Jval;
            Ind20u_BA += Jval;
            ExchInd20u_BA_termsp[A][b] = Kval;
            ExchInd20u_BA += Kval;
            Indu_BA_termsp[A][b] = Jval + Kval;
            Indu_BA += Jval + Kval;
        }
    
    } 

    double Ind20u = Ind20u_AB + Ind20u_BA;
    fprintf(outfile,"    Ind20,u (A<-B)      = %18.12lf H\n",Ind20u_AB);
    fprintf(outfile,"    Ind20,u (A->B)      = %18.12lf H\n",Ind20u_BA);
    fprintf(outfile,"    Ind20,u             = %18.12lf H\n",Ind20u);
    fflush(outfile);

    double ExchInd20u = ExchInd20u_AB + ExchInd20u_BA;
    fprintf(outfile,"    Exch-Ind20,u (A<-B) = %18.12lf H\n",ExchInd20u_AB);
    fprintf(outfile,"    Exch-Ind20,u (B<-A) = %18.12lf H\n",ExchInd20u_BA);
    fprintf(outfile,"    Exch-Ind20,u        = %18.12lf H\n",ExchInd20u);
    fprintf(outfile,"\n");
    fflush(outfile);

    double Ind = Ind20u + ExchInd20u;
    boost::shared_ptr<Matrix> Ind_AB_terms = Indu_AB_terms;
    boost::shared_ptr<Matrix> Ind_BA_terms = Indu_BA_terms;

    if (ind_resp_) {

        fprintf(outfile, "  COUPLED INDUCTION (You asked for it!):\n\n");

        // ==> Coupled Targets <== //

        boost::shared_ptr<Matrix> Ind20r_AB_terms(new Matrix("Ind20 [A<-B] (a x B)", na, nB));
        boost::shared_ptr<Matrix> Ind20r_BA_terms(new Matrix("Ind20 [B<-A] (A x b)", nA, nb));
        double** Ind20r_AB_termsp = Ind20r_AB_terms->pointer();
        double** Ind20r_BA_termsp = Ind20r_BA_terms->pointer();

        double Ind20r_AB = 0.0; 
        double Ind20r_BA = 0.0;

        boost::shared_ptr<Matrix> ExchInd20r_AB_terms(new Matrix("ExchInd20 [A<-B] (a x B)", na, nB));
        boost::shared_ptr<Matrix> ExchInd20r_BA_terms(new Matrix("ExchInd20 [B<-A] (A x b)", nA, nb));
        double** ExchInd20r_AB_termsp = ExchInd20r_AB_terms->pointer();
        double** ExchInd20r_BA_termsp = ExchInd20r_BA_terms->pointer();

        double ExchInd20r_AB = 0.0; 
        double ExchInd20r_BA = 0.0;

        boost::shared_ptr<Matrix> Indr_AB_terms(new Matrix("Ind [A<-B] (a x B)", na, nB));
        boost::shared_ptr<Matrix> Indr_BA_terms(new Matrix("Ind [B<-A] (A x b)", nA, nb));
        double** Indr_AB_termsp = Indr_AB_terms->pointer();
        double** Indr_BA_termsp = Indr_BA_terms->pointer();

        double Indr_AB = 0.0; 
        double Indr_BA = 0.0;

        // => JK Object <= //

        boost::shared_ptr<JK> jk = JK::build_JK();

        // TODO: Account for 2-index overhead in memory
        int nso = Cocc_A_->nrow();
        long int jk_memory = (long int)memory_;
        jk_memory -= 24 * nso * nso;
        jk_memory -=  4 * na * nso;
        jk_memory -=  4 * nb * nso;
        if (jk_memory < 0L) {
            throw PSIEXCEPTION("Too little static memory for ASAPT::induction");
        }
        jk->set_memory((unsigned long int )jk_memory);
        jk->set_do_J(true);
        jk->set_do_K(true);
        jk->initialize();
        jk->print_header();

        // ==> Master Loop over perturbing atoms <== //
    
        int nC = std::max(nA,nB);

        fseek(WBarf,0L,SEEK_SET);
        fseek(WAbsf,0L,SEEK_SET);
        
        for (int C = 0; C < nC; C++) {
            
            if (C < nB) fread(wBp[0],sizeof(double),na*nr,WBarf); 
            if (C < nA) fread(wAp[0],sizeof(double),nb*ns,WAbsf); 

            fprintf(outfile,"    Responses for (A <- Atom B = %3d) and (B <- Atom A = %3d)\n\n",
                    (C < nB ? C : nB - 1), (C < nA ? C : nA - 1));

            std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Matrix> > x_sol = compute_x(jk,wB,wA);
            xA = x_sol.first;
            xB = x_sol.second;
            xA->scale(-1.0);
            xB->scale(-1.0);

            if (C < nB) {
                // Backtransform the amplitude to LO
                boost::shared_ptr<Matrix> x2A = Matrix::doublet(Uocc_A_,xA,true,false);
                double** x2Ap = x2A->pointer();

                // Zip up the Ind20 contributions
                for (int a = 0; a < na; a++) {
                    double Jval = 2.0 * C_DDOT(nr,x2Ap[a],1,wBTp[a],1);
                    double Kval = 2.0 * C_DDOT(nr,x2Ap[a],1,uBTp[a],1);
                    Ind20r_AB_termsp[a][C] = Jval;
                    Ind20r_AB += Jval;
                    ExchInd20r_AB_termsp[a][C] = Kval;
                    ExchInd20r_AB += Kval;
                    Indr_AB_termsp[a][C] = Jval + Kval;
                    Indr_AB += Jval + Kval;
                }
            }

            if (C < nA) { 
                // Backtransform the amplitude to LO
                boost::shared_ptr<Matrix> x2B = Matrix::doublet(Uocc_B_,xB,true,false);
                double** x2Bp = x2B->pointer();

                // Zip up the Ind20 contributions
                for (int b = 0; b < nb; b++) {
                    double Jval = 2.0 * C_DDOT(ns,x2Bp[b],1,wATp[b],1);
                    double Kval = 2.0 * C_DDOT(ns,x2Bp[b],1,uATp[b],1);
                    Ind20r_BA_termsp[C][b] = Jval;
                    Ind20r_BA += Jval;
                    ExchInd20r_BA_termsp[C][b] = Kval;
                    ExchInd20r_BA += Kval;
                    Indr_BA_termsp[C][b] = Jval + Kval;
                    Indr_BA += Jval + Kval;
                }
            }
        }

        double Ind20r = Ind20r_AB + Ind20r_BA;
        fprintf(outfile,"    Ind20,r (A<-B)      = %18.12lf H\n",Ind20r_AB);
        fprintf(outfile,"    Ind20,r (A->B)      = %18.12lf H\n",Ind20r_BA);
        fprintf(outfile,"    Ind20,r             = %18.12lf H\n",Ind20r);
        fflush(outfile);

        double ExchInd20r = ExchInd20r_AB + ExchInd20r_BA;
        fprintf(outfile,"    Exch-Ind20,r (A<-B) = %18.12lf H\n",ExchInd20r_AB);
        fprintf(outfile,"    Exch-Ind20,r (B<-A) = %18.12lf H\n",ExchInd20r_BA);
        fprintf(outfile,"    Exch-Ind20,r        = %18.12lf H\n",ExchInd20r);
        fprintf(outfile,"\n");
        fflush(outfile);

        Ind = Ind20r + ExchInd20r;
        Ind_AB_terms = Indr_AB_terms;
        Ind_BA_terms = Indr_BA_terms;
    }

    // => Induction scaling <= //

    if (ind_scale_) {
        double dHF = 0.0;
        if (energies_["HF"] != 0.0) {
            dHF = energies_["HF"] - energies_["Elst10,r"] - energies_["Exch10"] - energies_["Ind20,r"] - energies_["Exch-Ind20,r"];
        }
        double IndSAPT0 = energies_["Ind20,r"] + energies_["Exch-Ind20,r"] + dHF;
        double scale = IndSAPT0 / Ind;
        Ind_AB_terms->scale(scale);
        Ind_BA_terms->scale(scale);
        fprintf(outfile,"    Scaling ASAPT Induction by %11.3E to match SAPT0\n\n", scale);
    }
    
    vis_->vars()["IndAB_aB"] = Ind_AB_terms;
    vis_->vars()["IndBA_Ab"] = Ind_BA_terms;
}
void ASAPT::disp()
{
    fprintf(outfile, "  DISPERSION:\n\n");

    // => Sizing <= //

    int nn = primary_->nbf();

    int naa = Caocc_A_->colspi()[0];
    int nab = Caocc_B_->colspi()[0];

    int na  = Locc_A_->colspi()[0]; 
    int nb  = Locc_B_->colspi()[0]; 

    int nfa = na - naa;
    int nfb = nb - nab;

    int nA = 0;
    std::vector<int> cA;
    for (int A = 0; A < monomer_A_->natom(); A++) {
        if (monomer_A_->Z(A) != 0.0) {
            nA++;
            cA.push_back(A);
        }
    }

    int nB = 0;
    std::vector<int> cB;
    for (int B = 0; B < monomer_B_->natom(); B++) {
        if (monomer_B_->Z(B) != 0.0) {
            nB++;
            cB.push_back(B);
        }
    }

    int nr = Cavir_A_->colspi()[0];
    int ns = Cavir_B_->colspi()[0];
    int nQ = mp2fit_->nbf();
    size_t naQ = naa * (size_t) nQ;
    size_t nbQ = nab * (size_t) nQ;

    int nT = 1;
    #ifdef _OPENMP
        nT = omp_get_max_threads();
    #endif

    // => Stashed Variables <= //

    boost::shared_ptr<Matrix> S   = vars_["S"];
    boost::shared_ptr<Matrix> U_A = vars_["U_A"];
    boost::shared_ptr<Matrix> L_A = vars_["L_A"];
    boost::shared_ptr<Matrix> D_A = vars_["D_A"];
    boost::shared_ptr<Matrix> P_A = vars_["P_A"];
    boost::shared_ptr<Matrix> V_A = vars_["V_A"];
    boost::shared_ptr<Matrix> J_A = vars_["J_A"];
    boost::shared_ptr<Matrix> K_A = vars_["K_A"];
    boost::shared_ptr<Matrix> U_B = vars_["U_B"];
    boost::shared_ptr<Matrix> L_B = vars_["L_B"];
    boost::shared_ptr<Matrix> D_B = vars_["D_B"];
    boost::shared_ptr<Matrix> P_B = vars_["P_B"];
    boost::shared_ptr<Matrix> V_B = vars_["V_B"];
    boost::shared_ptr<Matrix> J_B = vars_["J_B"];
    boost::shared_ptr<Matrix> K_B = vars_["K_B"];
    boost::shared_ptr<Matrix> K_O = vars_["K_O"];

    boost::shared_ptr<Matrix> Q2A = Q_A_;
    boost::shared_ptr<Matrix> Q2B = Q_B_;
    double** Q2Ap = Q2A->pointer();
    double** Q2Bp = Q2B->pointer();

    // => Auxiliary C matrices <= //

    boost::shared_ptr<Matrix> Cr1 = Matrix::triplet(D_B,S,Cavir_A_);
    Cr1->scale(-1.0);
    Cr1->add(Cavir_A_);
    boost::shared_ptr<Matrix> Cs1 = Matrix::triplet(D_A,S,Cavir_B_);
    Cs1->scale(-1.0);
    Cs1->add(Cavir_B_);
    boost::shared_ptr<Matrix> Ca2 = Matrix::triplet(D_B,S,Caocc_A_);
    boost::shared_ptr<Matrix> Cb2 = Matrix::triplet(D_A,S,Caocc_B_);
    boost::shared_ptr<Matrix> Cr3 = Matrix::triplet(D_B,S,Cavir_A_);
    boost::shared_ptr<Matrix> CrX = Matrix::triplet(Matrix::triplet(D_A,S,D_B),S,Cavir_A_);
    Cr3->subtract(CrX);
    Cr3->scale(2.0);
    boost::shared_ptr<Matrix> Cs3 = Matrix::triplet(D_A,S,Cavir_B_);
    boost::shared_ptr<Matrix> CsX = Matrix::triplet(Matrix::triplet(D_B,S,D_A),S,Cavir_B_);
    Cs3->subtract(CsX);
    Cs3->scale(2.0);
    boost::shared_ptr<Matrix> Ca4 = Matrix::triplet(Matrix::triplet(D_A,S,D_B),S,Caocc_A_);
    Ca4->scale(-2.0);
    boost::shared_ptr<Matrix> Cb4 = Matrix::triplet(Matrix::triplet(D_B,S,D_A),S,Caocc_B_);
    Cb4->scale(-2.0);

    // => Auxiliary V matrices <= //

    boost::shared_ptr<Matrix> Jbr = Matrix::triplet(Caocc_B_,J_A,Cavir_A_,true,false,false);
    Jbr->scale(2.0);
    boost::shared_ptr<Matrix> Kbr = Matrix::triplet(Caocc_B_,K_A,Cavir_A_,true,false,false);
    Kbr->scale(-1.0);

    boost::shared_ptr<Matrix> Jas = Matrix::triplet(Caocc_A_,J_B,Cavir_B_,true,false,false);
    Jas->scale(2.0);
    boost::shared_ptr<Matrix> Kas = Matrix::triplet(Caocc_A_,K_B,Cavir_B_,true,false,false);
    Kas->scale(-1.0);

    boost::shared_ptr<Matrix> KOas = Matrix::triplet(Caocc_A_,K_O,Cavir_B_,true,false,false);
    KOas->scale(1.0);
    boost::shared_ptr<Matrix> KObr = Matrix::triplet(Caocc_B_,K_O,Cavir_A_,true,true,false);
    KObr->scale(1.0);

    boost::shared_ptr<Matrix> JBas = Matrix::triplet(Matrix::triplet(Caocc_A_,S,D_B,true,false,false),J_A,Cavir_B_);
    JBas->scale(-2.0);
    boost::shared_ptr<Matrix> JAbr = Matrix::triplet(Matrix::triplet(Caocc_B_,S,D_A,true,false,false),J_B,Cavir_A_);
    JAbr->scale(-2.0);

    boost::shared_ptr<Matrix> Jbs = Matrix::triplet(Caocc_B_,J_A,Cavir_B_,true,false,false);
    Jbs->scale(4.0);
    boost::shared_ptr<Matrix> Jar = Matrix::triplet(Caocc_A_,J_B,Cavir_A_,true,false,false);
    Jar->scale(4.0);

    boost::shared_ptr<Matrix> JAas = Matrix::triplet(Matrix::triplet(Caocc_A_,J_B,D_A,true,false,false),S,Cavir_B_);
    JAas->scale(-2.0);
    boost::shared_ptr<Matrix> JBbr = Matrix::triplet(Matrix::triplet(Caocc_B_,J_A,D_B,true,false,false),S,Cavir_A_);
    JBbr->scale(-2.0);

    // Get your signs right Hesselmann!
    boost::shared_ptr<Matrix> Vbs = Matrix::triplet(Caocc_B_,V_A,Cavir_B_,true,false,false);
    Vbs->scale(2.0);
    boost::shared_ptr<Matrix> Var = Matrix::triplet(Caocc_A_,V_B,Cavir_A_,true,false,false);
    Var->scale(2.0);
    boost::shared_ptr<Matrix> VBas = Matrix::triplet(Matrix::triplet(Caocc_A_,S,D_B,true,false,false),V_A,Cavir_B_);
    VBas->scale(-1.0);
    boost::shared_ptr<Matrix> VAbr = Matrix::triplet(Matrix::triplet(Caocc_B_,S,D_A,true,false,false),V_B,Cavir_A_);
    VAbr->scale(-1.0);
    boost::shared_ptr<Matrix> VRas = Matrix::triplet(Matrix::triplet(Caocc_A_,V_B,P_A,true,false,false),S,Cavir_B_);
    VRas->scale(1.0);
    boost::shared_ptr<Matrix> VSbr = Matrix::triplet(Matrix::triplet(Caocc_B_,V_A,P_B,true,false,false),S,Cavir_A_);
    VSbr->scale(1.0);

    boost::shared_ptr<Matrix> Sas = Matrix::triplet(Caocc_A_,S,Cavir_B_,true,false,false);
    boost::shared_ptr<Matrix> Sbr = Matrix::triplet(Caocc_B_,S,Cavir_A_,true,false,false);

    boost::shared_ptr<Matrix> Qbr(Jbr->clone());
    Qbr->zero();
    Qbr->add(Jbr);
    Qbr->add(Kbr);
    Qbr->add(KObr);
    Qbr->add(JAbr);
    Qbr->add(JBbr);
    Qbr->add(VAbr);
    Qbr->add(VSbr);

    boost::shared_ptr<Matrix> Qas(Jas->clone());
    Qas->zero();
    Qas->add(Jas);
    Qas->add(Kas);
    Qas->add(KOas);
    Qas->add(JAas);
    Qas->add(JBas);
    Qas->add(VBas);
    Qas->add(VRas);

    boost::shared_ptr<Matrix> SBar = Matrix::triplet(Matrix::triplet(Caocc_A_,S,D_B,true,false,false),S,Cavir_A_);
    boost::shared_ptr<Matrix> SAbs = Matrix::triplet(Matrix::triplet(Caocc_B_,S,D_A,true,false,false),S,Cavir_B_);

    boost::shared_ptr<Matrix> Qar(Jar->clone());
    Qar->zero();
    Qar->add(Jar);
    Qar->add(Var);

    boost::shared_ptr<Matrix> Qbs(Jbs->clone());
    Qbs->zero();
    Qbs->add(Jbs);
    Qbs->add(Vbs);

    Jbr.reset();
    Kbr.reset();
    Jas.reset();
    Kas.reset();
    KOas.reset();
    KObr.reset();
    JBas.reset();
    JAbr.reset();
    Jbs.reset();
    Jar.reset();
    JAas.reset();
    JBbr.reset();
    Vbs.reset();
    Var.reset();
    VBas.reset();
    VAbr.reset();
    VRas.reset();
    VSbr.reset();

    S.reset();
    L_A.reset();
    D_A.reset();
    P_A.reset();
    V_A.reset();
    J_A.reset();
    K_A.reset();
    L_B.reset();
    D_B.reset();
    P_B.reset();
    V_B.reset();
    J_B.reset();
    K_B.reset();
    K_O.reset();

    vars_.clear();

    // => Memory <= //

    // => Integrals from the THCE <= //

    boost::shared_ptr<DFERI> df = DFERI::build(primary_,mp2fit_,Process::environment.options);
    df->clear();

    std::vector<boost::shared_ptr<Matrix> > Cs;
    Cs.push_back(Caocc_A_);
    Cs.push_back(Cavir_A_);
    Cs.push_back(Caocc_B_);
    Cs.push_back(Cavir_B_);
    Cs.push_back(Cr1);
    Cs.push_back(Cs1);
    Cs.push_back(Ca2);
    Cs.push_back(Cb2);
    Cs.push_back(Cr3);
    Cs.push_back(Cs3);
    Cs.push_back(Ca4);
    Cs.push_back(Cb4);
    boost::shared_ptr<Matrix> Call = Matrix::horzcat(Cs);
    Cs.clear();

    df->set_C(Call);
    df->set_memory(memory_ - Call->nrow() * Call->ncol());

    int offset = 0;
    df->add_space("a",offset,offset+Caocc_A_->colspi()[0]); offset += Caocc_A_->colspi()[0];
    df->add_space("r",offset,offset+Cavir_A_->colspi()[0]); offset += Cavir_A_->colspi()[0];
    df->add_space("b",offset,offset+Caocc_B_->colspi()[0]); offset += Caocc_B_->colspi()[0];
    df->add_space("s",offset,offset+Cavir_B_->colspi()[0]); offset += Cavir_B_->colspi()[0];
    df->add_space("r1",offset,offset+Cr1->colspi()[0]); offset += Cr1->colspi()[0];
    df->add_space("s1",offset,offset+Cs1->colspi()[0]); offset += Cs1->colspi()[0];
    df->add_space("a2",offset,offset+Ca2->colspi()[0]); offset += Ca2->colspi()[0];
    df->add_space("b2",offset,offset+Cb2->colspi()[0]); offset += Cb2->colspi()[0];
    df->add_space("r3",offset,offset+Cr3->colspi()[0]); offset += Cr3->colspi()[0];
    df->add_space("s3",offset,offset+Cs3->colspi()[0]); offset += Cs3->colspi()[0];
    df->add_space("a4",offset,offset+Ca4->colspi()[0]); offset += Ca4->colspi()[0];
    df->add_space("b4",offset,offset+Cb4->colspi()[0]); offset += Cb4->colspi()[0];

    // Disk stuff is all transposed for ab exposure, but transforms down to a or b first for speed

    df->add_pair_space("Aar", "a",  "r",  -1.0/2.0, true);
    df->add_pair_space("Abs", "b",  "s",  -1.0/2.0, true);
    df->add_pair_space("Bas", "a",  "s1", -1.0/2.0, true);
    df->add_pair_space("Bbr", "b",  "r1", -1.0/2.0, true);
    df->add_pair_space("Cas", "a2", "s",  -1.0/2.0, true);
    df->add_pair_space("Cbr", "b2", "r",  -1.0/2.0, true);
    df->add_pair_space("Dar", "a",  "r3", -1.0/2.0, true);
    df->add_pair_space("Dbs", "b",  "s3", -1.0/2.0, true);
    df->add_pair_space("Ear", "a4", "r",  -1.0/2.0, true);
    df->add_pair_space("Ebs", "b4", "s",  -1.0/2.0, true);

    Cr1.reset();
    Cs1.reset();
    Ca2.reset();
    Cb2.reset();
    Cr3.reset();
    Cs3.reset();
    Ca4.reset();
    Cb4.reset();
    Call.reset();

    df->print_header();
    df->compute();

    std::map<std::string, boost::shared_ptr<Tensor> >& ints = df->ints();

    boost::shared_ptr<Tensor> AarT = ints["Aar"];
    boost::shared_ptr<Tensor> AbsT = ints["Abs"];
    boost::shared_ptr<Tensor> BasT = ints["Bas"];
    boost::shared_ptr<Tensor> BbrT = ints["Bbr"];
    boost::shared_ptr<Tensor> CasT = ints["Cas"];
    boost::shared_ptr<Tensor> CbrT = ints["Cbr"];
    boost::shared_ptr<Tensor> DarT = ints["Dar"];
    boost::shared_ptr<Tensor> DbsT = ints["Dbs"];
    boost::shared_ptr<Tensor> EarT = ints["Ear"];
    boost::shared_ptr<Tensor> EbsT = ints["Ebs"];

    df.reset();

    // => Blocking <= //

    long int overhead = 0L;
    overhead += 5L * nT * na * nb;
    overhead += 2L * na * ns + 2L * nb * nr + 2L * na * nr + 2L * nb * ns;
    long int rem = memory_ - overhead;

    if (rem < 0L) {
        throw PSIEXCEPTION("Too little static memory for DFTSAPT::mp2_terms");
    }

    long int cost_r = 2L * naa * nQ + 2L * nab * nQ;
    long int max_r = rem / (2L * cost_r);
    long int max_s = max_r;
    max_r = (max_r > nr ? nr : max_r);
    max_s = (max_s > ns ? ns : max_s);
    if (max_s < 1L || max_s < 1L) {
        throw PSIEXCEPTION("Too little dynamic memory for DFTSAPT::mp2_terms");
    }

    // => Tensor Slices <= //

    boost::shared_ptr<Matrix> Aar(new Matrix("Aar",max_r*naa,nQ));
    boost::shared_ptr<Matrix> Abs(new Matrix("Abs",max_s*nab,nQ));
    boost::shared_ptr<Matrix> Bas(new Matrix("Bas",max_s*naa,nQ));
    boost::shared_ptr<Matrix> Bbr(new Matrix("Bbr",max_r*nab,nQ));
    boost::shared_ptr<Matrix> Cas(new Matrix("Cas",max_s*naa,nQ));
    boost::shared_ptr<Matrix> Cbr(new Matrix("Cbr",max_r*nab,nQ));
    boost::shared_ptr<Matrix> Dar(new Matrix("Dar",max_r*naa,nQ));
    boost::shared_ptr<Matrix> Dbs(new Matrix("Dbs",max_s*nab,nQ));

    // => Thread Work Arrays <= //

    std::vector<boost::shared_ptr<Matrix> > Tab;
    std::vector<boost::shared_ptr<Matrix> > Vab;
    std::vector<boost::shared_ptr<Matrix> > T2ab;
    std::vector<boost::shared_ptr<Matrix> > V2ab;
    std::vector<boost::shared_ptr<Matrix> > Iab;
    for (int t = 0; t < nT; t++) {
        Tab.push_back(boost::shared_ptr<Matrix>(new Matrix("Tab",naa,nab)));
        Vab.push_back(boost::shared_ptr<Matrix>(new Matrix("Vab",naa,nab)));
        T2ab.push_back(boost::shared_ptr<Matrix>(new Matrix("T2ab",na,nb)));
        V2ab.push_back(boost::shared_ptr<Matrix>(new Matrix("V2ab",na,nb)));
        Iab.push_back(boost::shared_ptr<Matrix>(new Matrix("Iab",naa,nb)));
    }

    // => Pointers <= //

    double** Aarp = Aar->pointer();
    double** Absp = Abs->pointer();
    double** Basp = Bas->pointer();
    double** Bbrp = Bbr->pointer();
    double** Casp = Cas->pointer();
    double** Cbrp = Cbr->pointer();
    double** Darp = Dar->pointer();
    double** Dbsp = Dbs->pointer();

    double** Sasp = Sas->pointer();
    double** Sbrp = Sbr->pointer();
    double** SBarp = SBar->pointer();
    double** SAbsp = SAbs->pointer();

    double** Qasp = Qas->pointer();
    double** Qbrp = Qbr->pointer();
    double** Qarp = Qar->pointer();
    double** Qbsp = Qbs->pointer();

    double*  eap  = eps_aocc_A_->pointer();
    double*  ebp  = eps_aocc_B_->pointer();
    double*  erp  = eps_avir_A_->pointer();
    double*  esp  = eps_avir_B_->pointer();

    // => File Pointers <= //

    FILE* Aarf = AarT->file_pointer();
    FILE* Absf = AbsT->file_pointer();
    FILE* Basf = BasT->file_pointer();
    FILE* Bbrf = BbrT->file_pointer();
    FILE* Casf = CasT->file_pointer();
    FILE* Cbrf = CbrT->file_pointer();
    FILE* Darf = DarT->file_pointer();
    FILE* Dbsf = DbsT->file_pointer();
    FILE* Earf = EarT->file_pointer();
    FILE* Ebsf = EbsT->file_pointer();

    // => Slice D + E -> D <= //

    fseek(Darf,0L,SEEK_SET);
    fseek(Earf,0L,SEEK_SET);
    for (int rstart = 0; rstart < nr; rstart += max_r) {
        int nrblock = (rstart + max_r >= nr ? nr - rstart : max_r);
        fread(Darp[0],sizeof(double),nrblock*naQ,Darf);
        fread(Aarp[0],sizeof(double),nrblock*naQ,Earf);
        C_DAXPY(nrblock*naQ,1.0,Aarp[0],1,Darp[0],1);
        fseek(Darf,sizeof(double)*rstart*naQ,SEEK_SET);
        fwrite(Darp[0],sizeof(double),nrblock*naQ,Darf);
    }

    fseek(Dbsf,0L,SEEK_SET);
    fseek(Ebsf,0L,SEEK_SET);
    for (int sstart = 0; sstart < ns; sstart += max_s) {
        int nsblock = (sstart + max_s >= ns ? ns - sstart : max_s);
        fread(Dbsp[0],sizeof(double),nsblock*nbQ,Dbsf);
        fread(Absp[0],sizeof(double),nsblock*nbQ,Ebsf);
        C_DAXPY(nsblock*nbQ,1.0,Absp[0],1,Dbsp[0],1);
        fseek(Dbsf,sizeof(double)*sstart*nbQ,SEEK_SET);
        fwrite(Dbsp[0],sizeof(double),nsblock*nbQ,Dbsf);
    }

    // => Targets <= //

    double Disp20 = 0.0;
    double ExchDisp20 = 0.0;

    // => Local Targets <= //

    std::vector<boost::shared_ptr<Matrix> > E_disp20_threads;
    std::vector<boost::shared_ptr<Matrix> > E_exch_disp20_threads;
    for (int t = 0; t < nT; t++) {
        E_disp20_threads.push_back(boost::shared_ptr<Matrix>(new Matrix("E_disp20",na,nb)));
        E_exch_disp20_threads.push_back(boost::shared_ptr<Matrix>(new Matrix("E_exch_disp20",na,nb)));
    }

    // => MO => LO Transform <= //

    double** UAp = Uocc_A_->pointer();
    double** UBp = Uocc_B_->pointer();

    // ==> Master Loop <== //

    fseek(Aarf,0L,SEEK_SET);
    fseek(Bbrf,0L,SEEK_SET);
    fseek(Cbrf,0L,SEEK_SET);
    fseek(Darf,0L,SEEK_SET);
    for (int rstart = 0; rstart < nr; rstart += max_r) {
        int nrblock = (rstart + max_r >= nr ? nr - rstart : max_r);

        fread(Aarp[0],sizeof(double),nrblock*naQ,Aarf);
        fread(Bbrp[0],sizeof(double),nrblock*nbQ,Bbrf);
        fread(Cbrp[0],sizeof(double),nrblock*nbQ,Cbrf);
        fread(Darp[0],sizeof(double),nrblock*naQ,Darf);

        fseek(Absf,0L,SEEK_SET);
        fseek(Basf,0L,SEEK_SET);
        fseek(Casf,0L,SEEK_SET);
        fseek(Dbsf,0L,SEEK_SET);
        for (int sstart = 0; sstart < ns; sstart += max_s) {
            int nsblock = (sstart + max_s >= ns ? ns - sstart : max_s);

            fread(Absp[0],sizeof(double),nsblock*nbQ,Absf);
            fread(Basp[0],sizeof(double),nsblock*naQ,Basf);
            fread(Casp[0],sizeof(double),nsblock*naQ,Casf);
            fread(Dbsp[0],sizeof(double),nsblock*nbQ,Dbsf);

            long int nrs = nrblock * nsblock;

            #pragma omp parallel for schedule(dynamic) reduction(+: Disp20, ExchDisp20)
            for (long int rs = 0L; rs < nrs; rs++) {
                int r = rs / nsblock;
                int s = rs % nsblock;

                int thread = 0;
                #ifdef _OPENMP
                    thread = omp_get_thread_num();
                #endif

                double** E_disp20Tp = E_disp20_threads[thread]->pointer();
                double** E_exch_disp20Tp = E_exch_disp20_threads[thread]->pointer();
                
                double** Tabp  = Tab[thread]->pointer();
                double** Vabp  = Vab[thread]->pointer();
                double** T2abp = T2ab[thread]->pointer();
                double** V2abp = V2ab[thread]->pointer();
                double** Iabp  = Iab[thread]->pointer();

                // => Amplitudes, Disp20 <= //

                C_DGEMM('N','T',naa,nab,nQ,1.0,Aarp[(r)*naa],nQ,Absp[(s)*nab],nQ,0.0,Vabp[0],nab);
                for (int a = 0; a < naa; a++) {
                    for (int b = 0; b < nab; b++) {
                        Tabp[a][b] = Vabp[a][b] / (eap[a] + ebp[b] - erp[r + rstart] - esp[s + sstart]);
                    }
                }

                C_DGEMM('N','N',naa,nb,nab,1.0,Tabp[0],nab,UBp[nfb],nb,0.0,Iabp[0],nb);
                C_DGEMM('T','N',na,nb,naa,1.0,UAp[nfa],na,Iabp[0],nb,0.0,T2abp[0],nb);
                C_DGEMM('N','N',naa,nb,nab,1.0,Vabp[0],nab,UBp[nfb],nb,0.0,Iabp[0],nb);
                C_DGEMM('T','N',na,nb,naa,1.0,UAp[nfa],na,Iabp[0],nb,0.0,V2abp[0],nb);

                for (int a = 0; a < na; a++) {
                    for (int b = 0; b < nb; b++) {
                        E_disp20Tp[a][b] += 4.0 * T2abp[a][b] * V2abp[a][b];
                        Disp20 += 4.0 * T2abp[a][b] * V2abp[a][b];
                    }
                }

                // => Exch-Disp20 <= //

                // > Q1-Q3 < //

                C_DGEMM('N','T',naa,nab,nQ,1.0,Basp[(s)*naa],nQ,Bbrp[(r)*nab],nQ,0.0,Vabp[0],nab);
                C_DGEMM('N','T',naa,nab,nQ,1.0,Casp[(s)*naa],nQ,Cbrp[(r)*nab],nQ,1.0,Vabp[0],nab);
                C_DGEMM('N','T',naa,nab,nQ,1.0,Aarp[(r)*naa],nQ,Dbsp[(s)*nab],nQ,1.0,Vabp[0],nab);
                C_DGEMM('N','T',naa,nab,nQ,1.0,Darp[(r)*naa],nQ,Absp[(s)*nab],nQ,1.0,Vabp[0],nab);

                // > V,J,K < //

                C_DGER(naa,nab,1.0,&Sasp[0][s + sstart], ns,&Qbrp[0][r + rstart], nr,Vabp[0],nab);
                C_DGER(naa,nab,1.0,&Qasp[0][s + sstart], ns,&Sbrp[0][r + rstart], nr,Vabp[0],nab);
                C_DGER(naa,nab,1.0,&Qarp[0][r + rstart], nr,&SAbsp[0][s + sstart],ns,Vabp[0],nab);
                C_DGER(naa,nab,1.0,&SBarp[0][r + rstart],nr,&Qbsp[0][s + sstart], ns,Vabp[0],nab);

                C_DGEMM('N','N',naa,nb,nab,1.0,Vabp[0],nab,UBp[nfb],nb,0.0,Iabp[0],nb);
                C_DGEMM('T','N',na,nb,naa,1.0,UAp[nfa],na,Iabp[0],nb,0.0,V2abp[0],nb);

                for (int a = 0; a < na; a++) {
                    for (int b = 0; b < nb; b++) {
                        E_exch_disp20Tp[a][b] -= 2.0 * T2abp[a][b] * V2abp[a][b];
                        ExchDisp20 -= 2.0 * T2abp[a][b] * V2abp[a][b];
                    }
                }
            }
        }
    }

    boost::shared_ptr<Matrix> E_disp20(new Matrix("E_disp20", na, nb));
    boost::shared_ptr<Matrix> E_exch_disp20(new Matrix("E_exch_disp20", na, nb));
    double** E_disp20p = E_disp20->pointer();
    double** E_exch_disp20p = E_exch_disp20->pointer();

    for (int t = 0; t < nT; t++) {
        E_disp20->add(E_disp20_threads[t]);
        E_exch_disp20->add(E_exch_disp20_threads[t]);
    }

    boost::shared_ptr<Matrix> E_disp(new Matrix("E_disp (a x b)", na, nb));
    double** E_dispp = E_disp->pointer();

    for (int a = 0; a < na; a++) {
        for (int b = 0; b < nb; b++) {
            E_dispp[a][b] = E_disp20p[a][b] +
                            E_exch_disp20p[a][b];
        }
    }

    //E_disp20->print();
    //E_exch_disp20->print();
    //E_disp->print();

    energies_["Disp20"] = Disp20;
    energies_["Exch-Disp20"] = ExchDisp20;
    fprintf(outfile,"    Disp20              = %18.12lf H\n",Disp20);
    fprintf(outfile,"    Exch-Disp20         = %18.12lf H\n",ExchDisp20);
    fprintf(outfile,"\n");
    fflush(outfile);

    vis_->vars()["Disp_ab"] = E_disp;
}
void ASAPT::analyze()
{
    vis_->analyze();
}

}}
