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
#include <libmints/mints.h>
#include <libmints/sieve.h>
#include <libmints/local.h>
#include <libfock/jk.h>
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
double ASAPT::compute_energy()
{
    energies_["HF"] = E_dimer_ - E_monomer_A_ - E_monomer_B_; // TODO: get dHF loaded correctly

    print_header();

    fock_terms();

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
    fprintf(outfile, "  To the pessimist, the glass is half empty.\n");
    fprintf(outfile, "  To the optimist, the glass is half full.\n");
    fprintf(outfile, "  To the engineer, the glass is twice as big as it needs to be.\n");
    fprintf(outfile, "\n");
}
void ASAPT::localize()
{
    fprintf(outfile," LOCALIZATION:\n\n");

    fprintf(outfile,"  Local Orbitals for Monomer A:\n\n");

    boost::shared_ptr<Localizer> localA = Localizer::build(primary_, Cocc_A_, Process::environment.options);
    localA->localize();
    Locc_A_ = localA->L();
    Uocc_A_ = localA->U();;

    fprintf(outfile,"  Local Orbitals for Monomer B:\n\n");

    boost::shared_ptr<Localizer> localB = Localizer::build(primary_, Cocc_B_, Process::environment.options);
    localB->localize();
    Locc_B_ = localB->L();
    Uocc_B_ = localB->U();;

    fflush(outfile);
}
void ASAPT::populate()
{
    fprintf(outfile,"  POPULATION:\n\n");

    // => Sizing <= //

    int na = Caocc_A_->colspi()[0];
    int nb = Caocc_B_->colspi()[0];

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

    // => S matrix <= //

    boost::shared_ptr<Matrix> S = vars_["S"];

    // => Raw population metrics <= //

    boost::shared_ptr<Matrix> QA(new Matrix("QA (atoms x a)",dimer_->natom(), na)); 
    boost::shared_ptr<Matrix> QB(new Matrix("QB (atoms x b)",dimer_->natom(), nb)); 
    double** QAp = QA->pointer();
    double** QBp = QB->pointer();

    boost::shared_ptr<Matrix> L_A = Locc_A_;
    boost::shared_ptr<Matrix> L_B = Locc_B_;

    if (population_type_ == "LOWDIN") {

        boost::shared_ptr<Matrix> S12(S->clone());
        S12->copy(S);
        S12->power(0.5);
        boost::shared_ptr<Matrix> S12A = doublet(S12,L_A);
        boost::shared_ptr<Matrix> S12B = doublet(S12,L_B);
        double** S12Ap = S12A->pointer();
        double** S12Bp = S12B->pointer();
        
        for (int a = 0; a < na; a++) {
            for (int M = 0; M < primary_->nshell(); M++) {
                int aM = primary_->shell(M).ncenter();
                int nM = primary_->shell(M).nfunction();
                int oM = primary_->shell(M).function_index();
                for (int m = 0; m < nM; m++) {
                    QAp[aM][a] += (S12Ap[m+oM][a]*S12Ap[m+oM][a]);
                }
            }
        }

        for (int b = 0; b < nb; b++) {
            for (int M = 0; M < primary_->nshell(); M++) {
                int aM = primary_->shell(M).ncenter();
                int nM = primary_->shell(M).nfunction();
                int oM = primary_->shell(M).function_index();
                for (int m = 0; m < nM; m++) {
                    QBp[aM][b] += (S12Bp[m+oM][b]*S12Bp[m+oM][b]);
                }
            }
        }

    } else if (population_type_ == "MULLIKEN") {
        
        boost::shared_ptr<Matrix> SA = doublet(S,L_A);
        boost::shared_ptr<Matrix> SB = doublet(S,L_B);
        double** SAp = SA->pointer();
        double** SBp = SB->pointer();
        double** LAp = L_A->pointer();
        double** LBp = L_B->pointer();

        for (int a = 0; a < na; a++) {
            for (int M = 0; M < primary_->nshell(); M++) {
                int aM = primary_->shell(M).ncenter();
                int nM = primary_->shell(M).nfunction();
                int oM = primary_->shell(M).function_index();
                for (int m = 0; m < nM; m++) {
                    QAp[aM][a] += (SAp[m+oM][a]*LAp[m+oM][a]);
                }
            }
        }

        for (int b = 0; b < nb; b++) {
            for (int M = 0; M < primary_->nshell(); M++) {
                int aM = primary_->shell(M).ncenter();
                int nM = primary_->shell(M).nfunction();
                int oM = primary_->shell(M).function_index();
                for (int m = 0; m < nM; m++) {
                    QBp[aM][b] += (SBp[m+oM][b]*LBp[m+oM][b]);
                }
            }
        }

    }  

    boost::shared_ptr<Matrix> Q2A(new Matrix("QA (A x a)", nA, na));
    boost::shared_ptr<Matrix> Q2B(new Matrix("QB (B x b)", nB, nb));
    double** Q2Ap = Q2A->pointer();
    double** Q2Bp = Q2B->pointer();

    for (int A = 0; A < nA; A++) {
        for (int a = 0; a < na; a++) {
            Q2Ap[A][a] = QAp[cA[A]][a];
        }
    }

    for (int B = 0; B < nB; B++) {
        for (int b = 0; b < nb; b++) {
            Q2Bp[B][b] = QBp[cB[B]][b];
        }
    }

    std::vector<double> normA;
    std::vector<double> normB;
    double maxA = 0.0;   
    double maxB = 0.0;
 
    for (int a = 0; a < na; a++) {
        double val = 0.0;
        for (int A = 0; A < nA; A++) {
            val += Q2Ap[A][a];
        }
        normA.push_back(1.0 - val);
        maxA = (maxA > fabs(1.0 - val) ? maxA : fabs(1.0 - val));
        C_DSCAL(nA,1.0/val,&Q2Ap[0][a],na);
    }
    for (int b = 0; b < nb; b++) {
        double val = 0.0;
        for (int B = 0; B < nB; B++) {
            val += Q2Bp[B][b];
        }
        normB.push_back(1.0 - val);
        maxB = (maxB > fabs(1.0 - val) ? maxB : fabs(1.0 - val));
        C_DSCAL(nB,1.0/val,&Q2Bp[0][b],nb);
    }
    
    if (print_ >= 0) {
        fprintf(outfile,"    Population type        = %11s\n", population_type_.c_str());
        fprintf(outfile,"    Max ghost population A = %11.3E\n",maxA);
        fprintf(outfile,"    Max ghost population B = %11.3E\n",maxB);
        fprintf(outfile,"\n");
    }
    if (print_ >= 1) {
        fprintf(outfile,"    Ghost populations for Monomoer A:\n");
        for (int a = 0; a < na; a++) {
            fprintf(outfile,"    %4d %11.3E\n", a+1, normA[a]);
        } 
        fprintf(outfile,"\n");

        fprintf(outfile,"    Ghost populations for Monomer B:\n");
        for (int b = 0; b < nb; b++) {
            fprintf(outfile,"    %4d %11.3E\n", b+1, normA[b]);
        } 
        fprintf(outfile,"\n");
    }
    if (print_ >= 2) {
        Q2A->print(); 
        Q2B->print(); 
    }

    fflush(outfile);

    Q_A_ = Q2A; 
    Q_B_ = Q2B; 

    R_A_ = doublet(Q_A_,Uocc_A_,false,true);
    R_B_ = doublet(Q_B_,Uocc_B_,false,true);
}
void ASAPT::elst()
{
    fprintf(outfile,"  ELECTROSTATICS:\n\n");

    // ==> Sizing <== //

    int nn = primary_->nbf();
    int na = Caocc_A_->colspi()[0];
    int nb = Caocc_B_->colspi()[0];

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

    // ==> DF ERI Setup (JKFIT Type, in Local Basis) <== //    

    std::vector<int> mAlist;
    std::vector<int> mBlist;
    std::vector<int> nlist;
    mAlist.push_back(0);
    mBlist.push_back(1);

    dimer_->set_orientation_fixed(true);
    dimer_->set_com_fixed(true);
    boost::shared_ptr<Molecule> mA = dimer_->extract_subsets(mAlist,nlist);
    boost::shared_ptr<Molecule> mB = dimer_->extract_subsets(mBlist,nlist);

    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    boost::shared_ptr<BasisSet> jkfitA = BasisSet::construct(parser, mA, "DF_BASIS_ELST");
    boost::shared_ptr<BasisSet> jkfitB = BasisSet::construct(parser, mB, "DF_BASIS_ELST");
    int nQA = jkfitA->nbf();
    int nQB = jkfitB->nbf();

    //mA->print();
    //mB->print();
    
    boost::shared_ptr<DFERI> dfA = DFERI::build(primary_,jkfitA,Process::environment.options);
    dfA->clear();
    std::vector<boost::shared_ptr<Matrix> > CsA;
    CsA.push_back(Cocc_A_);
    CsA.push_back(Cocc_B_);
    CsA.push_back(Cvir_B_);
    boost::shared_ptr<Matrix> CallA = Matrix::horzcat(CsA);
    CsA.clear();
    dfA->set_C(CallA);
    dfA->set_memory(memory_);
    int offsetA = 0;
    dfA->add_space("a",offsetA,offsetA + na); offsetA += na;
    dfA->add_space("b",offsetA,offsetA + nb); offsetA += nb;
    dfA->add_space("s",offsetA,offsetA + ns); offsetA += ns;
    dfA->add_pair_space("Aaa2", "a", "a", -1.0);
    dfA->add_pair_space("Vbs2", "b", "s", 0.0);
    dfA->set_keep_raw_integrals(true);
    fprintf(outfile,"  ==> Local DF for Monomer A <==\n\n");
    dfA->print_header();
    dfA->compute();
    std::map<std::string, boost::shared_ptr<Tensor> >& intsA = dfA->ints();
    boost::shared_ptr<Tensor> AaaT = intsA["Aaa2"];
    boost::shared_ptr<Tensor> VbsT = intsA["Vbs2_temp"]; // (P|bs) striping
    dfA.reset();

    boost::shared_ptr<DFERI> dfB = DFERI::build(primary_,jkfitB,Process::environment.options);
    dfB->clear();
    std::vector<boost::shared_ptr<Matrix> > CsB;
    CsB.push_back(Cocc_B_);
    CsB.push_back(Cocc_A_);
    CsB.push_back(Cvir_A_);
    boost::shared_ptr<Matrix> CallB = Matrix::horzcat(CsB);
    CsB.clear();
    dfB->set_C(CallB);
    dfB->set_memory(memory_);
    int offsetB = 0;
    dfB->add_space("b",offsetB,offsetB + nb); offsetB += nb;
    dfB->add_space("a",offsetB,offsetB + na); offsetB += na;
    dfB->add_space("r",offsetB,offsetB + nr); offsetB += nr;
    dfB->add_pair_space("Abb2", "b", "b", -1.0);
    dfB->add_pair_space("Var2", "a", "r", 0.0);
    dfB->set_keep_raw_integrals(true);
    fprintf(outfile,"  ==> Local DF for Monomer B <==\n\n");
    dfB->print_header();
    dfB->compute();
    std::map<std::string, boost::shared_ptr<Tensor> >& intsB = dfB->ints();
    boost::shared_ptr<Tensor> AbbT = intsB["Abb2"];
    boost::shared_ptr<Tensor> VarT = intsB["Var2_temp"]; // (Q|ar) striping
    dfB.reset();

    // ==> DF Nuclear Potential Setup (JKFIT Type, in Local Basis) <== //

    boost::shared_ptr<Matrix> VAQ(new Matrix("VAQ",nA,nQB));
    boost::shared_ptr<Matrix> VBQ(new Matrix("VBQ",nB,nQA));
    double** VAQp = VAQ->pointer();
    double** VBQp = VBQ->pointer();

    boost::shared_ptr<Matrix> Zxyz(new Matrix("Zxyz",1,4));
    double** Zxyzp = Zxyz->pointer();

    boost::shared_ptr<IntegralFactory> VfactA(new IntegralFactory(jkfitB,BasisSet::zero_ao_basis_set()));
    boost::shared_ptr<PotentialInt> VintA(static_cast<PotentialInt*>(VfactA->ao_potential()));
    VintA->set_charge_field(Zxyz);
    boost::shared_ptr<Matrix> VtempA(new Matrix("Vtemp",nQB,1));
    double** VtempAp = VtempA->pointer();
    for (int A = 0; A < nA; A++) {
        VtempA->zero();
        Zxyzp[0][0] = mA->Z(A);
        Zxyzp[0][1] = mA->x(A);
        Zxyzp[0][2] = mA->y(A);
        Zxyzp[0][3] = mA->z(A); 
        VintA->compute(VtempA);
        ::memcpy(VAQp[A],VtempAp[0],sizeof(double)*nQB);
    }

    boost::shared_ptr<IntegralFactory> VfactB(new IntegralFactory(jkfitA,BasisSet::zero_ao_basis_set()));
    boost::shared_ptr<PotentialInt> VintB(static_cast<PotentialInt*>(VfactB->ao_potential()));
    VintB->set_charge_field(Zxyz);
    boost::shared_ptr<Matrix> VtempB(new Matrix("Vtemp",nQA,1));
    double** VtempBp = VtempB->pointer();
    for (int B = 0; B < nB; B++) {
        VtempB->zero();
        Zxyzp[0][0] = mB->Z(B);
        Zxyzp[0][1] = mB->x(B);
        Zxyzp[0][2] = mB->y(B);
        Zxyzp[0][3] = mB->z(B); 
        VintB->compute(VtempB);
        ::memcpy(VBQp[B],VtempBp[0],sizeof(double)*nQA);
    }

    boost::shared_ptr<Matrix> AaQ(new Matrix("AaQ",na,nQA));
    boost::shared_ptr<Matrix> AbQ(new Matrix("AbQ",nb,nQB));
    double** AaQp = AaQ->pointer();
    double** AbQp = AbQ->pointer();
    FILE* Aaaf = AaaT->file_pointer();
    FILE* Abbf = AbbT->file_pointer();

    for (int a = 0; a < na; a++) {
        fseek(Aaaf,(a*na+a)*(size_t)nQA*sizeof(double),SEEK_SET);
        fread(AaQp[a],sizeof(double),nQA,Aaaf);
    }

    for (int b = 0; b < nb; b++) {
        fseek(Abbf,(b*nb+b)*(size_t)nQB*sizeof(double),SEEK_SET);
        fread(AbQp[b],sizeof(double),nQB,Abbf);
    }

    boost::shared_ptr<Matrix> dA(new Matrix("dA",nQA,1));
    double** dAp = dA->pointer();
    for (int a = 0; a < na; a++) {
        C_DAXPY(nQA,1.0,AaQp[a],1,dAp[0],1);
    }

    boost::shared_ptr<Matrix> dB(new Matrix("dB",nQB,1));
    double** dBp = dB->pointer();
    for (int b = 0; b < nb; b++) {
        C_DAXPY(nQB,1.0,AbQp[b],1,dBp[0],1);
    }

    boost::shared_ptr<Matrix> J(new Matrix("J",nQA,nQB));
    double** Jp = J->pointer();
    boost::shared_ptr<IntegralFactory> ABfact(new IntegralFactory(jkfitA,BasisSet::zero_ao_basis_set(),jkfitB,BasisSet::zero_ao_basis_set()));
    boost::shared_ptr<TwoBodyAOInt> JABint(ABfact->eri());
    for (int P = 0; P < jkfitA->nshell(); P++) {
        for (int Q = 0; Q < jkfitB->nshell(); Q++) {
            JABint->compute_shell(P,0,Q,0);
            const double* buffer = JABint->buffer();
            int nP = jkfitA->shell(P).nfunction(); 
            int nQ = jkfitB->shell(Q).nfunction(); 
            int oP = jkfitA->shell(P).function_index(); 
            int oQ = jkfitB->shell(Q).function_index(); 
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Jp[p+oP][q+oQ] = (*buffer++);
                }
            }
        }
    }

    std::vector<int> assA;
    for (int P = 0; P < jkfitA->nshell(); P++) {
        int cP = jkfitA->shell(P).ncenter();
        int nP = jkfitA->shell(P).nfunction(); 
        for (int p = 0; p < nP; p++) {
            assA.push_back(cP);    
        }
    }
    assA.push_back(nA);

    std::vector<int> assB;
    for (int P = 0; P < jkfitB->nshell(); P++) {
        int cP = jkfitB->shell(P).ncenter();
        int nP = jkfitB->shell(P).nfunction(); 
        for (int p = 0; p < nP; p++) {
            assB.push_back(cP);    
        }
    }
    assB.push_back(nB);

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

    // => a <-> b <- //

    for (int A = 0; A < nQA; A++) {
        for (int B = 0; B < nQB; B++) {
            double val = 4.0 * dAp[A][0] * Jp[A][B] * dBp[B][0]; 
            Elst10_terms[2] += val;
            Elst_atomsp[assA[A]][assB[B]] += val;
        }
    }

    // => A <-> b <- //

    for (int A = 0; A < nA; A++) {
        for (int B = 0; B < nQB; B++) {
            double val = 2.0 * VAQp[A][B] * dBp[B][0]; 
            Elst10_terms[1] += val;
            Elst_atomsp[A][assB[B]] += val;
        }
    }

    // => a <-> B <- //

    for (int A = 0; A < nQA; A++) {
        for (int B = 0; B < nB; B++) {
            double val = 2.0 * dAp[A][0] * VBQp[B][A]; 
            Elst10_terms[0] += val;
            Elst_atomsp[assA[A]][B] += val;
        }
    }

    local_vars_["Elst"] = Elst_atoms;

    for (int k = 0; k < Elst10_terms.size(); k++) {
        Elst10 += Elst10_terms[k];
    }
    if (debug_) {
        for (int k = 0; k < Elst10_terms.size(); k++) {
            fprintf(outfile,"    Elst10,r (%1d)        = %18.12lf H\n",k+1,Elst10_terms[k]);
        }
    }
    energies_["Elst10,r"] = Elst10;
    fprintf(outfile,"    Elst10,r (L-DF)     = %18.12lf H\n",Elst10);
    fprintf(outfile,"\n");
    fflush(outfile);

    // ==> Setup Atomic Electrostatic Fields for Induction <== //

    boost::shared_ptr<Tensor> WBarT = DiskTensor::build("WBar", "nB", nB, "na", na, "nr", nr, false, false);
    FILE* WBarf = WBarT->file_pointer();
    boost::shared_ptr<Tensor> WAbsT = DiskTensor::build("WAbs", "nA", nA, "nb", nb, "ns", ns, false, false);
    FILE* WAbsf = WAbsT->file_pointer();

    // => Nuclear Part (PITA) <= //
    
    boost::shared_ptr<IntegralFactory> Vfact(new IntegralFactory(primary_));
    boost::shared_ptr<PotentialInt> Vint(static_cast<PotentialInt*>(Vfact->ao_potential()));
    Vint->set_charge_field(Zxyz);
    boost::shared_ptr<Matrix> Vtemp(new Matrix("Vtemp",nn,nn));

    for (int A = 0; A < nA; A++) {
        Vtemp->zero();
        Zxyzp[0][0] = mA->Z(A);
        Zxyzp[0][1] = mA->x(A);
        Zxyzp[0][2] = mA->y(A);
        Zxyzp[0][3] = mA->z(A); 
        Vint->compute(Vtemp);
        boost::shared_ptr<Matrix> Vbs = triplet(Cocc_B_,Vtemp,Cvir_B_,true,false,false);
        double** Vbsp = Vbs->pointer();
        fwrite(Vbsp[0],sizeof(double),nb*ns,WAbsf);
    }

    for (int B = 0; B < nB; B++) {
        Vtemp->zero();
        Zxyzp[0][0] = mB->Z(B);
        Zxyzp[0][1] = mB->x(B);
        Zxyzp[0][2] = mB->y(B);
        Zxyzp[0][3] = mB->z(B); 
        Vint->compute(Vtemp);
        boost::shared_ptr<Matrix> Var = triplet(Cocc_A_,Vtemp,Cvir_A_,true,false,false);
        double** Varp = Var->pointer();
        fwrite(Varp[0],sizeof(double),na*nr,WBarf);
    }

    // => Electronic Part <= //

    FILE* Vbsf = VbsT->file_pointer();
    fseek(Vbsf,0L,SEEK_SET);
    boost::shared_ptr<Matrix> Jbs(new Matrix("Jbs",nb,ns));
    boost::shared_ptr<Matrix> Vbs(new Matrix("Vbs",nb,ns));
    double** Vbsp = Vbs->pointer();
    double** Jbsp = Jbs->pointer();
    int centerA = 0;
    for (int P = 0; P < nQA; P++) {
        fread(Vbsp[0],sizeof(double),nb*ns,Vbsf);
        C_DAXPY(nb*ns,dAp[P][0],Vbsp[0],1,Jbsp[0],1);
        if (assA[P+1] != centerA) {
            fseek(WAbsf,centerA*nb*ns*sizeof(double),SEEK_SET);
            fread(Vbsp[0],sizeof(double),nb*ns,WAbsf);
            Jbs->scale(2.0);
            Vbs->add(Jbs); 
            fseek(WAbsf,centerA*nb*ns*sizeof(double),SEEK_SET);
            fwrite(Vbsp[0],sizeof(double),nb*ns,WAbsf);
            Jbs->zero();
            centerA++;
        }
    }

    FILE* Varf = VarT->file_pointer();
    fseek(Varf,0L,SEEK_SET);
    boost::shared_ptr<Matrix> Jar(new Matrix("Jar",na,nr));
    boost::shared_ptr<Matrix> Var(new Matrix("Var",na,nr));
    double** Varp = Var->pointer();
    double** Jarp = Jar->pointer();
    int centerB = 0;
    for (int P = 0; P < nQB; P++) {
        fread(Varp[0],sizeof(double),na*nr,Varf);
        C_DAXPY(na*nr,dBp[P][0],Varp[0],1,Jarp[0],1);
        if (assB[P+1] != centerB) {
            fseek(WBarf,centerB*na*nr*sizeof(double),SEEK_SET);
            fread(Varp[0],sizeof(double),na*nr,WBarf);
            Jar->scale(2.0);
            Var->add(Jar); 
            fseek(WBarf,centerB*na*nr*sizeof(double),SEEK_SET);
            fwrite(Varp[0],sizeof(double),na*nr,WBarf);
            Jar->zero();
            centerB++;
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
    int na = Caocc_A_->colspi()[0];
    int nb = Caocc_B_->colspi()[0];

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
    Cs.push_back(Cavir_A_);
    Cs.push_back(L_B);
    Cs.push_back(Cavir_B_);
    boost::shared_ptr<Matrix> Call = Matrix::horzcat(Cs);
    Cs.clear();

    df->set_C(Call);
    df->set_memory(memory_);

    int offset = 0;
    df->add_space("a",offset,offset+Caocc_A_->colspi()[0]); offset += Caocc_A_->colspi()[0];
    df->add_space("r",offset,offset+Cavir_A_->colspi()[0]); offset += Cavir_A_->colspi()[0];
    df->add_space("b",offset,offset+Caocc_B_->colspi()[0]); offset += Caocc_B_->colspi()[0];
    df->add_space("s",offset,offset+Cavir_B_->colspi()[0]); offset += Cavir_B_->colspi()[0];

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

    boost::shared_ptr<Matrix> WAbs = triplet(Locc_B_,W_A,Cvir_B_,true,false,false);
    boost::shared_ptr<Matrix> WBar = triplet(Locc_A_,W_B,Cvir_A_,true,false,false);
    double** WBarp = WBar->pointer();
    double** WAbsp = WAbs->pointer();
    
    W_A.reset();
    W_B.reset();

    // ==> Exchange S^2 Computation <== //

    double Exch10_2 = 0.0;
    std::vector<double> Exch10_2_terms;
    Exch10_2_terms.resize(3);

    boost::shared_ptr<Matrix> Sab = triplet(L_A,S,L_B,true,false,false);
    boost::shared_ptr<Matrix> Sba = triplet(L_B,S,L_A,true,false,false);
    boost::shared_ptr<Matrix> Sas = triplet(L_A,S,Cvir_B_,true,false,false);
    boost::shared_ptr<Matrix> Sbr = triplet(L_B,S,Cvir_A_,true,false,false);
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

    //E_exch->print();

    boost::shared_ptr<Matrix> Exch_atoms(new Matrix("Exch (A x B)", nA, nB));
    double** Exch_atomsp = Exch_atoms->pointer();

    double** Q2Ap = Q_A_->pointer();
    double** Q2Bp = Q_B_->pointer();

    for (int b = 0; b < nb; b++) {
        for (int B = 0; B < nB; B++) {
            for (int a = 0; a < na; a++) {
                for (int A = 0; A < nA; A++) {
                    Exch_atomsp[A][B] += Q2Ap[A][a] * Q2Bp[B][b] * E_exchp[a][b];
                }
            }
        }
    }

    //Exch_atoms->print();
    local_vars_["Exch"] = Exch_atoms;

    for (int k = 0; k < Exch10_2_terms.size(); k++) {
        Exch10_2 += Exch10_2_terms[k];
    }
    if (debug_) {
        for (int k = 0; k < Exch10_2_terms.size(); k++) {
            fprintf(outfile,"    Exch10(S^2) (%1d)     = %18.12lf H\n",k+1,Exch10_2_terms[k]);
        }
    }
    energies_["Exch10(S^2)"] = Exch10_2;
    fprintf(outfile,"    Exch10(S^2)         = %18.12lf H\n",Exch10_2);
    fprintf(outfile, "\n");
    fflush(outfile);
}
void ASAPT::ind()
{
    fprintf(outfile, "  INDUCTION:\n\n");

    // => Sizing <= //

    int nn = primary_->nbf();

    int na = Caocc_A_->colspi()[0];
    int nb = Caocc_B_->colspi()[0];
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
    size_t naQ = na * (size_t) nQ;
    size_t nbQ = nb * (size_t) nQ;

    int nT = 1;
    #ifdef _OPENMP
        nT = omp_get_max_threads();
    #endif

    // ==> Stack Variables <== //

    double** RAp = R_A_->pointer();
    double** RBp = R_B_->pointer();
    
    double*  eap = eps_occ_A_->pointer();
    double*  ebp = eps_occ_B_->pointer();
    double*  erp = eps_vir_A_->pointer();
    double*  esp = eps_vir_B_->pointer();

    FILE* WAbsf = tensors_["WAbs"]->file_pointer();
    FILE* WBarf = tensors_["WBar"]->file_pointer();
    
    boost::shared_ptr<Matrix> S   = vars_["S"];
    boost::shared_ptr<Matrix> V_A = vars_["V_A"];
    boost::shared_ptr<Matrix> J_A = vars_["J_A"];
    boost::shared_ptr<Matrix> V_B = vars_["V_B"];
    boost::shared_ptr<Matrix> J_B = vars_["J_B"];
    boost::shared_ptr<Matrix> L_A = Locc_A_;
    boost::shared_ptr<Matrix> L_B = Locc_B_;

    // ==> Targets <== //

    boost::shared_ptr<Matrix> IndAB_terms(new Matrix("Ind [A<-B] (a x B)", na, nB));
    boost::shared_ptr<Matrix> IndBA_terms(new Matrix("Ind [B<-A] (b x A)", nb, nA));
    double** IndAB_termsp = IndAB_terms->pointer();
    double** IndBA_termsp = IndBA_terms->pointer();

    double Ind20_AB = 0.0; 
    double Ind20_BA = 0.0;

    // ==> MO Amplitudes/Sources <== //

    boost::shared_ptr<Matrix> xA(new Matrix("xA",na,nr));
    boost::shared_ptr<Matrix> xB(new Matrix("xB",nb,ns));
    double** xAp = xA->pointer(); 
    double** xBp = xB->pointer(); 

    boost::shared_ptr<Matrix> wB(new Matrix("wB",na,nr));
    boost::shared_ptr<Matrix> wA(new Matrix("wA",nb,ns));
    double** wBp = wB->pointer(); 
    double** wAp = wA->pointer(); 
    
    // ==> Total ESPs <== //

    boost::shared_ptr<Matrix> W_A(J_A->clone());
    W_A->copy(J_A);
    W_A->scale(2.0);
    W_A->add(V_A);

    boost::shared_ptr<Matrix> W_B(J_B->clone());
    W_B->copy(J_B);
    W_B->scale(2.0);
    W_B->add(V_B);

    boost::shared_ptr<Matrix> wAT = triplet(Locc_B_,W_A,Cvir_B_,true,false,false);
    boost::shared_ptr<Matrix> wBT = triplet(Locc_A_,W_B,Cvir_A_,true,false,false);
    double** wATp = wAT->pointer();
    double** wBTp = wBT->pointer();
    
    W_A.reset();
    W_B.reset();

    // ==> TODO: Total Exchange ESPs <== //

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
        boost::shared_ptr<Matrix> x2A = doublet(Uocc_A_,xA,true,false);
        double** x2Ap = x2A->pointer();

        // Zip up the Ind20 contributions
        for (int a = 0; a < na; a++) {
            double val = 2.0 * C_DDOT(nr,x2Ap[a],1,wBTp[a],1);
            IndAB_termsp[a][B] = val;
            Ind20_AB += val;
        }
    
        // TODO: ExchInd20

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
        boost::shared_ptr<Matrix> x2B = doublet(Uocc_B_,xB,true,false);
        double** x2Bp = x2B->pointer();

        // Zip up the Ind20 contributions
        for (int b = 0; b < nb; b++) {
            double val = 2.0 * C_DDOT(ns,x2Bp[b],1,wATp[b],1);
            IndBA_termsp[b][A] = val;
            Ind20_BA += val;
        }
    
        // TODO: ExchInd20

    } 

    // ==> LO -> A transform <== //
   
    boost::shared_ptr<Matrix> IndAB_atoms = doublet(Q_A_,IndAB_terms,false,false);
    boost::shared_ptr<Matrix> IndBA_atoms = doublet(IndBA_terms,Q_B_,true,true);

    IndAB_atoms->set_name("Ind [B<-A] (A x B)");
    IndBA_atoms->set_name("Ind [A<-B] (A x B)");

    local_vars_["IndAB"] = IndAB_atoms;
    local_vars_["IndBA"] = IndBA_atoms;

    double Ind20 = Ind20_AB + Ind20_BA;
    energies_["Ind20 (A<-B)"] = Ind20_AB;
    energies_["Ind20 (B->A)"] = Ind20_BA;
    energies_["Ind20"] = Ind20;
    fprintf(outfile,"    Ind20 (A<-B)       = %18.12lf H\n",Ind20_AB);
    fprintf(outfile,"    Ind20 (A->B)       = %18.12lf H\n",Ind20_BA);
    fprintf(outfile,"    Ind20              = %18.12lf H\n",Ind20);
    fprintf(outfile,"\n");
    fflush(outfile);
}
void ASAPT::disp()
{
    fprintf(outfile, "  DISPERSION:\n\n");

    // => Sizing <= //

    int nn = primary_->nbf();

    int na = Caocc_A_->colspi()[0];
    int nb = Caocc_B_->colspi()[0];
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
    size_t naQ = na * (size_t) nQ;
    size_t nbQ = nb * (size_t) nQ;

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

    boost::shared_ptr<Matrix> Cr1 = triplet(D_B,S,Cavir_A_);
    Cr1->scale(-1.0);
    Cr1->add(Cavir_A_);
    boost::shared_ptr<Matrix> Cs1 = triplet(D_A,S,Cavir_B_);
    Cs1->scale(-1.0);
    Cs1->add(Cavir_B_);
    boost::shared_ptr<Matrix> Ca2 = triplet(D_B,S,Caocc_A_);
    boost::shared_ptr<Matrix> Cb2 = triplet(D_A,S,Caocc_B_);
    boost::shared_ptr<Matrix> Cr3 = triplet(D_B,S,Cavir_A_);
    boost::shared_ptr<Matrix> CrX = triplet(triplet(D_A,S,D_B),S,Cavir_A_);
    Cr3->subtract(CrX);
    Cr3->scale(2.0);
    boost::shared_ptr<Matrix> Cs3 = triplet(D_A,S,Cavir_B_);
    boost::shared_ptr<Matrix> CsX = triplet(triplet(D_B,S,D_A),S,Cavir_B_);
    Cs3->subtract(CsX);
    Cs3->scale(2.0);
    boost::shared_ptr<Matrix> Ca4 = triplet(triplet(D_A,S,D_B),S,Caocc_A_);
    Ca4->scale(-2.0);
    boost::shared_ptr<Matrix> Cb4 = triplet(triplet(D_B,S,D_A),S,Caocc_B_);
    Cb4->scale(-2.0);

    // => Auxiliary V matrices <= //

    boost::shared_ptr<Matrix> Jbr = triplet(Caocc_B_,J_A,Cavir_A_,true,false,false);
    Jbr->scale(2.0);
    boost::shared_ptr<Matrix> Kbr = triplet(Caocc_B_,K_A,Cavir_A_,true,false,false);
    Kbr->scale(-1.0);

    boost::shared_ptr<Matrix> Jas = triplet(Caocc_A_,J_B,Cavir_B_,true,false,false);
    Jas->scale(2.0);
    boost::shared_ptr<Matrix> Kas = triplet(Caocc_A_,K_B,Cavir_B_,true,false,false);
    Kas->scale(-1.0);

    boost::shared_ptr<Matrix> KOas = triplet(Caocc_A_,K_O,Cavir_B_,true,false,false);
    KOas->scale(1.0);
    boost::shared_ptr<Matrix> KObr = triplet(Caocc_B_,K_O,Cavir_A_,true,true,false);
    KObr->scale(1.0);

    boost::shared_ptr<Matrix> JBas = triplet(triplet(Caocc_A_,S,D_B,true,false,false),J_A,Cavir_B_);
    JBas->scale(-2.0);
    boost::shared_ptr<Matrix> JAbr = triplet(triplet(Caocc_B_,S,D_A,true,false,false),J_B,Cavir_A_);
    JAbr->scale(-2.0);

    boost::shared_ptr<Matrix> Jbs = triplet(Caocc_B_,J_A,Cavir_B_,true,false,false);
    Jbs->scale(4.0);
    boost::shared_ptr<Matrix> Jar = triplet(Caocc_A_,J_B,Cavir_A_,true,false,false);
    Jar->scale(4.0);

    boost::shared_ptr<Matrix> JAas = triplet(triplet(Caocc_A_,J_B,D_A,true,false,false),S,Cavir_B_);
    JAas->scale(-2.0);
    boost::shared_ptr<Matrix> JBbr = triplet(triplet(Caocc_B_,J_A,D_B,true,false,false),S,Cavir_A_);
    JBbr->scale(-2.0);

    // Get your signs right Hesselmann!
    boost::shared_ptr<Matrix> Vbs = triplet(Caocc_B_,V_A,Cavir_B_,true,false,false);
    Vbs->scale(2.0);
    boost::shared_ptr<Matrix> Var = triplet(Caocc_A_,V_B,Cavir_A_,true,false,false);
    Var->scale(2.0);
    boost::shared_ptr<Matrix> VBas = triplet(triplet(Caocc_A_,S,D_B,true,false,false),V_A,Cavir_B_);
    VBas->scale(-1.0);
    boost::shared_ptr<Matrix> VAbr = triplet(triplet(Caocc_B_,S,D_A,true,false,false),V_B,Cavir_A_);
    VAbr->scale(-1.0);
    boost::shared_ptr<Matrix> VRas = triplet(triplet(Caocc_A_,V_B,P_A,true,false,false),S,Cavir_B_);
    VRas->scale(1.0);
    boost::shared_ptr<Matrix> VSbr = triplet(triplet(Caocc_B_,V_A,P_B,true,false,false),S,Cavir_A_);
    VSbr->scale(1.0);

    boost::shared_ptr<Matrix> Sas = triplet(Caocc_A_,S,Cavir_B_,true,false,false);
    boost::shared_ptr<Matrix> Sbr = triplet(Caocc_B_,S,Cavir_A_,true,false,false);

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

    boost::shared_ptr<Matrix> SBar = triplet(triplet(Caocc_A_,S,D_B,true,false,false),S,Cavir_A_);
    boost::shared_ptr<Matrix> SAbs = triplet(triplet(Caocc_B_,S,D_A,true,false,false),S,Cavir_B_);

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

    long int cost_r = 2L * na * nQ + 2L * nb * nQ;
    long int max_r = rem / (2L * cost_r);
    long int max_s = max_r;
    max_r = (max_r > nr ? nr : max_r);
    max_s = (max_s > ns ? ns : max_s);
    if (max_s < 1L || max_s < 1L) {
        throw PSIEXCEPTION("Too little dynamic memory for DFTSAPT::mp2_terms");
    }

    // => Tensor Slices <= //

    boost::shared_ptr<Matrix> Aar(new Matrix("Aar",max_r*na,nQ));
    boost::shared_ptr<Matrix> Abs(new Matrix("Abs",max_s*nb,nQ));
    boost::shared_ptr<Matrix> Bas(new Matrix("Bas",max_s*na,nQ));
    boost::shared_ptr<Matrix> Bbr(new Matrix("Bbr",max_r*nb,nQ));
    boost::shared_ptr<Matrix> Cas(new Matrix("Cas",max_s*na,nQ));
    boost::shared_ptr<Matrix> Cbr(new Matrix("Cbr",max_r*nb,nQ));
    boost::shared_ptr<Matrix> Dar(new Matrix("Dar",max_r*na,nQ));
    boost::shared_ptr<Matrix> Dbs(new Matrix("Dbs",max_s*nb,nQ));

    // => Thread Work Arrays <= //

    std::vector<boost::shared_ptr<Matrix> > Tab;
    std::vector<boost::shared_ptr<Matrix> > Vab;
    std::vector<boost::shared_ptr<Matrix> > T2ab;
    std::vector<boost::shared_ptr<Matrix> > V2ab;
    std::vector<boost::shared_ptr<Matrix> > Iab;
    for (int t = 0; t < nT; t++) {
        Tab.push_back(boost::shared_ptr<Matrix>(new Matrix("Tab",na,nb)));
        Vab.push_back(boost::shared_ptr<Matrix>(new Matrix("Vab",na,nb)));
        T2ab.push_back(boost::shared_ptr<Matrix>(new Matrix("T2ab",na,nb)));
        V2ab.push_back(boost::shared_ptr<Matrix>(new Matrix("V2ab",na,nb)));
        Iab.push_back(boost::shared_ptr<Matrix>(new Matrix("Iab",na,nb)));
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

    boost::shared_ptr<Matrix> E_disp20(new Matrix("E_disp20", na, nb));
    boost::shared_ptr<Matrix> E_exch_disp20(new Matrix("E_exch_disp20", na, nb));
    double** E_disp20p = E_disp20->pointer();
    double** E_exch_disp20p = E_exch_disp20->pointer();

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

                double** Tabp  = Tab[thread]->pointer();
                double** Vabp  = Vab[thread]->pointer();
                double** T2abp = T2ab[thread]->pointer();
                double** V2abp = V2ab[thread]->pointer();
                double** Iabp  = Iab[thread]->pointer();

                // => Amplitudes, Disp20 <= //

                C_DGEMM('N','T',na,nb,nQ,1.0,Aarp[(r)*na],nQ,Absp[(s)*nb],nQ,0.0,Vabp[0],nb);
                for (int a = 0; a < na; a++) {
                    for (int b = 0; b < nb; b++) {
                        Tabp[a][b] = Vabp[a][b] / (eap[a] + ebp[b] - erp[r + rstart] - esp[s + sstart]);
                    }
                }

                C_DGEMM('N','N',na,nb,nb,1.0,Tabp[0],nb,UBp[0],nb,0.0,Iabp[0],nb);
                C_DGEMM('T','N',na,nb,na,1.0,UAp[0],na,Iabp[0],nb,0.0,T2abp[0],nb);
                C_DGEMM('N','N',na,nb,nb,1.0,Vabp[0],nb,UBp[0],nb,0.0,Iabp[0],nb);
                C_DGEMM('T','N',na,nb,na,1.0,UAp[0],na,Iabp[0],nb,0.0,V2abp[0],nb);

                for (int a = 0; a < na; a++) {
                    for (int b = 0; b < nb; b++) {
                        E_disp20p[a][b] += 4.0 * T2abp[a][b] * V2abp[a][b];
                        Disp20 += 4.0 * T2abp[a][b] * V2abp[a][b];
                    }
                }

                // => Exch-Disp20 <= //

                // > Q1-Q3 < //

                C_DGEMM('N','T',na,nb,nQ,1.0,Basp[(s)*na],nQ,Bbrp[(r)*nb],nQ,0.0,Vabp[0],nb);
                C_DGEMM('N','T',na,nb,nQ,1.0,Casp[(s)*na],nQ,Cbrp[(r)*nb],nQ,1.0,Vabp[0],nb);
                C_DGEMM('N','T',na,nb,nQ,1.0,Aarp[(r)*na],nQ,Dbsp[(s)*nb],nQ,1.0,Vabp[0],nb);
                C_DGEMM('N','T',na,nb,nQ,1.0,Darp[(r)*na],nQ,Absp[(s)*nb],nQ,1.0,Vabp[0],nb);

                // > V,J,K < //

                C_DGER(na,nb,1.0,&Sasp[0][s + sstart], ns,&Qbrp[0][r + rstart], nr,Vabp[0],nb);
                C_DGER(na,nb,1.0,&Qasp[0][s + sstart], ns,&Sbrp[0][r + rstart], nr,Vabp[0],nb);
                C_DGER(na,nb,1.0,&Qarp[0][r + rstart], nr,&SAbsp[0][s + sstart],ns,Vabp[0],nb);
                C_DGER(na,nb,1.0,&SBarp[0][r + rstart],nr,&Qbsp[0][s + sstart], ns,Vabp[0],nb);

                C_DGEMM('N','N',na,nb,nb,1.0,Vabp[0],nb,UBp[0],nb,0.0,Iabp[0],nb);
                C_DGEMM('T','N',na,nb,na,1.0,UAp[0],na,Iabp[0],nb,0.0,V2abp[0],nb);

                for (int a = 0; a < na; a++) {
                    for (int b = 0; b < nb; b++) {
                        E_exch_disp20p[a][b] -= 2.0 * T2abp[a][b] * V2abp[a][b];
                        ExchDisp20 -= 2.0 * T2abp[a][b] * V2abp[a][b];
                    }
                }
            }
        }
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

    boost::shared_ptr<Matrix> Disp_atoms = triplet(Q2A,E_disp,Q2B,false,false,true);
    Disp_atoms->set_name("Disp (A x B)");

    //Disp_atoms->print();
    local_vars_["Disp"] = Disp_atoms;

    energies_["Disp20"] = Disp20;
    energies_["Exch-Disp20"] = ExchDisp20;
    fprintf(outfile,"    Disp20              = %18.12lf H\n",Disp20);
    fprintf(outfile,"    Exch-Disp20         = %18.12lf H\n",ExchDisp20);
    fprintf(outfile,"\n");
    fflush(outfile);
}
void ASAPT::analyze()
{
    std::vector<std::string> keys;
    keys.push_back("Elst"); 
    keys.push_back("Exch"); 
    keys.push_back("IndAB"); 
    keys.push_back("IndBA"); 
    keys.push_back("Disp"); 

    boost::shared_ptr<Matrix> Tot(local_vars_["Elst"]->clone());
    Tot->zero();
    Tot->set_name("Total SAPT0 (A x B)");
    for (int i = 0; i < keys.size(); i++) {
        Tot->add(local_vars_[keys[i]]);
    }
    local_vars_["Total"] = Tot;

    std::map<std::string, boost::shared_ptr<Matrix> > mA;
    std::map<std::string, boost::shared_ptr<Matrix> > mB;

    for (int i = 0; i < keys.size(); i++) {
        std::string key = keys[i];
        boost::shared_ptr<Matrix> T = local_vars_[key];
        int nA = T->nrow();
        int nB = T->ncol();
        boost::shared_ptr<Matrix> TA(new Matrix(key,nA,1));
        boost::shared_ptr<Matrix> TB(new Matrix(key,nB,1));
        double** Tp = T->pointer();
        double** TAp = TA->pointer();
        double** TBp = TB->pointer();
        for (int A = 0; A < nA; A++) {
            for (int B = 0; B < nB; B++) {
                TAp[A][0] += Tp[A][B];
                TBp[B][0] += Tp[A][B];
            }
        }
        mA[key] = TA;
        mB[key] = TB;
    } 
    
    /**
    fprintf(outfile, "  L-SAPT0 Analysis [H]:\n\n");

    fprintf(outfile, "   Diatomic Partition:\n\n");
    for (int i = 0; i < keys.size(); i++) {
        local_vars_[keys[i]]->print();
    }
     
    fprintf(outfile, "   Monomer A Partition:\n\n");
    for (int i = 0; i < keys.size(); i++) {
        mA[keys[i]]->print();
    }

    fprintf(outfile, "   Monomer B Partition:\n\n");
    for (int i = 0; i < keys.size(); i++) {
        mB[keys[i]]->print();
    }
    **/

    monomer_A_->save_xyz("mA.xyz",false);
    monomer_B_->save_xyz("mB.xyz",false);

    for (int i = 0; i < keys.size(); i++) {
        std::string key = keys[i];
        boost::shared_ptr<Matrix> TA = mA[key]; 
        boost::shared_ptr<Matrix> TB = mB[key]; 
        double** TAp = TA->pointer();
        double** TBp = TB->pointer();
        int nA = TA->nrow();
        int nB = TB->nrow();

        std::stringstream ssA;
        ssA << key << "_A.dat";
        FILE* fhA = fopen(ssA.str().c_str(),"w");
        for (int A = 0; A < nA; A++) {
            fprintf(fhA,"%24.16E\n", TAp[A][0]);
        }
        fclose(fhA);

        std::stringstream ssB;
        ssB << key << "_B.dat";
        FILE* fhB = fopen(ssB.str().c_str(),"w");
        for (int B = 0; B < nB; B++) {
            fprintf(fhB,"%24.16E\n", TBp[B][0]);
        }
        fclose(fhB);
    }
}

}}
