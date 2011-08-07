
/*
 *  dfmp2.cc
 *  
 *
 *  Created by M.Saitow on 11/07/21.
 *  Copyright 2011 M.Saitow. All rights reserved.
 *
 */

#include "psi4-dec.h"
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libmints/wavefunction.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.hpp>
#include <libplugin/plugin.h>
#include <liboptions/liboptions.h>
#include <libqt/qt.h>

#define TIME_DF_MP2 1

INIT_PLUGIN

using namespace boost;

namespace psi{ namespace dfmp2{
	
extern void formInvSqrtJ(double **&J_mhalf, shared_ptr<BasisSet> basis,
                         shared_ptr<BasisSet> ribasis, shared_ptr<BasisSet> zero);
	
extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "DF-MP2") {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        /*- Whether to compute the SCS energy -*/
        options.add_bool("DO_SCS", true);
        /*- Whether to compute the SCS-N energy -*/
	options.add_bool("DO_SCS-N", true);
	/*- The name of the orbital basis set -*/
        options.add_str("BASIS", "");
	/*- The name of the auxilliary basis set -*/
	options.add_str("RI_BASIS", "");
	/*- The opposite-spin scale factor for the SCS energy -*/
	options.add_double("SCALE_OS", 6.0/5.0);
	/*- The same-spin scale factor for the SCS energy -*/
        options.add_double("SCALE_SS", 1.0/3.0);
    }
		
    return true;
}
	
	
extern "C" PsiReturnType
df-mp2(Options &options)
{ 
  shared_ptr<Molecule> molecule = Process::environment.molecule();
  shared_ptr<Wavefunction> wfn = Process::environment.reference_wavefunction();
  if (wfn == NULL)
    throw PSIEXCEPTION("Could not find MO coefficients. Run scf first.");
		
    shared_ptr<MatrixFactory> factory(new MatrixFactory);
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser);
		
    std::string ri_basis = options.get_str("RI_BASIS");
    if(ri_basis == ""){
			throw InputException(std::string("Keyword not specified"),
													 std::string("RI_BASIS"), __FILE__, __LINE__);
    }
    shared_ptr<BasisSet> ribasis = BasisSet::construct(parser, molecule, ri_basis);
		
    std::string orbital_basis = options.get_str("BASIS");
    if(orbital_basis == ""){
			throw InputException(std::string("Keyword not specified"),
													 std::string("BASIS"), __FILE__, __LINE__);
    }
    fprintf(outfile, "Using the %s basis set for the orbitals, with the %s RI basis\n\n",
						orbital_basis.c_str(), ri_basis.c_str());
    shared_ptr<BasisSet> basis = BasisSet::construct(parser, molecule, orbital_basis);
		
    int nbf[] = { basis->nbf() };
    factory->init_with(1, nbf, nbf);
		
    int print = options.get_int("PRINT");
    if(print > 5) ribasis->print();
		
    shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());
		
    bool doSCSN = options.get_bool("DO_SCS-N");
    bool doSCS  = options.get_bool("DO_SCS");
    double SCSScaleOS = options.get_double("SCALE_OS");
    double SCSScaleSS = options.get_double("SCALE_SS");
    if(doSCS){
			fprintf(outfile,"\tSpin-Component Scaled RI-MP2 requested\n"
							"\tOpposite-spin scaled by %10.4lf\n"
							"\tSame-spin scaled by     %10.4lf\n", SCSScaleOS, SCSScaleSS);
    }
    // The integrals code below does not use symmetry, so we need to accumulate the
    // orbital info for each irrep
    int nirreps = Process::environment.reference_wavefunction()->nirrep();
    int *clsdpi = wfn->doccpi();
    int *orbspi = wfn->nmopi();
    int *frzcpi = wfn->frzcpi();
    int *frzvpi = wfn->frzvpi();
    int ndocc = 0;
    int nvirt = 0;
    int nfocc = 0;
    int nfvir = 0;
    int norbs = 0;
    int nact_docc = 0;
    int nact_virt = 0;
    for(int h=0; h < nirreps; ++h){
			nfocc     += frzcpi[h];
			nfvir     += frzvpi[h];
			ndocc     += clsdpi[h];
			nact_docc += clsdpi[h] - frzcpi[h];
			nvirt     += orbspi[h] - clsdpi[h];
			nact_virt += orbspi[h] - frzvpi[h] - clsdpi[h];
			norbs     += orbspi[h];
    }
		
    fprintf(outfile, "\n\t\t==============================================\n");
    fprintf(outfile, "\t\t #ORBITALS #RI  FOCC DOCC AOCC AVIR VIRT FVIR \n");
    fprintf(outfile, "\t\t----------------------------------------------\n");
    fprintf(outfile, "\t\t  %5d  %5d  %4d %4d %4d %4d %4d %4d\n",
            norbs,ribasis->nbf(),nfocc,ndocc,nact_docc,nact_virt,nvirt,nfvir);
    fprintf(outfile, "\t\t==============================================\n");
		
    // Read in MO coefficients
    shared_ptr<SimpleMatrix> C_so(factory->
																	create_simple_matrix("MO coefficients (SO basis)"));
		
    shared_ptr<Vector> orbital_energies = wfn->epsilon_a();
    shared_ptr<Matrix> C = wfn->Ca();
		
    double** Co   = block_matrix(norbs, nact_docc);
    double** Cv   = block_matrix(norbs, nact_virt);
    double** half = block_matrix(nact_docc, norbs);
    double* epsilon_act_docc = new double[nact_docc];
    double* epsilon_act_virt = new double[nact_virt];
    int*    docc_sym = new int[nact_docc];
    int*    virt_sym = new int[nact_virt];
    int offset = 0;
    int act_docc_count  = 0;
    int act_virt_count  = 0;
    for(int h=0; h<nirreps; ++h){
			// Skip over the frozen core orbitals in this irrep
			offset += frzcpi[h];
			// Copy over the info for active occupied orbitals
			for(int i=0; i<clsdpi[h]-frzcpi[h]; ++i){
				for (int mu=0; mu<norbs; ++mu){
					Co[mu][act_docc_count] = C->get(0, mu, offset);
				}
				epsilon_act_docc[act_docc_count] = orbital_energies->get(0, offset);
				docc_sym[act_docc_count] = h;
				++act_docc_count;
				++offset;
			}
			// Copy over the info for active virtual orbitals
			for(int a=0; a<orbspi[h]-clsdpi[h]-frzvpi[h]; ++a){
				for (int mu=0; mu<norbs; ++mu){
					Cv[mu][act_virt_count] = C->get(0, mu, offset);
				}
				epsilon_act_virt[act_virt_count] = orbital_energies->get(0, offset);
				virt_sym[act_virt_count] = h;
				++offset;
				++act_virt_count;
			}
			// Skip over the frozen virtual orbitals in this irrep
			offset += frzvpi[h];
    }
		
    if(print > 5){
			fprintf(outfile, "Co:\n");
			print_mat(Co, norbs, nact_docc, outfile);
			fprintf(outfile, "Cv:\n");
			print_mat(Cv, norbs, nact_virt, outfile);
    }
		
    double **J_mhalf;
    formInvSqrtJ(J_mhalf, basis, ribasis, zero);
		
    double **mo_p_ia = block_matrix(ribasis->nbf(),nact_docc*nact_virt);
		
    // find out the max number of P's in a P shell
    int maxPshell = 0;
    for (int Pshell=0; Pshell < ribasis->nshell(); ++Pshell) {
			int numPshell = ribasis->shell(Pshell)->nfunction();
			maxPshell = numPshell > maxPshell ? numPshell : maxPshell;
    }
    double*** temp = new double**[maxPshell];
    for (int P=0; P<maxPshell; P++) temp[P] = block_matrix(norbs, norbs);
		
#ifdef TIME_DF_MP2
    timer_on("Form mo_p_ia");
#endif
		
    shared_ptr<IntegralFactory>
		rifactory(new IntegralFactory(ribasis, zero, basis, basis));
    shared_ptr<TwoBodyAOInt> eri(rifactory->eri());
    const double *buffer = eri->buffer();
		
    for(int Pshell=0; Pshell < ribasis->nshell(); ++Pshell){
			int numPshell = ribasis->shell(Pshell)->nfunction();
			for(int P=0; P<numPshell; ++P){
				zero_mat(temp[P], norbs, norbs);
			}
			for(int MU=0; MU < basis->nshell(); ++MU){
				int nummu = basis->shell(MU)->nfunction();
				for(int NU=0; NU <= MU; ++NU){
					int numnu = basis->shell(NU)->nfunction();
					eri->compute_shell(Pshell, 0, MU, NU);
					for(int P=0, index=0; P < numPshell; ++P){
						for(int mu=0; mu < nummu; ++mu) {
							int omu = basis->shell(MU)->function_index() + mu;
							for(int nu=0; nu < numnu; ++nu, ++index){
								int onu = basis->shell(NU)->function_index() + nu;
								// (oP | omu onu) integral
								temp[P][omu][onu] = buffer[index];
								// (oP | onu omu) integral
								temp[P][onu][omu] = buffer[index];
							}
						}
					} // end loop over P in Pshell
				} // end loop over NU shell
			} // end loop over MU shell
			// now we've gone through all P, mu, nu for a given Pshell
			// transform the integrals for all P in the given P shell
			for(int P=0, index=0; P < numPshell; ++P){
				int oP = ribasis->shell(Pshell)->function_index() + P;
				// Do transform
				C_DGEMM('T', 'N', nact_docc, norbs, norbs, 1.0, Co[0],
								nact_docc, temp[P][0], norbs, 0.0, half[0], norbs);
				C_DGEMM('N', 'N', nact_docc, nact_virt, norbs, 1.0, half[0],
								norbs, Cv[0], nact_virt, 0.0, mo_p_ia[oP], nact_virt);
			}
    } // end loop over P shells; done with forming MO basis (P|ia)'s
		
    for(int P=0; P<maxPshell; P++) free_block(temp[P]);
		
#ifdef TIME_DF_MP2
    timer_off("Form mo_p_ia");
    timer_on("Form B_ia^P");
#endif
		
    // mo_p_ia has integrals
    // B_ia^P = Sum_Q (i a | Q) (J^-1/2)_QP
    double **B_ia_p = block_matrix(nact_docc * nact_virt, ribasis->nbf());
		
    C_DGEMM('T','N',nact_docc*nact_virt,ribasis->nbf(),ribasis->nbf(),
						1.0, mo_p_ia[0], nact_docc*nact_virt, J_mhalf[0], ribasis->nbf(),
						0.0, B_ia_p[0], ribasis->nbf());
		
    free_block(mo_p_ia);
    free_block(J_mhalf);
		
#ifdef TIME_DF_MP2
    timer_off("Form B_ia^P");
    timer_on("Compute EMP2");
#endif
		
    double *I = init_array(nact_virt * nact_virt);
    double emp2 = 0.0;
    double os_mp2 = 0.0, ss_mp2 = 0.0;
		
    // loop over i>=j pairs
    for(int i=0; i < nact_docc; ++i){
			for(int j=0; j <= i; ++j){
				int ijsym = docc_sym[i] ^ docc_sym[j];
				// get the integrals for this pair of occupied orbitals
#ifdef TIME_DF_MP2
				timer_on("Construct I");
#endif
				int ia, jb;
				for(int a=0,ab=0; a < nact_virt; ++a){
					ia = nact_virt*i + a;
					for(int b=0; b < nact_virt; ++b, ++ab){
						jb = nact_virt*j + b;
						int absym = virt_sym[a] ^ virt_sym[b];
						if (ijsym == absym)
							I[ab] = C_DDOT(ribasis->nbf(), B_ia_p[ia], 1, B_ia_p[jb], 1);
						else
							I[ab] = 0.0;
						// I[ab] = (ia|jb) for the fixed i,j
						// note (I[ab] = (ia|jb)) != (I[ba] = (ib|ja))
					}
				} // end loop over a
#ifdef TIME_DF_MP2
				timer_off("Construct I");
#endif
				double iajb, ibja, tval, denom;
				int ab, ba;
				for(int a=0,ab=0; a < nact_virt; ++a){
					for(int b=0; b < nact_virt; ++b,ab++){
						int ba = b * nact_virt + a;
						iajb = I[ab];
						ibja = I[ba];
						denom = 1.0 /
						(epsilon_act_docc[i] + epsilon_act_docc[j] -
						 epsilon_act_virt[a] - epsilon_act_virt[b]);
						
						tval = ((i==j) ? 0.5 : 1.0) * (iajb * iajb + ibja * ibja) * denom;
						os_mp2 += tval;
						ss_mp2 += tval - ((i==j) ? 1.0 : 2.0)*(iajb*ibja)*denom;
					}
				}
			} // end loop over j<=i
    } // end loop over i
    emp2 = os_mp2 + ss_mp2;
		
#ifdef TIME_DF_MP2
    timer_off("Compute EMP2");
    timer_done();
#endif
		
    free(I);
    free_block(B_ia_p);
    
    double escf = Process::environment.globals["CURRENT ENERGY"];
    fprintf(outfile,"\tRI-MP2 correlation energy         = %20.15f\n",emp2);
    fprintf(outfile,"      * RI-MP2 total energy               = %20.15f\n\n",
						escf + emp2);
    fprintf(outfile,"\tOpposite-Spin correlation energy  = %20.15f\n",os_mp2);
    fprintf(outfile,"\tSame-Spin correlation energy      = %20.15f\n\n",ss_mp2);
    fprintf(outfile,"      * SCS-RI-MP2 total energy           = %20.15f\n\n",
						escf + SCSScaleOS * os_mp2 + SCSScaleSS * ss_mp2);
    if(doSCSN){
			fprintf(outfile,"      * SCSN-RI-MP2 total energy          = %20.15f\n\n",
							escf + 1.76*ss_mp2);
    }
		
    return Success;
	}
	
}} // end namespaces
