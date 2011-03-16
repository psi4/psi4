#include "dcft.h"
#include <libdpd/dpd.h>
#include <libtrans/integraltransform.h>
#include <libdiis/diismanager.h>
#include "defines.h"

using namespace boost;

namespace psi{ namespace dcft{

/**
 * Computes the DCFT density matrix and energy
 */
void
DCFTSolver::compute()
{
    bool scfDone    = false;
    bool lambdaDone = false;
    bool densityConverged = false;
    bool energyConverged = false;
    double oldEnergy;
    scf_guess();
    mp2_guess();

    int cycle = 0;
    fprintf(outfile, "\n\n\t*=================================================================================*\n"
                     "\t* Cycle  RMS [F, Kappa]   RMS Lambda Error   delta E        Total Energy     DIIS *\n"
                     "\t*---------------------------------------------------------------------------------*\n");
    if(_options.get_str("ALGORITHM") == "TWOSTEP"){
        // This is the two-step update - in each macro iteration, update the orbitals first, then update lambda
        // to self-consistency, until converged.  When lambda is converged and only one scf cycle is needed to reach
        // the desired cutoff, we're done
        SharedMatrix tmp = shared_ptr<Matrix>(new Matrix("temp", _nIrreps, _soPI, _soPI));
        // Set up the DIIS manager
        dpdbuf4 Laa, Lab, Lbb;
        dpd_buf4_init(&Laa, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
        dpd_buf4_init(&Lab, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
        dpd_buf4_init(&Lbb, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
        DIISManager scfDiisManager(_maxDiis, "DCFT DIIS Orbitals",DIISManager::LargestError,DIISManager::InCore);
        scfDiisManager.set_error_vector_size(2, DIISEntry::Matrix, _aScfError.get(),
                                                DIISEntry::Matrix, _bScfError.get());
        scfDiisManager.set_vector_size(2, DIISEntry::Matrix, _Fa.get(),
                                          DIISEntry::Matrix, _Fb.get());
        DIISManager lambdaDiisManager(_maxDiis, "DCFT DIIS Lambdas",DIISManager::LargestError,DIISManager::InCore);
        lambdaDiisManager.set_error_vector_size(3, DIISEntry::DPDBuf4, &Laa,
                                                   DIISEntry::DPDBuf4, &Lab,
                                                   DIISEntry::DPDBuf4, &Lbb);
        lambdaDiisManager.set_vector_size(3, DIISEntry::DPDBuf4, &Laa,
                                             DIISEntry::DPDBuf4, &Lab,
                                             DIISEntry::DPDBuf4, &Lbb);
        dpd_buf4_close(&Laa);
        dpd_buf4_close(&Lab);
        dpd_buf4_close(&Lbb);
        _oldCa->copy(_Ca);
        _oldCb->copy(_Cb);
        // The macro-iterations
        while((!scfDone || !lambdaDone) && cycle++ < _maxNumIterations){
            // The lambda iterations
            int nLambdaIterations = 0;
            lambdaDiisManager.reset_subspace();
            fprintf(outfile, "\t                          *** Macro Iteration %d ***\n"
                             "\tCumulant Iterations\n",cycle);
            lambdaDone = false;
            while(!lambdaDone && nLambdaIterations++ < _options.get_int("LAMBDA_MAXITER")){
                std::string diisString;
                build_tensors();
                build_intermediates();
                _lambdaConvergence = compute_lambda_residual();
                if(_lambdaConvergence < _diisStartThresh){
                    //Store the DIIS vectors
                    dpdbuf4 Laa, Lab, Lbb, Raa, Rab, Rbb;
                    dpd_buf4_init(&Raa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                  ID("[O,O]"), ID("[V,V]"), 0, "R <OO|VV>");
                    dpd_buf4_init(&Rab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                                  ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
                    dpd_buf4_init(&Rbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                                  ID("[o,o]"), ID("[v,v]"), 0, "R <oo|vv>");
                    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
                    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
                    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
                    if(lambdaDiisManager.add_entry(6, &Raa, &Rab, &Rbb, &Laa, &Lab, &Lbb)){
                        diisString += "S";
                    }
                    if(lambdaDiisManager.subspace_size() == _maxDiis && _maxDiis > 0){
                        diisString += "/E";
                        lambdaDiisManager.extrapolate(3, &Laa, &Lab, &Lbb);
                        lambdaDiisManager.reset_subspace();
                    }else{
                        update_lambda_from_residual();
                    }
                    dpd_buf4_close(&Raa);
                    dpd_buf4_close(&Rab);
                    dpd_buf4_close(&Rbb);
                    dpd_buf4_close(&Laa);
                    dpd_buf4_close(&Lab);
                    dpd_buf4_close(&Lbb);
                }else{
                    update_lambda_from_residual();
                }
                oldEnergy = _newTotalEnergy;
                compute_energy();
                lambdaDone = _lambdaConvergence < _lambdaThreshold &&
                                fabs(_newTotalEnergy - oldEnergy) < _lambdaThreshold;
                fprintf(outfile, "\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n",
                        nLambdaIterations, _scfConvergence, _lambdaConvergence, _newTotalEnergy - oldEnergy,
                        _newTotalEnergy, diisString.c_str());
                fflush(outfile);
            }
            build_tau();
            // Update the orbitals
            int nSCFCycles = 0;
            densityConverged = false;
            scfDiisManager.reset_subspace();
            fprintf(outfile, "\tOrbital Updates\n");
            while((!densityConverged || _scfConvergence > _scfThreshold)
                    && nSCFCycles++ < _options.get_int("SCF_MAXITER")){
                std::string diisString;
                _Fa->copy(_soH);
                _Fb->copy(_soH);
                // This will build the new Fock matrix from the SO integrals
                process_so_ints();
                // The SCF energy has to be evaluated before adding Tau and orthonormalizing F
                oldEnergy = _newTotalEnergy;
                compute_scf_energy();
                _Fa->add(_aGTau);
                _Fb->add(_bGTau);
                _scfConvergence = compute_scf_error_vector();
                if(_scfConvergence < _diisStartThresh){
                    if(scfDiisManager.add_entry(4, _aScfError.get(), _bScfError.get(), _Fa.get(), _Fb.get()))
                        diisString += "S";
                }
                if(scfDiisManager.subspace_size() > _minDiisVecs){
                    diisString += "/E";
                    scfDiisManager.extrapolate(2, _Fa.get(), _Fb.get());
                }
                _Fa->transform(_sHalfInv);
                _Fa->diagonalize(tmp, _aEvals);
                _oldCa->copy(_Ca);
                _Ca->gemm(false, false, 1.0, _sHalfInv, tmp, 0.0);
                _Fb->transform(_sHalfInv);
                _Fb->diagonalize(tmp, _bEvals);
                _oldCb->copy(_Cb);
                _Cb->gemm(false, false, 1.0, _sHalfInv, tmp, 0.0);
                correct_mo_phases(false);
                find_occupation(_aEvals);
                densityConverged = update_scf_density() < _scfThreshold;
                compute_energy();
                fprintf(outfile, "\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n",
                    nSCFCycles, _scfConvergence, _lambdaConvergence, _newTotalEnergy - oldEnergy,
                    _newTotalEnergy, diisString.c_str());
                fflush(outfile);
            }
            write_orbitals_to_checkpoint();
            scfDone = nSCFCycles == 1;
            transform_integrals();
        }
    }else{
        // This is the simultaneous orbital/lambda update algorithm
        SharedMatrix tmp = shared_ptr<Matrix>(new Matrix("temp", _nIrreps, _soPI, _soPI));
        // Set up the DIIS manager
        DIISManager diisManager(_maxDiis, "DCFT DIIS vectors");
        dpdbuf4 Laa, Lab, Lbb;
        dpd_buf4_init(&Laa, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
        dpd_buf4_init(&Lab, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
        dpd_buf4_init(&Lbb, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
        diisManager.set_error_vector_size(5, DIISEntry::Matrix, _aScfError.get(),
                                             DIISEntry::Matrix, _bScfError.get(),
                                             DIISEntry::DPDBuf4, &Laa,
                                             DIISEntry::DPDBuf4, &Lab,
                                             DIISEntry::DPDBuf4, &Lbb);
        diisManager.set_vector_size(5, DIISEntry::Matrix, _Fa.get(),
                                       DIISEntry::Matrix, _Fb.get(),
                                       DIISEntry::DPDBuf4, &Laa,
                                       DIISEntry::DPDBuf4, &Lab,
                                       DIISEntry::DPDBuf4, &Lbb);
        dpd_buf4_close(&Laa);
        dpd_buf4_close(&Lab);
        dpd_buf4_close(&Lbb);
        while((!scfDone || !lambdaDone || !densityConverged || !energyConverged)
                && cycle++ < _maxNumIterations){
            std::string diisString;
            oldEnergy = _newTotalEnergy;
            build_tau();
            _Fa->copy(_soH);
            _Fb->copy(_soH);
            // This will build the new Fock matrix from the SO integrals
            process_so_ints();
            // The SCF energy has to be evaluated before adding Tau and orthonormalizing F
            compute_scf_energy();
            _Fa->add(_aGTau);
            _Fb->add(_bGTau);
            _scfConvergence = compute_scf_error_vector();
            scfDone = _scfConvergence < _scfThreshold;
            build_intermediates();
            _lambdaConvergence = compute_lambda_residual();
            lambdaDone = _lambdaConvergence < _lambdaThreshold;
            update_lambda_from_residual();
            compute_energy();
            energyConverged = fabs(oldEnergy - _newTotalEnergy) < _lambdaConvergence;
            if(_scfConvergence < _diisStartThresh && _lambdaConvergence < _diisStartThresh){
                //Store the DIIS vectors
                dpdbuf4 Laa, Lab, Lbb, Raa, Rab, Rbb;
                dpd_buf4_init(&Raa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                              ID("[O,O]"), ID("[V,V]"), 0, "R <OO|VV>");
                dpd_buf4_init(&Rab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                              ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
                dpd_buf4_init(&Rbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                              ID("[o,o]"), ID("[v,v]"), 0, "R <oo|vv>");
                dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                              ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
                dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                              ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
                dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                              ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
                if(diisManager.add_entry(10, _aScfError.get(), _bScfError.get(), &Raa, &Rab, &Rbb,
                                           _Fa.get(), _Fb.get(), &Laa, &Lab, &Lbb)){
                    diisString += "S";
                }
                if(diisManager.subspace_size() > _minDiisVecs){
                    diisString += "/E";
                    diisManager.extrapolate(5, _Fa.get(), _Fb.get(), &Laa, &Lab, &Lbb);
                }
                dpd_buf4_close(&Raa);
                dpd_buf4_close(&Rab);
                dpd_buf4_close(&Rbb);
                dpd_buf4_close(&Laa);
                dpd_buf4_close(&Lab);
                dpd_buf4_close(&Lbb);
            }
            _Fa->transform(_sHalfInv);
            _Fa->diagonalize(tmp, _aEvals);
            _oldCa->copy(_Ca);
            _Ca->gemm(false, false, 1.0, _sHalfInv, tmp, 0.0);
            _Fb->transform(_sHalfInv);
            _Fb->diagonalize(tmp, _bEvals);
            _oldCb->copy(_Cb);
            _Cb->gemm(false, false, 1.0, _sHalfInv, tmp, 0.0);
            if(!correct_mo_phases(false)){
                fprintf(outfile,"\t\tThere was a problem correcting the MO phases.\n"
                                "\t\tIf this does not converge, try ALGORITHM=TWOSTEP\n");
            }
            write_orbitals_to_checkpoint();
            transform_integrals();
            find_occupation(_aEvals);
            densityConverged = update_scf_density() < _scfThreshold;
            // If we've performed enough lambda updates since the last orbitals
            // update, reset the counter so another SCF update is performed
            fprintf(outfile, "\t* %-3d   %12.3e      %12.3e   %12.3e  %21.15f  %-3s *\n",
                    cycle, _scfConvergence, _lambdaConvergence, _newTotalEnergy - oldEnergy,
                    _newTotalEnergy, diisString.c_str());
            fflush(outfile);
        }
    }
    if(!scfDone || !lambdaDone || !densityConverged)
        throw ConvergenceError<int>("DCFT", _maxNumIterations, _lambdaThreshold,
                               _lambdaConvergence, __FILE__, __LINE__);

    Process::environment.globals["CURRENT ENERGY"] = _newTotalEnergy;
    Process::environment.globals["DCFT ENERGY"] = _newTotalEnergy;
    Process::environment.globals["DCFT SCF ENERGY"] = _scfEnergy;

    fprintf(outfile, "\t*=================================================================================*\n");
    fprintf(outfile, "\n\t*DCFT SCF Energy            = %20.15f\n", _scfEnergy);
    fprintf(outfile, "\t*DCFT Total Energy          = %20.15f\n", _newTotalEnergy);

    if(!_options.get_bool("RELAX_ORBITALS")){
        fprintf(outfile, "Warning!  The orbitals were not relaxed\n");
    }

    print_opdm();


    if(_options.get_bool("COMPUTE_TPDM")) dump_density();
    mulliken_charges();
    check_n_representability();

    if(!_options.get_bool("RELAX_ORBITALS") && _options.get_bool("IGNORE_TAU")){
        _psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
        /*
         * Comout the CEPA-0 correlation energy
         * E = 1/4 L_IJAB <IJ||AB>
         *        +L_IjAb <Ij|Ab>
         *    +1/4 L_ijab <ij||ab>
         */
        dpdbuf4 I, L;
        // Alpha - Alpha
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
        dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
        double eAA = 0.25 * dpd_buf4_dot(&L, &I);
        dpd_buf4_close(&I);
        dpd_buf4_close(&L);

        // Alpha - Beta
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
        dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
        double eAB = dpd_buf4_dot(&L, &I);
        dpd_buf4_close(&I);
        dpd_buf4_close(&L);

        // Beta - Beta
        dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
        dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
        double eBB = 0.25 * dpd_buf4_dot(&L, &I);
        dpd_buf4_close(&I);
        dpd_buf4_close(&L);
        fprintf(outfile, "\t!CEPA0 Total Energy         = %20.15f\n",
                _scfEnergy + eAA + eAB + eBB);
        _psio->close(PSIF_LIBTRANS_DPD, 1);
    }
}

}} // Namespaces

