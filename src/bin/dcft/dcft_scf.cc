#include <map>
#include "dcft.h"
#include "libiwl/iwl.hpp"
#include "libdpd/dpd.h"
#include "libqt/qt.h"
#include "libmints/matrix.h"

namespace psi{ namespace dcft{

/**
 * Performs the SCF procedure to update the orbitals, and transforms the integrals
 * to the updated basis
 */
void
DCFTSolver::perform_scf()
{
    scf_guess();
    bool converged = false;
    _nScfIterations = 0;
    Matrix tmp("temp", _nIrreps, _soPI, _soPI);
    fprintf(outfile, "\t*===================================================================*\n"
                     "\t* Cycle   RMS Density Change   delta E            SCF Energy        *\n"
                     "\t*-------------------------------------------------------------------*\n");
    while(!converged && _nScfIterations++ < _scfMaxIter){
        double oldEnergy = _scfEnergy;
        _Fa.copy(_soH);
        _Fb.copy(_soH);
        // This will build the new Fock matrix from the SO integrals
        process_so_ints();
        compute_scf_energy();
        // Having computed the energy using the standard Fock matrix, add the external
        // potential before diagonalizing
        _Fa.add(_aGTau);
        _Fb.add(_bGTau);
        _Fa.transform(_sHalfInv);
        _Fb.transform(_sHalfInv);
        _Fa.diagonalize(tmp, _aEvals);
        _oldCa.copy(_Ca);
        _Ca.gemm(false, false, 1.0, _sHalfInv, tmp, 0.0);
        _Fb.diagonalize(tmp, _bEvals);
        _oldCb.copy(_Cb);
        _Cb.gemm(false, false, 1.0, _sHalfInv, tmp, 0.0);
        find_occupation(_aEvals);
        double rmsDeltaP = update_scf_density();
        converged = rmsDeltaP < pow(10.0,-_options.get_int("SCF_CONV"));
        fprintf(outfile, "\t*  %-3d      %10.2e       %10.2e      %20.16f  *\n",
                _nScfIterations, rmsDeltaP, _scfEnergy - oldEnergy, _scfEnergy);
    }
    fprintf(outfile, "\t*===================================================================*\n\n");
    
    if(!converged){
        throw ConvergenceError<int>("DCFT SCF guess", _scfMaxIter, _scfThreshold, 
                               _scfConvergence, __FILE__, __LINE__);
    }

    // Dump information to checkpoint
    _nMo = _nSo;
    _chkpt->wt_nmo(_nMo);
    int *temp = init_int_array(_nIrreps);
    _chkpt->wt_orbspi(_soPI);
    for(int h = 0; h < _nIrreps; ++h) temp[h] = _nBOccPI[h];
    _chkpt->wt_clsdpi(temp);
    for(int h = 0; h < _nIrreps; ++h) temp[h] = _nAOccPI[h] - _nBOccPI[h];
    _chkpt->wt_openpi(temp);
    // TODO implement frozen core
    for(int h = 0; h < _nIrreps; ++h) temp[h] = 0;
    _chkpt->wt_frzcpi(temp);
    _chkpt->wt_frzvpi(temp);
    delete [] temp;
    double *aEvals = _aEvals.to_block_vector();
    _chkpt->wt_alpha_evals(aEvals);
    delete [] aEvals;
    double *bEvals = _bEvals.to_block_vector();
    _chkpt->wt_beta_evals(bEvals);
    delete [] bEvals;
    _chkpt->wt_escf(_scfEnergy);
    _chkpt->wt_eref(_scfEnergy);
    write_orbitals_to_checkpoint();
}


/**
 * Dumps the current Fock matrix eigenvectors to the checkpoint file
 */
void
DCFTSolver::write_orbitals_to_checkpoint()
{
    double **aEvecs = _Ca.to_block_matrix();
    _chkpt->wt_alpha_scf(aEvecs);
    free_block(aEvecs);
    double **bEvecs = _Cb.to_block_matrix();
    _chkpt->wt_beta_scf(bEvecs);
    free_block(bEvecs);
}


/**
 * Checks to make sure that the phase, and ordering, of the MOs is consistent
 * with the previous set.
 *
 * @return Whether the phase correction was successful
 */
bool
DCFTSolver::correct_mo_phases(bool dieOnError)
{
    Matrix temp("temp", _nIrreps, _soPI, _soPI);
    Matrix overlap("Old - New Overlap", _nIrreps, _soPI, _soPI);

    temp.gemm(true, false, 1.0, _oldCa, _aoS, 0.0);
    overlap.gemm(false, false, 1.0, temp, _Ca, 0.0);
    temp.copy(_Ca);
    std::map<int, int> mosUsed;
    int offset = 0;
    for(int h = 0; h < _nIrreps; ++h){
        for(int oldMO = 0; oldMO < _soPI[h]; ++oldMO){
            int bestMO = 0;
            double maximumProjection = 0.0;
            double prefactor = 0.0;
            for(int newMO = 0; newMO < _soPI[h]; ++newMO){
                double val = overlap.get(h, oldMO, newMO);
                if(fabs(val) > maximumProjection){
                    maximumProjection = fabs(val);
                    bestMO = newMO;
                    prefactor = val < 0.0 ? -1.0 : 1.0;
                }
            }
            // Now we've found the MO to use, check it's not been used already then
            // copy it over.
            if(mosUsed[bestMO + offset]++){
                if(dieOnError){
                    overlap.print();
                    _oldCa.print();
                    temp.print();
                    throw SanityCheckError("Duplicate MOs used in phase check", __FILE__, __LINE__);
                }else{
                    // Copy Ca back from temp
                    _Ca.copy(temp);
                    return false;
                }
            }
            for(int so = 0; so < _soPI[h]; ++so){
                _Ca.set(h, so, oldMO, prefactor * temp.get(h, so, bestMO));
            }
        }
        offset += _soPI[h];
    }

    temp.gemm(true, false, 1.0, _oldCb, _aoS, 0.0);
    overlap.gemm(false, false, 1.0, temp, _Cb, 0.0);
    temp.copy(_Cb);
    mosUsed.clear();
    offset = 0;
    for(int h = 0; h < _nIrreps; ++h){
        for(int oldMO = 0; oldMO < _soPI[h]; ++oldMO){
            int bestMO = 0;
            double bestOverlap = 0.0;
            double prefactor = 0.0;
            for(int newMO = 0; newMO < _soPI[h]; ++newMO){
                double val = overlap.get(h, oldMO, newMO);
                if(fabs(val) > bestOverlap){
                    bestOverlap = fabs(val);
                    bestMO = newMO;
                    prefactor = val < 0.0 ? -1.0 : 1.0;
                }
            }
            // Now we've found the MO to use, check it's not been used already then
            // copy it over.
            if(mosUsed[bestMO + offset]++){
                if(dieOnError){
                    overlap.print();
                    _oldCb.print();
                    temp.print();
                    throw SanityCheckError("Duplicate MOs used in phase check", __FILE__, __LINE__);
                }else{
                    // Copy Cb back from temp
                    _Cb.copy(temp);
                    return false;
                }
            }
            for(int so = 0; so < _soPI[h]; ++so){
                _Cb.set(h, so, oldMO, prefactor * temp.get(h, so, bestMO));
            }
        }
        offset += _soPI[h];
    }
    return true;
}

/**
 * Figures out the orbital symmetries and prints them, ordered by energy and type
 */
void
DCFTSolver::print_orbital_energies()
{
    std::vector<std::pair<double, int> > aPairs;
    std::vector<std::pair<double, int> > bPairs;
    for (int h = 0; h < _nIrreps; ++h) {
        for (int i=0; i < _soPI[h]; ++i){
            aPairs.push_back(make_pair(_aEvals.get(h, i), h));
            bPairs.push_back(make_pair(_bEvals.get(h, i), h));
        }
    }
    sort(aPairs.begin(), aPairs.end());
    sort(bPairs.begin(), bPairs.end());

    int *aIrrepCount = init_int_array(_nIrreps);
    int *bIrrepCount = init_int_array(_nIrreps);
    char **irrepLabels = _chkpt->rd_irr_labs();

    fprintf(outfile, "\n\tOrbital energies (a.u.):\n\t\tAlpha occupied orbitals\n\t\t");
    for (int i = 0, count = 0; i < _nAOcc; ++i, ++count) {
        int irrep = aPairs[i].second;
        fprintf(outfile, "%4d%-4s%11.6f  ", ++aIrrepCount[irrep], irrepLabels[irrep], aPairs[i].first);
        if (count % 4 == 3 && i != _nAOcc)
            fprintf(outfile, "\n\t\t");
    }
    fprintf(outfile, "\n\n\t\tBeta occupied orbitals\n\t\t");
    for (int i = 0, count = 0; i < _nBOcc; ++i, ++count) {
        int irrep = bPairs[i].second;
        fprintf(outfile, "%4d%-4s%11.6f  ", ++bIrrepCount[irrep], irrepLabels[irrep], bPairs[i].first);
        if (count % 4 == 3 && i != _nBOcc)
            fprintf(outfile, "\n\t\t");
    }
    fprintf(outfile, "\n\n\t\tAlpha virtual orbitals\n\t\t");
    for (int i = _nAOcc, count = 0; i < _nMo; ++i, ++count) {
        int irrep = aPairs[i].second;
        fprintf(outfile, "%4d%-4s%11.6f  ", ++aIrrepCount[irrep], irrepLabels[irrep], aPairs[i].first);
        if (count % 4 == 3 && i != _nMo)
            fprintf(outfile, "\n\t\t");
    }
    fprintf(outfile, "\n\n\t\tBeta virtual orbitals\n\t\t");
    for (int i = _nBOcc, count = 0; i < _nMo; ++i, ++count) {
        int irrep = bPairs[i].second;
        fprintf(outfile, "%4d%-4s%11.6f  ", ++bIrrepCount[irrep], irrepLabels[irrep], bPairs[i].first);
        if (count % 4 == 3 && i != _nMo)
            fprintf(outfile, "\n\t\t");
    }
    fprintf(outfile, "\n\n");

    fprintf(outfile, "\n\tIrrep              ");
    for(int h = 0; h < _nIrreps; ++h) fprintf(outfile, "%4s ", irrepLabels[h]);
    fprintf(outfile, "\n\t-------------------");
    for(int h = 0; h < _nIrreps; ++h) fprintf(outfile, "-----");
    fprintf(outfile, "\n\t#Symmetry Orbitals ");
    for(int h = 0; h < _nIrreps; ++h) fprintf(outfile, "%4d ", _soPI[h]);
    fprintf(outfile, "\n\t#Alpha Occupied    ");
    for(int h = 0; h < _nIrreps; ++h) fprintf(outfile, "%4d ", _nAOccPI[h]);
    fprintf(outfile, "\n\t#Beta Occupied     ");
    for(int h = 0; h < _nIrreps; ++h) fprintf(outfile, "%4d ", _nBOccPI[h]);
    fprintf(outfile, "\n\t#Alpha Virtual     ");
    for(int h = 0; h < _nIrreps; ++h) fprintf(outfile, "%4d ", _nAVirPI[h]);
    fprintf(outfile, "\n\t#Beta Virtual      ");
    for(int h = 0; h < _nIrreps; ++h) fprintf(outfile, "%4d ", _nBVirPI[h]);
    fprintf(outfile, "\n\t-------------------");
    for(int h = 0; h < _nIrreps; ++h) fprintf(outfile, "-----");
    fprintf(outfile, "\n\n");

    if(_print > 2){
        _Ca.print();
        _Cb.print();
    }
    for (int h = 0; h < _nIrreps; ++h)
        delete [] irrepLabels[h];
    delete[] irrepLabels;
    delete[] aIrrepCount;
    delete[] bIrrepCount;
}
/**
 * Computes the initial Hartree-Fock orbitals by either reading them from the
 * checkpoint file or computing a core Hamiltonian guess.
 */
void
DCFTSolver::scf_guess()
{
    Matrix T("SO basis kinetic energy integrals", _nIrreps, _soPI, _soPI);
    Matrix V("SO basis potential energy integrals", _nIrreps, _soPI, _soPI);
    double *ints = init_array(_nTriSo);

    IWL::read_one(_psio.get(), PSIF_OEI, PSIF_SO_T, ints, _nTriSo, 0, 0, outfile);
    T.set(ints);
    IWL::read_one(_psio.get(), PSIF_OEI, PSIF_SO_V, ints, _nTriSo, 0, 0, outfile);
    V.set(ints);
    delete [] ints;

    _soH.add(T);
    _soH.add(V);

    Matrix tmp("temp", _nIrreps, _soPI, _soPI);

    _Fa.copy(_soH);
    _Fa.transform(_sHalfInv);
    _Fa.diagonalize(tmp, _aEvals);
    _Ca.gemm(false, false, 1.0, _sHalfInv, tmp, 0.0);
    _Cb.copy(_Ca);
    find_occupation(_aEvals);
    update_scf_density();
}


/**
 * Computes the SCF energy from the latest Fock and density matrices.
 */
void
DCFTSolver::compute_scf_energy()
{
    _scfEnergy = _eNuc;
    _scfEnergy += 0.5 * _aKappa.vector_dot(_soH);
    _scfEnergy += 0.5 * _aKappa.vector_dot(_Fa);
    _scfEnergy += 0.5 * _bKappa.vector_dot(_soH);
    _scfEnergy += 0.5 * _bKappa.vector_dot(_Fb);
}


/**
 * Computes the SCF error vector by transforming the Fock matrices to the
 * MO basis and computing [F, Kappa], and the RMS value of this quantity.
 * @return RMS error
 */
double
DCFTSolver::compute_scf_error_vector()
{
    size_t nElements = 0;
    double sumOfSquares = 0.0;
    Matrix tmp("tmp", _nIrreps, _soPI, _soPI);
    Matrix moF("MO Basis Fock Matrix", _nIrreps, _soPI, _soPI);
    tmp.gemm(true, false, 1.0, _Ca, _Fa, 0.0);
    moF.gemm(false, false, 1.0, tmp, _Ca, 0.0);
    for(int h = 0; h < _nIrreps; ++h){
        for(int i = 0; i < _nAOccPI[h]; ++i){
            for(int a = 0; a < _nAVirPI[h]; ++a){
                double kappaF = moF.get(h, i, a + _nAOccPI[h]);
                _aScfError.set(h, i, a, -kappaF);
                ++nElements;
                sumOfSquares += pow(-kappaF, 2.0);
            }
        }
    }
    tmp.gemm(true, false, 1.0, _Cb, _Fb, 0.0);
    moF.gemm(false, false, 1.0, tmp, _Cb, 0.0);
    for(int h = 0; h < _nIrreps; ++h){
        for(int i = 0; i < _nBOccPI[h]; ++i){
            for(int a = 0; a < _nBVirPI[h]; ++a){
                double kappaF = moF.get(h, i, a + _nBOccPI[h]);
                _aScfError.set(h, i, a, -kappaF);
                ++nElements;
                sumOfSquares += pow(-kappaF, 2.0);
            }
        }
    }
    return sqrt(sumOfSquares / nElements);
}


/**
 * Uses the MO coefficients to form the SCF density matrices in the SO basis.
 * @return RMS density change
 */
double
DCFTSolver::update_scf_density()
{
    size_t nElements = 0;
    double sumOfSquares = 0.0;
    Matrix old(_aKappa);
    for (int h = 0; h < _nIrreps; ++h) {
        for (int mu = 0; mu < _soPI[h]; ++mu) {
            for (int nu = 0; nu < _soPI[h]; ++nu) {
                double val = 0.0;
                for (int i = 0; i < _nAOccPI[h]; ++i)
                    val += _Ca.get(h, mu, i) * _Ca.get(h, nu, i);
                _aKappa.set(h, mu, nu, val);
                ++nElements;
                sumOfSquares += pow(val - old.get(h, mu, nu), 2.0);
            }
        }
    }
    old.copy(_bKappa);
    for (int h = 0; h < _nIrreps; ++h) {
        for (int mu = 0; mu < _soPI[h]; ++mu) {
            for (int nu = 0; nu < _soPI[h]; ++nu) {
                double val = 0.0;
                for (int i = 0; i < _nBOccPI[h]; ++i)
                    val += _Cb.get(h, mu, i) * _Cb.get(h, nu, i);
                _bKappa.set(h, mu, nu, val);
                ++nElements;
                sumOfSquares += pow(val - old.get(h, mu, nu), 2.0);
            }
        }
    }
    // We're not converged until the RMS error vector *and* the RMS density
    // changes are below the threshold
    return sqrt(sumOfSquares / nElements);
}


/**
 * Builds the G matrix for Fock matrix, the external potential (tau contracted with the integrals)
 * and other tensors, if requested, out-of-core using the SO integrals.
 * All quantities are built simultaneously to reduce I/O.
 */
void
DCFTSolver::process_so_ints()
{

    IWL *iwl = new IWL(_psio.get(), PSIF_SO_TEI, _intTolerance, 1, 1);

    Label *lblptr = iwl->labels();
    Value *valptr = iwl->values();

    double *Da = init_array(_nTriSo);
    double *Db = init_array(_nTriSo);
    double *Ta = init_array(_nTriSo);
    double *Tb = init_array(_nTriSo);
    double *Ga = init_array(_nTriSo);
    double *Gb = init_array(_nTriSo);
    double *Va = init_array(_nTriSo);
    double *Vb = init_array(_nTriSo);
    int soOffset = 0;
    for(int h = 0; h < _nIrreps; ++h){
        for(int mu = 0; mu < _soPI[h]; ++ mu){
            for(int nu = 0; nu <= mu; ++ nu){
                int muNu = INDEX((nu+soOffset), (mu+soOffset));
                Da[muNu] = _aKappa.get(h, mu, nu);
                Db[muNu] = _bKappa.get(h, mu, nu);
                Ta[muNu] = _aTau[h][mu][nu];
                Tb[muNu] = _bTau[h][mu][nu];
            }
        }
        soOffset += _soPI[h];
    }

    double value;
    int Gp, Gq, Gr, Gs, Gpr, Grp, Gps, Gsp, Gsq, Gqs, Gqr, Grq;
    int pq, rs, pr, rp, ps, sp, qr, rq, qs, sq;
    int pqArr, qpArr, rsArr, srArr, qrArr, rqArr;
    int qsArr, sqArr, psArr, spArr, prArr, rpArr;
    int labelIndex, p, q, r, s;
    dpdbuf4 *tau1_AO, *tau2_AO;

    bool buildTensors = false;
    bool lastBuffer;
    do{
        lastBuffer = iwl->last_buffer();
        for(int index = 0; index < iwl->buffer_count(); ++index){
            labelIndex = 4*index;
            p = abs((int) lblptr[labelIndex++]);
            q = (int) lblptr[labelIndex++];
            r = (int) lblptr[labelIndex++];
            s = (int) lblptr[labelIndex++];
            value = (double) valptr[index];
            if(buildTensors){
                Gp = tau1_AO->params->psym[p];
                Gq = tau1_AO->params->psym[q];
                Gr = tau1_AO->params->psym[r];
                Gs = tau1_AO->params->psym[s];

                Gpr = Grp = Gp^Gr;
                Gps = Gsp = Gp^Gs;
                Gqr = Grq = Gq^Gr;
                Gqs = Gsq = Gq^Gs;

                pq = tau1_AO->params->rowidx[p][q];
                rs = tau1_AO->params->rowidx[r][s];

                pr = tau1_AO->params->rowidx[p][r];
                rp = tau1_AO->params->rowidx[r][p];
                ps = tau1_AO->params->rowidx[p][s];
                sp = tau1_AO->params->rowidx[s][p];
                qr = tau1_AO->params->rowidx[q][r];
                rq = tau1_AO->params->rowidx[r][q];
                qs = tau1_AO->params->rowidx[q][s];
                sq = tau1_AO->params->rowidx[s][q];
            }

            qpArr = pqArr = INDEX(p, q);
            srArr = rsArr = INDEX(r, s);
            prArr = rpArr = INDEX(p, r);
            qsArr = sqArr = INDEX(q, s);
            spArr = psArr = INDEX(p, s);
            qrArr = rqArr = INDEX(q, r);

            /* (pq|rs) */
            Ga[rsArr] += (Da[pqArr] + Db[pqArr]) * value;
            Gb[rsArr] += (Da[pqArr] + Db[pqArr]) * value;
            Va[rsArr] += (Ta[pqArr] + Tb[pqArr]) * value;
            Vb[rsArr] += (Ta[pqArr] + Tb[pqArr]) * value;
            if(q >= r){
                Ga[qrArr] -= Da[psArr] * value;
                Gb[qrArr] -= Db[psArr] * value;
                Va[qrArr] -= Ta[psArr] * value;
                Vb[qrArr] -= Tb[psArr] * value;
            }
            if(buildTensors && tau1_AO->params->coltot[Gpr])
                C_DAXPY(tau1_AO->params->coltot[Gpr], value, tau1_AO->matrix[Gpr][qs], 1,
                tau2_AO->matrix[Gpr][pr], 1);

            if(p!=q && r!=s && pqArr!=rsArr){
                /* (pq|sr) */
                if(s >= r){
                    Ga[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                    Gb[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                    Va[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                    Vb[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                }
                if(q >= s){
                    Ga[qsArr] -= Da[prArr] * value;
                    Gb[qsArr] -= Db[prArr] * value;
                    Va[qsArr] -= Ta[prArr] * value;
                    Vb[qsArr] -= Tb[prArr] * value;
                }
                if(buildTensors && tau1_AO->params->coltot[Gps])
                    C_DAXPY(tau1_AO->params->coltot[Gps], value, tau1_AO->matrix[Gps][qr], 1,
                    tau2_AO->matrix[Gps][ps], 1);

                /* (qp|rs) */
                if(r >= s){
                    Ga[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Gb[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Va[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                    Vb[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                }
                if(p >= r){
                    Ga[prArr] -= Da[qsArr] * value;
                    Gb[prArr] -= Db[qsArr] * value;
                    Va[prArr] -= Ta[qsArr] * value;
                    Vb[prArr] -= Tb[qsArr] * value;
                }
                if(buildTensors && tau1_AO->params->coltot[Gqr])
                    C_DAXPY(tau1_AO->params->coltot[Gqr], value, tau1_AO->matrix[Gqr][ps], 1,
                    tau2_AO->matrix[Gqr][qr], 1);

                /* (qp|sr) */
                if(s >= r){
                    Ga[srArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Gb[srArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Va[srArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                    Vb[srArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                }
                if(p >= s){
                    Ga[psArr] -= Da[qrArr] * value;
                    Gb[psArr] -= Db[qrArr] * value;
                    Va[psArr] -= Ta[qrArr] * value;
                    Vb[psArr] -= Tb[qrArr] * value;
                }
                if(buildTensors && tau1_AO->params->coltot[Gqs])
                    C_DAXPY(tau1_AO->params->coltot[Gqs], value, tau1_AO->matrix[Gqs][pr], 1,
                    tau2_AO->matrix[Gqs][qs], 1);

                /* (rs|pq) */
                if(p >= q){
                    Ga[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Va[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                    Vb[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                }
                if(s >= p){
                    Ga[spArr] -= Da[rqArr] * value;
                    Gb[spArr] -= Db[rqArr] * value;
                    Va[spArr] -= Ta[rqArr] * value;
                    Vb[spArr] -= Tb[rqArr] * value;
                }
                if(buildTensors && tau1_AO->params->coltot[Grp])
                    C_DAXPY(tau1_AO->params->coltot[Grp], value, tau1_AO->matrix[Grp][sq], 1,
                    tau2_AO->matrix[Grp][rp], 1);

                /* (sr|pq) */
                if(p >= q){
                    Ga[pqArr] += (Da[srArr] + Db[srArr]) * value;
                    Gb[pqArr] += (Da[srArr] + Db[srArr]) * value;
                    Va[pqArr] += (Ta[srArr] + Tb[srArr]) * value;
                    Vb[pqArr] += (Ta[srArr] + Tb[srArr]) * value;
                }
                if(r >= p){
                    Ga[rpArr] -= Da[sqArr] * value;
                    Gb[rpArr] -= Db[sqArr] * value;
                    Va[rpArr] -= Ta[sqArr] * value;
                    Vb[rpArr] -= Tb[sqArr] * value;
                }
                if(buildTensors && tau1_AO->params->coltot[Gsp])
                    C_DAXPY(tau1_AO->params->coltot[Gsp], value, tau1_AO->matrix[Gsp][rq], 1,
                    tau2_AO->matrix[Gsp][sp], 1);

                /* (rs|qp) */
                if(q >= p){
                    Ga[qpArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[qpArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Va[qpArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                    Vb[qpArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                }
                if(s >= q){
                    Ga[sqArr] -= Da[rpArr] * value;
                    Gb[sqArr] -= Db[rpArr] * value;
                    Va[sqArr] -= Ta[rpArr] * value;
                    Vb[sqArr] -= Tb[rpArr] * value;
                }
                if(buildTensors && tau1_AO->params->coltot[Grq])
                    C_DAXPY(tau1_AO->params->coltot[Grq], value, tau1_AO->matrix[Grq][sp], 1,
                    tau2_AO->matrix[Grq][rq], 1);

                /* (sr|qp) */
                if(q >= p){
                    Ga[qpArr] += (Da[srArr] + Db[srArr]) * value;
                    Gb[qpArr] += (Da[srArr] + Db[srArr]) * value;
                    Va[qpArr] += (Ta[srArr] + Tb[srArr]) * value;
                    Vb[qpArr] += (Ta[srArr] + Tb[srArr]) * value;
                }
                if(r >= q){
                    Ga[rqArr] -= Da[spArr] * value;
                    Gb[rqArr] -= Db[spArr] * value;
                    Va[rqArr] -= Ta[spArr] * value;
                    Vb[rqArr] -= Tb[spArr] * value;
                }
                if(buildTensors && tau1_AO->params->coltot[Gsq])
                    C_DAXPY(tau1_AO->params->coltot[Gsq], value, tau1_AO->matrix[Gsq][rp], 1,
                    tau2_AO->matrix[Gsq][sq],1 );
            }else if(p!=q && r!=s && pqArr==rsArr){
                /* (pq|sr) */
                if(s >= r){
                    Ga[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                    Gb[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                    Va[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                    Vb[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                }
                if(q >= s){
                    Ga[qsArr] -= Da[prArr] * value;
                    Gb[qsArr] -= Db[prArr] * value;
                    Va[qsArr] -= Ta[prArr] * value;
                    Vb[qsArr] -= Tb[prArr] * value;
                }
                if(buildTensors && tau1_AO->params->coltot[Gps])
                    C_DAXPY(tau1_AO->params->coltot[Gps], value, tau1_AO->matrix[Gps][qr], 1,
                    tau2_AO->matrix[Gps][ps], 1);
                /* (qp|rs) */
                if(r >= s){
                    Ga[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Gb[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Va[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                    Vb[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                }
                if(p >= r){
                    Ga[prArr] -= Da[qsArr] * value;
                    Gb[prArr] -= Db[qsArr] * value;
                    Va[prArr] -= Ta[qsArr] * value;
                    Vb[prArr] -= Tb[qsArr] * value;
                }
                if(buildTensors && tau1_AO->params->coltot[Gqr])
                    C_DAXPY(tau1_AO->params->coltot[Gqr], value, tau1_AO->matrix[Gqr][ps], 1,
                    tau2_AO->matrix[Gqr][qr], 1);

                /* (qp|sr) */
                if(s >= r){
                    Ga[srArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Gb[srArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Va[srArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                    Vb[srArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                }
                if(p >= s){
                    Ga[psArr] -= Da[qrArr] * value;
                    Gb[psArr] -= Db[qrArr] * value;
                    Va[psArr] -= Ta[qrArr] * value;
                    Vb[psArr] -= Tb[qrArr] * value;
                }
                if(buildTensors && tau1_AO->params->coltot[Gqs])
                    C_DAXPY(tau1_AO->params->coltot[Gqs], value, tau1_AO->matrix[Gqs][pr], 1,
                    tau2_AO->matrix[Gqs][qs], 1);
            }else if(p!=q && r==s){
                /* (qp|rs) */
                if(r >= s){
                    Ga[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Gb[rsArr] += (Da[qpArr] + Db[qpArr]) * value;
                    Va[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                    Vb[rsArr] += (Ta[qpArr] + Tb[qpArr]) * value;
                }
                if(p >= r){
                    Ga[prArr] -= Da[qsArr] * value;
                    Gb[prArr] -= Db[qsArr] * value;
                    Va[prArr] -= Ta[qsArr] * value;
                    Vb[prArr] -= Tb[qsArr] * value;
                }
                if(buildTensors && tau1_AO->params->coltot[Gqr])
                    C_DAXPY(tau1_AO->params->coltot[Gqr], value, tau1_AO->matrix[Gqr][ps], 1,
                    tau2_AO->matrix[Gqr][qr], 1);

                /* (rs|pq) */
                if(p >= q){
                    Ga[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Va[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                    Vb[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                }
                if(s >= p){
                    Ga[spArr] -= Da[rqArr] * value;
                    Gb[spArr] -= Db[rqArr] * value;
                    Va[spArr] -= Ta[rqArr] * value;
                    Vb[spArr] -= Tb[rqArr] * value;
                }
                if(buildTensors && tau1_AO->params->coltot[Grp])
                    C_DAXPY(tau1_AO->params->coltot[Grp], value, tau1_AO->matrix[Grp][sq], 1,
                    tau2_AO->matrix[Grp][rp], 1);

                /* (rs|qp) */
                if(q >= p){
                    Ga[qpArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[qpArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Va[qpArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                    Vb[qpArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                }
                if(s >= q){
                    Ga[sqArr] -= Da[rpArr] * value;
                    Gb[sqArr] -= Db[rpArr] * value;
                    Va[sqArr] -= Ta[rpArr] * value;
                    Vb[sqArr] -= Tb[rpArr] * value;
                }
                if(buildTensors && tau1_AO->params->coltot[Grq])
                    C_DAXPY(tau1_AO->params->coltot[Grq], value, tau1_AO->matrix[Grq][sp], 1,
                    tau2_AO->matrix[Grq][rq], 1);
            }else if(p==q && r!=s){
                /* (pq|sr) */
                if(s >= r){
                    Ga[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                    Gb[srArr] += (Da[pqArr] + Db[pqArr]) * value;
                    Va[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                    Vb[srArr] += (Ta[pqArr] + Tb[pqArr]) * value;
                }
                if(q >= s){
                    Ga[qsArr] -= Da[prArr] * value;
                    Gb[qsArr] -= Db[prArr] * value;
                    Va[qsArr] -= Ta[prArr] * value;
                    Vb[qsArr] -= Tb[prArr] * value;
                }
                if(buildTensors && tau1_AO->params->coltot[Gps])
                    C_DAXPY(tau1_AO->params->coltot[Gps], value, tau1_AO->matrix[Gps][qr], 1,
                    tau2_AO->matrix[Gps][ps], 1);

                /* (rs|pq) */
                if(p >= q){
                    Ga[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Va[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                    Vb[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                }
                if(s >= p){
                    Ga[spArr] -= Da[rqArr] * value;
                    Gb[spArr] -= Db[rqArr] * value;
                    Va[spArr] -= Ta[rqArr] * value;
                    Vb[spArr] -= Tb[rqArr] * value;
                }
                if(buildTensors && tau1_AO->params->coltot[Grp])
                    C_DAXPY(tau1_AO->params->coltot[Grp], value, tau1_AO->matrix[Grp][sq], 1,
                    tau2_AO->matrix[Grp][rp], 1);

                /* (sr|pq) */
                if(p >= q){
                    Ga[pqArr] += (Da[srArr] + Db[srArr]) * value;
                    Gb[pqArr] += (Da[srArr] + Db[srArr]) * value;
                    Va[pqArr] += (Ta[srArr] + Tb[srArr]) * value;
                    Vb[pqArr] += (Ta[srArr] + Tb[srArr]) * value;
                }
                if(r >= p){
                    Ga[rpArr] -= Da[sqArr] * value;
                    Gb[rpArr] -= Db[sqArr] * value;
                    Va[rpArr] -= Ta[sqArr] * value;
                    Vb[rpArr] -= Tb[sqArr] * value;
                }
                if(buildTensors && tau1_AO->params->coltot[Gsp])
                    C_DAXPY(tau1_AO->params->coltot[Gsp], value, tau1_AO->matrix[Gsp][rq], 1,
                    tau2_AO->matrix[Gsp][sp], 1);
            }else if(p==q && r==s && pqArr!=rsArr){
                /* (rs|pq) */
                if(p >= q){
                    Ga[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Gb[pqArr] += (Da[rsArr] + Db[rsArr]) * value;
                    Va[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                    Vb[pqArr] += (Ta[rsArr] + Tb[rsArr]) * value;
                }
                if(s >= p){
                    Ga[spArr] -= Da[rqArr] * value;
                    Gb[spArr] -= Db[rqArr] * value;
                    Va[spArr] -= Ta[rqArr] * value;
                    Vb[spArr] -= Tb[rqArr] * value;
                }
                if(buildTensors && tau1_AO->params->coltot[Grp])
                    C_DAXPY(tau1_AO->params->coltot[Grp], value, tau1_AO->matrix[Grp][sq], 1,
                    tau2_AO->matrix[Grp][rp], 1);
            }
        } /* end loop through current buffer */
        if(!lastBuffer) iwl->fetch();
    }while(!lastBuffer);
    iwl->set_keep_flag(1);
    delete iwl;

    // Build the Fock matrices from the H and G matrices
    soOffset = 0;
    for(int h = 0; h < _nIrreps; ++h){
        for(int mu = 0; mu < _soPI[h]; ++mu){
            for(int nu = 0; nu <= mu; ++nu){
                int muNu = INDEX((nu+soOffset), (mu+soOffset));
                double aVal   = Ga[muNu];
                double bVal   = Gb[muNu];
                double aGTVal = Va[muNu];
                double bGTVal = Vb[muNu];
                _Fa.add(h, mu, nu, aVal);
                _Fb.add(h, mu, nu, bVal);
                _aGTau.set(h, mu, nu, aGTVal);
                _bGTau.set(h, mu, nu, bGTVal);
                if(mu != nu){
                    _Fa.add(h, nu, mu, aVal);
                    _Fb.add(h, nu, mu, bVal);
                    _aGTau.set(h, nu, mu, aGTVal);
                    _bGTau.set(h, nu, mu, bGTVal);
                }
            }
        }
        soOffset += _soPI[h];
    }

    delete [] Ta;
    delete [] Tb;
    delete [] Va;
    delete [] Vb;
    delete [] Da;
    delete [] Db;
    delete [] Ga;
    delete [] Gb;
}

/**
 * Uses the orbital energies to determine the occupation according to
 * the aufbau principle
 */
void
DCFTSolver::find_occupation(Vector & evals, bool forcePrint)
{
    std::vector<std::pair<double, int> > pairs;
    for (int h=0; h<evals.nirreps(); ++h) {
        for (int i=0; i<evals.dimpi()[h]; ++i)
            pairs.push_back(make_pair(evals.get(h, i), h));
    }
    sort(pairs.begin(),pairs.end());

    if(_inputDocc){
        for(int h = 0; h < _nIrreps; ++h)
            _nBOccPI[h] = _options["DOCC"][h].to_integer();
    }else{
        memset(_nBOccPI, 0, sizeof(int) * _nIrreps);
        for (int i=0; i < _nBOcc; ++i)
            _nBOccPI[pairs[i].second]++;
    }
    if(_inputSocc){
        for(int h = 0; h < _nIrreps; ++h)
            _nAOccPI[h] = _nBOccPI[h] + _options["SOCC"][h].to_integer();
    }else{
        for(int h = 0; h < _nIrreps; ++h)
            _nAOccPI[h] = _nBOccPI[h];
        for (int i=_nBOcc; i < _nAOcc; ++i)
            _nAOccPI[pairs[i].second]++;
    }

    if(_print > 1 || forcePrint){
        fprintf(outfile, "\t\t\t\tDOCC: [");
        for (int h = 0; h < evals.nirreps(); ++h){
            fprintf(outfile, "%3d ", _nBOccPI[h]);
        }
        fprintf(outfile, "]\n");
        fprintf(outfile, "\t\t\t\tSOCC: [");
        for (int h = 0; h < evals.nirreps(); ++h){
            fprintf(outfile, "%3d ", _nAOccPI[h] - _nBOccPI[h]);
        }
        fprintf(outfile, "]\n");
    }
}

}} // Namespaces


