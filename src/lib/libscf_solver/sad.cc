/*
 *  sad.cc
 *
 * Routines for the high-meantenance SAD guess
 * and daul-basis projections
 *
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libparallel/parallel.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include <libmints/mints.h>

#include "hf.h"
#include "sad.h"

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace scf {

SADGuess::SADGuess(boost::shared_ptr<BasisSet> basis, int nalpha, int nbeta, Options& options) :
    basis_(basis), nalpha_(nalpha), nbeta_(nbeta), options_(options)
{
    common_init();
}
SADGuess::~SADGuess()
{
}
void SADGuess::common_init()
{
    molecule_ = basis_->molecule();

    boost::shared_ptr<IntegralFactory> ints(new IntegralFactory(basis_));
    boost::shared_ptr<PetiteList> petite(new PetiteList(basis_,ints));
    AO2SO_ =  petite->aotoso();

    print_ = options_.get_int("SAD_PRINT");
    debug_ = options_.get_int("DEBUG");
}
void SADGuess::compute_guess()
{
    form_D();
    form_C();
}
void SADGuess::form_D()
{
    // Build Neutral D in AO basis (block diagonal)
    SharedMatrix DAO = form_D_AO();

    // Transform Neutral D from AO to SO basis
    Da_ = SharedMatrix(new Matrix("Da SAD",AO2SO_->colspi(),AO2SO_->colspi()));

    double* temp = new double[AO2SO_->rowspi()[0] * (ULI) AO2SO_->max_ncol()];
    for (int h = 0; h < Da_->nirrep(); h++) {
        int nao = AO2SO_->rowspi()[h];
        int nso = AO2SO_->colspi()[h];
        if (!nao || !nso) continue;

        double** DAOp = DAO->pointer();
        double** DSOp = Da_->pointer(h);
        double** Up = AO2SO_->pointer(h);

        C_DGEMM('N','N',nao,nso,nao,1.0,DAOp[0],nao,Up[0],nso,0.0,temp,nso);
        C_DGEMM('T','N',nso,nso,nao,1.0,Up[0],nso,temp,nso,0.0,DSOp[0],nso);
    }
    delete[] temp;

    // Scale Da to true electron count
    double npair = 0.0;
    for (int A = 0; A < molecule_->natom(); A++) {
        npair += 0.5 * molecule_->Z(A);
    }
    Da_->scale(((double) nalpha_) / npair);

    // Build/Scale Db if needed
    if (nalpha_ == nbeta_) {
        Db_ = Da_;
    } else {
        Db_ = SharedMatrix(Da_->clone());
        Db_->set_name("Db SAD");
        Db_->scale(((double) nbeta_) / ((double) nalpha_));
    }

    if (debug_) {
        Da_->print();
        Db_->print();
    }
}
void SADGuess::form_C()
{
    Ca_ = Da_->partial_cholesky_factorize(options_.get_double("SAD_CHOL_TOLERANCE"));
    Ca_->set_name("Ca SAD");
    if (nalpha_ == nbeta_) {

        Cb_ = Ca_;
    } else {
        Cb_ = SharedMatrix(Ca_->clone());
        Cb_->set_name("Cb SAD");
        Cb_->scale(sqrt(((double)nbeta_)/((double)nalpha_)));
    }

    if (debug_) {
        Ca_->print();
        Cb_->print();
    }
}
SharedMatrix SADGuess::form_D_AO()
{
    std::vector<boost::shared_ptr<BasisSet> > atomic_bases;

    if (print_ > 6) {
        fprintf(outfile,"\n  Constructing atomic basis sets\n  Molecule:\n");
        molecule_->print();
    }

    //Build the atomic basis sets for libmints use in UHF
    for (int A = 0; A<molecule_->natom(); A++) {
        atomic_bases.push_back(basis_->atomic_basis_set(A));
        if (print_>6) {
            fprintf(outfile,"  SAD: Atomic Basis Set %d\n", A);
            atomic_bases[A]->molecule()->print();
            fprintf(outfile,"\n");
            atomic_bases[A]->print(outfile);
            fprintf(outfile,"\n");
        }
    }

    //Spin occupations per atom, to be determined by Hund's Rules
    //or user input
    int* nalpha = init_int_array(molecule_->natom());
    int* nbeta = init_int_array(molecule_->natom());
    int* nelec = init_int_array(molecule_->natom());
    int* nhigh = init_int_array(molecule_->natom());
    int tot_elec = 0;

    //Ground state high spin occupency array, atoms 0 to 36 (see Giffith's Quantum Mechanics, pp. 217)
    const int reference_S[] = {0,1,0,1,0,1,2,3,2,1,0,1,0,1,2,3,2,1,0,1,0,1,2,3,6,5,4,3,2,1,0,1,2,3,2,1,0};
    const int MAX_Z = 36;

    if (print_>1)
        fprintf(outfile,"  Determining Atomic Occupations\n");
    for (int A = 0; A<molecule_->natom(); A++) {
        int Z = molecule_->Z(A);
        if (Z>MAX_Z) {
            throw std::domain_error(" Only Atoms up to 36 (Kr) are currently supported with SAD Guess");
        }
        nhigh[A] = reference_S[Z];
        nelec[A] = Z;
        tot_elec+= nelec[A];
        nbeta[A] = (nelec[A]-nhigh[A])/2;
        nalpha[A] = nelec[A]-nbeta[A];
        if (print_>1)
            fprintf(outfile,"  Atom %d, Z = %d, nelec = %d, nhigh = %d, nalpha = %d, nbeta = %d\n",A,Z,nelec[A],nhigh[A],nalpha[A],nbeta[A]);
    }
    fflush(outfile);

    // Determine redundant atoms
    int* unique_indices = init_int_array(molecule_->natom()); // All atoms to representative unique atom
    int* atomic_indices = init_int_array(molecule_->natom()); // unique atom to first representative atom
    int* offset_indices = init_int_array(molecule_->natom()); // unique atom index to rank
    int nunique = 0;
    for (int l = 0; l < molecule_->natom(); l++) {
        unique_indices[l] = l;
        atomic_indices[l] = l;
    }
    for (int l = 0; l < molecule_->natom() - 1; l++) {
        for (int m = l + 1; m < molecule_->natom(); m++) {
            if (unique_indices[m] != m)
                continue; //Already assigned
            if (molecule_->Z(l) != molecule_->Z(m))
                continue;
            if (nalpha[l] != nalpha[m])
                continue;
            if (nbeta[l] != nbeta[m])
                continue;
            if (nhigh[l] != nhigh[m])
                continue;
            if (nelec[l] != nelec[m])
                continue;
            if (atomic_bases[l]->nbf() != atomic_bases[m]->nbf())
                continue;
            if (atomic_bases[l]->nshell() != atomic_bases[m]->nshell())
                continue;
            if (atomic_bases[l]->nprimitive() != atomic_bases[m]->nprimitive())
                continue;
            if (atomic_bases[l]->max_am() != atomic_bases[m]->max_am())
                continue;
            if (atomic_bases[l]->max_nprimitive() != atomic_bases[m]->max_nprimitive())
                continue;
            if (atomic_bases[l]->has_puream() !=  atomic_bases[m]->has_puream())
                continue;

            // Semi-Rigorous match obtained
            unique_indices[m] = l;
        }
    }
    for (int l = 0; l < molecule_->natom(); l++) {
        if (unique_indices[l] == l) {
            atomic_indices[nunique] = l;
            offset_indices[l] = nunique;
            nunique++;
        }
    }

    //Atomic D matrices within the atom specific AO basis
    double*** atomic_D = (double***)malloc(nunique*sizeof(double**));
    for (int A = 0; A<nunique; A++) {
        atomic_D[A] = block_matrix(atomic_bases[atomic_indices[A]]->nbf(),atomic_bases[atomic_indices[A]]->nbf());
    }

    if (print_ > 1)
        fprintf(outfile,"\n  Performing Atomic UHF Computations:\n");
    for (int A = 0; A<nunique; A++) {
        int index = atomic_indices[A];
        if (print_ > 1)
            fprintf(outfile,"\n  UHF Computation for Unique Atom %d which is Atom %d:",A, index);
        getUHFAtomicDensity(atomic_bases[index],nelec[index],nhigh[index],atomic_D[A]);
    }
    if (print_)
        fprintf(outfile,"\n");

    fflush(outfile);

    //Add atomic_D into D (scale by 1/2, we like effective pairs)
    SharedMatrix DAO = SharedMatrix(new Matrix("D_SAD (AO)", basis_->nbf(), basis_->nbf()));
    for (int A = 0, offset = 0; A < molecule_->natom(); A++) {
        int norbs = atomic_bases[A]->nbf();
        int back_index = unique_indices[A];
        for (int m = 0; m<norbs; m++)
            for (int n = 0; n<norbs; n++)
                DAO->set(0,m+offset,n+offset,0.5*atomic_D[offset_indices[back_index]][m][n]);
        offset+=norbs;
    }

    for (int A = 0; A<nunique; A++)
        free_block(atomic_D[A]);
    free(atomic_D);

    free(nelec);
    free(nhigh);
    free(nalpha);
    free(nbeta);

    free(atomic_indices);
    free(unique_indices);
    free(offset_indices);

    if (debug_) {
        DAO->print();
    }

    return DAO;
}
void SADGuess::getUHFAtomicDensity(boost::shared_ptr<BasisSet> bas, int nelec, int nhigh, double** D)
{
    boost::shared_ptr<Molecule> mol = bas->molecule();

    int nbeta = (nelec-nhigh)/2;
    int nalpha = nelec-nbeta;
    int natom = mol->natom();
    int norbs = bas->nbf();

    if (print_>1) {
        fprintf(outfile,"\n");
        bas->print(outfile);
        fprintf(outfile,"  Occupation: nalpha = %d, nbeta = %d, norbs = %d\n",nalpha,nbeta,norbs);
        fprintf(outfile,"\n  Atom:\n");
        mol->print();
    }

    if (natom != 1) {
        throw std::domain_error("SAD Atomic UHF has been given a molecule, not an atom");
    }

    double** Dold = block_matrix(norbs,norbs);
    double **Shalf = block_matrix(norbs, norbs);
    double** Ca = block_matrix(norbs,norbs);
    double** Cb = block_matrix(norbs,norbs);
    double** Da = block_matrix(norbs,norbs);
    double** Db = block_matrix(norbs,norbs);
    double** Fa = block_matrix(norbs,norbs);
    double** Fb = block_matrix(norbs,norbs);
    double** Fa_old = block_matrix(norbs,norbs);
    double** Fb_old = block_matrix(norbs,norbs);
    double** Ga = block_matrix(norbs,norbs);
    double** Gb = block_matrix(norbs,norbs);

    IntegralFactory integral(bas, bas, bas, bas);
    MatrixFactory mat;
    mat.init_with(1,&norbs,&norbs);
    OneBodyAOInt *S_ints = integral.ao_overlap();
    OneBodyAOInt *T_ints = integral.ao_kinetic();
    OneBodyAOInt *V_ints = integral.ao_potential();
    TwoBodyAOInt *TEI = integral.eri();

    //Compute Shalf;
    //Fill S
    SharedMatrix S_UHF(mat.create_matrix("S_UHF"));
    S_ints->compute(S_UHF);
    double** S = S_UHF->to_block_matrix();

    if (print_>6) {
        fprintf(outfile,"  S:\n");
        print_mat(S,norbs,norbs,outfile);
    }
    // S^{-1/2}

    // First, diagonalize S
    // the C_DSYEV call replaces the original matrix J with its eigenvectors
    double* eigval = init_array(norbs);
    int lwork = norbs * 3;
    double* work = init_array(lwork);
    int stat = C_DSYEV('v','u',norbs,S[0],norbs,eigval, work,lwork);
    if (stat != 0) {
        fprintf(outfile, "C_DSYEV failed\n");
        exit(PSI_RETURN_FAILURE);
    }
    free(work);

    double **S_copy = block_matrix(norbs, norbs);
    C_DCOPY(norbs*norbs,S[0],1,S_copy[0],1);

    // Now form S^{-1/2} = U(T)*s^{-1/2}*U,
    // where s^{-1/2} is the diagonal matrix of the inverse square roots
    // of the eigenvalues, and U is the matrix of eigenvectors of S
    for (int i=0; i<norbs; i++) {
        if (eigval[i] < 1.0E-10)
            eigval[i] = 0.0;
        else
            eigval[i] = 1.0 / sqrt(eigval[i]);

        // scale one set of eigenvectors by the diagonal elements s^{-1/2}
        C_DSCAL(norbs, eigval[i], S[i], 1);
    }
    free(eigval);

    // Smhalf = S_copy(T) * S
    C_DGEMM('t','n',norbs,norbs,norbs,1.0,
            S_copy[0],norbs,S[0],norbs,0.0,Shalf[0],norbs);

    free_block(S);
    free_block(S_copy);

    if (print_>6) {
        fprintf(outfile,"  S^-1/2:\n");
        print_mat(Shalf,norbs,norbs,outfile);
    }

    //Compute H
    SharedMatrix T_UHF(mat.create_matrix("T_UHF"));
    T_ints->compute(T_UHF);
    SharedMatrix V_UHF(mat.create_matrix("V_UHF"));
    V_ints->compute(V_UHF);
    SharedMatrix H_UHF(mat.create_matrix("H_UHF"));
    H_UHF->zero();
    H_UHF->add(T_UHF);
    H_UHF->add(V_UHF);
    double** H = H_UHF->to_block_matrix();

    delete S_ints;
    delete T_ints;
    delete V_ints;

    if (print_>6) {
        fprintf(outfile,"  H:\n");
        print_mat(H,norbs,norbs,outfile);
    }

    //Compute initial Ca and Da from core guess
    atomicUHFHelperFormCandD(nalpha,norbs,Shalf,H,Ca,Da);
    //Compute initial Cb and Db from core guess
    atomicUHFHelperFormCandD(nbeta,norbs,Shalf,H,Cb,Db);
    //Compute intial D
    C_DCOPY(norbs*norbs,Da[0],1,D[0],1);
    C_DAXPY(norbs*norbs,1.0,Db[0],1,D[0],1);
    if (print_>6) {
        fprintf(outfile,"  Ca:\n");
        print_mat(Ca,norbs,norbs,outfile);

        fprintf(outfile,"  Cb:\n");
        print_mat(Cb,norbs,norbs,outfile);

        fprintf(outfile,"  Da:\n");
        print_mat(Da,norbs,norbs,outfile);

        fprintf(outfile,"  Db:\n");
        print_mat(Db,norbs,norbs,outfile);

        fprintf(outfile,"  D:\n");
        print_mat(D,norbs,norbs,outfile);
    }

    //Compute inital E for reference
    double E = C_DDOT(norbs*norbs,D[0],1,H[0],1);
    E += C_DDOT(norbs*norbs,Da[0],1,Fa[0],1);
    E += C_DDOT(norbs*norbs,Db[0],1,Fb[0],1);
    E *= 0.5;

    const double* buffer = TEI->buffer();

    double E_tol = options_.get_double("SAD_E_CONVERGENCE");
    double D_tol = options_.get_double("SAD_D_CONVERGENCE");
    int maxiter = options_.get_int("SAD_MAXITER");
    int f_mixing_iteration = options_.get_int("SAD_F_MIX_START");

    double E_old;
    int iteration = 0;

    bool converged = false;
    if (print_>1) {
        fprintf(outfile, "\n  Initial Atomic UHF Energy:    %14.10f\n\n",E);
        fprintf(outfile, "                                         Total Energy            Delta E              Density RMS\n\n");
        fflush(outfile);
    }
    do {

        iteration++;

        //Copy the old values over for error analysis
        E_old = E;
        //I'm only going to use the total for now, could be expanded later
        C_DCOPY(norbs*norbs,D[0],1,Dold[0],1);
        //And old Fock matrices for level shift
        C_DCOPY(norbs*norbs,Fa[0],1,Fa_old[0],1);
        C_DCOPY(norbs*norbs,Fb[0],1,Fb_old[0],1);

        //Form Ga and Gb via integral direct
        memset((void*) Ga[0], '\0',norbs*norbs*sizeof(double));
        memset((void*) Gb[0], '\0',norbs*norbs*sizeof(double));

        //At the moment this is 8-fold slower than it could be, we'll see if it is signficant
        for (int MU = 0; MU < bas->nshell(); MU++) {
        int numMU = bas->shell(MU).nfunction();
        for (int NU = 0; NU < bas->nshell(); NU++) {
        int numNU = bas->shell(NU).nfunction();
        for (int LA = 0; LA < bas->nshell(); LA++) {
        int numLA = bas->shell(LA).nfunction();
        for (int SI = 0; SI < bas->nshell(); SI++) {
        int numSI = bas->shell(SI).nfunction();
        TEI->compute_shell(MU,NU,LA,SI);
        for (int m = 0, index = 0; m < numMU; m++) {
        int omu = bas->shell(MU).function_index() + m;
        for (int n = 0; n < numNU; n++) {
        int onu = bas->shell(NU).function_index() + n;
        for (int l = 0; l < numLA; l++) {
        int ola = bas->shell(LA).function_index() + l;
        for (int s = 0; s < numSI; s++, index++) {
        int osi = bas->shell(SI).function_index() + s;
             //fprintf(outfile,"  Integral (%d, %d| %d, %d) = %14.10f\n",omu,onu,ola,osi,buffer[index]);
             Ga[omu][onu] += D[ola][osi]*buffer[index];
             //Ga[ola][osi] += D[omu][onu]*buffer[index];
             Ga[omu][osi] -= Da[onu][ola]*buffer[index];
             Gb[omu][onu] += D[ola][osi]*buffer[index];
             //Gb[ola][osi] += D[omu][onu]*buffer[index];
             Gb[omu][osi] -= Db[onu][ola]*buffer[index];
        } } } } } } } }

        //Form Fa and Fb
        C_DCOPY(norbs*norbs,H[0],1,Fa[0],1);
        C_DAXPY(norbs*norbs,1.0,Ga[0],1,Fa[0],1);
        C_DCOPY(norbs*norbs,H[0],1,Fb[0],1);
        C_DAXPY(norbs*norbs,1.0,Gb[0],1,Fb[0],1);

        //Compute E
        E = C_DDOT(norbs*norbs,D[0],1,H[0],1);
        E += C_DDOT(norbs*norbs,Da[0],1,Fa[0],1);
        E += C_DDOT(norbs*norbs,Db[0],1,Fb[0],1);
        E *= 0.5;

        //Perform any required convergence stabilization
        //F-mixing (should help with oscillation)
        //20% old, 80% new Fock matrix for now
        if (iteration >= f_mixing_iteration) {
            C_DSCAL(norbs*norbs,0.8,Fa[0],1);
            C_DSCAL(norbs*norbs,0.8,Fb[0],1);
            C_DAXPY(norbs*norbs,0.2,Fa_old[0],1,Fa[0],1);
            C_DAXPY(norbs*norbs,0.2,Fb_old[0],1,Fb[0],1);
        }

        //Diagonalize Fa and Fb to from Ca and Cb and Da and Db
        atomicUHFHelperFormCandD(nalpha,norbs,Shalf,Fa,Ca,Da);
        atomicUHFHelperFormCandD(nbeta,norbs,Shalf,Fb,Cb,Db);

        //Form D
        C_DCOPY(norbs*norbs,Da[0],1,D[0],1);
        C_DAXPY(norbs*norbs,1.0,Db[0],1,D[0],1);

        //Form delta D and Drms
        C_DAXPY(norbs*norbs,-1.0,D[0],1,Dold[0],1);
        double Drms = sqrt(1.0/(1.0*norbs*norbs)*C_DDOT(norbs*norbs,Dold[0],1,Dold[0],1));

        double deltaE = fabs(E-E_old);

        if (print_>6) {
            fprintf(outfile,"  Fa:\n");
            print_mat(Fa,norbs,norbs,outfile);

            fprintf(outfile,"  Fb:\n");
            print_mat(Fb,norbs,norbs,outfile);

            fprintf(outfile,"  Ga:\n");
            print_mat(Ga,norbs,norbs,outfile);

            fprintf(outfile,"  Gb:\n");
            print_mat(Gb,norbs,norbs,outfile);

            fprintf(outfile,"  Ca:\n");
            print_mat(Ca,norbs,norbs,outfile);

            fprintf(outfile,"  Cb:\n");
            print_mat(Cb,norbs,norbs,outfile);

            fprintf(outfile,"  Da:\n");
            print_mat(Da,norbs,norbs,outfile);

            fprintf(outfile,"  Db:\n");
            print_mat(Db,norbs,norbs,outfile);

            fprintf(outfile,"  D:\n");
            print_mat(D,norbs,norbs,outfile);
        }
        if (print_>1)
            fprintf(outfile, "  @Atomic UHF iteration %3d energy: %20.14f    %20.14f %20.14f\n", iteration, E, E-E_old, Drms);
        if (iteration > 1 && deltaE < E_tol && Drms < D_tol)
            converged = true;

        if (iteration > maxiter) {
            fprintf(outfile, "\n WARNING: Atomic UHF is not converging! Try casting from a smaller basis or call Rob at CCMST.\n");
            break;
        }

        //Check convergence
    } while (!converged);
    if (converged && print_ > 1)
        fprintf(outfile, "  @Atomic UHF Final Energy for atom %s: %20.14f\n", mol->symbol(0).c_str(),E);

    delete TEI;
    free_block(Dold);
    free_block(Ca);
    free_block(Cb);
    free_block(Da);
    free_block(Db);
    free_block(Fa);
    free_block(Fb);
    free_block(Fa_old);
    free_block(Fb_old);
    free_block(Ga);
    free_block(Gb);
    free_block(H);
    free_block(Shalf);
}
void SADGuess::atomicUHFHelperFormCandD(int nelec, int norbs,double** Shalf, double**F, double** C, double** D)
{
    //Forms C in the AO basis for SAD Guesses
    double **Temp = block_matrix(norbs,norbs);
    double **Fp = block_matrix(norbs,norbs);
    double **Cp = block_matrix(norbs,norbs);

    //Form F' = X'FX = XFX for symmetric orthogonalization
    C_DGEMM('N','N',norbs,norbs,norbs,1.0,Shalf[0],norbs,F[0],norbs,0.0,Temp[0],norbs);
    C_DGEMM('N','N',norbs,norbs,norbs,1.0,Temp[0],norbs,Shalf[0],norbs,0.0,Fp[0],norbs);

    //Form C' = eig(F')
    double *eigvals = init_array(norbs);
    sq_rsp(norbs, norbs, Fp,  eigvals, 1, Cp, 1.0e-14);
    free(eigvals);

    //Form C = XC'
    C_DGEMM('N','N',norbs,norbs,norbs,1.0,Shalf[0],norbs,Cp[0],norbs,0.0,C[0],norbs);

    //Form D = Cocc*Cocc'
    C_DGEMM('N','T',norbs,norbs,nelec,1.0,C[0],norbs,C[0],norbs,0.0,D[0],norbs);

    free_block(Temp);
    free_block(Cp);
    free_block(Fp);
}

void HF::compute_SAD_guess()
{
    boost::shared_ptr<SADGuess> guess(new SADGuess(basisset_,nalpha_,nbeta_,options_));
    guess->compute_guess();

    Da_->copy(guess->Da());
    Db_->copy(guess->Db());

    for (int h = 0; h < Da_->nirrep(); h++) {

        int nso = guess->Ca()->rowspi()[h];
        int nmo = guess->Ca()->colspi()[h];
        if (nmo > X_->colspi()[h])
            nmo = X_->colspi()[h];

        sad_nocc_[h] = nmo;

        if (!nso || !nmo) continue;

        double** Cap = Ca_->pointer(h);
        double** Cbp = Cb_->pointer(h);
        double** Ca2p = guess->Ca()->pointer(h);
        double** Cb2p = guess->Cb()->pointer(h);

        for (int i = 0; i < nso; i++) {
            ::memcpy((void*) Cap[i], (void*) Ca2p[i], nmo*sizeof(double));
            ::memcpy((void*) Cbp[i], (void*) Cb2p[i], nmo*sizeof(double));
        }
    }

    int temp_nocc;
    for (int h = 0 ; h < Da_->nirrep(); h++) {
        temp_nocc = sad_nocc_[h];
        sad_nocc_[h] = doccpi_[h];
        nalphapi_[h] = temp_nocc;
        nbetapi_[h]  = temp_nocc;
        doccpi_[h]   = temp_nocc;
        soccpi_[h]   = 0;
    }

    fflush(outfile);

    E_ = 0.0; // This is the -1th iteration
}
SharedMatrix HF::BasisProjection(SharedMatrix C_A, int* noccpi, boost::shared_ptr<BasisSet> old_basis, boost::shared_ptr<BasisSet> new_basis)
{

    //Based on Werner's method from Mol. Phys. 102, 21-22, 2311
    boost::shared_ptr<IntegralFactory> newfactory(new IntegralFactory(new_basis,new_basis,new_basis,new_basis));
    boost::shared_ptr<IntegralFactory> hybfactory(new IntegralFactory(old_basis,new_basis,old_basis,new_basis));
    boost::shared_ptr<OneBodySOInt> intBB(newfactory->so_overlap());
    boost::shared_ptr<OneBodySOInt> intAB(hybfactory->so_overlap());

    boost::shared_ptr<PetiteList> pet(new PetiteList(new_basis, newfactory));
    SharedMatrix AO2USO(pet->aotoso());

    SharedMatrix SAB(new Matrix("S_AB", C_A->nirrep(), C_A->rowspi(), AO2USO->colspi()));
    SharedMatrix SBB(new Matrix("S_BB", C_A->nirrep(), AO2USO->colspi(), AO2USO->colspi()));

    intAB->compute(SAB);
    intBB->compute(SBB);

    //SAB->print();
    //SBB->print();

    newfactory.reset();
    hybfactory.reset();
    intAB.reset();
    intBB.reset();
    pet.reset();

    // Constrained to the same symmetry at the moment, we can relax this soon
    SharedMatrix C_B(new Matrix("C_B", C_A->nirrep(),AO2USO->colspi(), noccpi));

    // Block over irreps (soon united irreps)
    for (int h = 0; h < C_A->nirrep(); h++) {

        int nocc = noccpi[h];
        int na = C_A->rowspi()[h];
        int nb = AO2USO->colspi()[h];

        if (nocc == 0 || na == 0 || nb == 0) continue;

        double** Ca = C_A->pointer(h);
        double** Cb = C_B->pointer(h);
        double** Sab = SAB->pointer(h);
        double** Sbb = SBB->pointer(h);

        int CholError = C_DPOTRF('L',nb,Sbb[0],nb);
        if (CholError !=0 )
            throw std::domain_error("S_BB Matrix Cholesky failed!");

        //Inversion (in place)
        int IError = C_DPOTRI('L',nb,Sbb[0],nb);
        if (IError !=0 )
            throw std::domain_error("S_BB Inversion Failed!");

        //LAPACK is smart and all, only uses half of the thing
        for (int m = 0; m<nb; m++)
            for (int n = 0; n<m; n++)
                Sbb[m][n] = Sbb[n][m];

        //Form T
        double** Temp1 = block_matrix(nb,nocc);
        C_DGEMM('T','N',nb,nocc,na,1.0,Sab[0],nb,Ca[0],nocc,0.0,Temp1[0],nocc);

        //fprintf(outfile," Temp1:\n");
        //print_mat(Temp1,nb,nocc,outfile);

        double** Temp2 = block_matrix(nb,nocc);
        C_DGEMM('N','N',nb,nocc,nb,1.0,Sbb[0],nb,Temp1[0],nocc,0.0,Temp2[0],nocc);

        //fprintf(outfile," Temp2:\n");
        //print_mat(Temp2,nb,nocc,outfile);

        double** Temp3 = block_matrix(na,nocc);
        C_DGEMM('N','N',na,nocc,nb,1.0,Sab[0],nb,Temp2[0],nocc,0.0,Temp3[0],nocc);

        //fprintf(outfile," Temp3:\n");
        //print_mat(Temp3,na,nocc,outfile);

        double** T = block_matrix(nocc,nocc);
        C_DGEMM('T','N',nocc,nocc,na,1.0,Ca[0],nocc,Temp3[0],nocc,0.0,T[0],nocc);

        //fprintf(outfile," T:\n");
        //print_mat(T,nocc,nocc,outfile);

        //Find T^-1/2
        // First, diagonalize T
        // the C_DSYEV call replaces the original matrix T with its eigenvectors
        double* eigval = init_array(nocc);
        int lwork = nocc * 3;
        double* work = init_array(lwork);
        int stat = C_DSYEV('v','u',nocc,T[0],nocc,eigval, work,lwork);
        if (stat != 0) {
            fprintf(outfile, "C_DSYEV failed\n");
            exit(PSI_RETURN_FAILURE);
        }
        free(work);

        // Now T contains the eigenvectors of the original T
        // Copy T to T_copy
        double **T_mhalf = block_matrix(nocc, nocc);
        double **T_copy = block_matrix(nocc, nocc);
        C_DCOPY(nocc*nocc,T[0],1,T_copy[0],1);

        // Now form T^{-1/2} = U(T)*t^{-1/2}*U,
        // where t^{-1/2} is the diagonal matrix of the inverse square roots
        // of the eigenvalues, and U is the matrix of eigenvectors of T
        for (int i=0; i<nocc; i++) {
            if (eigval[i] < 1.0E-10)
                eigval[i] = 0.0;
            else
                eigval[i] = 1.0 / sqrt(eigval[i]);

            // scale one set of eigenvectors by the diagonal elements t^{-1/2}
            C_DSCAL(nocc, eigval[i], T[i], 1);
        }
        free(eigval);

        // T_mhalf = T_copy(T) * T
        C_DGEMM('t','n',nocc,nocc,nocc,1.0,
                T_copy[0],nocc,T[0],nocc,0.0,T_mhalf[0],nocc);

        //Form CB
        C_DGEMM('N','N',nb,nocc,nocc,1.0,Temp2[0],nocc,T_mhalf[0],nocc,0.0,Cb[0],nocc);

        free_block(Temp1);
        free_block(Temp2);
        free_block(Temp3);
        free_block(T);
        free_block(T_copy);
        free_block(T_mhalf);

    }
    return C_B;
}

}}
