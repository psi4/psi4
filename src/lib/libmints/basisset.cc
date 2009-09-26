#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.hpp>
#include <psifiles.h>

#include "basisset.h"
#include "integral.h"
#include "symmetry.h"
#include "factory.h"

using namespace psi;

BasisSet::BasisSet() :
    max_nprimitives_(0), shell_first_basis_function_(NULL), shell_first_ao_(NULL), shell_center_(NULL), 
    max_stability_index_(0), uso2ao_(NULL), simple_mat_uso2ao_(NULL), uso2bf_(NULL), 
    simple_mat_uso2bf_(NULL)
{ }

BasisSet::BasisSet(Chkpt* chkpt, std::string basiskey) :
    max_nprimitives_(0), shell_first_basis_function_(NULL), shell_first_ao_(NULL), shell_center_(NULL), 
    max_stability_index_(0), uso2ao_(NULL), simple_mat_uso2ao_(NULL), uso2bf_(NULL), 
    simple_mat_uso2bf_(NULL)
{
    // This requirement holds no matter what.
    puream_ = chkpt->rd_puream(basiskey.c_str()) ? true : false;

    // Initialize molecule, retrieves number of center and geometry
    molecule_ = new Molecule;
    molecule_->init_with_chkpt(chkpt);

    // Initialize the shells
    initialize_shells(chkpt, basiskey);
}

BasisSet::BasisSet(Chkpt* chkpt, std::string genbas_filename, std::string genbas_basis) :
    max_nprimitives_(0), shell_first_basis_function_(NULL), shell_first_ao_(NULL), shell_center_(NULL), 
    max_stability_index_(0), uso2ao_(NULL), simple_mat_uso2ao_(NULL), uso2bf_(NULL), 
    simple_mat_uso2bf_(NULL)
{
    // This requirement holds no matter what.
    puream_ = chkpt->rd_puream() ? true : false;

    // Initialize molecule, retrieves number of centers and geometry
    molecule_ = new Molecule;
    molecule_->init_with_chkpt(chkpt);

    if (!genbas_filename.empty() && !genbas_basis.empty()) {
        fprintf(outfile, "  Initializing BasisSet object with data in: %s (basis: %s)\n", 
            genbas_filename.c_str(), genbas_basis.c_str());
        initialize_shells_via_genbas(genbas_filename, genbas_basis);
    } else {
        fprintf(outfile, "  BasisSet: When using GENBAS you must provide basis set name.\n");
        exit(EXIT_FAILURE);
    }
}

BasisSet::~BasisSet()
{
    if (shell_first_basis_function_)
        Chkpt::free(shell_first_basis_function_);
    if (shell_first_ao_)
        Chkpt::free(shell_first_ao_);
    if (shell_center_)
        Chkpt::free(shell_center_);
    if (uso2ao_)
        Chkpt::free(uso2ao_);
    if (uso2bf_)
        Chkpt::free(uso2bf_);
    if (molecule_)
        delete molecule_;
    if (simple_mat_uso2ao_)
        delete simple_mat_uso2ao_;
    if (simple_mat_uso2bf_)
        delete simple_mat_uso2bf_;
    if (sotransform_)
        delete sotransform_;
    
    for  (int i=0; i < nshells_; ++i)
        delete shells_[i];
    delete[] shells_;
}

void BasisSet::initialize_shells(Chkpt *chkpt, std::string& basiskey)
{
    // Initialize some data from checkpoint.
    nshells_      = chkpt->rd_nshell(basiskey.c_str());
    nprimitives_  = chkpt->rd_nprim(basiskey.c_str());
    nao_          = chkpt->rd_nao(basiskey.c_str());

    // Psi3 only allows either all Cartesian or all Spherical harmonic
    nbf_          = chkpt->rd_nso(basiskey.c_str());
    max_am_       = chkpt->rd_max_am(basiskey.c_str());
    uso2ao_       = chkpt->rd_usotao(basiskey.c_str());
    uso2bf_       = chkpt->rd_usotbf(basiskey.c_str());

    simple_mat_uso2ao_ = new SimpleMatrix("Unique SO to AO transformation matrix", nbf_, nao_);
    simple_mat_uso2ao_->set(uso2ao_);
    // simple_mat_uso2ao_.print();
    
    simple_mat_uso2bf_ = new SimpleMatrix("Unique SO to BF transformation matrix", nbf_, nbf_);
    simple_mat_uso2bf_->set(uso2bf_);
    // simple_mat_uso2bf_.print();
    
    // Allocate memory for the shells
    shells_       = new GaussianShell*[nshells_];
    
    // Retrieve angular momentum of each shell (1=s, 2=p, ...)
    int *shell_am = chkpt->rd_stype(basiskey.c_str());
    
    // Retrieve number of primitives per shell
    int *shell_num_prims = chkpt->rd_snumg(basiskey.c_str());
    
    // Retrieve exponents of primitive Gaussians
    double *exponents = chkpt->rd_exps(basiskey.c_str());
    fprintf(outfile, "Exponents:\n");
    for (int i=0; i<nprimitives_; i++) {
      fprintf(outfile, "%f\n", exponents[i]);
    }
 
    // Retrieve coefficients of primitive Gaussian
    double **ccoeffs = chkpt->rd_contr_full(basiskey.c_str());
    
    // Retrieve pointer to first primitive in shell
    int *shell_fprim = chkpt->rd_sprim(basiskey.c_str());
    
    // Retrieve pointer to first basis function in shell
    shell_first_basis_function_ = chkpt->rd_sloc_new(basiskey.c_str());
    
    // Retrieve pointer to first AO in shell
    shell_first_ao_ = chkpt->rd_sloc(basiskey.c_str());
    
    // Retrieve location of shells (which atom it's centered on)
    shell_center_ = chkpt->rd_snuc(basiskey.c_str());
        
    // Initialize SphericalTransform
    for (int i=0; i<=max_am_; ++i) {
        sphericaltransforms_.push_back(SphericalTransform(i));
    }
    
    // Initialize SOTransform
    sotransform_ = new SOTransform;
    sotransform_->init(nshells_);
    
    int *so2symblk = new int[nbf_];
    int *so2index  = new int[nbf_];
    int *sopi = chkpt->rd_sopi(basiskey.c_str());
    int nirreps = chkpt->rd_nirreps();
    
    // Create so2symblk and so2index
    int ij = 0; int offset = 0;
    for (int h=0; h<nirreps; ++h) {
        for (int i=0; i<sopi[h]; ++i) {
            so2symblk[ij] = h;
            so2index[ij] = ij-offset;
            
            ij++;
        }
        offset += sopi[h];
    }
    
    // Currently all basis sets are treated as segmented contractions
    // even though GaussianShell is generalized (well not really).
    int ncontr = 1;
    int ao_start = 0;
    int puream_start = 0;
    
    for (int i=0; i<nshells_; ++i) {
        int *am = new int[ncontr];
        am[0] = shell_am[i] - 1;
        int fprim = shell_fprim[i] - 1;
        int nprims = shell_num_prims[i];
        Vector3 center = molecule_->xyz(shell_center_[i] - 1);
        double **cc = new double*[nprims];
        for (int p=0; p<nprims; ++p) {
            cc[p] = new double[ncontr];
            cc[p][0] = ccoeffs[fprim+p][am[0]];
        }
        
        // Construct a new shell. GaussianShell copies the data to new memory
        shells_[i] = new GaussianShell(ncontr, nprims, &(exponents[fprim]), am, 
            puream_ ? GaussianShell::Pure : GaussianShell::Cartesian, cc, shell_center_[i]-1, center,
            puream_start);
            
        if (nprims > max_nprimitives_)
            max_nprimitives_ = nprims;
            
        for (int p=0; p<nprims; p++) {
            delete[] cc[p];
        }
        delete[] cc;
        
        // OK, for a given number of AO functions in a shell INT_NCART(am)
        // beginning at column ao_start go through all rows finding where this
        // AO function contributes to an SO.
        for (int ao = 0; ao < INT_NCART(am[0]); ++ao) {
            int aooffset = ao_start + ao;
            for (int so = 0; so < nbf_; ++so) {
                if (fabs(uso2ao_[so][aooffset]) >= 1.0e-14)
                    sotransform_->add_transform(i, so2symblk[so], so2index[so], uso2ao_[so][aooffset], ao, so);
            }
        }        
        
        // Shift the ao starting index over to the next shell
        ao_start += INT_NCART(am[0]);
        puream_start += INT_NFUNC(puream_, am[0]);
        delete[] am;
    }
        
    delete[] so2symblk;
    delete[] so2index;
    Chkpt::free(sopi);
    Chkpt::free(ccoeffs);
    Chkpt::free(exponents);
    Chkpt::free(shell_am);
    Chkpt::free(shell_num_prims);
    Chkpt::free(shell_fprim);
}

void BasisSet::initialize_shells_via_genbas(std::string& genbas_filename, std::string& genbas_basis)
{
    char search_string[80];
    char line[120];
    int len, j, tval, k, l;
    int totalAM;
    int *num_primitives, *num_contractedto;
    Vector3 center;
    int puream_start = 0;
    int ao_start = 0;
    
    // Initial maximum angular momentum
    max_am_ = 0;
    // Initial maximum number of primitives
    max_nprimitives_ = 0;
    
    // Temporary storage for shells, will be moved into shells_ once we know the total number.
    std::vector<GaussianShell* > shells;
    
    // Attempt to open user specific genbas file
    FILE *genbas = fopen(genbas_filename.c_str(), "r");
    if (genbas == NULL) {
        fprintf(outfile, " Unable to open %s.\n", genbas_filename.c_str());
        exit(EXIT_FAILURE);
    }
    
    // Loop through each atom and load in basis set from genbas
    for (int i=0; i < molecule_->natom(); ++i) {
        rewind(genbas);
        
        // Retrieve center
        center = molecule_->xyz(i);
        
        // Construct search string
        sprintf(search_string, "%s:%s", molecule_->label(i).c_str(), genbas_basis.c_str());
        
        // Check for empty file.
        if (fgets(line, 120, genbas) == NULL)
        {
            fprintf(stderr, "%s is empty!\n", genbas_filename.c_str());
            exit(1);
        }
        len = strlen(search_string);

        // Search the file for "search_string"
        while (strncasecmp(line, search_string, len) != 0)
        {
            if (fgets(line, 120, genbas) == NULL)
            {
                fprintf(stderr, "End of file encountered! (%s)\n", search_string);
                exit(1);
            }
        }
        
        // If we make it here we must have found it.
        // get the comment line
        fgets(line, 120, genbas);

        // Get some useful information
        fscanf(genbas, "%d", &totalAM);

        // Save maximum angular momentum, needed to create sphericaltransforms_
        if (totalAM - 1 > max_am_)
            max_am_ = totalAM - 1;
            
        // Create the arrays num_primitives, num_contractedto
        num_primitives = new int[totalAM];
        num_contractedto = new int[totalAM];

        // Now just loop on the numbers that define the angular momentum, I'm just assuming it
        // goes in order
        for (j = 0; j < totalAM; j++)
            fscanf(genbas, "%d", &tval);
        
        // Now start reading in the number of actual functions for each angular momentum
        for (j = 0; j < totalAM; j++)
            fscanf(genbas, "%d", &(num_contractedto[j]));

        // Now start reading in the number of primitives for each angular momentum
        for (j = 0; j < totalAM; j++)
            fscanf(genbas, "%d", &(num_primitives[j]));

        // Read in the primitives and contractions for each AM
        for (j = 0; j < totalAM; j++) {
            double *exponents = new double[num_primitives[j]];
            double **contractions = init_matrix(num_primitives[j], num_contractedto[j]);

            // Read in the exponents
            for (k = 0; k < num_primitives[j]; k++) 
                fscanf(genbas, "%lf", &(exponents[k]));

            // Read in the coefficients
            for (k = 0; k < num_primitives[j]; k++) {
                for (l = 0; l < num_contractedto[j]; l++) {
                    fscanf(genbas, "%lf", &(contractions[k][l]));
                }
            }

            // Go through and create the needed Gaussian and Contracted Gaussians
            for (l = 0; l < num_contractedto[j]; l++) {
                k = 0;
                while (contractions[k][l] == 0.0) k++; // Find the first non-zero contraction

                // Count how many coefficients there are in this contraction.
                int kount = 0;
                for (int m=k; m<num_primitives[j]; ++m)
                    if (fabs(contractions[m][l]) > 0.0) 
                        kount++;
                
                // Okay, we know how many there are allocate memory to hold information for GaussianShell creation.
                int ncontr = 1;
                int *am = new int[ncontr];
                am[0] = j;
                double **cc = new double*[kount];
                double *e = new double[kount];
                for (int p=0; p<kount; ++p) {
                    cc[p] = new double[ncontr];
                    cc[p][0] = contractions[k+p][l];
                    e[p] = exponents[k+p];
                }
                
                GaussianShell* pshell = new GaussianShell(ncontr, kount, e, am, puream_ ? GaussianShell::Pure : GaussianShell::Cartesian,
                                         cc, i, center, puream_start, GaussianShell::Unnormalized);
                shells.push_back(pshell);
                
                if (kount > max_nprimitives_)
                    max_nprimitives_ = kount;

                for (int p=0; p<kount; ++p) {
                    delete[] cc[p];
                }
                delete[] cc;
                delete[] e;
                ao_start += INT_NCART(am[0]);
                puream_start += INT_NFUNC(puream_, am[0]);
                delete[] am;
            }

            delete[] exponents;
            free_matrix(contractions, num_primitives[j]);
        }

        // Clean up from this atom
        delete[] num_primitives;
        delete[] num_contractedto;
    }
    
    // Close the file
    fclose(genbas);

    // Initialize SphericalTransform
    for (int i=0; i<=max_am_; ++i) {
        sphericaltransforms_.push_back(SphericalTransform(i));
    }

    // Counters
    nao_ = ao_start;
    nbf_ = puream_start;
    
    // Okay, we now know how many shells we have, allocate permanent home for it.
    nshells_ = shells.size();
    shells_ = new GaussianShell*[nshells_];
    
    // Reading basis sets from GENBAS means that symmetry must be ignored. But we
    // must have an sotransform_ object created. Create an identity matrix that will
    // server as our transform.
    uso2ao_ = Chkpt::matrix<double>(nbf_, nao_);
    for (int i=0; i<nbf_; ++i) {
        uso2ao_[i][i] = 1.0;
    }
    
    // Certain information could not be created unless we knew the number of shells.
    // Now that we do know create all the necessary information for integral computations.
    // Retrieve location of shells (which atom it's centered on)
    shell_center_ = new int[nshells_];

    // Initialize SOTransform
    sotransform_ = new SOTransform;
    sotransform_->init(nshells_);
    
    // Reset ao_start
    ao_start = 0;
    for (int i=0; i<nshells_; ++i) {
        shells_[i] = shells[i];
        
        // Save which center each shell is on.
        shell_center_[i] = shells[i]->ncenter();
        
        // OK, for a given number of AO functions in a shell INT_NCART(am)
        // beginning at column ao_start go through all rows finding where this
        // AO function contributes to an SO.
        for (int ao = 0; ao < INT_NCART(shells[i]->am(0)); ++ao) {
            int aooffset = ao_start + ao;
            for (int so = 0; so < nbf_; ++so) {
                if (fabs(uso2ao_[so][aooffset]) >= 1.0e-14)
                    // 0 = always A1, so always equals index
                    sotransform_->add_transform(i, 0, aooffset, uso2ao_[so][aooffset], ao, so);
            }
        }
        ao_start += INT_NCART(shells[i]->am(0));
    }
}

void BasisSet::print(FILE *out) const
{
    fprintf(out, "  Basis Set\n");
    fprintf(out, "    Number of shells: %d\n", nshell());
    fprintf(out, "    Number of basis function: %d\n", nbf());
    fprintf(out, "    Number of Cartesian functions: %d\n", nao());
    fprintf(out, "    Spherical Harmonics?: %s\n", has_puream() ? "true" : "false");
    fprintf(out, "    Max angular momentum: %d\n\n", max_am());
    
    fprintf(out, "    Shells:\n\n");
    for (int s=0; s<nshell(); ++s)
        shells_[s]->print(out);
}

GaussianShell* BasisSet::shell(int si) const
{
    #ifdef DEBUG
    assert(si < nshell());
    #endif
    return shells_[si];
}

BasisSet* BasisSet::zero_basis_set()
{
    BasisSet *new_basis = new BasisSet();

    // Setup all the parameters needed for a zero basis set
    new_basis->shell_first_basis_function_ = NULL;
    new_basis->shell_first_ao_ = NULL;
    new_basis->shell_center_ = new int[1];
    new_basis->shell_center_[0] = 0;
    new_basis->uso2ao_ = NULL;
    new_basis->max_nprimitives_ = 1;
    new_basis->max_stability_index_ = 0;
    new_basis->max_am_ = 0;

    new_basis->puream_ = false;

    // Add "ghost" atom to the molecule for this basis
    new_basis->molecule_ = new Molecule;
    new_basis->molecule_->add_atom(0, 0.0, 0.0, 0.0);
    Vector3 center = new_basis->molecule_->xyz(0);

    new_basis->nshells_ = 1;
    new_basis->nprimitives_ = 1;
    new_basis->nao_ = 1;
    new_basis->nbf_ = 1;
    new_basis->uso2ao_ = Chkpt::matrix<double>(1, 1);
    new_basis->uso2ao_[0][0] = 1.0;
    new_basis->uso2bf_ = Chkpt::matrix<double>(1, 1);
    new_basis->uso2bf_[0][0] = 1.0;
 
    // Create shell array
    new_basis->shells_ = new GaussianShell*[1];

    // Spherical and SO-transforms are expected even if not used.
    new_basis->sphericaltransforms_.push_back(SphericalTransform(0));
    new_basis->sotransform_ = new SOTransform;
    new_basis->sotransform_->init(1);
    
    // Create out basis set arrays
    int *am = new int[1];
    double *e = new double[1];
    double **c = new double*[1];
    c[0] = new double[1];

    // null basis set
    am[0] = 0;
    e[0] = 0.0;
    c[0][0] = 1.0;

    // Add the null-s-function
    new_basis->shells_[0] = new GaussianShell(1, 1, e, am, GaussianShell::Cartesian, c, 0, center, 0, GaussianShell::Normalized);

    // Delete the basis set arrays, GaussianShell makes internal copies of everything
    delete[] am;
    delete[] e;
    delete[] c[0];
    delete[] c;

    // Add s-function SO transform.
    new_basis->sotransform_->add_transform(0, 0, 0, 1.0, 0, 0);

    return new_basis;
}

