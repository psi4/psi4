#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libmints/view.h>
#include <libpsio/psio.hpp>
#include <libiwl/iwl.hpp>
#include <libtrans/integraltransform.h>
#include <libtrans/mospace.h>
#include <libdpd/dpd.h>
#include <libciomr/libciomr.h>
#include <libfock/apps.h>
#include <libqt/qt.h>
#include <vector>

#include <psifiles.h>

#include <fstream>
#include <algorithm>

#include <boost/python.hpp>
#include <boost/python/dict.hpp>

namespace psi{ namespace mrcc {

namespace {

void write_oei_to_disk(FILE* fort55, SharedMatrix moH)
{
    // Walk through moH and save the non-zero values
    int offset = 0;
    for (int h=0; h<moH->nirrep(); ++h) {
        for (int m=0; m<moH->rowdim(h); ++m) {
            for (int n=0; n<=m; ++n) {
                if (fabs(moH->get(h, m, n)) > 1.0e-12) {
                    fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", moH->get(h, m, n), m+offset+1, n+offset+1, 0, 0);
                }
            }
        }
        offset += moH->rowdim(h);
    }
}

void write_tei_to_disk(FILE* fort55, int nirrep, dpdbuf4& K, double ints_tolerance)
{
    for(int h = 0; h < nirrep; ++h){
        dpd_buf4_mat_irrep_init(&K, h);
        dpd_buf4_mat_irrep_rd(&K, h);
        for(int pq = 0; pq < K.params->rowtot[h]; ++pq){
            int p = K.params->roworb[h][pq][0];
            int q = K.params->roworb[h][pq][1];
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                int r = K.params->colorb[h][rs][0];
                int s = K.params->colorb[h][rs][1];

                if (fabs(K.matrix[h][pq][rs]) > ints_tolerance)
                    fprintf(fort55, "%28.20E%4d%4d%4d%4d\n",
                            K.matrix[h][pq][rs], p+1, q+1, r+1, s+1);
            }
        }
        dpd_buf4_mat_irrep_close(&K, h);
    }
}

void print_dim(const std::string& name, const Dimension& dim)
{
    fprintf(outfile, "        %15s [ ", name.c_str());
    for (int h=0; h<dim.n(); ++h) {
        fprintf(outfile, "%3d", dim[h]);
        if (h != dim.n()-1)
            fprintf(outfile, ",");
    }
    fprintf(outfile, "]\n");
}

class DPDFillerFunctor {
private:
    dpdfile4 *file_;
    dpdparams4 *params_;
    int **bucket_map_;
    int **bucket_offset_;
    bool symmetrize_;
    bool have_bra_ket_sym_;
public:
    DPDFillerFunctor(dpdfile4 *file, int **bucket_map,
                     int **bucket_offset, bool symmetrize, bool have_bra_ket_sym):
        file_(file), bucket_map_(bucket_map),
        bucket_offset_(bucket_offset), symmetrize_(symmetrize),
        have_bra_ket_sym_(have_bra_ket_sym)
    {
        params_ = file_->params;
    }
    void operator()(int this_bucket, int p, int q, int r, int s, double value)
    {
        if(symmetrize_){
            // Symmetrize the quantity (used in density matrix processing)
            if(p!=q) value *= 0.5;
            if(r!=s) value *= 0.5;
        }

        bool bra_ket_different = !(p==r && q==s);

        /* Get the orbital symmetries */
        int p_sym = params_->psym[p];
        int q_sym = params_->qsym[q];
        int r_sym = params_->rsym[r];
        int s_sym = params_->ssym[s];
        int pq_sym = p_sym^q_sym;
        int rs_sym = r_sym^s_sym;

        /* The allowed (Mulliken) permutations are very simple in this case */
        if(bucket_map_[p][q] == this_bucket) {

            /* Get the row and column indices and assign the value */
            int pq = params_->rowidx[p][q];
            int rs = params_->colidx[r][s];
            int offset = bucket_offset_[this_bucket][pq_sym];

            if((pq-offset >= params_->rowtot[pq_sym]) || (rs >= params_->coltot[rs_sym]))
                error("MP Params_make: pq, rs", p,q,r,s,pq,rs,pq_sym,rs_sym);
            file_->matrix[pq_sym][pq-offset][rs] += value;
        }

        /*
         * We also add in the bra-ket transposed value, as a result of the matrix
         * storage, but we need to make sure we don't duplicate "diagonal" values.
         * We don't do this if the quantity does not have bra-ket symmetry, like
         * in the Alpha-Beta TPDM.
         */
        if(bucket_map_[r][s] == this_bucket && bra_ket_different && have_bra_ket_sym_) {
            int rs = params_->rowidx[r][s];
            int pq = params_->colidx[p][q];
            int offset = bucket_offset_[this_bucket][rs_sym];
            if((rs-offset >= params_->rowtot[rs_sym])||(pq >= params_->coltot[pq_sym]))
                error("MP Params_make: rs, pq", p,q,r,s,rs,pq,rs_sym,pq_sym);
            file_->matrix[rs_sym][rs-offset][pq] += value;
        }
    }

private:
    void error(const char *message, int p, int q, int r, int s,
               int pq, int rs, int pq_sym, int rs_sym)
    {

        fprintf(outfile, "\n\tDPD Parameter Error in %s\n", message);
        fprintf(outfile,"\t-------------------------------------------------\n");
        fprintf(outfile,"\t    p      q      r      s  [   pq]  [   rs] pq_symm rs_symm\n");
        fprintf(outfile,"\t%5d  %5d  %5d  %5d  [%5d]  [%5d]   %1d   %1d\n", p,q,r,s,
                pq,rs,pq_sym,rs_sym);
        throw PsiException("DPD idx failure.", __FILE__, __LINE__);
    }
};

class MRCCRestrictedReader
{
    enum { line_length = 45 };

    FILE* ccdensities_;
    const double tolerance_;
    char *batch_;

    SharedMatrix one_particle_;

    off_t opdm_start_;

    int* abs_mo_to_rel_;
    int* abs_mo_to_irrep_;

public:
    MRCCRestrictedReader(FILE* ccdensities, const double tolerance, SharedMatrix one_particle)
        : ccdensities_(ccdensities), tolerance_(tolerance), batch_(0), one_particle_(one_particle)
    {
        batch_ = new char[line_length*1000+1];

        const Dimension& nmopi = one_particle_->rowspi();
        int nmo = nmopi.sum();
        abs_mo_to_rel_ = new int[nmo];
        abs_mo_to_irrep_ = new int[nmo];

        int count=0;
        for (int h=0; h<nmopi.n(); ++h) {
            for (int i=0; i<nmopi[h]; ++i) {
                abs_mo_to_rel_[count] = i;
                abs_mo_to_irrep_[count] = h;
                count++;
            }
        }
    }
    ~MRCCRestrictedReader() {
        delete[] abs_mo_to_irrep_;
        delete[] abs_mo_to_rel_;
        delete[] batch_;
    }

    template<typename Filler>
    void operator()(Filler& filler, int bucket) {
        // ensure we're at the beginning.
        fseek(ccdensities_, 0, SEEK_CUR);

        // each line in CCDENSITIES is 45 characters long.
        // read in a batch of 1000 lines; add one for '\0'.
        char *batch = new char[line_length*1000+1];

        double value;
        int p, q, r, s;

        size_t readin=0;
        size_t loops=0;
        while ((readin = fread(batch, line_length, 1000, ccdensities_))  ) {
            off_t offset = 0;

            for (size_t i=0; i<readin; ++i) {
                if (sscanf(batch+offset, "%lE %d %d %d %d\n", &value, &p, &q, &r, &s) != 5) {
                    std::string line(batch+offset, line_length);
                    fprintf(stderr, "Malformed line: %s\n", line.c_str());
                    throw PSIEXCEPTION("MRCC interface: Unable to interpret line.");
                }

                // The density arrives in Dirac notation, so we reorder to Mulliken
                // It's also normalized differently to Psi's, by a factor of 2.
                if (r != 0 && s != 0) {
                    if (p >= r && q >= s)
                        if (fabs(value) > tolerance_)
                            filler(bucket, p - 1, r - 1, q - 1, s - 1, value * 0.5);
                }
                else
                    one_particle_->set(abs_mo_to_irrep_[p-1], abs_mo_to_rel_[p-1], abs_mo_to_rel_[q-1], value);

                offset += line_length;
                loops++;
            }
        }
    }
};

class DPDBucketFiller
{
    dpdfile4 *I_;

    psio_address next_;

    int nbucket_;
    int **bucket_map_;
    int **bucket_offset_;
    int **bucket_row_dim_;
    int **bucket_size_;
public:
    DPDBucketFiller(dpdfile4* I, size_t memory_limit)
        : I_(I), next_(PSIO_ZERO)
    {
        memory_limit /= sizeof(double);
        int nirrep = I_->params->nirreps;
        int nump = 0, numq = 0;
        for(int h=0; h < nirrep; ++h){
            nump += I_->params->ppi[h];
            numq += I_->params->qpi[h];
        }
        bucket_map_ = init_int_matrix(nump, numq);

        /* Room for one bucket to begin with */
        bucket_offset_ = (int **) malloc(sizeof(int *));
        bucket_offset_[0] = init_int_array(nirrep);
        bucket_row_dim_ = (int **) malloc(sizeof(int *));
        bucket_row_dim_[0] = init_int_array(nirrep);
        bucket_size_ = (int **) malloc(sizeof(int *));
        bucket_size_[0] = init_int_array(nirrep);

        /* Figure out how many passes we need and where each p,q goes */
        nbucket_ = 1;
        for(int h = 0; h < nirrep; ++h){
            size_t row_length = (size_t) I_->params->coltot[h^(I_->my_irrep)];
            for(int row=0; row < I_->params->rowtot[h]; ++row) {
                if((memory_limit - row_length) >= 0){  // <-- This is aways true (unsigned - unsigned >= 0)
                    memory_limit -= row_length;
                    bucket_row_dim_[nbucket_-1][h]++;
                    bucket_size_[nbucket_-1][h] += row_length;
                }
                else {
                    nbucket_++;
                    memory_limit = memory_limit - row_length;
                    /* Make room for another bucket */
                    bucket_offset_ = (int **) realloc((void *) bucket_offset_,
                                                 nbucket_ * sizeof(int *));
                    bucket_offset_[nbucket_-1] = init_int_array(nirrep);
                    bucket_offset_[nbucket_-1][h] = row;

                    bucket_row_dim_ = (int **) realloc((void *) bucket_row_dim_,
                                                 nbucket_ * sizeof(int *));
                    bucket_row_dim_[nbucket_-1] = init_int_array(nirrep);
                    bucket_row_dim_[nbucket_-1][h] = 1;

                    bucket_size_ = (int **) realloc((void *) bucket_size_,
                                                    nbucket_ * sizeof(int *));
                    bucket_size_[nbucket_-1] = init_int_array(nirrep);
                    bucket_size_[nbucket_-1][h] = row_length;
                }
                int p = I_->params->roworb[h][row][0];
                int q = I_->params->roworb[h][row][1];
                bucket_map_[p][q] = nbucket_ - 1;
            }
        }
    }

    ~DPDBucketFiller() {
        free_int_matrix(bucket_map_);

        for(int n=0; n < nbucket_; ++n) {
            free(bucket_offset_[n]);
            free(bucket_row_dim_[n]);
            free(bucket_size_[n]);
        }
        free(bucket_offset_);
        free(bucket_row_dim_);
        free(bucket_size_);
    }

    DPDFillerFunctor dpd_filler(bool symmetrize, bool have_bra_ket_sym) {
        return DPDFillerFunctor(I_, bucket_map_, bucket_offset_, symmetrize, have_bra_ket_sym);
    }

    template<typename IntegralProcessor>
    void operator()(IntegralProcessor& mrccreader) {
        DPDFillerFunctor filler = dpd_filler(false, false);
        next_ = PSIO_ZERO;
        for(int n=0; n < nbucket_; ++n) { /* nbuckets = number of passes */
            /* Prepare target matrix */
            for(int h=0; h < I_->params->nirreps; h++) {
                I_->matrix[h] = block_matrix(bucket_row_dim_[n][h], I_->params->coltot[h]);
            }

            mrccreader(filler, n);

            for(int h=0; h < I_->params->nirreps; ++h) {
                if(bucket_size_[n][h])
                    _default_psio_lib_->write(I_->filenum, I_->label, (char *) I_->matrix[h][0],
                                              bucket_size_[n][h]*((long int) sizeof(double)), next_, &next_);
                free_block(I_->matrix[h]);
            }
        } /* end loop over buckets/passes */
    }
};

/**
  * Loads a RHF TPDM from CCDENSITIES into an IWL buffer.
  * \param ccdensities File handle to process
  * \param tolerance Any values below this are neglected
  * \param active_mopi Dimension object of the MOs per irrep, needed to form OPDM matrix
  * \param ints IntegralTransform object needed to determine DPD ID numbers.
  */
void load_restricted(FILE* ccdensities, double tolerance, const Dimension& active_mopi, IntegralTransform& ints)
{
    // => Sizing <= //

    int debug = Process::environment.options.get_int("DEBUG");

    boost::shared_ptr<Wavefunction> ref = Process::environment.reference_wavefunction();
    Dimension focc = ref->frzcpi();
    Dimension docc = ref->doccpi();
    Dimension nmopi = ref->nmopi();
    Dimension fvir = ref->frzvpi();
    Dimension aocc = docc - focc;
    Dimension avir = nmopi - fvir - docc;

    // => Read results from MRCC <= //
    
    dpdfile4 I;
    _default_psio_lib_->open(PSIF_TPDM_PRESORT, PSIO_OPEN_NEW);

    // Just in case the buffer exists, delete it now
    dpdbuf4 Ibuf;
    dpd_buf4_init(&Ibuf, PSIF_TPDM_PRESORT, 0, ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[A>=A]+"),
                  ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[A>=A]+"), 0, "MO TPDM (AA|AA)");
    dpd_buf4_scm(&Ibuf, 0.0);
    dpd_buf4_close(&Ibuf);


    dpd_file4_init(&I, PSIF_TPDM_PRESORT, 0,
                   ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[A>=A]+"), "MO TPDM (AA|AA)");

    SharedMatrix one_particle(new Matrix("MO-basis OPDM", active_mopi, active_mopi));

    DPDBucketFiller bucket(&I, Process::environment.get_memory());
    // While MRCCRestricedReader is processing the two particle values it will
    // go ahead and handle the one particle
    MRCCRestrictedReader mrccreader(ccdensities, tolerance, one_particle);

    // Read the density matrices into DPD buffers
    bucket(mrccreader);

    // Close DPD file
    dpd_file4_close(&I);
    

    /****START OF HACK****/
    fprintf(outfile, "    Beginning integral transformation.\n");
    // This transforms everything (OEI and TEI)
    ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(ints.get_dpd_id());
    fprintf(outfile, "    Transformation complete.\n\n");
    /****END HACK****/

    /*
     * Form the coupled cluster Lagrangian, X
     * [eq. 39, Scheiner et al., JCP, 87 5361 (1987)]
     */

    // One-electron contribution: Xpq <- h_pr D_rq
    SharedMatrix H(new Matrix(PSIF_MO_FZC, nmopi, nmopi));
    // TODO make sure the density is frozen appropriately with frozen core
    H->load(_default_psio_lib_,PSIF_OEI);
    SharedMatrix X(new Matrix("X (1e contribution)", nmopi, nmopi));
    X->gemm(false, false, 1.0, one_particle, H, 0.0);

    // Two-electron contribution: Xpq <- 2 (pr|st) G_qrst
    dpdbuf4 G;
    dpdbuf4 D;    
    dpdfile2 X2;


    _default_psio_lib_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    _default_psio_lib_->tocprint(PSIF_LIBTRANS_DPD);
    dpd_buf4_init(&G,PSIF_LIBTRANS_DPD, 0, ints.DPD_ID("[A,A]"), ints.DPD_ID("[A,A]"),
                  ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
    dpd_buf4_init(&D,PSIF_TPDM_PRESORT, 0, ints.DPD_ID("[A,A]"), ints.DPD_ID("[A,A]"),
                  ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[A>=A]+"), 0, "MO TPDM (AA|AA)");
    dpd_file2_init(&X2, PSIF_LIBTRANS_DPD, 0, ints.DPD_ID('A'), ints.DPD_ID('A'),
                   "X (2e contribution)");

    // Check energy
    double enuc = Process::environment.molecule()->nuclear_repulsion_energy();
    double E1e = one_particle->vector_dot(H);
    double E2e = dpd_buf4_dot(&G, &D);
    fprintf(outfile, "\tEnergies recomputed from MRCC's density matrices:\n");
    fprintf(outfile, "\t\tOne-electron energy = %16.10f\n",E1e);
    fprintf(outfile, "\t\tTwo-electron energy = %16.10f\n",E2e);
    fprintf(outfile, "\t\tTotal energy        = %16.10f\n",enuc + E1e + E2e);

    dpd_contract442(&G, &D, &X2, 0, 0, 2.0, 0.0);
    SharedMatrix X2mat(new Matrix(&X2));
X->print();
X2mat->print();
    X->add(X2mat);
    X->set_name("Full X");
X->print();

    // Symmetrize X, to form the lagrangian
    SharedMatrix Lag(X->clone());
    Lag->add(X->transpose());
    Lag->set_name("Coupled cluster Lagrangian");
    Lag->print();

    // Dump the lagrangian to disk, then we're done with it!

    // Orbital lagrangian: Lxy = Xxy - Xyx
    SharedMatrix Lxy(X->clone());
    Lxy->set_name("Full orbital Lagrangian");
    Lxy->subtract(X->transpose());
    Lxy->print();

    dpd_buf4_close(&D);
    dpd_buf4_close(&G);

    // At this point, the OPDM is the response OPDM
    if (debug) {
        one_particle->print();
    } 

    // Form the energy-weighted OPDM 

    // Form the Lagrangian
    //SharedMatrix Lia(new Matrix("Lia", naocc, navir));
    //for (int h = 0; h < Lia->nirrep(); h++) {
    //    int ni = naocc[h];
    //    int na = navir[h];
    //    if (!ni || !na) continue;
    //    
    //    double** Liap = Lia->pointer(h);
    //    double** Wp = W->pointer(h);
    //    for (int i = 0; i < ni; i++) {
    //        for (int a = 0; a < na; a++) {
    //            Liap[i][a] = Wp[i + focc[h]][a + docc[h]]
    //        }
    //    }
    //}
    
    // Form the orbital response contributions to the relaxed OPDM
    View ia_view(one_particle, aocc, avir, focc, docc);
    SharedMatrix Pia = ia_view();
    Pia->set_name("Pia (MRCC OPDM ov Block)");

    if (debug) {
        Pia->print();
    }

    // Construct a RCPHF Object 
    boost::shared_ptr<RCPHF> cphf(new RCPHF());
    cphf->preiterations();

    // TODO: Add pre-CPHF A-matrix correction
    boost::shared_ptr<JK> jk = cphf->jk();

    // Task and solve orbital Z-Vector Equations
    std::map<std::string, SharedMatrix>& b = cphf->b();  
    b["Orbital Z-Vector"] = Pia;
    
    cphf->compute_energy();
        
    std::map<std::string, SharedMatrix>& x = cphf->x();  
    SharedMatrix Xia = x["Orbital Z-Vector"]; 

    if (debug) {
        Xia->print();
    }

    // TODO: Add post-CPHF A-matrix correction

    for (int h = 0; h < Xia->nirrep(); h++) {
        int naocc = Xia->rowspi()[h];
        int navir = Xia->colspi()[h];
        if (!naocc || !navir) continue;
        double** Piap = one_particle->pointer(h);
        double** Xiap = Xia->pointer(h);
        for (int i = 0; i < naocc; i++) {
            C_DCOPY(navir,Xiap[i],1,&Piap[focc[h]+i][docc[h]],1);
            C_DCOPY(navir,Xiap[i],1,&Piap[docc[h]][focc[h]+i],nmopi[h]);
        }
    }

    cphf->postiterations();

    if (debug) {
        one_particle->print();
    } 

    one_particle->save(_default_psio_lib_, PSIF_MO_OPDM, Matrix::Full);

    one_particle->print();

    _default_psio_lib_->close(PSIF_TPDM_PRESORT, 1);
    _default_psio_lib_->close(PSIF_LIBTRANS_DPD, 1);
}

} // namespace anonymous

PsiReturnType mrcc_load_ccdensities(Options& options, const boost::python::dict& level)
{
    tstart();

    fprintf(outfile, "  PSI4 interface to MRCC:\n");

    // Ensure the dict provided has everything we need.
    if (!level.has_key("method") ||
        !level.has_key("order") ||
        !level.has_key("fullname"))
        throw PSIEXCEPTION("MRCC interface: Provided dictionary is incomplete.");

    int method  = boost::python::extract<int>(level.get("method"));
    int exlevel = boost::python::extract<int>(level.get("order"));
    std::string fullname = boost::python::extract<std::string>(level.get("fullname"));
    bool pertcc = exlevel > 0 ? false : true;
    exlevel = abs(exlevel);

    fprintf(outfile, "    Loading gradient data for %s.\n\n", fullname.c_str());

    // Currently, only do RHF case
    if (options.get_str("REFERENCE") != "RHF")
        throw PSIEXCEPTION("MRCC: Gradient interface only coded for RHF reference.");

    // Check the reference.
    bool restricted = true;

    if (options.get_str("REFERENCE") == "UHF")
        restricted = false;

    if (pertcc && options.get_str("REFERENCE") == "ROHF")
        throw PSIEXCEPTION("Perturbative methods not implemented for ROHF references.");

    // Use libtrans to initialize DPD
    boost::shared_ptr<Wavefunction> wave = Process::environment.reference_wavefunction();
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);
    IntegralTransform ints(wave, spaces, restricted ? IntegralTransform::Restricted : IntegralTransform::Unrestricted);

    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(ints.get_dpd_id());

    // Obtain a single handle to the CCDENSITIES file
    FILE* ccdensities = fopen("CCDENSITIES", "r");
    if (ccdensities == NULL)
        throw PSIEXCEPTION("MRCC interface: Unable to open CCDENSITIES. Did MRCC finish successfully?");

    const Dimension active_mopi = wave->nmopi() - wave->frzcpi() - wave->frzvpi();
    if (restricted) {
        load_restricted(ccdensities, options.get_double("INTS_TOLERANCE"),
                        active_mopi, ints);
    }
//    else
//        load_unrestricted_tpdm(ccdensities);

    tstop();
    return Success;
}

PsiReturnType mrcc_generate_input(Options& options, const boost::python::dict& level)
{
    tstart();

    fprintf(outfile, "  PSI4 interface to MRCC:\n");

    // Ensure the dict provided has everything we need.
    if (!level.has_key("method") ||
        !level.has_key("order") ||
        !level.has_key("fullname"))
        throw PSIEXCEPTION("MRCC interface: Provided dictionary is incomplete.");

    int method  = boost::python::extract<int>(level.get("method"));
    int exlevel = boost::python::extract<int>(level.get("order"));
    std::string fullname = boost::python::extract<std::string>(level.get("fullname"));
    bool pertcc = exlevel > 0 ? false : true;
    exlevel = abs(exlevel);

    fprintf(outfile, "    Generating inputs for %s.\n\n", fullname.c_str());

    fprintf(outfile, "    Automatically determined settings:\n");
    fprintf(outfile, "        method %d\n        exlevel %d\n        fullname %s\n\n",
            method, exlevel, fullname.c_str());

    boost::shared_ptr<Wavefunction> wave     = Process::environment.reference_wavefunction();
    boost::shared_ptr<Molecule>     molecule = wave->molecule();

    // Orbitals spaces
    Dimension docc        = wave->doccpi();
    Dimension frzcpi      = wave->frzcpi();
    Dimension frzvpi      = wave->frzvpi();
    Dimension active_docc = docc - frzcpi;
    Dimension active_socc = wave->soccpi();
    Dimension active_mopi = wave->nmopi() - frzcpi - frzvpi;

    int nbf = active_mopi.sum();
    int nirrep = wave->nirrep();
    int nelectron = 2 * active_docc.sum() + active_socc.sum();

    fprintf(outfile, "    Orbital Information:\n\n");
    print_dim("Frozen Core", frzcpi);
    print_dim("Active DOCC", active_docc);
    print_dim("SOCC", active_socc);
    print_dim("Frozen Virtual", frzvpi);

    fprintf(outfile, "\n");
    print_dim("Total MOs", active_mopi);

    fprintf(outfile, "\n");
    FILE* fort55 = fopen("fort.55", "w");
    fprintf(fort55, "%22d%22d\n", nbf, nelectron);

    // Print out orbital symmetries
    int count=0;
    for (int h=0; h<active_mopi.n(); ++h) {
        for (int n=0; n<active_mopi[h]; ++n) {
            fprintf(fort55, "%22d", h+1);  // 1 based irrep ordering

            count++;
            if (count % 3 == 0) // We only want 3 per line
                fprintf(fort55, "\n");
        }
    }

    if (count % 3) // integrals start on their own line
        fprintf(fort55, "\n");

    fprintf(outfile, "    Beginning integral transformation.\n");

    // Define the orbital space of the MO integrals we need.
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);

    // Check the reference.
    bool restricted = true;

    if (options.get_str("REFERENCE") == "UHF")
        restricted = false;

    if (pertcc && options.get_str("REFERENCE") == "ROHF")
        throw PSIEXCEPTION("Perturbative methods not implemented for ROHF references.");

    // Create integral transformation object
    IntegralTransform ints(wave, spaces, restricted ? IntegralTransform::Restricted : IntegralTransform::Unrestricted);

    // This transforms everything (OEI and TEI)
    ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);

    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(ints.get_dpd_id());

    fprintf(outfile, "    Transformation complete.\n\n");
    fprintf(outfile, "  Generating fort.55 integral file...");

    double ints_tolerance = options.get_double("INTS_TOLERANCE");

    _default_psio_lib_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    dpdbuf4 K;

    fprintf(fort55, " 150000\n");

    // RHF
    if (restricted) {

        // We want only the permutationally unique integrals, hence [A>=A]+, see libtrans documenation for details
        dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0,
                      ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[A>=A]+"), // In memory
                      ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[A>=A]+"), // On disk
                      0,
                      "MO Ints (AA|AA)");
        write_tei_to_disk(fort55, nirrep, K, ints_tolerance);
        dpd_buf4_close(&K);

        // Load in frozen core operator, in the event of FREEZE_CORE = FALSE this is the MO OEI
        SharedMatrix moH(new Matrix(PSIF_MO_FZC, wave->nmopi(), wave->nmopi()));
        moH->load(_default_psio_lib_, PSIF_OEI);
        View vmoH(moH, active_mopi, active_mopi, wave->frzcpi(), wave->frzcpi());
        moH = vmoH();
        write_oei_to_disk(fort55, moH);

        // Print nuclear repulsion energy.
        // Eventually needs to be changed to frozen core energy + nuclear repulsion energy
        fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", ints.get_frozen_core_energy() + molecule->nuclear_repulsion_energy(), 0, 0, 0, 0);
    }
    else {

        // We want only the permutationally unique integrals, hence [A>=A]+, see libtrans documenation for details

        // Load up alpha alpha
        dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0,
                      ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[A>=A]+"), // In memory
                      ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[A>=A]+"), // On disk
                      0,
                      "MO Ints (AA|AA)");
        write_tei_to_disk(fort55, nirrep, K, ints_tolerance);
        dpd_buf4_close(&K);

        // Write out separator
        fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", 0.0, 0, 0, 0, 0);

        // Load up beta beta
        dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0,
                      ints.DPD_ID("[a>=a]+"), ints.DPD_ID("[a>=a]+"), // In memory
                      ints.DPD_ID("[a>=a]+"), ints.DPD_ID("[a>=a]+"), // On disk
                      0,
                      "MO Ints (aa|aa)");
        write_tei_to_disk(fort55, nirrep, K, ints_tolerance);
        dpd_buf4_close(&K);

        // Write out separator
        fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", 0.0, 0, 0, 0, 0);

        // Load up alpha beta
        dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0,
                      ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[a>=a]+"), // In memory
                      ints.DPD_ID("[A>=A]+"), ints.DPD_ID("[a>=a]+"), // On disk
                      0,
                      "MO Ints (AA|aa)");
        write_tei_to_disk(fort55, nirrep, K, ints_tolerance);
        dpd_buf4_close(&K);

        // Write out separator
        fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", 0.0, 0, 0, 0, 0);

        // Load in alpha frozen core operator, in the event of FREEZE_CORE = FALSE this is the MO OEI
        SharedMatrix moH(new Matrix(PSIF_MO_A_FZC, wave->nmopi(), wave->nmopi()));
        moH->load(_default_psio_lib_, PSIF_OEI);
        View vmoH(moH, active_mopi, active_mopi, wave->frzcpi(), wave->frzcpi());
        moH = vmoH();
        write_oei_to_disk(fort55, moH);

        // Write out separator
        fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", 0.0, 0, 0, 0, 0);

        // Load in beta frozen core operator, in the event of FREEZE_CORE = FALSE this is the MO OEI
        SharedMatrix moHb(new Matrix(PSIF_MO_B_FZC, wave->nmopi(), wave->nmopi()));
        moHb->load(_default_psio_lib_, PSIF_OEI);
        View vmoHb(moHb, active_mopi, active_mopi, wave->frzcpi(), wave->frzcpi());
        moHb = vmoHb();
        write_oei_to_disk(fort55, moHb);

        // Write out separator
        fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", 0.0, 0, 0, 0, 0);

        // Print nuclear repulsion energy.
        // Eventually needs to be changed to frozen core energy + nuclear repulsion energy
        fprintf(fort55, "%28.20E%4d%4d%4d%4d\n", ints.get_frozen_core_energy() + molecule->nuclear_repulsion_energy(), 0, 0, 0, 0);
    }
    _default_psio_lib_->close(PSIF_LIBTRANS_DPD, 1);

    fclose(fort55);

    fprintf(outfile, "done.\n  Generating fort.56 input file...");

    // Determine energy convergence to pass to MRCC
    double user_e = fabs(Process::environment.options.get_double("E_CONVERGENCE"));
    int e_conv = 0;

    if (user_e >= 1.0)
        throw InputException("Your energy convergence criterea is absurb.", "E_CONVERGENCE", user_e, __FILE__, __LINE__);

    while (user_e <= 0.1) {
        user_e *= 10.0;
        e_conv++;
    }

    if (e_conv == 0)
        throw PSIEXCEPTION("e_conv became 0. This isn't right.");

    // The following two are expert options.
    if (options["MRCC_METHOD"].has_changed())
        method = options.get_int("MRCC_METHOD");
    if (options["MRCC_LEVEL"].has_changed())
        exlevel = options.get_int("MRCC_LEVEL");

    int nsing = 1;
    int closed_shell = 1;
    int spatial_orbitals = 1;
    int ndoub = 0;
    if (pertcc && (options.get_str("REFERENCE") == "ROHF" || options.get_str("REFERENCE") == "UHF")) {
        nsing = 0;
        closed_shell = 0;
        spatial_orbitals = 0;
        ndoub = 1;
    }
    if (options["MRCC_NUM_SINGLET_ROOTS"].has_changed())
        nsing = options.get_int("MRCC_NUM_SINGLET_ROOTS");

    int symm = 0;
    for (int h=0; h<nirrep; ++h)
        for (int n=0; n<active_socc[h]; ++n)
            symm ^= h;
    symm += 1; // stupid 1 based fortran

    FILE* fort56 = fopen("fort.56", "w");
    fprintf(fort56, "%6d%6d%6d%6d%6d      0     0%6d     0%6d%6d     1%6d      0      0%6d     0     0    0.00    0%6lu\n",
            exlevel,                                         // # 1
            nsing,                                           // # 2
            options.get_int("MRCC_NUM_TRIPLET_ROOTS"),       // # 3
            options.get_int("MRCC_RESTART"),                 // # 4
            method,                                          // # 5
            symm,                                            // # 8
            closed_shell,                                    // #10
            spatial_orbitals,                                // #11
            ndoub,                                           // #13
            e_conv,                                          // #16
            Process::environment.get_memory() / 1000 / 1000  // #21
        );

    fprintf(fort56, "ex.lev,nsing,ntrip, rest,CC/CI, dens,conver, symm, diag,   CS,spatial, HF,ndoub, nacto, nactv,  tol,maxex, sacc,   freq,dboc,  mem\n");

    for (int h=0; h < nirrep; ++h) {
        for (int nmo = 0; nmo < active_mopi[h]; ++nmo) {
            int val = 0;
            if (nmo < active_docc[h])
                val = 2;
            else if (nmo < (active_docc[h] + active_socc[h]))
                val = 1;

            fprintf(fort56, "%3d", val);
        }
    }
    fprintf(fort56, "\n");
    fclose(fort56);

    fprintf(outfile, "done.\n");
    fflush(outfile);

    tstop();

    return Success;
}

}} // End Namespaces

