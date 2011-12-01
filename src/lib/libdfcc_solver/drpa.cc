#include "drpa.h"

#include <time.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libmints/mints.h>
#include <utility> 

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace boost;
using namespace psi;
using namespace std;

namespace psi { namespace dfcc {

dRPA::dRPA(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt)
  : CC(options, psio, chkpt)
{
  print_header();
  psio_->open(DFCC_INT_FILE,PSIO_OPEN_NEW);
  df_integrals();
}

dRPA::~dRPA()
{
  psio_->close(DFCC_INT_FILE,1);
}

double dRPA::compute_energy()
{
  double energy;

  if (options_.get_str("RPA_ALGORITHM") == "DF")
    energy = df_compute_energy();
  if (options_.get_str("RPA_ALGORITHM") == "CD")
    energy = cd_compute_energy();

  return (energy);
}

double dRPA::df_compute_energy()
{
  time_t start = time(NULL);
  time_t stop;

  if (options_.get_bool("DIIS"))
    diis_ = boost::shared_ptr<DFCCDIIS>(new DFCCDIIS(DFCC_DIIS_FILE,naocc_*naocc_*
      navir_*navir_,options_.get_int("MAX_DIIS_VECS"),psio_));

  B_p_IA_ = block_matrix(naocc_*navir_,ndf_);
  Th_p_IA_ = block_matrix(naocc_*navir_,ndf_);
  tIAJB_ = init_array((long int) naocc_*navir_*naocc_*navir_);
  t2IAJB_ = init_array((long int) naocc_*navir_*naocc_*navir_);
  xIAJB_ = init_array((long int) naocc_*navir_*naocc_*navir_);

  psio_->read_entry(DFCC_INT_FILE,"OV DF Integrals",(char *) &(B_p_IA_[0][0]),
    sizeof(double)*naocc_*navir_*ndf_);

  C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,ndf_,1.0,B_p_IA_[0],ndf_,
    B_p_IA_[0],ndf_,0.0,tIAJB_,naocc_*navir_);

  apply_denom();
  C_DCOPY((long int) naocc_*navir_*naocc_*navir_,tIAJB_,1,t2IAJB_,1);

  C_DGEMM('N','N',naocc_*navir_,ndf_,naocc_*navir_,1.0,tIAJB_,naocc_*navir_,
    B_p_IA_[0],ndf_,0.0,Th_p_IA_[0],ndf_); 

  double e_mp2j = df_energy();

  fprintf(outfile,"  Iter       Energy (H)          dE (H)             RMS (H)     Time (s)\n");
  fflush(outfile);

  int iter = 1;
  int done = 0;
  double e_new;
  double e_old = e_mp2j;
  double rms = 1.0;

  do {
    C_DGEMM('N','N',naocc_*navir_,ndf_,naocc_*navir_,1.0,tIAJB_,naocc_*navir_,
      B_p_IA_[0],ndf_,0.0,Th_p_IA_[0],ndf_); 

    C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,ndf_,1.0,B_p_IA_[0],ndf_,
      B_p_IA_[0],ndf_,0.0,tIAJB_,naocc_*navir_);
    C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,ndf_,1.0,B_p_IA_[0],ndf_,
      Th_p_IA_[0],ndf_,1.0,tIAJB_,naocc_*navir_);
    C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,ndf_,1.0,Th_p_IA_[0],ndf_,
      B_p_IA_[0],ndf_,1.0,tIAJB_,naocc_*navir_);
    C_DGEMM('N','T',naocc_*navir_,naocc_*navir_,ndf_,1.0,Th_p_IA_[0],ndf_,
      Th_p_IA_[0],ndf_,1.0,tIAJB_,naocc_*navir_);
  
    apply_denom();
    
    C_DGEMM('N','N',naocc_*navir_,ndf_,naocc_*navir_,1.0,tIAJB_,naocc_*navir_,
      B_p_IA_[0],ndf_,0.0,Th_p_IA_[0],ndf_);
    
    e_new = df_energy();

    rms = df_store_error_vecs();

    stop = time(NULL);
    fprintf(outfile,"  %4d %16.8lf %17.9lf %17.9lf %12ld",iter,e_new,
      e_old-e_new,rms,stop-start);
    fflush(outfile);

    if (options_.get_int("MIN_DIIS_VECS") <= iter &&
        options_.get_bool("DIIS")) {
      diis_->get_new_vector(tIAJB_,xIAJB_);
      fprintf(outfile,"  DIIS\n");
    }
    else {
      fprintf(outfile,"\n");
    }
    fflush(outfile);

    C_DCOPY((long int) naocc_*navir_*naocc_*navir_,tIAJB_,1,t2IAJB_,1);
    iter++;

    if (iter > options_.get_int("MAXITER")) done = 1;
    if (fabs(e_old-e_new) < options_.get_double("E_CONVERGE")
      && rms < options_.get_double("T_CONVERGE")) done = 1;

    e_old = e_new;
  }
  while(!done);

  double scale = options_.get_double("RPA_ALPHA");

  fprintf(outfile,"\n");
  fprintf(outfile,"  Reference Energy                    %18.10lf [H]\n",Eref_);
  fprintf(outfile,"  Total CD-dRPA Energy                %18.10lf [H]\n",Eref_ + e_old);
  fprintf(outfile,"  DF-MP2J Energy                      %18.10lf [H]\n",e_mp2j);
  fprintf(outfile,"  DF-dRPA Energy                      %18.10lf [H]\n",e_old);
  fprintf(outfile,"  Delta DF-dRPA/MP2J Energy           %18.10lf [H]\n",e_old - e_mp2j);
  fprintf(outfile,"  Delta DF-dRPA/MP2J Scale            %18.10lf [-]\n",scale); 
  fprintf(outfile,"  Scaled Delta DF-dRPA/MP2J Energy    %18.10lf [H]\n",scale*(e_old - e_mp2j));

  Process::environment.globals["RPA TOTAL ENERGY"]        = Eref_ + e_old;
  Process::environment.globals["RPA MP2J ENERGY"]         = e_mp2j;
  Process::environment.globals["RPA DRPA ENERGY"]         = e_old;
  Process::environment.globals["RPA DELTA ENERGY"]        = e_old - e_mp2j;
  Process::environment.globals["RPA SCALED DELTA ENERGY"] = scale*(e_old - e_mp2j);
  Process::environment.globals["CURRENT ENERGY"]          = Eref_ + e_old;
 
  fprintf(outfile,"\n");
  fprintf(outfile,"  Reference Energy            %18.10lf\n",Eref_);
  fprintf(outfile,"  Correlation Energy          %18.10lf\n",e_old);
  fprintf(outfile,"  Total DF-dRPA Energy        %18.10lf\n\n",Eref_+e_old);

  free_block(B_p_IA_);
  free_block(Th_p_IA_);
  free(tIAJB_);
  free(xIAJB_);
  free(t2IAJB_);

  return(0.0);
}

double dRPA::cd_compute_energy()
{
    // ==> Initialization <== // 
    fprintf(outfile,"  CD-dRPA Algorithm Selected.\n\n");

    double PLUS_EPSILON = options_.get_double("RPA_PLUS_EPSILON");

    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif

    time_t start = time(NULL);
    time_t stop;

    int naux = ndf_;
    int nocc = naocc_;
    int nvir = navir_;
    ULI nov = nocc*(ULI)nvir;
 
    if (options_.get_bool("DIIS"))
      diis_ = boost::shared_ptr<DFCCDIIS>(new DFCCDIIS(DFCC_DIIS_FILE,nov*naux,
        options_.get_int("MAX_DIIS_VECS"),psio_));
  
    if (debug_) {
        evals_aocc_->print();
        evals_avir_->print();
    }

    // => DF Integrals <= //  
    SharedMatrix Qia(new Matrix("(Q|ia) Integrals", nov,naux));
    double** Qiap = Qia->pointer();
 
    // TODO run out of core with (Q|ia) striping
    psio_->read_entry(DFCC_INT_FILE,"OV DF Integrals",(char *)
        &(Qiap[0][0]),sizeof(double)*nov*naux);

    if (debug_)
        Qia->print();

    // => Initial ZiaQ tensor <= //
    SharedMatrix ZiaQ (new Matrix("Z_ia^Q Tensor", nov, naux)); 
    double** ZiaQp = ZiaQ->pointer();
    C_DCOPY(nov*naux,Qiap[0],1,ZiaQp[0],1);

    if (debug_)
        ZiaQ->print();

    // Restripe for the preferred (Qia) order
    Qia.reset();
    Qia = SharedMatrix(new Matrix("(Q|ia) Integrals", naux, nov));
    Qiap = Qia->pointer();

    for (int Q = 0; Q < naux; Q++) {
        C_DCOPY(nov, &ZiaQp[0][Q], naux, Qiap[Q], 1);
    }   

    if (debug_)
        Qia->print();

    // => Iteration Control <= //
    int iter = 0;
    int done = 0;
    double e_new;
    double e_mp2j = 0.0;
    double e_old = e_mp2j;
    double rms = 1.0;

    double delta = options_.get_double("RPA_DELTA"); 
 
    // => Iterations <= //
    fprintf(outfile,"  Iter       Energy (H)          dE (H)             RMS (H)     Time (s) Vectors\n");
    do {
 
        // => Validation <= //
        SharedMatrix Texact; 
        if (debug_) {
            Texact = SharedMatrix(new Matrix("-T exact", nov, nov));
            double** Texactp = Texact->pointer();
            C_DGEMM('N','T',nov,nov,naux,1.0,ZiaQp[0],naux,ZiaQp[0],naux,0.0,Texactp[0],nov);

            for (int ia = 0; ia < nov; ia++) {
                int i = ia / nvir;
                int a = ia % nvir;
                for (int jb = 0; jb < nov; jb++) {
                    int j = jb / nvir;
                    int b = jb % nvir;
                    Texactp[ia][jb] /= evals_avirp_[a] + evals_avirp_[b] - evals_aoccp_[i] - evals_aoccp_[j];
                }
            }    

            Texact->print();
        }

        // => Superdiagonal <= // 
        boost::shared_ptr<Vector> t_iaia(new Vector("-t_ia^ia",nov));
        double* t_iaiap = t_iaia->pointer();
        std::vector<std::pair<double, int> > super;
        for (int i = 0, ia = 0; i < nocc; i++) {
            for (int a = 0; a < nvir; a++, ia++) {
                t_iaiap[ia] = C_DDOT(naux,ZiaQp[ia],1,ZiaQp[ia],1) / \
                    (2.0 * (evals_avirp_[a] - evals_aoccp_[i]));
                super.push_back(make_pair(t_iaiap[ia], ia));
            }
        }
      
        if (debug_) 
            t_iaia->print();
 
        // => Sort <= //
        std::sort(super.begin(), super.end(), greater<std::pair<double, int> >() );
        boost::shared_ptr<IntVector> order(new IntVector("ia' Order", nov));
        int* orderp = order->pointer();
        for (int ia = 0; ia < nov; ia++) {
            orderp[ia] = super[ia].second;
        }

        if (debug_) {
            order->print(outfile);
            fflush(outfile);
        }   

        // => More Validation <= //
        SharedMatrix Tsort;
        if (debug_) {
            Tsort = SharedMatrix(new Matrix("-t_ia^jb Unsorted",nov,nov));
            double** Tsortp = Tsort->pointer();
            double** Texactp = Texact->pointer();

            for (int ia = 0; ia < nov; ia++) {
                for (int jb = 0; jb < nov; jb++) {
                    Tsortp[ia][jb] = Texactp[orderp[ia]][orderp[jb]];
                }
            }

            Tsort->print();

            C_DPOTRF('L',nov,Tsortp[0],nov);
            //Tsort->zero_lower();
            Tsort->set_name("-t_ia_jb Cholesky Unsorted Exact");
            Tsort->print();
        } 

        // => Eigenstructure Validation <= //
        if (debug_) {
            SharedMatrix Temp(  new Matrix("-T Temp",nov,nov));
            SharedMatrix Eigmat(new Matrix("-T Eigenvectors",nov,nov));
            boost::shared_ptr<Vector> Eigval(new Vector("-T Eigenvalues" ,nov));

            Eigmat->copy(Texact);
            Texact->diagonalize(Eigmat.get(),Eigval.get());

            Eigval->print();
        }

        // => Cholesky <= //
        std::vector<double*> tau;
        double* temp = new double[nov];
        int nP = 0;

        while (nP < nov) {
            nP++;
            int P = nP - 1;
            int ia = orderp[P];
            int i = ia / nvir;
            int a = ia % nvir;

            double* tau_ia = new double[nov];
            memset(static_cast<void*>(tau_ia), '\0', nov*sizeof(double));
            tau.push_back(tau_ia);

            // => L_ii type element <= //
            double diag = C_DDOT(naux,ZiaQp[ia],1,ZiaQp[ia],1) / \
                (2.0 * (evals_avirp_[a] - evals_aoccp_[i]));

            //fprintf(outfile," P = %d, ia = %d, t = %14.10f\n",P,ia,diag);

            for (int R = 0; R < P; R++) {
                diag -= tau[R][P] * tau[R][P];
            }

            // Throw if the user is feeling conservative and the diagonal is imaginary
            if (diag < PLUS_EPSILON && !options_.get_bool("RPA_RISKY")) 
                throw PSIEXCEPTION("dRPA T amplitudes are not numerically positive definite");

            // Else, ignore the row if imaginary
            if (diag > PLUS_EPSILON) {
                diag = sqrt(diag);

                // => L_ij type element <= //
                C_DGEMV('N',nov,naux,1.0,ZiaQp[0],naux,ZiaQp[ia],1,0.0,temp,1);
                for (int j = 0, jb = 0; j < nocc; j++) {
                    for (int b = 0; b < nvir; b++, jb++) {
                        temp[jb] /= (evals_avirp_[a] + evals_avirp_[b] - \
                            evals_aoccp_[i] - evals_aoccp_[j]);
                    }
                }         

                // Sort
                for (int jb = P+1; jb < nov; jb++) {
                    tau_ia[jb] = temp[orderp[jb]];
                }

                // TODO OpenMP that guy
                for (int R = 0; R < P; R++) {
                    C_DAXPY(nov - P - 1,-tau[R][P], &tau[R][P+1], 1, &tau_ia[P+1], 1); 
                }

                C_DSCAL(nov,1.0/diag,tau_ia,1);
                
                // Must not forget the diagonal 
                //fprintf(outfile," P = %d, diagonal = %14.10f\n", P, diag);
                tau_ia[P] = diag;
            }
 
            // => Convergence Check <= //
            
            // Assume OK
            bool OK = true;

            // Conventional wisdom
            if (diag * diag >= delta) OK = false;

            // Case where degeneracies occur
            //fprintf(outfile,"nP = %3d\n", nP);
            for (int R = P + 1; R < nov && OK; R++) {
                double errR = t_iaiap[orderp[R]];
                for (int Q = 0; Q < nP; Q++) {
                    errR -= tau[Q][R] * tau[Q][R];
                }
                //fprintf(outfile,"  R = %d, errR = %14.10f\n", R, errR);
                if (fabs(errR) > delta) OK = false; 
            }
        
            if (OK) break;
        }

        if (debug_) {
            SharedMatrix Tau_u(new Matrix("Tau Unsorted", nP, nov));
            double** Tau_up = Tau_u->pointer();
            for (int P = 0; P < nP; P++) {
                C_DCOPY(nov,tau[P],1,Tau_up[P],1);
            }
            Tau_u->print(); 
        }

        // Make a contiguous block, and backsort
        SharedMatrix Tau(new Matrix("Tau_P^ia", nP, nov));
        double** Taup = Tau->pointer();

        // Painful unsort
        #pragma omp parallel for num_threads(nthread)
        for (int P = 0; P < nP; P++) {
            for (int ia = 0; ia < nov; ia++) {
                Taup[P][orderp[ia]] = tau[P][ia];
            }
        }

        if (debug_) {
            Tau->print();
            SharedMatrix Tapp(new Matrix("-T Approximate",nov,nov));
            double** Tappp = Tapp->pointer();
            C_DGEMM('T','N',nov,nov,nP,1.0,Taup[0],nov,Taup[0],nov,0.0,Tappp[0],nov);
            Tapp->print();
            fflush(outfile);
        }
       
        if (debug_) {
            SharedMatrix Tapp(new Matrix("-T Approximate (from Exact Cholesky)",nov,nov));
            double** Tappp = Tapp->pointer();
            double** Tsortp = Tsort->pointer();
            C_DGEMM('T','N',nov,nov,nP,1.0,Tsortp[0],nov,Tsortp[0],nov,0.0,Tappp[0],nov);
            Tapp->print();
            fflush(outfile);
            
        }
 
        // Clear the memory off
        delete[] temp;
        for (int P = 0; P < nP; P++)
            delete[] tau[P]; 
   
        // => Energy Evaluation <= //
        SharedMatrix A(new Matrix("A_PQ", nP, naux));
        double** Ap = A->pointer();

        C_DGEMM('N','T',nP,naux,nov,1.0,Taup[0],nov,Qiap[0],nov,0.0,Ap[0],naux);

        double E = 0.0;
        for (ULI PR = 0L; PR < naux*(ULI)nP; PR++) {
            E += Ap[0][PR] * Ap[0][PR];    
        } 
        E *= -2.0; 

        // => Store the MP2J energy <= //
        if (iter == 0) e_mp2j = E;  
        e_new = E; 

        // => X_PQ <= //
        SharedMatrix X(new Matrix("X_PQ", nP, naux));
        double** Xp = X->pointer();
        
        C_DGEMM('N','T',nP,naux,nov,1.0,Taup[0],nov,Qiap[0],nov,0.0,Xp[0],naux);
       
        if (debug_) 
            X->print();
 
        // => Y_ia_Q <= //
        ZiaQ->set_name("Y_ia_Q");
        C_DGEMM('T','N',nov,naux,nP,-1.0,Taup[0],nov,Xp[0],naux,0.0,ZiaQp[0],naux);       
        X.reset();

        if (debug_)
            ZiaQ->print();
 
        // => DIIS Y_ia^Q <= //
        if (options_.get_bool("DIIS")) {
            if (iter == 0) {
              psio_->write_entry(DFCC_INT_FILE,"Old Y_ia_Q",(char *)
                  &(ZiaQp[0][0]),sizeof(double)*nov*naux);
            }
            else {
                diis_->store_current_vector((char *) &(ZiaQp[0][0]));
  
                double *temp = init_array(naux);
                psio_address next_YiaQ = PSIO_ZERO;
                for (int ia=0; ia < nov; ia++) {
                    psio_->read(DFCC_INT_FILE,"Old Y_ia_Q",(char *)
                        &(temp[0]),sizeof(double)*naux,next_YiaQ,&next_YiaQ);
                    C_DAXPY(naux,-1.0,temp,1,ZiaQp[ia],1);
                }
                diis_->store_error_vector((char *) &(ZiaQp[0][0]));
                free(temp);

                diis_->increment_vectors();

                rms = C_DDOT(nov*naux,ZiaQp[0],1,ZiaQp[0],1);
                rms /= (double) nov*naux;
                rms = sqrt(rms);

                char *veclabel = diis_->get_last_vec_label();
                psio_->read_entry(DFCC_DIIS_FILE,veclabel,(char *)
                    &(ZiaQp[0][0]),sizeof(double)*nov*naux);
                free(veclabel);

                psio_->write_entry(DFCC_INT_FILE,"Old Y_ia_Q",(char *)
                    &(ZiaQp[0][0]),sizeof(double)*nov*naux);
            }
            if (options_.get_int("MIN_DIIS_VECS") <= iter) {
              diis_->get_new_vector(ZiaQp,naux);
            }
        }
        else {
            if (iter == 0) {
              psio_->write_entry(DFCC_INT_FILE,"Old Y_ia_Q",(char *)
                  &(ZiaQp[0][0]),sizeof(double)*nov*naux);
            }
            else {
                rms = 0.0;
                double *temp = init_array(naux);
                psio_address next_YiaQ = PSIO_ZERO;
                for (int ia=0; ia < nov; ia++) {
                    psio_->read(DFCC_INT_FILE,"Old Y_ia_Q",(char *)
                        &(temp[0]),sizeof(double)*naux,next_YiaQ,&next_YiaQ);
                    C_DAXPY(naux,-1.0,ZiaQp[ia],1,temp,1);
                    rms += C_DDOT(naux,temp,1,temp,1);
                }
                free(temp);

                rms /= (double) nov*naux;
                rms = sqrt(rms);

                psio_->write_entry(DFCC_INT_FILE,"Old Y_ia_Q",(char *)
                    &(ZiaQp[0][0]),sizeof(double)*nov*naux);
            }
        }
 
        // => Z_ia_Q <= //
        ZiaQ->set_name("Z_ia_Q");
        for (int Q = 0; Q < naux; Q++) {
            C_DAXPY(nov,1.0,Qiap[Q],1,&ZiaQp[0][Q],naux);
        }
        
        if (debug_) 
            ZiaQ->print();
        
        // => Convergence Check <= // 
        stop = time(NULL);
        if (iter == 0)
            fprintf(outfile,
                "  %4d %16.8lf %17.9lf                   %12ld %6d",
                iter,e_new,e_old-e_new,stop-start, nP);
        else
            fprintf(outfile,"  %4d %16.8lf %17.9lf %17.9lf %12ld %6d",
                iter,e_new,e_old-e_new,rms,stop-start, nP);

        if (options_.get_int("MIN_DIIS_VECS") <= iter && 
            options_.get_bool("DIIS"))
            fprintf(outfile," DIIS\n");
        else 
            fprintf(outfile,"\n");
 
        fflush(outfile);
        iter++;
    
        if (iter > options_.get_int("MAXITER")) done = 1;
        if (fabs(e_old-e_new) < options_.get_double("E_CONVERGE"))
          done = 1;
        if (rms < options_.get_double("T_CONVERGE")) done = 1;
    
        e_old = e_new;
    }
    while(!done);
  
    double scale = options_.get_double("RPA_ALPHA");

    fprintf(outfile,"\n");
    fprintf(outfile,"  Reference Energy                    %18.10lf [H]\n",Eref_);
    fprintf(outfile,"  Total CD-dRPA Energy                %18.10lf [H]\n",Eref_ + e_old);
    fprintf(outfile,"  CD-MP2J Energy                      %18.10lf [H]\n",e_mp2j);
    fprintf(outfile,"  CD-dRPA Energy                      %18.10lf [H]\n",e_old);
    fprintf(outfile,"  Delta CD-dRPA/MP2J Energy           %18.10lf [H]\n",e_old - e_mp2j);
    fprintf(outfile,"  Delta CD-dRPA/MP2J Scale            %18.10lf [-]\n",scale); 
    fprintf(outfile,"  Scaled Delta CD-dRPA/MP2J Energy    %18.10lf [H]\n",scale*(e_old - e_mp2j));

    Process::environment.globals["RPA TOTAL ENERGY"]        = Eref_ + e_old;
    Process::environment.globals["RPA MP2J ENERGY"]         = e_mp2j;
    Process::environment.globals["RPA DRPA ENERGY"]         = e_old;
    Process::environment.globals["RPA DELTA ENERGY"]        = e_old - e_mp2j;
    Process::environment.globals["RPA SCALED DELTA ENERGY"] = scale*(e_old - e_mp2j);
    Process::environment.globals["CURRENT ENERGY"]          = Eref_ + e_old;
 
    return Eref_ + e_old;
}

void dRPA::print_header()
{
    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*                        dRPA                          *\n");
    fprintf(outfile, "\t*           Direct Random Phase Approximation          *\n");
    fprintf(outfile, "\t*                with all sorts of shit                *\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*            Rob Parrish and Ed Hohenstein             *\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\n");
    CC::print_header();

}

void dRPA::apply_denom()
{
  for (int i=0,ia=0; i<naocc_; i++) {
  for (int a=0; a<navir_; a++,ia++) {
    for (int j=0,jb=0; j<naocc_; j++) {
    for (int b=0; b<navir_; b++,jb++) {
      long int iajb = (long int) ia*naocc_*navir_ + jb;
      tIAJB_[iajb] /= evals_aoccp_[i] + evals_aoccp_[j] - evals_avirp_[a]
        - evals_avirp_[b];
    }}
  }}
}

double dRPA::df_energy()
{
  // Because Ed's a pussy
  return(2.0*C_DDOT(naocc_*navir_*ndf_,Th_p_IA_[0],1,B_p_IA_[0],1));
}

double dRPA::df_store_error_vecs()
{
  double rms;

  C_DAXPY((long int) naocc_*navir_*naocc_*navir_,-1.0,tIAJB_,1,t2IAJB_,1);

  if (options_.get_bool("DIIS"))
    diis_->store_vectors(tIAJB_,t2IAJB_);

  rms = C_DDOT((long int) naocc_*navir_*naocc_*navir_,t2IAJB_,1,t2IAJB_,1);
  rms /= (double) naocc_*navir_*naocc_*navir_;

  return(sqrt(rms));
}

}}

