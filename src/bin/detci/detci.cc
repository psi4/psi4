/*! \defgroup DETCI detci: The Determinant CI code */

/*! \file
    \ingroup DETCI
    \brief Determinant-based CI program

   DETCI

   DETERMINANT CI Program, incorporating Abelian point-group symmetry

   C. David Sherrill
   Center for Computational Quantum Chemistry
   University of Georgia
   August 1994

   Updated 3/95 to do frozen core and virtuals correctly
   Updated 5/95 to do RAS CI's again
   Updated 2/96 to clean up code and rename DETCI

*/

//#include <sstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sys/types.h>
#include <unistd.h>
#include <cstring>
#include <psifiles.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libmints/wavefunction.h>
#include <libmints/molecule.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <libqt/slaterdset.h>
#include <masses.h>
#include "structs.h"
#include "globals.h"
#include "ci_tol.h"
#include <pthread.h>
#include "tpool.h"
#include <iostream>
#include "odometer.h"
#include "slaterd.h"
#include "civect.h"
#include "ciwave.h"

namespace psi { namespace detci {

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

unsigned char ***Occs;
struct stringwr **alplist;
struct stringwr **betlist;

extern void eivout_t(double **evecs, double *evals, int rows, int cols,
   FILE *outfile);
extern void read_integrals(void);
extern void tf_onel_ints(int printflg, FILE *outfile);
extern void zapt_shift(double *TEI, int nirreps, int nmo, int *doccpi,
   int *soccpi, int *orbspi, int *frzdoccpi, int *reorder);
extern void form_gmat(int printflg, FILE *outfile);
extern void get_mo_info(Options &);
extern void print_vec(unsigned int nprint, int *Iacode, int *Ibcode,
   int *Iaidx, int *Ibidx, double *coeff,
   struct olsen_graph *AlphaG, struct olsen_graph *BetaG,
   struct stringwr **alplist, struct stringwr **betlist,
   FILE *outfile);
extern void print_config(int nbf, int num_alp_el, int num_bet_el,
   struct stringwr *stralp, struct stringwr *strbet,
   int num_fzc_orbs, char *outstring);
extern void init_stringwr_temps(int nel, int num_ci_orbs, int nsym);
extern void free_stringwr_temps(int nsym);
extern void str_abs2rel(int absidx, int *relidx, int *listnum,
   struct olsen_graph *Graph);
extern int str_rel2abs(int relidx, int listnum, struct olsen_graph *Graph);
extern void H0block_init(unsigned int size);
extern void H0block_fill(struct stringwr **alplist,
   struct stringwr **betlist);
extern void H0block_free(void);
extern void H0block_print(void);
extern void H0block_setup(int num_blocks, int *Ia_code, int *Ib_code);
extern void H0block_pairup(int guess);
extern void H0block_spin_cpl_chk(void);
extern void H0block_filter_setup(void);
extern void sem_test(double **A, int N, int M, int L, double **evecs,
   double *evals, double **b, double conv_e, double conv_rms,
   int maxiter, double offst, int *vu, int maxnvect, FILE *outfile);
extern void form_ov(struct stringwr **alplist);
extern void write_energy(int nroots, double *evals, double offset);

void init_io(void);
void close_io(void);
void cleanup(void);
void title(void);
void quote(void);
void init_ioff(void);
void diag_h(struct stringwr **strlista, struct stringwr **strlistb);
void mpn(struct stringwr **strlista, struct stringwr **strlistb);
void form_opdm(void);
void form_tpdm(void);
extern void get_parameters(Options &);
extern void print_parameters(void);
extern void set_ras_parms(void);
extern void print_ras_parms(void);
extern void form_strings(void);
extern void mitrush_iter(CIvect &Hd,
   struct stringwr **alplist, struct stringwr **betlist,
   int nroots, double *evals, double conv_rms, double conv_e, double enuc,
   double efzc, int maxiter, int maxnvect, FILE *outfile,
   int print_lvl);
extern void sem_iter(CIvect &Hd, struct stringwr **alplist, struct stringwr
   **betlist, double *evals, double conv_e,
   double conv_rms, double enuc, double efzc,
   int nroots, int maxiter, int maxnvect, FILE *outfile, int print_lvl);
extern void mpn_generator(CIvect &Hd, struct stringwr **alplist,
   struct stringwr **betlist);
extern void opdm(struct stringwr **alplist, struct stringwr **betlist,
   int transdens, int dipmom,
   int Inroots, int Iroot, int Inunits, int Ifirstunit,
   int Jnroots, int Jroot, int Jnunits, int Jfirstunit,
   int targetfile, int writeflag, int printflag);
extern void tpdm(struct stringwr **alplist, struct stringwr **betlist,
   int Inroots, int Inunits, int Ifirstunit,
   int Jnroots, int Jnunits, int Jfirstunit,
   int targetfile, int writeflag, int printflag);
extern void compute_cc(void);
extern void calc_mrpt(void);

PsiReturnType detci(Options &options);

}} // namespace psi::detci


namespace psi { namespace detci {

CIWavefunction::CIWavefunction(boost::shared_ptr<Wavefunction> reference_wavefunction, Options &options)
    : Wavefunction(options, _default_psio_lib_)
{
    set_reference_wavefunction(reference_wavefunction);
    init();
}

CIWavefunction::~CIWavefunction()
{

}

void CIWavefunction::init()
{
    // TODO-CDS:
    // The CC codes destroy the checkpoint object created by Wavefunction.
    // We'd like to be able to do the same here.  Need to convert everything
    // such that we don't explicitly need checkpoint
    // Destroy it. Otherwise we will see a "file already open" error.
    chkpt_.reset();

    // We're copying this stuff over like it's done in ccenergy, but
    // these Wavefunction member data are not actually used by the code
    // yet (Sept 2011 CDS).  Copying these only in case they're needed
    // by someone else who uses the CIWavefunction and expects them to
    // be available.

    // CDS-TODO: Note, some of these data should be updated to reflect what
    // the CI wavefunction is doing, not what the reference wavefunction
    // is doing, in case these are actually used elsewhere.
    nso_        = reference_wavefunction_->nso();
    nirrep_     = reference_wavefunction_->nirrep();
    nmo_        = reference_wavefunction_->nmo();
    for(int h = 0; h < nirrep_; ++h){
        soccpi_[h] = reference_wavefunction_->soccpi()[h];
        doccpi_[h] = reference_wavefunction_->doccpi()[h];
        frzcpi_[h] = reference_wavefunction_->frzcpi()[h];
        frzvpi_[h] = reference_wavefunction_->frzvpi()[h];
        nmopi_[h]  = reference_wavefunction_->nmopi()[h];
        nsopi_[h]  = reference_wavefunction_->nsopi()[h];
    }
}

double CIWavefunction::compute_energy()
{
    energy_ = 0.0;
    PsiReturnType ci_return;
    if ((ci_return = psi::detci::detci(options_)) == Success) {
        // Get the total energy
        energy_ = Process::environment.globals["CURRENT ENERGY"];
    }

    return energy_;
}

PsiReturnType detci(Options &options)
{
   Parameters.print_lvl = 1;
   Parameters.have_special_conv = 0;
   int fci_norb_check = 0;
   int i = 0;

   boost::shared_ptr<PSIO> psio = PSIO::shared_object();

   init_io();                   /* parse cmd line and open input and output */
   get_parameters(options);     /* get running params (convergence, etc)    */
   init_ioff();                 /* set up the ioff array                    */
   title();                     /* print program identification             */
   get_mo_info(options);        /* read DOCC, SOCC, frozen, nmo, etc        */
   set_ras_parms();             /* set fermi levels and the like            */

   if (Parameters.print_lvl) {
     print_parameters();       /* print running parameters                 */
     print_ras_parms();
   }
   fflush(outfile);

   form_strings();              /* form the alpha/beta strings              */
   if (Parameters.nthreads > 1)
     tpool_init(&thread_pool, Parameters.nthreads, CalcInfo.num_alp_str, 0);
                                /* initialize thread pool */
   init_time_new(detci_time);             /* initialize timing routines */
   fflush(outfile);

   if (Parameters.istop) {      /* Print size of space, other stuff, only   */
     close_io();
     cleanup();
     Process::environment.globals["CURRENT ENERGY"] = 0.0;
     Process::environment.globals["CURRENT CORRELATION ENERGY"] = 0.0;
     Process::environment.globals["CI TOTAL ENERGY"] = 0.0;
     Process::environment.globals["CI CORRELATION ENERGY"] = 0.0;

     return Success;
   }


   read_integrals();            /* get the 1 and 2 elec MO integrals        */

   if(Parameters.zaptn)         /* Shift SCF eigenvalues for ZAPTn          */
      zapt_shift(CalcInfo.twoel_ints, CalcInfo.nirreps, CalcInfo.nmo,
         CalcInfo.docc, CalcInfo.socc, CalcInfo.orbs_per_irr,
         CalcInfo.frozen_docc, CalcInfo.reorder);

   if (Parameters.bendazzoli) /* form the Bendazzoli OV arrays            */
      form_ov(alplist);
                                /* lump together one-electron contributions */
   tf_onel_ints((Parameters.print_lvl>3), outfile);

                                /* form the RAS g matrix (eq 28-29)         */
   form_gmat((Parameters.print_lvl>3), outfile);
   fflush(outfile);

   if (Parameters.mpn)
     mpn(alplist, betlist);
   else if (Parameters.cc)
     compute_cc();
   else
     diag_h(alplist, betlist);

   if (Parameters.opdm || Parameters.transdens) form_opdm();
   if (Parameters.tpdm) form_tpdm();
   if (Parameters.print_lvl) print_time_new(detci_time);
   if (Parameters.nthreads > 1) tpool_destroy(thread_pool, 1);
   if (Parameters.print_lvl > 0) quote();
   close_io();
   cleanup();
   return Success;
}


/*
** init_io(): Figure out command-line arguments and start up I/O for
** input and output files.
**
*/
//void init_io(int argc, char *argv[])
void init_io(void)
{
   /*
   int i, num_unparsed;
   char *argv_unparsed[100];
   */

   /*
   Parameters.write_energy = 0;

   for (i=1,num_unparsed=0; i<argc; i++) {
      if (strcmp(argv[i], "--quiet") == 0) {
         Parameters.print_lvl = 0;
      }
      else if (strcmp(argv[i], "-e") == 0) {
         Parameters.write_energy = 1;
      }
      else if (strcmp(argv[i], "-c") == 0) {
         Parameters.have_special_conv = 1;
         if (i+1 >= argc) {
            fprintf(stderr, "detci: -c flag requires an argument\n");
            exit(1);
         }
         if (sscanf(argv[i+1], "%lf", &(Parameters.special_conv)) != 1) {
            fprintf(stderr, "detci: trouble reading argument to -c flag\n");
            exit(1);
         }
      i++;
      }
      else {
        argv_unparsed[num_unparsed++] = argv[i];
      }
   }
   */

   /* initialize input and output files.  We want to pass the list
    * of arguments starting after the last one we have parsed, so
    * we do pointer arithmetic on argv
    */
   /* init_in_out(argc-parsed,argv+parsed); */

    // errcod = psi_start(&infile,&outfile,&psi_file_prefix,num_unparsed,argv_unparsed,0);

   if (Parameters.print_lvl) tstart();

   /*
    * this stuff is now inside psi_start
   ip_set_uppercase(1);
   ip_initialize(infile, outfile);
   ip_cwk_clear();
   ip_cwk_add(":DEFAULT");
   */

   /*
   ip_cwk_add(":DETCI");
   psio_init(); psio_ipv1_config();
   */

}



/*
** init_ioff(): Set up the ioff array for quick indexing
**
*/
void init_ioff(void)
{
   int i;

   /* set offsets for ij-type canonical ordering */
   ioff = (int *) malloc (IOFF_MAX * sizeof(int)) ;
   ioff[0] = 0;
   for (i = 1; i < IOFF_MAX ; i++) {
      ioff[i] = ioff[i-1] + i;
      }
}



/*
** close_io(): Function closes down I/O and exits
*/
void close_io(void)
{
   // psio_done();
   if (Parameters.print_lvl) tstop();
   //errcod = psi_stop(infile,outfile,psi_file_prefix);
}


/*
** cleanup(): Free any allocated memory that wasn't already freed elsewhere
*/
void cleanup(void)
{
 
}

/*
** title(): Function prints a program identification
*/
void title(void)
{
  if (Parameters.print_lvl) {
   fprintf(outfile,"\n");
   fprintf(outfile,"*******************************************************\n");
   fprintf(outfile,"                       D E T C I  \n");
   fprintf(outfile,"\n");
   fprintf(outfile,"                   C. David Sherrill\n") ;
   fprintf(outfile,"                   Matt L. Leininger\n") ;
   fprintf(outfile,"                     18 June 1999\n") ;
   fprintf(outfile,"*******************************************************\n");
   fprintf(outfile,"\n\n\n");
   }
  else {
   fprintf(outfile,
   "\nD E T C I : C. David Sherrill and Matt L. Leininger, 18 June 1999\n");
   }
  fflush(outfile);
}




/*
** diag_h(): Function diagonalizes the hamiltonian
**
** Parameters:
**    alplist = list of alpha strings
**    betlist = list of beta strings
**
** Returns: none
*/
void diag_h(struct stringwr **alplist, struct stringwr **betlist)
{
   int nroots, i, j, size;
   double conv_rms, conv_e, *evals, **evecs, nucrep, efzc, tval;
   int *tptr;
   double *cbuf;
   char e_label[PSIO_KEYLEN]; /* 80... */

   nroots = Parameters.num_roots;

   conv_rms = Parameters.   convergence;
   conv_e = Parameters.energy_convergence;

   if (Parameters.have_special_conv) {
     tval = sqrt(conv_rms) * 10.0;
     if (tval > conv_e) conv_e = tval;
     if (Parameters.special_conv > conv_rms)
       conv_rms = Parameters.special_conv;
   }
   size = CIblks.vectlen;
   if (Parameters.nprint > size) Parameters.nprint = size;
   nucrep = CalcInfo.enuc;
   efzc = CalcInfo.efzc;
   tmp_ras_array = init_array(1024);

   H0block.size = 0;
   H0block.osize = 0;

   if (Parameters.bendazzoli)
      fprintf(outfile, "\nBendazzoli algorithm selected for sigma3\n");

   /* Direct Method --- use RSP diagonalization routine */
   if (Parameters.diag_method == METHOD_RSP) {

      CIvect Cvec(CIblks.vectlen, CIblks.num_blocks, 1, Parameters.Ms0,
         CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
         CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
         CalcInfo.nirreps, AlphaG->subgr_per_irrep, 1, 0, 0,
         CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);
      // shouldn't need to open I/O files for this fake CIvec, unit=0 

      double **H, **rsp_evecs;
      int Iarel, Ialist, Ibrel, Iblist;
      unsigned long int ii, jj;
      SlaterDeterminant I, J;
      int *mi_iac, *mi_ibc, *mi_iaidx, *mi_ibidx;
      double *mi_coeff;

      if (Parameters.print_lvl) {
         fprintf(outfile, "\nFind all roots with RSP\n") ;
         fprintf(outfile, "\n") ;
         fflush(outfile) ;
         }

      /* construct and print one block at a time for debugging */
      /*
      int ii2, jj2, blk, blk2, det1, det2;
      double **Hpart;

      for (blk = 0; blk < CIblks.num_blocks; blk++) {
        for (blk2 = 0; blk2 < CIblks.num_blocks; blk2++) {
          Hpart = init_matrix(CIblks.Ia_size[blk]*CIblks.Ib_size[blk],
                              CIblks.Ia_size[blk2]*CIblks.Ib_size[blk2]);
          for (ii=0,det1=0; ii<CIblks.Ia_size[blk]; ii++) {
            for (jj=0; jj<CIblks.Ib_size[blk]; jj++, det1++) {
              I.set(CalcInfo.num_alp_expl,alplist[CIblks.Ia_code[blk]][ii].occs,
                   CalcInfo.num_bet_expl,betlist[CIblks.Ib_code[blk]][jj].occs);
              for (ii2=0,det2=0; ii2<CIblks.Ia_size[blk2]; ii2++) {
                for (jj2=0; jj2<CIblks.Ib_size[blk2]; jj2++,det2++) {
                  J.set(CalcInfo.num_alp_expl,
                        alplist[CIblks.Ia_code[blk2]][ii2].occs,
                        CalcInfo.num_bet_expl,
                        betlist[CIblks.Ib_code[blk2]][jj2].occs);
                  Hpart[det1][det2] = matrix_element(&I,&J);
                }
              }
            }
          }
          if (Parameters.print_lvl > 4 && size < 200) {
            fprintf(outfile, "\nBlock %d %d of ", blk, blk2);
            fprintf(outfile, "Hamiltonian matrix:\n");
            print_mat(Hpart, CIblks.Ia_size[blk]*CIblks.Ib_size[blk],
                             CIblks.Ia_size[blk2]*CIblks.Ib_size[blk2],
                      outfile);
          }
          free_matrix(Hpart, CIblks.Ia_size[blk]*CIblks.Ib_size[blk]);
        }
      }
      */
      /* end block-at-a-time stuff */

      H = init_matrix(size, size);
      rsp_evecs = init_matrix(size, size) ;
      evals = init_array(size) ;
      for (ii=0; ii<size; ii++) {
         Cvec.det2strings(ii, &Ialist, &Iarel, &Iblist, &Ibrel);
         I.set(CalcInfo.num_alp_expl,
               alplist[Ialist][Iarel].occs, CalcInfo.num_bet_expl,
               betlist[Iblist][Ibrel].occs);
         /* introduce symmetry or other restrictions here */
         for (jj=0; jj<ii; jj++) {
            Cvec.det2strings(jj, &Ialist, &Iarel, &Iblist, &Ibrel);
            J.set(CalcInfo.num_alp_expl,
               alplist[Ialist][Iarel].occs, CalcInfo.num_bet_expl,
               betlist[Iblist][Ibrel].occs);
            H[ii][jj] = H[jj][ii] = matrix_element(&I, &J);
            }
         Cvec.det2strings(jj, &Ialist, &Iarel, &Iblist, &Ibrel);
         J.set(CalcInfo.num_alp_expl,
            alplist[Ialist][Iarel].occs, CalcInfo.num_bet_expl,
            betlist[Iblist][Ibrel].occs);
         H[jj][jj] = matrix_element(&J, &J) + CalcInfo.efzc;
         }

      if (Parameters.print_lvl > 4 && size < 200) {
         fprintf(outfile, "\nHamiltonian matrix:\n");
         print_mat(H, size, size, outfile);
         }

      sq_rsp(size, size, H, evals, 1, rsp_evecs, 1.0E-14);
      if (Parameters.print_lvl > 4) {
         eivout(rsp_evecs, evals, size, nroots, outfile) ;
         }
      evecs = init_matrix(nroots, size);
      for (ii=0; ii<nroots; ii++) {
         for (jj=0; jj<size; jj++) {
            evecs[ii][jj] = rsp_evecs[jj][ii];
            }
         }
      free_matrix(rsp_evecs, size);
      free_matrix(H, size);

      if (Parameters.print_lvl) {
         mi_iac = init_int_array(Parameters.nprint);
         mi_ibc = init_int_array(Parameters.nprint);
         mi_iaidx = init_int_array(Parameters.nprint);
         mi_ibidx = init_int_array(Parameters.nprint);
         mi_coeff = init_array(Parameters.nprint);

         for (i=0; i<nroots; i++) {
            fprintf(outfile,
               "\n\n* ROOT %2d CI total energy = %17.13lf\n",
               i+1,evals[i]+nucrep);
            Cvec.setarray(evecs[i], size);
            zero_arr(mi_coeff, Parameters.nprint);
            Cvec.max_abs_vals(Parameters.nprint, mi_iac, mi_ibc, mi_iaidx,
               mi_ibidx, mi_coeff, Parameters.neg_only);
            print_vec(Parameters.nprint, mi_iac, mi_ibc, mi_iaidx, mi_ibidx,
               mi_coeff, AlphaG, BetaG, alplist, betlist, outfile);
            }

         free(mi_iac);  free(mi_ibc);
         free(mi_iaidx);  free(mi_ibidx);
         free(mi_coeff);
         }

      if (Parameters.write_energy) write_energy(nroots, evals, nucrep);

      /* Dump the vector to a PSIO file
         Added by Edward valeev (June 2002) */
      if (Parameters.export_ci_vector) {
        StringSet alphastrings, betastrings;
        SlaterDetSet dets;
        SlaterDetVector vec;
        short int *fzc_occ;
        unsigned char *newocc;

        if (CalcInfo.num_fzc_orbs > 0) {
          fzc_occ = (short int *)
            malloc(CalcInfo.num_fzc_orbs*sizeof(short int));
          for (int l=0; l<CalcInfo.num_fzc_orbs; l++) {
            fzc_occ[l] = CalcInfo.order[l]; /* put it in Pitzer order */
          }
        }

        newocc = (unsigned char *)
          malloc(((AlphaG->num_el > BetaG->num_el) ?
            AlphaG->num_el : BetaG->num_el)*sizeof(unsigned char));

        stringset_init(&alphastrings,AlphaG->num_str,AlphaG->num_el,
                       CalcInfo.num_fzc_orbs, fzc_occ);
        int list_gr = 0;
        int offset = 0;
        for(int irrep=0; irrep<AlphaG->nirreps; irrep++) {
          for(int gr=0; gr<AlphaG->subgr_per_irrep; gr++,list_gr++) {
            int nlists_per_gr = AlphaG->sg[irrep][gr].num_strings;
            for(int l=0; l<nlists_per_gr; l++) {
              /* convert occs to Pitzer order */
              for (int n=0; n<AlphaG->num_el; n++) {
                newocc[n] = (unsigned char)
                  CalcInfo.order[alplist[list_gr][l].occs[n] +
                                CalcInfo.num_fzc_orbs];
              }
              stringset_add(&alphastrings,l+offset,newocc);
            }
            offset += nlists_per_gr;
          }
        }

        stringset_init(&betastrings,BetaG->num_str,BetaG->num_el,
                       CalcInfo.num_fzc_orbs, fzc_occ);
        list_gr = 0;
        offset = 0;
        for(int irrep=0; irrep<BetaG->nirreps; irrep++) {
          for(int gr=0; gr<BetaG->subgr_per_irrep; gr++,list_gr++) {
            int nlists_per_gr = BetaG->sg[irrep][gr].num_strings;
            for(int l=0; l<nlists_per_gr; l++) {
              /* convert occs to Pitzer order */
              for (int n=0; n<BetaG->num_el; n++) {
                newocc[n] = (unsigned char)
                  CalcInfo.order[betlist[list_gr][l].occs[n] +
                                CalcInfo.num_fzc_orbs];
              }
              stringset_add(&betastrings,l+offset,newocc);
            }
            offset += nlists_per_gr;
          }
        }
        free(newocc);
        if (CalcInfo.num_fzc_orbs > 0)
          free(fzc_occ);

        int Iarel, Ialist, Ibrel, Iblist;
        slaterdetset_init(&dets,size,&alphastrings,&betastrings);
        for (int ii=0; ii<size; ii++) {
          Cvec.det2strings(ii, &Ialist, &Iarel, &Iblist, &Ibrel);
          int irrep = Ialist/AlphaG->subgr_per_irrep;
          int gr = Ialist%AlphaG->subgr_per_irrep;
          int Ia = Iarel + AlphaG->list_offset[Ialist];
          irrep = Iblist/BetaG->subgr_per_irrep;
          gr = Iblist%BetaG->subgr_per_irrep;
          int Ib = Ibrel + BetaG->list_offset[Iblist];
          slaterdetset_add(&dets, ii, Ia, Ib);
        }

        slaterdetvector_init(&vec, &dets);
        slaterdetvector_set(&vec, evecs[0]);
        slaterdetvector_write(PSIF_CIVECT,"CI vector",&vec);

        slaterdetvector_delete(&vec);
        slaterdetset_delete(&dets);
        stringset_delete(&alphastrings);
        stringset_delete(&betastrings);
      }
    } /* end RSP section */

   /* RSP test of Davidson/Liu (SEM) diagonalization routine */
   else if (Parameters.diag_method == METHOD_RSPTEST_OF_SEM) {

      // in-core CIvectors, shouldn't need to open files
      CIvect Cvec(CIblks.vectlen, CIblks.num_blocks, 1, Parameters.Ms0,
         CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
         CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
         CalcInfo.nirreps, AlphaG->subgr_per_irrep, 1, 0, 0,
         CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);
      CIvect Hd(CIblks.vectlen, CIblks.num_blocks, 1, Parameters.Ms0,
         CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size, CIblks.Ib_size,
         CIblks.offset, CIblks.num_alp_codes, CIblks.num_bet_codes,
         CalcInfo.nirreps, AlphaG->subgr_per_irrep, 1, 0, 0,
         CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);

      double **H, **b;
      int Ia, Ib, Iarel, Ialist, Ibrel, Iblist, ij, k, l, tmpi, L;
      unsigned long int ii, jj;
      SlaterDeterminant I, J;
      int *mi_iac, *mi_ibc, *mi_iaidx, *mi_ibidx;
      double *mi_coeff;
      int sm_tridim;
      double *sm_evals, *sm_mat, **sm_evecs, tval;

      if (Parameters.print_lvl) {
         fprintf(outfile, "\nFind the roots by the SEM Test Method\n");
         fprintf(outfile, "(n.b. this is for debugging purposes only!)\n");
         fprintf(outfile, "Energy convergence = %3g\n", conv_e);
         fprintf(outfile, "RMS CI vector convergence = %3g\n\n", conv_rms);
         fflush(outfile) ;
         }

      H0block_init(size);

      /* get the diagonal elements of H into an array Hd */

      Hd.diag_mat_els(alplist, betlist, CalcInfo.onel_ints,
         CalcInfo.twoel_ints, efzc, CalcInfo.num_alp_expl,
         CalcInfo.num_bet_expl, CalcInfo.num_ci_orbs, Parameters.hd_ave);

      /* get the biggest elements and put in H0block */
      if (H0block.size) {
         Hd.max_abs_vals(H0block.size, H0block.alplist, H0block.betlist,
            H0block.alpidx, H0block.betidx, H0block.H00, Parameters.neg_only);
         }

    /* MLL added this line 5-21-98 */
    /*
      if (Parameters.hd_otf) rclose(Parameters.first_hd_tmp_unit,4);
    */

      H0block_setup(CIblks.num_blocks, CIblks.Ia_code, CIblks.Ib_code);
      if (Parameters.hd_ave) {
        H0block_spin_cpl_chk();
         if (H0block.osize - H0block.size) {
            fprintf(outfile,"H0block size reduced by %d to %d to ensure"
             "completion of spin-coupling sets\n",
             (H0block.osize - H0block.size), H0block.size);
            fflush(outfile);
            }
        }
      if (Parameters.Ms0) {
         H0block_pairup(0);
         if (H0block.osize - H0block.size) {
            fprintf(outfile,"H0block size reduced by %d to ensure pairing.\n",
               (H0block.osize - H0block.size));
            fflush(outfile);
            }
         }

      if (Parameters.print_lvl > 4 && Parameters.hd_otf == FALSE) {
         fprintf(outfile, "\nDiagonal elements of the Hamiltonian\n");
         Hd.print(outfile);
         }


      if (H0block.size) {
         H0block_fill(alplist, betlist);
         }

      if (Parameters.print_lvl > 2 && H0block.size) {
         H0block_print();
         }

      if (Parameters.print_lvl > 3 && H0block.size) {
         fprintf(outfile, "\n\nH0 Block:\n");
         print_mat(H0block.H0b, H0block.size, H0block.size, outfile);
         }

      H = init_matrix(size, size);
      evals = init_array(size) ;
      for (ii=0; ii<size; ii++) {
         Cvec.det2strings(ii, &Ialist, &Iarel, &Iblist, &Ibrel);
         I.set(CalcInfo.num_alp_expl,
               alplist[Ialist][Iarel].occs, CalcInfo.num_bet_expl,
               betlist[Iblist][Ibrel].occs);
         /* introduce symmetry or other restrictions here */
         for (jj=0; jj<ii; jj++) {
            Cvec.det2strings(jj, &Ialist, &Iarel, &Iblist, &Ibrel);
            J.set(CalcInfo.num_alp_expl,
               alplist[Ialist][Iarel].occs, CalcInfo.num_bet_expl,
               betlist[Iblist][Ibrel].occs);
            H[ii][jj] = H[jj][ii] = matrix_element(&I, &J);
            }
         Cvec.det2strings(jj, &Ialist, &Iarel, &Iblist, &Ibrel);
         J.set(CalcInfo.num_alp_expl,
            alplist[Ialist][Iarel].occs, CalcInfo.num_bet_expl,
            betlist[Iblist][Ibrel].occs);
         H[jj][jj] = matrix_element(&J, &J) + CalcInfo.efzc;
         }

      /* obtain a set of L orthonormal trial vectors, L > nroots */
      b = (double **) malloc (Parameters.maxnvect * sizeof(double *)) ;
      for (i=0; i<Parameters.maxnvect; i++) {
         if (i<Parameters.num_init_vecs) b[i] = init_array(size) ;
         else b[i] = NULL ;
         }

      evecs = init_matrix(Parameters.num_roots, size);

      L = H0block.size;
      sm_tridim = L * (L + 1) / 2 ;
      sm_mat = init_array(sm_tridim) ;
      sm_evals = init_array(L) ;
      sm_evecs = init_matrix(L, L) ;
      for (i=0, ij=0; i<L; i++)
         for (j=0; j<=i; j++, ij++)
            sm_mat[ij] = H0block.H0b[i][j] ;
      rsp(L, L, sm_tridim, sm_mat, sm_evals, 1, sm_evecs, 1E-14) ;

      /*
      if (Parameters.precon == PRECON_GEN_DAVIDSON) {
        for (i=0; i<H0block.size; i++) {
           H0block.H0b_eigvals[i] = sm_evals[i];
           for (j=0; j<H0block.size; i++)
              H0block.H0b_diag[i][j] = sm_evecs[i][j];
           }
        }
      */

      /* need to fill out sm_evecs into b (pad w/ 0's) */
      cbuf = *(Cvec.blockptr(0));
      Cvec.buf_unlock();
      for (i=0,k=0; i<L && k < Parameters.num_init_vecs; i++) {

         /* check sm_evecs[i] to see if it has the correct spin symmetry */
         for (j=0,tmpi=0; Parameters.Ms0 && j<L && !tmpi; j++) {
            l = H0block.pair[j];
            if (l == -1) {
               printf("(diag_h sem_test): unpaired H0block member!\n");
               exit(1);
               }
            tval = sm_evecs[l][i];
            if ((int) Parameters.S%2) tval = -tval;
            if (sm_evecs[j][i] - tval > 1.0E-12) tmpi=1;
            }
         if (tmpi) continue;

         for (j=0; j<L; j++) sm_evals[j] = sm_evecs[j][i];

         Cvec.buf_lock(b[k]);
         Cvec.init_vals(0, L, H0block.alplist, H0block.alpidx,
            H0block.betlist, H0block.betidx, H0block.blknum, sm_evals);
         Cvec.buf_unlock();
         k++;
         }
      Cvec.buf_lock(cbuf);
      free(sm_mat);
      free(sm_evals);
      free_matrix(sm_evecs, L);
      for (i=k; i<Parameters.num_init_vecs; i++) free(b[i]);
      L = k;
      if (L < Parameters.num_roots) {
         printf("(diag_h sem_test): Ooops! L < num_roots!\n");
         exit(1);
         }
      sem_test(H, size, Parameters.num_roots, L, evecs, evals, b, conv_e,
         conv_rms, Parameters.maxiter, (nucrep+CalcInfo.efzc), &i,
         Parameters.maxnvect, outfile);

      fprintf(outfile, "SEM used %d expansion vectors\n", i);

      if (Parameters.print_lvl > 4) {
         eivout(evecs, evals, size, nroots, outfile) ;
         }
      free_matrix(H, size);

      if (Parameters.print_lvl) {
         mi_iac = init_int_array(Parameters.nprint);
         mi_ibc = init_int_array(Parameters.nprint);
         mi_iaidx = init_int_array(Parameters.nprint);
         mi_ibidx = init_int_array(Parameters.nprint);
         mi_coeff = init_array(Parameters.nprint);

         for (i=0; i<nroots; i++) {
            fprintf(outfile,
               "\n\n* ROOT %2d CI total energy = %17.13lf\n",
               i+1,evals[i]+nucrep);
            Cvec.setarray(evecs[i], size);
            zero_arr(mi_coeff, Parameters.nprint);
            Cvec.max_abs_vals(Parameters.nprint, mi_iac, mi_ibc, mi_iaidx,
               mi_ibidx, mi_coeff, Parameters.neg_only);
            print_vec(Parameters.nprint, mi_iac, mi_ibc, mi_iaidx, mi_ibidx,
               mi_coeff, AlphaG, BetaG, alplist, betlist, outfile);
            }
         free(mi_iac);  free(mi_ibc);
         free(mi_iaidx);  free(mi_ibidx);
         free(mi_coeff);
         free_matrix(evecs, Parameters.num_roots);
         }

      if (Parameters.write_energy) write_energy(nroots, evals, nucrep);

      } /* end Test of Davidson/Liu section */


   /*
    * Davidson/Liu Simultaneous Expansion Method OR
    * Mitrushenkov's Olsen-modified Davidson Algorithm
    */

   else {

      /* prepare the H0 block */
      H0block_init(size);

      CIvect Hd(CIblks.vectlen, CIblks.num_blocks, Parameters.icore,
         Parameters.Ms0, CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size,
         CIblks.Ib_size, CIblks.offset, CIblks.num_alp_codes,
         CIblks.num_bet_codes, CalcInfo.nirreps, AlphaG->subgr_per_irrep, 1,
         Parameters.num_hd_tmp_units, Parameters.first_hd_tmp_unit,
         CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);

      bool open_old = false;
      if (Parameters.restart) open_old = true;
      Hd.init_io_files(open_old);

      /* get the diagonal elements of H into an array Hd */
      if (!Parameters.restart || (Parameters.restart && Parameters.hd_otf)) {
         if (Parameters.print_lvl > 1) {
            fprintf(outfile, "\nForming diagonal elements of H\n");
            fflush(outfile);
           }
         Hd.diag_mat_els(alplist, betlist, CalcInfo.onel_ints,
            CalcInfo.twoel_ints, efzc, CalcInfo.num_alp_expl,
            CalcInfo.num_bet_expl, CalcInfo.num_ci_orbs, Parameters.hd_ave);
         }
      else {
         Hd.read(0,0);
         }

      /* get the biggest elements and put in H0block */
      if (H0block.size) {

         if (Parameters.print_lvl > 1) {
            fprintf(outfile, "\nForming H0 block\n");
            fflush(outfile);
           }

         if (!Parameters.hd_otf)
           Hd.max_abs_vals(H0block.size+H0block.coupling_size,
              H0block.alplist, H0block.betlist, H0block.alpidx,
              H0block.betidx, H0block.H00, Parameters.neg_only);
         }

      //if (Parameters.hd_otf) rclose(Parameters.first_hd_tmp_unit,4);
      if (Parameters.hd_otf) psio_close(Parameters.first_hd_tmp_unit,1);

      H0block_setup(CIblks.num_blocks, CIblks.Ia_code, CIblks.Ib_code);
      if (Parameters.filter_guess) H0block_filter_setup();
      if (Parameters.hd_ave) {
        H0block_spin_cpl_chk();
         if ((H0block.osize - H0block.size) && Parameters.print_lvl > 1) {
            fprintf(outfile,"H0block size reduced by %d to ensure "
             "completion of spin-coupling sets\n",
             (H0block.osize - H0block.size));
            H0block.osize = H0block.size;
            }
         if ((H0block.oguess_size - H0block.guess_size) &&
             Parameters.print_lvl > 1) {
           fprintf(outfile,"H0block guess size reduced by %d to ensure "
             "completion of spin-coupling sets\n",
             (H0block.oguess_size - H0block.guess_size));
            H0block.oguess_size = H0block.guess_size;
           }
         if ((H0block.ocoupling_size - H0block.coupling_size) &&
             Parameters.print_lvl > 1) {
           fprintf(outfile,"H0block coupling size reduced by %d to ensure "
             "completion of spin-coupling sets\n",
             (H0block.ocoupling_size - H0block.coupling_size));
            H0block.ocoupling_size = H0block.coupling_size;
           }
        fflush(outfile);
        }
      if (Parameters.Ms0) {
         /* if (H0block.guess_size < H0block.size) */
         H0block_pairup(0); /* pairup h0block size */
         H0block_pairup(1); /* pairup guess_size */
         H0block_pairup(2); /* pairup coupling size */
         if ((H0block.osize - H0block.size) && Parameters.print_lvl > 1) {
           fprintf(outfile,"H0block size reduced by %d to ensure pairing"
               "and spin-coupling.\n", (H0block.osize - H0block.size));
            }
         if ((H0block.oguess_size - H0block.guess_size) &&
             Parameters.print_lvl > 1) {
           fprintf(outfile,"H0block guess size reduced by %d to "
               "ensure pairing and spin-coupling.\n",
               (H0block.oguess_size - H0block.guess_size));
            }
         if ((H0block.ocoupling_size - H0block.coupling_size) &&
             Parameters.print_lvl > 1) {
           fprintf(outfile,"H0block coupling size reduced by %d to "
               "ensure pairing and spin-coupling.\n",
               (H0block.ocoupling_size - H0block.coupling_size));
            }
         fflush(outfile);
         }

      Parameters.neg_only = 0; /* MLL 7-2-97 */
      if (Parameters.print_lvl > 4) {
         fprintf(outfile, "\nDiagonal elements of the Hamiltonian\n");
         Hd.print(outfile);
         }

      if (H0block.size) {
         H0block_fill(alplist, betlist);
         }

      if (Parameters.print_lvl > 2 && H0block.size) {
         H0block_print();
         }

      if (Parameters.print_lvl > 3 && H0block.size) {
         fprintf(outfile, "\n\nH0 Block:\n");
         print_mat(H0block.H0b, H0block.size, H0block.size, outfile);
         }

      /* Davidson/Liu Simultaneous Expansion Method */
      if (Parameters.diag_method == METHOD_DAVIDSON_LIU_SEM) {

         if (Parameters.print_lvl) {
            fprintf(outfile,
               "\nFind the roots by the Simultaneous Expansion Method ");
            fprintf(outfile, "(Block Davidson Method)\n");
            fprintf(outfile, "Energy convergence = %3g\n", conv_e);
            fprintf(outfile, "RMS CI vector convergence = %3g\n\n", conv_rms);
            fflush(outfile);
            }

         evals = init_array(nroots);

         sem_iter(Hd, alplist, betlist, evals, conv_e, conv_rms,
            nucrep, efzc, nroots, Parameters.maxiter,
            Parameters.maxnvect, outfile, Parameters.print_lvl);
         }

      /* Mitrushenkov's Olsen Method */
      else {

         if (Parameters.print_lvl) {
            if (Parameters.diag_method == METHOD_MITRUSHENKOV)
              fprintf(outfile,
                "\nFind the roots with Mitrushenkov's two vector algorithm\n");
            else if (Parameters.diag_method == METHOD_OLSEN)
              fprintf(outfile,
                "\nFind the roots with Olsen's single vector algorithm\n");
            fprintf(outfile, "Energy convergence = %3g\n", conv_e);
            fprintf(outfile, "RMS CI vector convergence = %3g\n", conv_rms);
            fflush(outfile);
            }

         evals = init_array(nroots);

         mitrush_iter(Hd, alplist, betlist, nroots, evals, conv_rms, conv_e,
            nucrep, efzc, Parameters.maxiter, Parameters.maxnvect, outfile,
            Parameters.print_lvl);

         H0block_free();
         }

      if (Parameters.write_energy) write_energy(nroots, evals, nucrep+efzc);

      } /* end the Davidson-Liu/Mitrushenkov-Olsen-Davidson section */

   /* write the CI energy to PSIF_CHKPT: later fix this to loop over roots */
   chkpt_init(PSIO_OPEN_OLD);
   tval = evals[Parameters.root]+efzc+nucrep;
   chkpt_wt_etot(tval);

   Process::environment.globals["CURRENT ENERGY"] = tval;
   Process::environment.globals["CURRENT CORRELATION ENERGY"] = tval - CalcInfo.escf;
   Process::environment.globals["CURRENT REFERENCE ENERGY"] = CalcInfo.escf;
   Process::environment.globals["CI TOTAL ENERGY"] = tval;
   // eref seems wrong for open shells so replace it with escf below
   // until I fix it ---CDS 11/5/11
   Process::environment.globals["CI CORRELATION ENERGY"] = tval - CalcInfo.escf;

   if (Parameters.fci) {
     Process::environment.globals["FCI TOTAL ENERGY"] = tval;
     Process::environment.globals["FCI CORRELATION ENERGY"] = tval - CalcInfo.escf;
   }
   else {
     if (Parameters.ex_lvl == 2) {
       Process::environment.globals["CISD TOTAL ENERGY"] = tval;
       Process::environment.globals["CISD CORRELATION ENERGY"] = tval - CalcInfo.escf;
     }
     else if (Parameters.ex_lvl == 3) {
       Process::environment.globals["CISDT TOTAL ENERGY"] = tval;
       Process::environment.globals["CISDT CORRELATION ENERGY"] = tval - CalcInfo.escf;
     }
     else if (Parameters.ex_lvl == 4) {
       Process::environment.globals["CISDTQ TOTAL ENERGY"] = tval;
       Process::environment.globals["CISDTQ CORRELATION ENERGY"] = tval - CalcInfo.escf;
     }
     else {
       /*- Process::environment.globals["CIn TOTAL ENERGY"] -*/
       /*- Process::environment.globals["CIn CORRELATION ENERGY"] -*/
       std::stringstream s;
       s << "CI" << Parameters.ex_lvl << " TOTAL ENERGY";
       Process::environment.globals[s.str()] = tval;
       s.str(std::string());
       s << "CI" << Parameters.ex_lvl << " CORRELATION ENERGY";
       Process::environment.globals[s.str()] = tval - CalcInfo.escf;
     }
   }

   for (i=0; i<nroots; i++) {
     sprintf(e_label,"Root %2d energy",i);
     tval = evals[i]+efzc+nucrep;
     chkpt_wt_e_labeled(e_label, tval);

     /*- Process::environment.globals["CI ROOT n TOTAL ENERGY"] -*/
     /*- Process::environment.globals["CI ROOT n CORRELATION ENERGY"] -*/
     std::stringstream s;
     s << "CI ROOT " << (i+1) << " TOTAL ENERGY";
     Process::environment.globals[s.str()] = tval;
     s.str(std::string());
     s << "CI ROOT " << (i+1) << " CORRELATION ENERGY";
     // eref seems wrong for open shells so replace it with escf below
     // until I fix it ---CDS 11/5/11
     Process::environment.globals[s.str()] = tval - CalcInfo.escf;
   }

   if (Parameters.average_num > 1) {
     tval = 0.0;
     for (i=0; i<Parameters.average_num; i++)
       tval += Parameters.average_weights[i] *
               (efzc+nucrep+evals[Parameters.average_states[i]]);
     chkpt_wt_e_labeled("State averaged energy",tval);
     Process::environment.globals["CI STATE-AVERAGED TOTAL ENERGY"] = tval;
     // eref seems wrong for open shells so replace it with escf below
     // until I fix it ---CDS 11/5/11
     Process::environment.globals["CI STATE-AVERAGED CORRELATION ENERGY"] =
       tval - CalcInfo.escf;
   }

   chkpt_close();

}



void H0block_fill(struct stringwr **alplist, struct stringwr **betlist)
{
   int i, j, size;
   int Ia, Ib, Ja, Jb;
   int Ialist, Iblist;
   SlaterDeterminant I, J;
   double *evals, **evecs;

   /* fill lower triangle */
   for (i=0; i<H0block.size; i++) {

      Ialist = H0block.alplist[i];
      Iblist = H0block.betlist[i];
      Ia = H0block.alpidx[i];
      Ib = H0block.betidx[i];
      I.set(CalcInfo.num_alp_expl,
          alplist[Ialist][Ia].occs, CalcInfo.num_bet_expl,
          betlist[Iblist][Ib].occs);
      for (j=0; j<=i; j++) {
         Ialist = H0block.alplist[j];
         Iblist = H0block.betlist[j];
         Ia = H0block.alpidx[j];
         Ib = H0block.betidx[j];
         J.set(CalcInfo.num_alp_expl,
            alplist[Ialist][Ia].occs, CalcInfo.num_bet_expl,
            betlist[Iblist][Ib].occs);

         /* pointers in next line avoids copying structures I and J */
         H0block.H0b[i][j] = matrix_element(&I, &J);
         if (i==j) H0block.H0b[i][i] += CalcInfo.efzc;
         /* fprintf(outfile," i = %d   j = %d\n",i,j); */
         }

      H0block.H00[i] = H0block.H0b[i][i];
      }

   /* fill upper triangle */
   fill_sym_matrix(H0block.H0b, H0block.size);

   /*
   evals = init_array(H0block.size);
   evecs = init_matrix(H0block.size, H0block.size);
   */
   evals = init_array(H0block.guess_size);
   evecs = init_matrix(H0block.guess_size, H0block.guess_size);
   fflush(outfile);
   if (Parameters.precon == PRECON_GEN_DAVIDSON)
     size = H0block.size;
   else
     size = H0block.guess_size;

   if (Parameters.print_lvl > 2) {
     fprintf(outfile,"H0block size = %d in H0block_fill\n",H0block.size);
     fprintf(outfile,
             "H0block guess size = %d in H0block_fill\n",H0block.guess_size);
     fprintf(outfile,
             "H0block coupling size = %d in H0block_fill\n",
             H0block.coupling_size);
     fprintf(outfile,"Diagonalizing H0block.H0b size %d in h0block_fill in"
                     " detci.cc ... ", size);
     fflush(outfile);
   }

   sq_rsp(size, size, H0block.H0b, H0block.H0b_eigvals, 1,
          H0block.H0b_diag, 1.0E-14);

   if (Parameters.print_lvl) {
      fprintf(outfile, "\n*** H0 Block Eigenvalue = %12.8lf\n",
             H0block.H0b_eigvals[0] + CalcInfo.enuc);
      fflush(outfile);
      }

   if (Parameters.print_lvl > 5 && size < 1000) {
      for (i=0; i<size; i++) H0block.H0b_eigvals[i] += CalcInfo.enuc;
      fprintf(outfile, "\nH0 Block Eigenvectors\n");
      eivout(H0block.H0b_diag, H0block.H0b_eigvals,
             size, size, outfile);
      fprintf(outfile, "\nH0b matrix\n");
      print_mat(H0block.H0b, size, size, outfile);
      }
}


void form_opdm(void)
{
  int i, j, natom;
  double *zvals, **geom;

  Process::environment.molecule()->print();

  //if (Parameters.dipmom) {
  //  chkpt_init(PSIO_OPEN_OLD);
  //  natom = chkpt_rd_natom();
  //  zvals = chkpt_rd_zvals();
  //  geom = chkpt_rd_geom();
  //  chkpt_close();

  //  fprintf(outfile, "   Cartesian Coordinates of Nuclear Centers (a.u.)\n\n");
  //  fprintf(outfile,
  //     "   Center           X                   Y                    Z\n");
  //  fprintf(outfile,
  //     "   ------   -----------------   -----------------   -----------------\n");

  //  for(i=0;i<natom;i++){
  //    fprintf(outfile,"   %4s ",atomic_labels[(int) zvals[i]]);
  //    for(j=0;j<3;j++)
  //      fprintf(outfile,"   %17.12lf",geom[i][j]);
  //    fprintf(outfile,"\n");
  //  }
  //  fprintf(outfile,"\n");
  //  fflush(outfile);

  //  free(zvals);
  //  free_block(geom);
  //}

  /* don't need Parameters.root since it writes all opdm's */
  if (Parameters.transdens) {
    opdm(alplist, betlist, 1, Parameters.dipmom,
      Parameters.num_roots, 0,
      Parameters.num_d_tmp_units, Parameters.first_d_tmp_unit,
      Parameters.num_roots, 0,
      Parameters.num_d_tmp_units, Parameters.first_d_tmp_unit,
      Parameters.opdm_file, Parameters.tdm_write, Parameters.tdm_print);
  }
  if (Parameters.opdm) {
    opdm(alplist, betlist, 0, Parameters.dipmom,
      Parameters.num_roots, 0,
      Parameters.num_d_tmp_units, Parameters.first_d_tmp_unit,
      Parameters.num_roots, 0,
      Parameters.num_d_tmp_units, Parameters.first_d_tmp_unit,
      Parameters.opdm_file, Parameters.opdm_write, Parameters.opdm_print);
  }

}


void form_tpdm(void)
{
  tpdm(alplist, betlist,
       Parameters.num_roots,
       Parameters.num_d_tmp_units, Parameters.first_d_tmp_unit,
       Parameters.num_roots,
       Parameters.num_d_tmp_units, Parameters.first_d_tmp_unit,
       Parameters.tpdm_file, Parameters.tpdm_write, Parameters.tpdm_print);
}

void quote(void)
{
   fprintf(outfile,"\t\t \"A good bug is a dead bug\" \n\n");
   fprintf(outfile,"\t\t\t - Starship Troopers\n\n");
   fprintf(outfile,"\t\t \"I didn't write FORTRAN.  That's the problem.\"\n\n");
   fprintf(outfile,"\t\t\t - Edward Valeev\n\n");
   fflush(outfile);
}

void H0block_coupling_calc(double E, struct stringwr **alplist, struct
                           stringwr **betlist)
{
   static int first_call = 1;
   int i, j, size, size2;
   double tval1, tval2, tval3;
   double *delta_2, *gamma_1, *gamma_2, *H_12, *delta_1;
   SlaterDeterminant I, J;
   int Ia, Ib, Ja, Jb;
   int Ialist, Iblist;
   double detH0;

   size = H0block.size;
   size2 = H0block.size + H0block.coupling_size;

   H_12 = init_array(H0block.coupling_size);
   delta_1 = init_array(H0block.size);
   delta_2 = init_array(H0block.coupling_size);
   gamma_1 = init_array(H0block.size);
   gamma_2 = init_array(H0block.coupling_size);

   if (Parameters.print_lvl > 5) {
      fprintf(outfile, "\nc0b in H0block_coupling_calc = \n");
      print_mat(&(H0block.c0b), 1, size2, outfile);
      fprintf(outfile, "\nc0bp in H0block_coupling_calc = \n");
      print_mat(&(H0block.c0bp), 1, size2, outfile);
      }

     /* copy to delta_1 */
     for (i=0; i<size; i++)
        delta_1[i] = H0block.c0bp[i];

     /* form delta_2 array  (D-E)^-1 r_2 */
     for (i=size; i<size2; i++) {
        tval1 = H0block.H00[i] - E;
        if (fabs(tval1) > HD_MIN)
          H0block.c0bp[i] = H0block.c0b[i]/tval1;
        else H0block.c0bp[i] = 0.0;
        delta_2[i-size] = H0block.c0bp[i];
        }
/*
     for (i=0; i<size2; i++)
        fprintf(outfile,"In Hcc H0block.c0bp[%d] = %lf\n", i, H0block.c0bp[i]);
*/

     zero_arr(gamma_2, size);
     /* Construct H_12 coupling block on-the-fly */
     for (i=0; i<size; i++) {
        Ialist = H0block.alplist[i];
        Iblist = H0block.betlist[i];
        Ia = H0block.alpidx[i];
        Ib = H0block.betidx[i];
        I.set(CalcInfo.num_alp_expl, alplist[Ialist][Ia].occs,
              CalcInfo.num_bet_expl, betlist[Iblist][Ib].occs);
        for (j=size; j<size2; j++) {
           Ialist = H0block.alplist[j];
           Iblist = H0block.betlist[j];
           Ia = H0block.alpidx[j];
           Ib = H0block.betidx[j];
           J.set(CalcInfo.num_alp_expl, alplist[Ialist][Ia].occs,
                 CalcInfo.num_bet_expl, betlist[Iblist][Ib].occs);
           H_12[j-size] = matrix_element(&I, &J);
           } /* end loop over j */

        dot_arr(H_12, delta_2, H0block.coupling_size, &tval2);
        gamma_1[i] = tval2;
        for (j=0; j<H0block.coupling_size; j++)
           gamma_2[j] += H_12[j] * delta_1[i];

        } /* end loop over i */


     /* Construct delta_1 = (H_11)^-1 gamma_1, delta_2 = (D_2-E)^-1 * gamma_2 */
     /* First delta_2 */
     for (i=size; i<size2; i++) {
        tval1 = H0block.H00[i] - E;
        if (fabs(tval1) > HD_MIN)
          delta_2[i-size] = gamma_2[i-size]/tval1;
        else delta_2[i-size] = 0.0;
        }

     /* Now delta_1 */

     /* form H0b-E and take its inverse */
     for (i=0; i<size; i++) {
        delta_1[i] = gamma_1[i];
        for (j=0; j<size; j++) {
           H0block.tmp1[i][j] = H0block.H0b[i][j];
           if (i==j) H0block.tmp1[i][i] -= E;
           }
        }

     if (Parameters.print_lvl > 4) {
        fprintf(outfile, "\n E = %lf\n", E);
        fprintf(outfile, " H0 - E\n");
        print_mat(H0block.tmp1, H0block.size, H0block.size, outfile);
        }

/*
       for (i=0; i<size; i++)
          fprintf(outfile,"gamma_1[%d] = %lf\n", i, gamma_1[i]);

       pople(H0block.tmp1, delta_1, size, 1, 1e-9, outfile,
             Parameters.print_lvl);
*/
       flin(H0block.tmp1, delta_1, size, 1, &tval1);

     /*
       detH0 = invert_matrix(H0block.tmp1, H0block.H0b_inv, size, outfile);
       mmult(H0block.H0b_inv,0,&(gamma_1),1,&(delta_1),1,size,size,1,0);
     */

    /*
       if (Parameters.update == UPDATE_OLSEN) {
         for (i=0; i<size; i++)
            for (j=0; j<size; j++) {
               H0block.tmp1[i][j] = H0block.H0b[i][j];
               if (i==j) H0block.tmp1[i][i] -= E;
               }
         pople(H0block.tmp1,H0block.s0bp,size,1,1e-9,outfile,
               Parameters.print_lvl);
         }
    */

      /* Construction of delta_1 and delta_2 completed */
      /* Now modify correction vectors in H0block structure */

      for (i=0; i<size; i++) H0block.c0bp[i] -= delta_1[i];
      for (i=size; i<size2; i++) H0block.c0bp[i] -= delta_2[i-size];

     /*
      for (i=0; i<size2; i++) {
         if (i>=H0block.coupling_size)
           H0block.c0bp[i] -= delta_2[i-size];
         else H0block.c0bp[i] -= delta_1[i];
         }
    */

    /*
      free(gamma_1);
      free(gamma_2);
      free(delta_2);
     */

}


/*
** mpn(): Function which sets up and generates the mpn series
**
** Parameters:
**    alplist = list of alpha strings
**    betlist = list of beta strings
**
** Returns: none
*/
void mpn(struct stringwr **alplist, struct stringwr **betlist)
{
  int i, j, irrep, cnt;
  struct stringwr *stralp, *strbet;
  int **fzc_orbs;
  double tval;

  H0block_init(CIblks.vectlen);

  CIvect Hd(CIblks.vectlen, CIblks.num_blocks, Parameters.icore,
         Parameters.Ms0, CIblks.Ia_code, CIblks.Ib_code, CIblks.Ia_size,
         CIblks.Ib_size, CIblks.offset, CIblks.num_alp_codes,
         CIblks.num_bet_codes, CalcInfo.nirreps, AlphaG->subgr_per_irrep, 1,
         Parameters.num_hd_tmp_units, Parameters.first_hd_tmp_unit,
         CIblks.first_iablk, CIblks.last_iablk, CIblks.decode);

  Hd.init_io_files(false);

  /* Compute E0 from orbital energies */
  stralp = alplist[CalcInfo.ref_alp_list] + CalcInfo.ref_alp_rel;
  strbet = betlist[CalcInfo.ref_bet_list] + CalcInfo.ref_bet_rel;

  fzc_orbs = init_int_matrix(CalcInfo.nirreps, CalcInfo.num_fzc_orbs);
  cnt = 0;
  for (irrep=0; irrep<CalcInfo.nirreps; irrep++)
     for (i=0; i<CalcInfo.frozen_docc[irrep]; i++)
        fzc_orbs[irrep][i] = cnt++;

  /* Loop over alp occs */
  //CalcInfo.e0 = CalcInfo.efzc;
  CalcInfo.e0 = 0.0;
  CalcInfo.e0_fzc = 0.0;
  for (i=0; i<CalcInfo.num_fzc_orbs; i++) {
     fprintf(outfile," orb_energy[%d] = %lf\n", i, CalcInfo.scfeigval[i]);
     tval = 2.0 * CalcInfo.scfeigval[i];
     CalcInfo.e0 += tval;
     CalcInfo.e0_fzc += tval;
     }

  if(Parameters.zaptn) {
    for (i=0; i<CalcInfo.num_alp_expl; i++) {
       j = (stralp->occs)[i] + CalcInfo.num_fzc_orbs;
       CalcInfo.e0 += CalcInfo.scfeigvala[j];
       }

    for (i=0; i<CalcInfo.num_bet_expl; i++) {
       j = (strbet->occs)[i] + CalcInfo.num_fzc_orbs;
       CalcInfo.e0 += CalcInfo.scfeigvalb[j];
       }
    } else {
    for (i=0; i<CalcInfo.num_alp_expl; i++) {
       j = (stralp->occs)[i] + CalcInfo.num_fzc_orbs;
       CalcInfo.e0 += CalcInfo.scfeigval[j];
       }

    for (i=0; i<CalcInfo.num_bet_expl; i++) {
       j = (strbet->occs)[i] + CalcInfo.num_fzc_orbs;
       CalcInfo.e0 += CalcInfo.scfeigval[j];
       }
    }

   /* prepare the H0 block */

   Hd.diag_mat_els(alplist, betlist, CalcInfo.onel_ints,
          CalcInfo.twoel_ints, CalcInfo.e0_fzc, CalcInfo.num_alp_expl,
          CalcInfo.num_bet_expl, CalcInfo.num_ci_orbs, Parameters.hd_ave);

   H0block_setup(CIblks.num_blocks, CIblks.Ia_code, CIblks.Ib_code);

   mpn_generator(Hd, alplist, betlist);
}

BIGINT strings2det(int alp_code, int alp_idx, int bet_code, int bet_idx) {

   int blknum;
   BIGINT addr;

   blknum = CIblks.decode[alp_code][bet_code];
   addr = CIblks.offset[blknum];
   addr += alp_idx * CIblks.Ib_size[blknum] + bet_idx;

   return(addr);

}


}} // namespace psi::detci


