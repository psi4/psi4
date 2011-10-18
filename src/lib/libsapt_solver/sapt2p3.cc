#include "sapt2p3.h"

namespace psi { namespace sapt {

SAPT2p3::SAPT2p3(Options& options, boost::shared_ptr<PSIO> psio, 
  boost::shared_ptr<Chkpt> chkpt) : SAPT2p(options, psio, chkpt)
{
  third_order_ = options_.get_bool("DO_THIRD_ORDER");
}

SAPT2p3::~SAPT2p3()
{
}

double SAPT2p3::compute_energy()
{
  print_header();

  timer_on("DF Integrals       ");
    df_integrals();
  timer_off("DF Integrals       ");
  timer_on("Omega Integrals    ");
    w_integrals();
  timer_off("Omega Integrals    ");
  timer_on("Amplitudes         "); 
    amplitudes();
  timer_off("Amplitudes         "); 
  timer_on("Elst10             ");
    elst10();
  timer_off("Elst10             ");
  timer_on("Exch10 S^2         ");
    exch10_s2();
  timer_off("Exch10 S^2         ");
  timer_on("Exch10             ");
    exch10();
  timer_off("Exch10             ");
  timer_on("Ind20,r            ");
    ind20r();
  timer_off("Ind20,r            ");
  timer_on("Exch-Ind20,r       ");
    exch_ind20r();
  timer_off("Exch-Ind20,r       ");
  timer_on("Disp20             ");
    disp20();
  timer_off("Disp20             ");
  timer_on("Exch-Disp20        ");
    exch_disp20();
  timer_off("Exch-Disp20        ");
  timer_on("Elst12             ");
    elst12();
  timer_off("Elst12             ");
  timer_on("Exch11             ");
    exch11();
  timer_off("Exch11             ");
  timer_on("Exch12             ");
    exch12();
  timer_off("Exch12             ");
  timer_on("Ind22              ");
    ind22();
  timer_off("Ind22              ");
  timer_on("Disp21             ");
    disp21();
  timer_off("Disp21             ");
  timer_on("Disp22 (SDQ)       ");
    disp22sdq();
  timer_off("Disp22 (SDQ)       ");
  timer_on("Disp22 (T)         ");
    disp22t();
  timer_off("Disp22 (T)         ");
  timer_on("Elst13             ");
    elst13();
  timer_off("Elst13             ");
  timer_on("Disp30             ");
    disp30();
  timer_off("Disp30             ");
  if (third_order_) {
    timer_on("ExchDisp30         ");
      exch_disp30();
    timer_off("ExchDisp30         ");
    timer_on("Ind30              ");
      ind30();
    timer_off("Ind30              ");
    timer_on("Ind30,r            ");
      ind30r();
    timer_off("Ind30,r            ");
    timer_on("Exch-Ind30         ");
      exch_ind30();
    timer_off("Exch-Ind30         ");
    timer_on("IndDisp30          ");
      ind_disp30();
    timer_off("IndDisp30          ");
    timer_on("ExchIndDisp30      ");
      exch_ind_disp30();
    timer_off("ExchIndDisp30      ");
  }

  print_results();

  return (e_sapt0_);
}

void SAPT2p3::print_header()
{
  if (third_order_)
    fprintf(outfile,"       SAPT2+3   \n");
  else
    fprintf(outfile,"      SAPT2+(3)  \n");
  fprintf(outfile,"    Ed Hohenstein\n") ;
  fprintf(outfile,"     6 June 2009\n") ;
  fprintf(outfile,"\n");
  fprintf(outfile,"      Orbital Information\n");
  fprintf(outfile,"  --------------------------\n");
  fprintf(outfile,"    NSO        = %9d\n",nso_);
  if (nsoA_ != nsoB_) {
    fprintf(outfile,"    NSO A      = %9d\n",nsoA_);
    fprintf(outfile,"    NSO B      = %9d\n",nsoB_);
  }
  fprintf(outfile,"    NMO        = %9d\n",nmo_);
  if (nmoA_ != nmoB_) {
    fprintf(outfile,"    NMO A      = %9d\n",nmoA_);
    fprintf(outfile,"    NMO B      = %9d\n",nmoB_);
  }
  fprintf(outfile,"    NRI        = %9d\n",ndf_);
  fprintf(outfile,"    NOCC A     = %9d\n",noccA_);
  fprintf(outfile,"    NOCC B     = %9d\n",noccB_);
  fprintf(outfile,"    FOCC A     = %9d\n",foccA_);
  fprintf(outfile,"    FOCC B     = %9d\n",foccB_);
  fprintf(outfile,"    NVIR A     = %9d\n",nvirA_);
  fprintf(outfile,"    NVIR B     = %9d\n",nvirB_);
  fprintf(outfile,"\n");
  fflush(outfile);
}

void SAPT2p3::print_results()
{
  e_sapt0_ = eHF_ + e_disp20_ + e_exch_disp20_;
  e_sapt2_ = e_sapt0_ + e_elst12_ + e_exch11_ + e_exch12_  + e_ind22_ 
    + e_exch_ind22_;
  e_sapt2p_ = e_sapt2_ + e_disp21_ + e_disp22sdq_ + e_disp22t_;
  e_sapt2pp3_ = e_sapt2p_ + e_elst13_ + e_disp30_;
  e_sapt2p3_ = e_sapt2pp3_ + e_exch_disp30_ + e_ind_disp30_ 
    + e_exch_ind_disp30_;

  e_exch_ind30r_ = e_ind30r_ * (e_exch_ind30_/e_ind30_);

  double dHF2 = eHF_ - (e_elst10_ + e_exch10_ + e_ind20_ + e_exch_ind20_);
  double dHF3 = eHF_ - (e_elst10_ + e_exch10_ + e_ind20_ + e_exch_ind20_
    + e_ind30r_ + e_exch_ind30r_);

  fprintf(outfile,"\n    SAPT Results  \n");
  fprintf(outfile,"  ------------------------------------------------------------------\n");
  fprintf(outfile,"    E_HF             %16.8lf mH %16.8lf kcal mol^-1\n",
    eHF_*1000.0,eHF_*627.5095);
  fprintf(outfile,"    Elst10,r         %16.8lf mH %16.8lf kcal mol^-1\n",
    e_elst10_*1000.0,e_elst10_*627.5095);
  fprintf(outfile,"    Elst12,r         %16.8lf mH %16.8lf kcal mol^-1\n",
    e_elst12_*1000.0,e_elst12_*627.5095);
  fprintf(outfile,"    Elst13,r         %16.8lf mH %16.8lf kcal mol^-1\n",
    e_elst13_*1000.0,e_elst13_*627.5095);
  fprintf(outfile,"    Exch10           %16.8lf mH %16.8lf kcal mol^-1\n",
    e_exch10_*1000.0,e_exch10_*627.5095);
  fprintf(outfile,"    Exch10(S^2)      %16.8lf mH %16.8lf kcal mol^-1\n",
    e_exch10_s2_*1000.0,e_exch10_s2_*627.5095);
  fprintf(outfile,"    Exch11(S^2)      %16.8lf mH %16.8lf kcal mol^-1\n",
    e_exch11_*1000.0,e_exch11_*627.5095);
  fprintf(outfile,"    Exch12(S^2)      %16.8lf mH %16.8lf kcal mol^-1\n",
    e_exch12_*1000.0,e_exch12_*627.5095);
  fprintf(outfile,"    Ind20,r          %16.8lf mH %16.8lf kcal mol^-1\n",
    e_ind20_*1000.0,e_ind20_*627.5095);
  if (third_order_)
    fprintf(outfile,"    Ind30,r          %16.8lf mH %16.8lf kcal mol^-1\n",
      e_ind30r_*1000.0,e_ind30r_*627.5095);
  fprintf(outfile,"    Ind22            %16.8lf mH %16.8lf kcal mol^-1\n",
    e_ind22_*1000.0,e_ind22_*627.5095);
  fprintf(outfile,"    Exch-Ind20,r     %16.8lf mH %16.8lf kcal mol^-1\n",
    e_exch_ind20_*1000.0,e_exch_ind20_*627.5095);
  if (third_order_)
    fprintf(outfile,"    Exch-Ind30,r     %16.8lf mH %16.8lf kcal mol^-1\n",
      e_exch_ind30r_*1000.0,e_exch_ind30r_*627.5095);
  fprintf(outfile,"    Exch-Ind22       %16.8lf mH %16.8lf kcal mol^-1\n",
    e_exch_ind22_*1000.0,e_exch_ind22_*627.5095);
  fprintf(outfile,"    delta HF,r (2)   %16.8lf mH %16.8lf kcal mol^-1\n",
    dHF2*1000.0,dHF2*627.5095);
  fprintf(outfile,"    delta HF,r (3)   %16.8lf mH %16.8lf kcal mol^-1\n",
    dHF3*1000.0,dHF3*627.5095);
  fprintf(outfile,"    Disp20           %16.8lf mH %16.8lf kcal mol^-1\n",
    e_disp20_*1000.0,e_disp20_*627.5095);
  fprintf(outfile,"    Disp30           %16.8lf mH %16.8lf kcal mol^-1\n",
    e_disp30_*1000.0,e_disp30_*627.5095);
  fprintf(outfile,"    Disp21           %16.8lf mH %16.8lf kcal mol^-1\n",
    e_disp21_*1000.0,e_disp21_*627.5095);
  fprintf(outfile,"    Disp22 (SDQ)     %16.8lf mH %16.8lf kcal mol^-1\n",
    e_disp22sdq_*1000.0,e_disp22sdq_*627.5095);
  fprintf(outfile,"    Disp22 (T)       %16.8lf mH %16.8lf kcal mol^-1\n",
    e_disp22t_*1000.0,e_disp22t_*627.5095);
  fprintf(outfile,"    Exch-Disp20      %16.8lf mH %16.8lf kcal mol^-1\n",
    e_exch_disp20_*1000.0,e_exch_disp20_*627.5095);
  if (third_order_) {
    fprintf(outfile,"    Exch-Disp30      %16.8lf mH %16.8lf kcal mol^-1\n",
      e_exch_disp30_*1000.0,e_exch_disp30_*627.5095);
    fprintf(outfile,"    Ind-Disp30       %16.8lf mH %16.8lf kcal mol^-1\n",
      e_ind_disp30_*1000.0,e_ind_disp30_*627.5095);
    fprintf(outfile,"    Exch-Ind-Disp30  %16.8lf mH %16.8lf kcal mol^-1\n",
      e_exch_ind_disp30_*1000.0,e_exch_ind_disp30_*627.5095);
  }
  fprintf(outfile,"\n    Total SAPT0      %16.8lf mH %16.8lf kcal mol^-1\n",
    e_sapt0_*1000.0,e_sapt0_*627.5095);
  fprintf(outfile,"    Total SAPT2      %16.8lf mH %16.8lf kcal mol^-1\n",
    e_sapt2_*1000.0,e_sapt2_*627.5095);
  fprintf(outfile,"    Total SAPT2+     %16.8lf mH %16.8lf kcal mol^-1\n",
    e_sapt2p_*1000.0,e_sapt2p_*627.5095);
  fprintf(outfile,"    Total SAPT2+(3)  %16.8lf mH %16.8lf kcal mol^-1\n",
    e_sapt2pp3_*1000.0,e_sapt2pp3_*627.5095);
  if (third_order_)
    fprintf(outfile,"    Total SAPT2+3    %16.8lf mH %16.8lf kcal mol^-1\n",
      e_sapt2p3_*1000.0,e_sapt2p3_*627.5095);

  double tot_elst = e_elst10_ + e_elst12_ + e_elst13_;
  double tot_exch = e_exch10_ + e_exch11_ + e_exch12_;
  double tot_ind = e_ind20_ + e_exch_ind20_ + dHF3 + e_ind22_ + e_exch_ind22_
    + e_ind30r_ + e_exch_ind30r_;
  double tot_ct = e_ind20_ + e_exch_ind20_ + e_ind22_ + e_exch_ind22_
    + e_ind30r_ + e_exch_ind30r_;
  double tot_disp = e_disp20_ + e_exch_disp20_ + e_disp21_ + e_disp22sdq_
    + e_disp22t_ + e_disp30_ + e_exch_disp30_;

  Process::environment.globals["SAPT ELST ENERGY"] = tot_elst;
  Process::environment.globals["SAPT EXCH ENERGY"] = tot_exch;
  Process::environment.globals["SAPT IND ENERGY"] = tot_ind;
  Process::environment.globals["SAPT CT ENERGY"] = tot_ct;
  Process::environment.globals["SAPT DISP ENERGY"] = tot_disp;
  Process::environment.globals["SAPT SAPT0 ENERGY"] = e_sapt0_;
  Process::environment.globals["SAPT SAPT2 ENERGY"] = e_sapt2_;
  Process::environment.globals["SAPT SAPT2+(3) ENERGY"] = e_sapt2pp3_;
  if (third_order_) {
    Process::environment.globals["SAPT SAPT2+(3) ENERGY"] = e_sapt2p3_;
    Process::environment.globals["SAPT ENERGY"] = e_sapt2p3_;
  }
  else
    Process::environment.globals["SAPT ENERGY"] = e_sapt2pp3_;
}

}}
