#ifndef SAPT0_H
#define SAPT0_H

#include "sapt.h"

namespace boost {

template<class T> class shared_ptr;

}

namespace psi { namespace sapt {

struct SAPTDFInts;
struct Iterator;

class SAPT0 : public SAPT {
private:
  virtual void print_header();
  virtual void print_results();

  void check_memory();

  void df_integrals();
  void df_integrals_aio();
  void w_integrals();

  void first_order_terms();
  void oo_df_integrals();

  SAPTDFInts set_A_AA();
  SAPTDFInts set_B_BB();
  SAPTDFInts set_A_AR();
  SAPTDFInts set_B_BS();
  SAPTDFInts set_A_AB();
  SAPTDFInts set_B_AB();
  SAPTDFInts set_A_RB();
  SAPTDFInts set_B_RB();
  SAPTDFInts set_A_AS();
  SAPTDFInts set_B_AS();
  SAPTDFInts set_C_AA();
  SAPTDFInts set_C_AR();
  SAPTDFInts set_C_RR();
  SAPTDFInts set_C_BB();
  SAPTDFInts set_C_BS();
  SAPTDFInts set_C_SS();

  SAPTDFInts set_act_A_AR();
  SAPTDFInts set_act_B_BS();
  SAPTDFInts set_act_C_AR();
  SAPTDFInts set_act_C_BS();

  SAPTDFInts set_H2_BS();
  SAPTDFInts set_H2_AS();
  SAPTDFInts set_H4_AR();
  SAPTDFInts set_H4_RB();
  SAPTDFInts set_Q2_AR();
  SAPTDFInts set_Q6_BS();
  SAPTDFInts set_Q12_AS();
  SAPTDFInts set_Q12_RB();
  SAPTDFInts set_Q13_BS();
  SAPTDFInts set_Q14_AR();

  Iterator get_iterator(long int, SAPTDFInts*, bool alloc=true);
  Iterator set_iterator(int, SAPTDFInts*, bool alloc=true);

  Iterator get_iterator(long int, SAPTDFInts*, SAPTDFInts*, bool alloc=true);
  Iterator set_iterator(int, SAPTDFInts*, SAPTDFInts*, bool alloc=true);

  void read_all(SAPTDFInts*);
  void read_block(Iterator *, SAPTDFInts *);
  void read_block(Iterator *, SAPTDFInts *, SAPTDFInts *);

  void ind20rA_B();
  void ind20rB_A();
  void ind20rA_B_aio();
  void ind20rB_A_aio();

  void get_denom();

  void theta_ar();
  void theta_bs();
  void test_theta();

  void arbs();
  void v1();
  void h1();
  double h2();
  void h3();
  double h4();
  void q1();
  double q2();
  void q3();
  void q5();
  double q6();
  void q7();
  void q10();
  void q11();
  void q12();
  double q13();
  double q14();

protected:
  bool no_response_;
  bool aio_cphf_;
  bool aio_dfints_;

  int maxiter_;
  double e_conv_;
  double d_conv_;

  double e_elst10_;
  double e_exch10_;
  double e_exch10_s2_;
  double e_ind20_;
  double e_exch_ind20_;
  double e_disp20_;
  double e_exch_disp20_;

  double e_disp20_ss_;
  double e_disp20_os_;
  double e_exch_disp20_ss_;
  double e_exch_disp20_os_;

  double e_sapt0_;
  double e_sapt0_scs_;

  double **wBAR_;
  double **wABS_;

public:
  SAPT0(Options& options, boost::shared_ptr<PSIO> psio, 
    boost::shared_ptr<Chkpt> chkpt);
  virtual ~SAPT0();

  virtual double compute_energy();

  void elst10();
  void exch10();
  void exch10_s2();
  void ind20();
  void ind20r();
  void exch_ind20A_B();
  void exch_ind20B_A();
  void disp20();
  void exch_disp20_n4();
  void exch_disp20_n5();

};

struct SAPTDFInts {

  bool dress_;
  bool dress_disk_;
  bool active_;

  int i_length_;
  int j_length_;
  int ij_length_;
  int i_start_;
  int j_start_;

  double **B_p_;
  double **B_d_;

  int filenum_;
  const char *label_;

  psio_address next_DF_;

  SAPTDFInts() { next_DF_ = PSIO_ZERO; B_p_ = NULL; B_d_ = NULL; };
  ~SAPTDFInts() {
    if (B_p_ != NULL) free_block(B_p_);
    if (B_d_ != NULL) free_block(B_d_); };
  void rewind() { next_DF_ = PSIO_ZERO; };
  void clear() { free_block(B_p_); B_p_ = NULL; next_DF_ = PSIO_ZERO; };
  void done() {
    free_block(B_p_); if (dress_) free_block(B_d_);
    B_p_ = NULL; B_d_ = NULL; };
};

struct Iterator {

  int num_blocks;
  int* block_size;

  int curr_block;
  long int curr_size;

  ~Iterator() { free(block_size); };
  void rewind() { curr_block = 1; curr_size = 0; };
};

}}

#endif
