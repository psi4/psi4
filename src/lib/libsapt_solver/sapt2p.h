#ifndef SAPT2p_H
#define SAPT2p_H

#include "sapt2.h"

namespace psi { namespace sapt {

class SAPT2p : public SAPT2 {
private:
  virtual void print_header();
  virtual void print_results();

protected:
  double e_disp21_;
  double e_disp22sdq_;
  double e_disp22t_;
  double e_est_disp22t_;
  double e_sapt2p_;


  void gARARxtARBS(int, const char *, const char, int, const char *, 
    const char *, const char *, int, int, int, int, int, int, int, 
    const char *);

  double disp21_1(int, const char *, const char *, int, int, int, int);
  double disp21_2(int, const char *, const char *, int, int);

  double disp211();
  double disp220s(int, const char *, const char *, int, const char *, 
    const char *, int, int, int);
  double disp220d_1(int, const char *, const char *, int, const char *,
    int, int, int);
  double disp220d_2(int, const char *, const char *, int, const char *,
    int, int, int, int, int, int, double *, double *, const char);
  double disp220q_1(int, const char *, const char *, const char *, int, int);
  double disp220q_2(int, const char *, const char *, const char *, int, 
    const char *, int, int, int);
  double disp220q_3(int, const char *, const char *, const char, int, 
    const char *, int, int, int, int, int, int);
  double disp220q_4(int, const char *, const char *, const char, int, 
    const char *, int, int, int, int, int, int);

  double disp220t(int, const char *, const char *, const char *, int, 
    const char *, int, const char *, int, int, int, int, int, int, double *, 
    double *, const char* = NULL);

  // CCD Dispersion Values

  double e_disp22t_ccd_;
  double e_est_disp22t_ccd_;
  
  // CCD Dispersion Parameters
  bool ccd_disp_;
  int ccd_maxiter_;
  int min_ccd_vecs_;
  int max_ccd_vecs_;
  double ccd_e_conv_;
  double ccd_t_conv_;

  // CCD Dispersion Methods
  void r_ccd_prep(char *, char *, char *, char *, char *, char *, char *, 
    char *, char *, char *, char *, char *, char *, char *, char *, char *, 
    char *, int, char *, int, char *, double *, double *, int, int, int, 
    int, int, int);
  double r_ccd_energy(char *, char *, int, int, int, int);
  double r_ccd_iterate(char *, char *, char *, char *, char *, char *, 
    char *, char *, char *, char *, char *, char *, double *, double *, 
    int, int, int, int, int, int);
  double r_ccd_amplitudes(char *, char *, char *, char *, char *, char *, 
    char *, char *, char *, char *, char *, double *, double *, int, int, 
    int, int, int, int);

  void s_ccd_prep(char *, char *, char *, char *, char *, char *, char *, 
    double *, int, int, int, int, int, int);
  double s_ccd_iterate(char *, char *, char *, char *, char *, char *, 
    char *, char *, char *, char *, char *, char *, char *, char *, char *,
    double *, int, int, int, double **, int);
  double s_ccd_amplitudes(char *, char *, char *, char *, char *, char *, 
    char *, char *, char *, char *, char *, char *, char *, char *, char *,
    double *, int, int, int, double **, int);

  void disp_s_prep(char *, char *, char *, char *, int, char *, char *, 
    char *, int, char *, double *, int, int, int, int, int, int);
  void natural_orbitalify_ccd();

  void ccd_prep(char *, char *, char *, char *, char *, char *, char *, 
    char *, char *, int, char *, char *, char *, double *, int, int, int, 
    int, double **, char *);
  double ccd_energy(char *, char *, int, int);
  void ccd_iterate(char *, char *, char *, char *, char *, char *, char *,
    char *, char *, char *, double *, int, int, int, double **, int);
  double ccd_amplitudes(char *, char *, char *, char *, char *, char *, 
    char *, char *, char *, char *, double *, int, int, int, double **, int);

  double **vvvv_ccd(char *, char *, char *, int, int, double **, int, int);

  double **read_IJKL(int, char *, int, int);
  void write_IJKL(double **, int, char *, int, int);
  
  double disp220tccd(int, char *, int, char *, char *, int, char *, int, char *,
    char *, double *, double *, int, int, int, int, int, int);

public:
  SAPT2p(Options& options, boost::shared_ptr<PSIO> psio, 
    boost::shared_ptr<Chkpt> chkpt);
  virtual ~SAPT2p();

  virtual double compute_energy();

  virtual void amplitudes();

  // PT Dispersion

  void disp21();
  void disp22sdq();
  void disp22t();

  // CCD Dispersion

  void disp2ccd();
  void disp22tccd();

};

/**
 * SAPTDIIS is a legacy helper for CCD
 **/
class SAPTDIIS {

private:
    int filenum_;
    char *vec_label_;
    char *err_label_;
    int max_diis_vecs_;

    int diis_file_;
    int vec_length_;

    int curr_vec_;
    int num_vecs_;

    char *get_err_label(int);
    char *get_vec_label(int);

protected:
    boost::shared_ptr<PSIO> psio_;

public:
    SAPTDIIS(int, char *, char *, int, int, boost::shared_ptr<PSIO>);
    ~SAPTDIIS();

    void store_vectors();
    void get_new_vector();
};

}}

#endif
