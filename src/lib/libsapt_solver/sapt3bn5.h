/*
 *  Header file for SAPT0 objects 
 *  Created by Rob Parrish on 07/21/2010
 *
 */

#ifndef SAPT3BN5_H
#define SAPT3BN5_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <psi4-dec.h>

#include "sapt3b.h"

using namespace psi;

namespace psi { namespace sapt {

class SAPT3BN5 : public SAPT3B {
private:
    double exch_s2(double **,double **,double **,double **,int,int,int,int,
      int);
    double exch_s3(double **,double **,double **,double **,double **,int,char*,
      int,char*,int,int,int,int,int);
    double exch_s4(double **,double **,double **,double **,double **,double **,
      int,char*,int,char*,int,int,int,int,int);
    double exch_s4_1(double **, double **, int, int, int);
    double exch_s4_2(double **, double **, double **, int, int, int, int);
    double exch_s4_3(double **, double **, double **, int, int, int, int);
    double exch_s4_4(double **, double **, double **, int, int, int, int);
    double exch_s4_5(double **, double **, int, char *, int, char *, int,
      int, int, int, int);
    double exch_s4_6(double **, double **, int, char *, int, char *, int,
      int, int, int, int);

    double ind110_1(double **,double **,double **,double **,int,int);
    double ind210_0(double **, double **, double **, double **, double **,
      double **, double **, int, char *, int, char *, int, int,int, int);
    double ind210_1(double **, double **, double **, double **, int, int);
    double ind210_2(double **, double **, double **, int, int);
    double ind111_1(double **, double **, int, char *, int, char *, int, int,
      int, int);

    double exch_ind200_s2_0(double **,double **,double **,double **,double **,
      double **,int,char*,int,char*,int,int,int,int,int);
    double exch_ind110_s2_0(double **, double **, double **, double **, 
      double **, double **, double **, double **, double **, double **,
      double **, double **, double **, double **, int, char *, int, char *, 
      int, char *, int, int, int, int, int, int);
    double exch_ind110_s2_1(double **, double **, double **, double **, 
      double **, int, int, int, int, int, int);
    double exch_ind110_s2_2(double **, double **, double **, double **, 
      double **, int, int, int, int, int, int);
    double exch_ind110_s2_3(double **, double **, double **, int, char *, int,
      char *, int, int, int, int, int);
    double exch_ind200_s3_0(double **, double **, double **, double **, 
      double **, double **, double **, double **, int, char *, int, char *, 
      int, int, int, int, int);
    double exch_ind200_s3_1(double **, double **, double **, double **, 
      double **, double **, int, int, int, int, int);
    double exch_ind200_s3_2(double **, double **, double **, double **, int,
      char *, int, char *, int, int, int, int, int);
    double exch_ind110_s3_0(double **, double **, double **, double **, 
      double **, double **, double **, double **, double **, double **,
      double **, double **, double **, double **, int, char *, int, char *, 
      int, char *, int, int, int, int, int, int);
    double exch_ind110_s3_1(double **, double **, double **, double **, 
      double **, double **, int, int, int, int, int, int);
    double exch_ind110_s3_2(double **, double **, double **, double **, 
      double **, double **, int, int, int, int, int, int);
    double exch_ind110_s3_3(double **, double **, double **, double **, 
      double **, int, char *, int, char *, int, int, int, int, int, int);

    double exch_disp200_s2_0(double **, double **, double **, double **, 
      char *, char *, int, char *, int, char *, char *, char, int, int, int, 
      int, int);
    double exch_disp200_s2_1(double **, char *, int, char *, int, int, int);
    double exch_disp200_s2_2(double **, double **, double **, double **, 
      char *, char, int, int, int, int, int);
    double exch_disp110_s2_0(double **, double **, double **, double **,
      double **, double **, double **, char *, char, char *, char, char *, 
      char *, int, char *, int, char *, int, int, int, int, int, int);
    double exch_disp110_s2_1(double **, char *, int, char *, int, int, int, 
      int);
    double exch_disp110_s2_2(double **, double **, double **, char *, char,
      int, int, int, int, int, int);
    double exch_disp200_s3_0(double **, double **, double **, double **, 
      double **, double **, double **, double **, char *, char *, int, char *,
      int, char *, char *, char, char, int, int, int, int, int, int);
    double exch_disp200_s3_1(double **, double **, double **, double **, 
      double **, double **, char *, int, char *, int, int, int, int, int, int);
    double exch_disp200_s3_2(double **, double **, double **, double **,
      double **, double **, int, char *, int, char *, char *, char, int, int,
      int, int, int, int);
    double exch_disp200_s3_3(double **, double **, double **, double **, 
      double **, double **, int, char *, int, char *, char *, char, int, int, 
      int, int, int, int);
    double exch_disp200_s3_4(double **, double **, double **, double **, 
      double **, double **, double **, char *, char, int, int, int, int, int,
      int);
    double exch_disp110_s3_0(double **, double **, double **, double **, 
      double **, double **, double **, double **, double **, double **, char *,
      char *, int, char *, int, char *, int, char *, char *, char, char *, 
      char, int, int, int, int, int, int);
    double exch_disp110_s3_1(double **, double **, double **, double **, 
      double **, double **, char *, int, char *, int, int, int, int, int, int);
    double exch_disp110_s3_2(double **, double **, double **, double **, 
      double **, double **, int, char *, int, char *, char *, char, int, int, 
      int, int, int, int);
    double exch_disp110_s3_3(double **, double **, double **, double **, 
      double **, double **, double **, char *, char, int, int, int, int, int, 
      int);
    double exch_disp110_s3_4(double **, double **, double **, double **, 
      double **, double **, double **, char *, char, int, int, int, int, int,
      int);

    double disp111_1(char *, char *, int, int);

    double ind_disp210_0(double **, char *, int, char *, char *, char *, char,
      double **, double **, int, int, int, int);
    double ind_disp210_1(double **, char *, int, char *, char *, int, int);
    double ind_disp210_2(char *, char, double **, double **, int, int, int, 
      int);

protected:
    virtual void print_header();

    void exch100_s2();
    void exch100_s3();
    void exch100_s4();
    void ind110();
    void ind210();
    void ind111();
    void exch_ind200_s2();
    void exch_ind110_s2();
    void exch_ind200_s3();
    void exch_ind110_s3();
    void exch_disp200_s2();
    void exch_disp110_s2();
    void exch_disp200_s3();
    void exch_disp110_s3();
    void disp111();
    void ind_disp210();

    virtual double print_results();

public:
    SAPT3BN5(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    virtual ~SAPT3BN5();

    virtual double compute_energy();

};

}}

#endif
