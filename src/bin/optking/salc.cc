/*! \file
    \ingroup OPTKING
    \brief SALC.CC : member functions for Salc_set
*/

#include "salc.h"

namespace psi { namespace optking {

Salc_set :: Salc_set(string key_in, int irrep_in, int na,
  int nsimp) throw(bad_intco_io)
{
  int i,j,b, errcod, nrow=0, a, simple_length;
  char *buffer;
  string lbl_one;
  vector<int> simples_one;
  vector<double> coeff_one;

  // set Salc_set scalar members
  irrep = irrep_in;
  keyword = key_in;
  natom = na;
  nsimples = nsimp;

  char *key_ip = c_string(keyword);
  buffer = new char[MAX_LINELENGTH];

  if (!ip_exist(key_ip,0))
    return;

  if (ip_count(key_ip,&nrow,0) != IPE_OK)
    throw bad_intco_io(key_in,1);

  for (i=0;i<nrow;++i) {
    a = 0;
    if (ip_count(key_ip,&a,1,i) != IPE_OK)
      throw bad_intco_io(key_in, i+1);
    if (a<2 || a>3)
      throw bad_intco_io(key_in, i+1);

    //get label
    buffer[0] = '\0';
    if (ip_string(key_ip,&buffer,2,i,0) != IPE_OK)
      throw bad_intco_io(key_in, i+1);

    // read simples
    simple_length = 0;
    if (ip_count(key_ip,&simple_length,2,i,1) != IPE_OK)
      throw bad_intco_io(key_in, i+1);

    simples_one.resize(simple_length);
    coeff_one.resize(simple_length);

    for (j=0;j<simple_length;++j)
      if (ip_data(key_ip,"%d",&(simples_one.at(j)),3,i,1,j) != IPE_OK)
        throw bad_intco_io(key_in, i+1);

    if (a==3) { //read coefficients
      if (ip_count(key_ip,&b,2,i,2) != IPE_OK) 
        throw bad_intco_io(key_in, i+1);
      if (b != simple_length)
        throw bad_intco_io(key_in, i+1);
      for (j=0;j<simple_length;++j)
        if (ip_data(key_ip,"%lf",&(coeff_one.at(j)),3,i,2,j) != IPE_OK)
          throw bad_intco_io(key_in, i+1);
    }
    else if (a==2) {
      for (j=0;j<simple_length;++j)
        coeff_one.at(j) = 1.0;
    }

    lbl_one = buffer;

    { // use block to make lsalc go out of scope and get destructed
     One_salc lsalc(lbl_one, simple_length, simples_one, coeff_one);
     salc.push_back(lsalc);
    }

    simples_one.clear();
    coeff_one.clear();
 }

 delete [] buffer;
 free(key_ip);
 printf("leaving salc_set constructor\n");
 return;
}

void One_salc::set_prefactor(void) {
  int j;
  double sum =0.0;
  for (j=0; j<length; ++j)
    sum += sqr(coeff.at(j));

  prefactor = 1.0/sqrt(sum);
}

// default arguments = 0
void One_salc::print(int print_flag) const {
  int i, col = 0;
  fprintf(outfile,"    (");
  fprintf(outfile,"\"%s\"", lbl.c_str());
  fprintf(outfile,"  (");
  for (col=0, i=0;i<length;++i, ++col) {
    fprintf(outfile," %d", simple.at(i));
    if ((col == 9) && (i!=length-1)) { fprintf(outfile,"\n    "); col=0; }
  }
  fprintf(outfile,")  (");
  for (col=0, i=0;i<length;++i, ++col) {
    fprintf(outfile,"%8.4f", coeff.at(i));
    if ((col == 9) && (i!=length-1)) { fprintf(outfile,"\n    "); col=0;}
  }
  fprintf(outfile,"))\n");
}

//default argument is 0
void Salc_set::print(int print_flag) const {
  if (salc.size()) {
    fprintf(outfile,"\n  %-s = (\n", keyword.c_str());
    for (int i=0; i<salc.size(); ++i)
      salc.at(i).print(print_flag);
    fprintf(outfile,"  )\n");
  }
  return;
}

void Salc_set::build_matrix(void) {
  int simple, i, j;
  double coeff, pre;

  if (matrix_present) free_block(matrix);
  block_matrix(salc.size(),nsimples);
  matrix_present = 1;

  for (i=0; i<salc.size(); ++i) {
    pre = salc.at(i).prefactor;
    for (j=0;j<salc.at(i).length;++j) {
      simple = salc.at(i).simple.at(j)-1;
      coeff = salc.at(i).coeff.at(j);
      matrix[i][simple] += pre * coeff;
      matrix[i][simple] += pre * coeff;
    }
  }
}

void Salc_set::print_matrix(void) const {
  if (matrix_present)
    print_mat(matrix,salc.size(),nsimples,outfile);
}

Salc_set::Salc_set(const Salc_set & ss) {
  printf("using salc_set copy constructor -not yet tested\n");
  keyword = ss.keyword;
  irrep = ss.irrep;
  natom = ss.natom;
  salc = ss.salc; //deep copy works?
  if (ss.matrix_present) {
    matrix = block_matrix(salc.size(),nsimples);
    for (int i=0; i<salc.size(); ++i)
      for (int j=0; j<nsimples; ++j)
        matrix[i][j] = ss.matrix[i][j];
    matrix_present = 1;
  }
  else
    matrix_present = 0;
}

Salc_set::~Salc_set() {
  printf("salc_set destructor called\n");
  if (matrix_present)
    free_block(matrix);
  salc.clear(); // calls One_salc destructors
}

inline char * c_string(const string s) {
  char * buf = new char [s.size()+1];
  for (int i=0; i<s.size(); ++i) buf[i] = s[i];
  buf[s.size()] = '\0';
  return buf;
}

}}
