/*! \file
    \ingroup OPTKING
    \brief SALC.CC : constructor functions for salc_set
*/

#include "salc.h"

namespace psi { namespace optking {

void salc_class::set_label(char *new_label) { strcpy(label,new_label); }

salc_set :: salc_set() {
  int i,j,k,num_type,a,b,count;
  int coeff_length, simple_length;
  int label_pos = 0;
  int simple_pos = 1;
  int coeff_pos = 2;
  double f;
  char *buffer, *keyword;

  keyword =  new char[MAX_LINELENGTH];
  buffer = new char[MAX_LINELENGTH];
  salc_array = new salc_class[MAX_SALCS];
  name = new char[MAX_LINELENGTH];
  strcpy(name,"Symmetry Adapted Internal Coordinates");

  count = -1;
  int count_type;
  num=0;
  for (k=0;k<2;++k) {
    if (k == 0)
      strcpy(const_cast<char *>(keyword),"SYMM");
    if (k == 1)
      strcpy(const_cast<char *>(keyword),"ASYMM");
    if (ip_exist(const_cast<char *>(keyword),0)) {
      num_type = 0;
      ip_count(const_cast<char *>(keyword),&num_type,0);
      //fprintf(outfile,"num_type: %d\n",num_type);
      num += num_type;
      count_type = -1;
      for (i=0;i<num_type;++i) {
        ++count_type;
        ++count;
        label_pos = 0;
        simple_pos = 1;
        coeff_pos = 2;
        a = 0;
        ip_count(const_cast<char *>(keyword),&a,1,count_type);
        //fprintf(outfile,"for intco %d, a=%d\n",i,a);
        if (a < 4)  {
          if ( ip_string(const_cast<char *>(keyword),&buffer,2,i,label_pos) != 0 ) {
            --simple_pos;
            --coeff_pos;
          }
          salc_array[count].set_label(buffer);
          //fprintf(outfile,"buffer: %s\n",buffer);
          simple_length = 0;
          ip_count(const_cast<char *>(keyword),&simple_length,2,i,simple_pos);
          salc_array[count].set_length(simple_length);
          for (j=0;j<simple_length;++j) {
            ip_data(const_cast<char *>(keyword),"%d",&a,3,i,simple_pos,j);
            salc_array[count].set_simple(j,a);
          }
          coeff_length = 0;
          ip_count(const_cast<char *>(keyword),&coeff_length,2,i,coeff_pos);
          if (coeff_length == 0) {
            for (j=0;j<simple_length;++j)
              salc_array[count].set_coeff(j,1.0);
          }
          else if (coeff_length == simple_length) {
            for (j=0;j<coeff_length;++j) {
              ip_data(const_cast<char *>(keyword),"%lf",&f,3,i,coeff_pos,j);
              salc_array[count].set_coeff(j,f);
            }
          }
          else {
            fprintf(outfile,"Coefficient length does not match salc length\n");
            exit(2);
          }
        }
        else {
          fprintf(outfile,"%s vector %d has too many parts\n",keyword,count+1);
          exit(2);
        }
      }
    }
  }

  /* computing normalization constants */
  //fprintf(outfile,"num: %d\n",num);

  double sum;
  for (i=0;i<num;++i) {
    sum = 0.0;
    for (j=0;j<salc_array[i].get_length();++j)
      sum += SQR(salc_array[i].get_coeff(j));
    salc_array[i].set_prefactor(1.0/sqrt(sum));
  }

  delete [] keyword;
  delete [] buffer;

  return;
}

salc_set :: salc_set(const char *keyword)
{
  int i,j,a,b;
  int coeff_length, simple_length;
  int label_pos = 0;
  int simple_pos = 1;
  int coeff_pos = 2;
  double f;
  char *buffer;

  buffer = (char*) malloc(sizeof(char)*MAX_LINELENGTH);
  salc_array = new salc_class[MAX_SALCS];
  name = new char[MAX_LINELENGTH];

  strcpy(name,keyword);

  num=0;
  ip_count(const_cast<char *>(keyword),&num,0);
  for (i=0;i<num;++i) {
    label_pos = 0;
    simple_pos = 1;
    coeff_pos = 2;
    a = 0;
    ip_count(const_cast<char *>(keyword),&a,1,i);
    if (a < 4)  {
      buffer[0] = '\0';
      if ( ip_string(const_cast<char *>(keyword),&buffer,2,i,label_pos) != 0 ) {
        --simple_pos;
        --coeff_pos;
      }
      salc_array[i].set_label(buffer);
      simple_length = 0;
      ip_count(const_cast<char *>(keyword),&simple_length,2,i,simple_pos);
      salc_array[i].set_length(simple_length);
      for (j=0;j<simple_length;++j) {
        ip_data(const_cast<char *>(keyword),"%d",&a,3,i,simple_pos,j);
        salc_array[i].set_simple(j,a);
      }
      coeff_length = 0;
      ip_count(const_cast<char *>(keyword),&coeff_length,2,i,coeff_pos);
      if (coeff_length == 0) {
        for (j=0;j<simple_length;++j)
          salc_array[i].set_coeff(j,1.0);
      }
      else if (coeff_length == simple_length) {
        for (j=0;j<coeff_length;++j) {
          ip_data(const_cast<char *>(keyword),"%lf",&f,3,i,coeff_pos,j);
          salc_array[i].set_coeff(j,f);
        }
      }
      else {
        fprintf(outfile,"Coefficient length does not match salc length\n");
        exit(2);
      }
    }
    else {
      fprintf(outfile,"%s vector %d has too many parts\n",keyword,i+1);
      exit(2);
    }
  }

  /* computing normalization constants */
  double sum;
  for (i=0;i<num;++i) {
    sum = 0.0;
    for (j=0;j<salc_array[i].get_length();++j)
      sum += SQR(salc_array[i].get_coeff(j));
    salc_array[i].set_prefactor(1.0/sqrt(sum));
  }

  free(buffer);

  return;
}

}} /* namespace psi::optking */

