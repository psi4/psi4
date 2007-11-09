
/*!
  \file slaterdset.h
  Edward Valeev, June 2002
*/

#ifndef _psi3_src_lib_libqt_slaterdset_h_
#define _psi3_src_lib_libqt_slaterdset_h_

#include <libpsio/psio.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
  String is a string
*/
typedef struct {
  int index;
  short int *occ;     /* Orbital indices in QT order */
                      /* CDS: I'm gonna use Pitzer order actually */
} String;

/*!
  StringSet is a set of strings
*/
typedef struct {
  String *strings;
  int size;
  int nelec;
  int nfzc;
  short int *fzc_occ;
} StringSet;

void stringset_init(StringSet *stringset, int size, int nelec, int nfzc,
  short int *frozen_occ);
void stringset_delete(StringSet *stringset);
void stringset_add(StringSet *stringset, int index, unsigned char *Occ);
void stringset_write(ULI unit, char *prefix, StringSet *sset);
void stringset_read(ULI unit, char *prefix, StringSet **sset);
void stringset_reindex(StringSet *stringset, short int* mo_map);

/*!
  SlaterDet is a Slater determinant
*/
typedef struct {
  int index;
  int alphastring, betastring;
} SlaterDet;

/*!
  SlaterDetSet is a set of Slater determinants
*/
typedef struct {
  int size;
  SlaterDet *dets;
  StringSet *alphastrings, *betastrings;
} SlaterDetSet;

void slaterdetset_init(SlaterDetSet *sdset, int size, StringSet *alphastrings, StringSet *betastrings);
void slaterdetset_delete(SlaterDetSet *sdset);
void slaterdetset_delete_full(SlaterDetSet *sdset);
void slaterdetset_add(SlaterDetSet *sdset, int index, int alphastring, int betastring);
void slaterdetset_write(ULI unit, char *prefix, SlaterDetSet *sdset);
void slaterdetset_read(ULI unit, char *prefix, SlaterDetSet **sdset);

/*!
  SlaterDetVector is a vector in the space of determinants
*/
typedef struct {
  int size;
  SlaterDetSet *sdset;
  double *coeffs;
} SlaterDetVector;

void slaterdetvector_init(SlaterDetVector *sdvector, SlaterDetSet *sdset);
void slaterdetvector_delete(SlaterDetVector *sdvector);
void slaterdetvector_delete_full(SlaterDetVector *sdvector);
void slaterdetvector_set(SlaterDetVector *sdvector, double *coeffs);

/*!
  slaterdetvector_write()

  Writes a vector in the space of Slater determinants along with the set itself to a PSIO file.

  \param
  ULI unit:  the unit number of the <em>uninitialized</em> PSIO file
  
  \param
  char *prefix:  the prefix to be used for PSIO keys of the file entries
  associated with the vector

  \param
  SlaterDetVector *vector:  the vector

*/
void slaterdetvector_write(ULI unit, char *prefix, SlaterDetVector *vector);


/*!
  slaterdetset_write_vect()

  Writes a vector in the space of Slater determinants to a PSIO file.  
  Similar to slaterdetvector_write() except that this one writes only 
  the vector, not the string and determinant information.  This does
  not actually depend on the presence of a SlaterDetVector object, so it 
  is called a SlaterDetSet function.

  \param
  ULI unit:  the unit number of the <em>uninitialized</em> PSIO file
  
  \param
  char *prefix:  the prefix to be used for PSIO keys of the file entries
  associated with the vector

  \param
  double *coeffs: the coefficient array

  \param 
  int size: the number of coefficients to write

  \param 
  int vectnum: an integer to be written in the TOC key.  Usually 
  corresponds to which root is being written.  Start numbering from 0.

*/
void slaterdetset_write_vect(ULI unit, char *prefix,
  double *coeffs, int size, int vectnum);


/*!
  slaterdetvector_read()

  Reads a vector in the space of Slater determinants from a PSIO file. 
  Complementary and analogous to slaterdetvector_write()

*/
void slaterdetvector_read(ULI unit, char *prefix, SlaterDetVector **vector);

/*
  slaterdetset_read_vect()
                                                                                
  Reads a vector in the space of Slater determinants from a PSIO file.
  Similar to slaterdetvector_read() except that this one reads only
  the vector, not the string and determinant information.  This does
  not actually depend on the presence of a SlaterDetVector object, so it
  is called a SlaterDetSet function.

  \param
  ULI unit:  the unit number of the <em>uninitialized</em> PSIO file
                                                                                
  \param
  char *prefix:  the prefix to be used for PSIO keys of the file entries
  associated with the vector
                                                                                
  \param
  double *coeffs:  the vector

  \param
  int size: the length of the coefficient array 
                                                                                
  \param
  int vectnum: which vector number to read (start from 0)

*/
void slaterdetset_read_vect(ULI unit, char *prefix, double *coeffs,
  int size, int vectnum);

 
#define STRINGSET_KEY_SIZE "StringSet Size"
#define STRINGSET_KEY_NELEC "StringSet Num. of Electrons"
#define STRINGSET_KEY_NFZC "StringSet Num. of Frozen DOCCs"
#define STRINGSET_KEY_FZC_OCC "StringSet Frz Core Occs"
#define STRINGSET_KEY_STRINGS "StringSet Strings"
#define SDSET_KEY_SIZE "SlaterDetSet Size"
#define SDSET_KEY_DETERMINANTS "SlaterDetSet Determinants"
#define SDSET_KEY_ALPHASTRINGS "Alpha Strings"
#define SDSET_KEY_BETASTRINGS "Beta Strings"
#define SDVECTOR_KEY_VECTOR "Vector"

#ifdef __cplusplus
}
#endif

#endif  /* _psi3_src_lib_libqt_slaterdset_h_ */

