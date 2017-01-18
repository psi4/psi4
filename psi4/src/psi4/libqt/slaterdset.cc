/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */


/*!
** \file
** \brief Utility functions for importing/exporting sets of Slater determinants
**   from the CI codes
** \ingroup QT
**
** Edward Valeev, June 2002
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "psi4/libpsio/psio.h"
#include "psi4/libciomr/libciomr.h"
#include "slaterdset.h"
#include "psi4/psi4-dec.h"
namespace psi {

#define PSIO_INIT if (!psio_state()) { \
    psio_init(); \
    need_to_init_psio = 1; \
  }

#define PSIO_OPEN(u,n) if (!psio_open_check(u)) { \
    psio_open((u),n); \
    unit_opened = 0; \
  }

#define PSIO_CLOSE(u) if (!unit_opened) \
    psio_close((u),1);

#define PSIO_DONE if (need_to_init_psio) \
    psio_done();


/*!
** stringset_init(): Initialize a set of alpha/beta strings
**
** \param sset       = pointer to StringSet (contains occs in Pitzer order)
** \param size       = number of strings
** \param nelec      = number of electrons
** \param ndrc       = number of dropped core orbitals
** \param drc_occ    = array of dropped occupied orbitals (Pitzer numbering!)
**
** Returns: none
** \ingroup QT
*/
void stringset_init(StringSet *sset, int size, int nelec, int ndrc,
  short int *drc_occ)
{
  int i;

  sset->size = size;
  sset->nelec = nelec;
  sset->ndrc = ndrc;
  sset->strings = (String *) malloc(size*sizeof(String));
  memset(sset->strings,0,size*sizeof(String));
  if (ndrc > 0) {
    sset->drc_occ = (short int *) malloc(ndrc * sizeof(short int));
    for (i=0; i<ndrc; i++) {
      sset->drc_occ[i] = drc_occ[i];
    }
  }
}


/*!
** stringset_delete(): Delete a StringSet
**
** \param sset = pointer to StringSet to delete
**
** Returns: none
**
** \ingroup QT
*/
void stringset_delete(StringSet *sset)
{
  if (sset->ndrc > 0) free(sset->drc_occ);
  sset->size = 0;
  sset->nelec = 0;
  sset->ndrc = 0;
  if (sset->strings) free(sset->strings);
  sset->strings = NULL;
}

/*!
** stringset_add(): Add a string (in Pitzer order, given by Occ) to
** the StringSet, writing to position index.
**
** \param sset  = StringSet to add to
** \param index = location in StringSet to add to
** \param Occ   = orbital occupations (Pitzer order) of the string to add
**
** Returns: none
** \ingroup QT
*/
void stringset_add(StringSet *sset, int index, unsigned char *Occ)
{
  int i;
  int nact = sset->nelec - sset->ndrc;
  String *s;

  if (index < sset->size && index >= 0) {
    s = sset->strings + index;
  }
  s->index = index;
  s->occ = (short int*) malloc(nact*sizeof(short int));
  for(i=0;i<nact;i++)
    s->occ[i] = Occ[i];
}

/*!
** stringset_reindex(): Remap orbital occupations from one ordering to
** another.
**
** \param sset   = pointer to StringSet
** \param mo_map = mapping array from original orbital order to new order
**
** Returns: none
** \ingroup QT
*/
void stringset_reindex(StringSet* sset, short int* mo_map)
{
  int s, mo, core;
  short int* occ;
  int nstrings = sset->size;
  int nact = sset->nelec - sset->ndrc;

  for (core=0; core<sset->ndrc; core++) {
    sset->drc_occ[core] = mo_map[sset->drc_occ[core]];
  }

  for(s=0; s<nstrings; s++) {
    occ = (sset->strings + s)->occ;
    for(mo=0; mo<nact; mo++)
      occ[mo] = mo_map[occ[mo]];
  }
}

/*!
** stringset_write(): Write a stringset to a PSI file
**
** \param unit   = file number to write to
** \param prefix = prefix string to come before libpsio entry keys
** \param sset   = pointer to StringSet to write
**
** Returns: none
** \ingroup QT
*/
void stringset_write(ULI unit, const char *prefix, StringSet *sset)
{
  int i, size, nact;
  int need_to_init_psio = 0;
  int unit_opened = 1;
  char *size_key, *nelec_key, *ndrc_key, *strings_key, *drc_occ_key;
  psio_address ptr;

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  size_key = (char *) malloc(strlen(prefix) + strlen(STRINGSET_KEY_SIZE) + 3);
  sprintf(size_key,":%s:%s",prefix,STRINGSET_KEY_SIZE);
  nelec_key = (char *) malloc(strlen(prefix) + strlen(STRINGSET_KEY_NELEC) + 3);
  sprintf(nelec_key,":%s:%s",prefix,STRINGSET_KEY_NELEC);
  ndrc_key = (char *) malloc(strlen(prefix) + strlen(STRINGSET_KEY_NDRC) + 3);
  sprintf(ndrc_key,":%s:%s",prefix,STRINGSET_KEY_NDRC);
  drc_occ_key = (char *) malloc(strlen(prefix) +
    strlen(STRINGSET_KEY_DRC_OCC) + 3);
  sprintf(drc_occ_key,":%s:%s",prefix,STRINGSET_KEY_DRC_OCC);
  strings_key = (char *) malloc(strlen(prefix) + strlen(STRINGSET_KEY_STRINGS) + 3);
  sprintf(strings_key,":%s:%s",prefix,STRINGSET_KEY_STRINGS);

  psio_write_entry( unit, size_key, (char *)&sset->size, sizeof(int));
  psio_write_entry( unit, nelec_key, (char *)&sset->nelec, sizeof(int));
  psio_write_entry( unit, ndrc_key, (char *)&sset->ndrc, sizeof(int));
  if (sset->ndrc) {
    psio_write_entry( unit, drc_occ_key, (char *)sset->drc_occ,
      sset->ndrc*sizeof(short int));
  }

  ptr = PSIO_ZERO;
  size = sset->size;
  nact = sset->nelec - sset->ndrc;
  for(i=0; i<size; i++) {
    psio_write( unit, strings_key, (char *) &(sset->strings[i].index),
      sizeof(int), ptr, &ptr);
    psio_write( unit, strings_key, (char *) sset->strings[i].occ,
      nact*sizeof(short int), ptr, &ptr);
  }

PSIO_CLOSE(unit)
PSIO_DONE

  free(size_key);
  free(nelec_key);
  free(ndrc_key);
  free(strings_key);
  free(drc_occ_key);
}


/*!
** stringset_read(): Read a StringSet from disk
**
** \param unit      = file number to read from
** \param prefix    = prefix string to come before libpsio entry keys
** \param stringset = double pointer to StringSet allocated by this function
**
** Returns: none
** \ingroup QT
*/
void stringset_read(ULI unit, const char *prefix, StringSet **stringset)
{
  int i, size, nelec, ndrc, nact;
  int need_to_init_psio = 0;
  int unit_opened = 1;
  char *size_key, *nelec_key, *ndrc_key, *drc_occ_key, *strings_key;
  short int *drc_occ;
  psio_address ptr;
  StringSet *sset = (StringSet *) malloc(sizeof(StringSet));

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  size_key = (char *) malloc( strlen(prefix) + strlen(STRINGSET_KEY_SIZE) + 3);
  sprintf(size_key,":%s:%s",prefix,STRINGSET_KEY_SIZE);
  nelec_key = (char *) malloc( strlen(prefix) + strlen(STRINGSET_KEY_NELEC)+ 3);
  sprintf(nelec_key,":%s:%s",prefix,STRINGSET_KEY_NELEC);
  ndrc_key = (char *) malloc( strlen(prefix) + strlen(STRINGSET_KEY_NDRC) + 3);
  sprintf(ndrc_key,":%s:%s",prefix,STRINGSET_KEY_NDRC);
  drc_occ_key = (char *) malloc(strlen(prefix) +
    strlen(STRINGSET_KEY_DRC_OCC) + 3);
  sprintf(drc_occ_key,":%s:%s",prefix,STRINGSET_KEY_DRC_OCC);
  strings_key = (char *) malloc( strlen(prefix) + strlen(STRINGSET_KEY_STRINGS)
    + 3);
  sprintf(strings_key,":%s:%s",prefix,STRINGSET_KEY_STRINGS);

  psio_read_entry( unit, size_key, (char *)&size, sizeof(int));
  psio_read_entry( unit, nelec_key, (char *)&nelec, sizeof(int));
  psio_read_entry( unit, ndrc_key, (char *)&ndrc, sizeof(int));
  if (ndrc > 0) {
    drc_occ = (short int *) malloc(ndrc*sizeof(short int));
    psio_read_entry( unit, drc_occ_key, (char *)drc_occ,
      ndrc*sizeof(short int));
  }
  else drc_occ = NULL;

  stringset_init(sset, size, nelec, ndrc, drc_occ);

  nact = nelec - ndrc;
  ptr = PSIO_ZERO;
  for(i=0; i<size; i++) {
    psio_read( unit, strings_key, (char *) &(sset->strings[i].index),
      sizeof(int), ptr, &ptr);
    sset->strings[i].occ = (short int*) malloc(nact*sizeof(short int));
    psio_read( unit, strings_key, (char *) sset->strings[i].occ,
      nact*sizeof(short int), ptr, &ptr);
  }

PSIO_CLOSE(unit)
PSIO_DONE

  free(size_key);
  free(nelec_key);
  free(ndrc_key);
  free(drc_occ_key);
  free(strings_key);
  if (ndrc > 0) free(drc_occ);
  *stringset = sset;
}


/*!
** slaterdetset_init(): Initialize a Slater Determinant Set
**
** \param sdset        = pointer to SlaterDetSet being initialized
** \param size         = number of SlaterDets to be held
** \param alphastrings = pointer to StringSet of alpha strings
** \param betastrings  = pointer to StringSet of beta strings
**
** Returns: none
** \ingroup QT
*/
void slaterdetset_init(SlaterDetSet *sdset, int size, StringSet *alphastrings,
  StringSet *betastrings)
{
  sdset->size = size;
  sdset->dets = (SlaterDet *) malloc(size*sizeof(SlaterDet));
  memset(sdset->dets,0,size*sizeof(SlaterDet));
  sdset->alphastrings = alphastrings;
  sdset->betastrings = betastrings;
}

/*!
** slaterdetset_delete(): Delete a Slater Determinant Set.
**
** Does not free the members alphastrings and betastrings.  See also:
**  slaterdetset_delete_full() which does this.
**
** \param sdset = pointer to SlaterDetSet to be de-allocated
**
** Returns: none
** \ingroup QT
*/
void slaterdetset_delete(SlaterDetSet *sdset)
{
  sdset->size = 0;
  if (sdset->dets) {
    free(sdset->dets);
    sdset->dets = NULL;
  }
  sdset->alphastrings = NULL;
  sdset->betastrings = NULL;
}

/*!
** slaterdetset_delete_full(): De-allocate a Slater Determinant Set.
**
** Frees memory including alpha and beta strings.  See
** slaterdetset_delete() for a similar version which does not free the
** alpha and beta strings.
**
** \param sdset = pointer to SlaterDetSet to be de-allocated
**
** Returns: none
** \ingroup QT
*/
void slaterdetset_delete_full(SlaterDetSet *sdset)
{
  sdset->size = 0;
  if (sdset->dets) {
    free(sdset->dets);
    sdset->dets = NULL;
  }
  if (sdset->alphastrings) {
    stringset_delete(sdset->alphastrings);
    sdset->alphastrings = NULL;
  }
  if (sdset->betastrings) {
    stringset_delete(sdset->betastrings);
    sdset->betastrings = NULL;
  }
}

/*!
** slaterdetset_add(): Add info for a particular Slater determinant to
** a SlaterDetSet.
**
** \param sdset       = pointer to SlaterDetSet to add to
** \param index       = location in the set to add to
** \param alphastring = alpha string ID for the new determinant
** \param betastring  = beta string ID for the new determinant
**
** Returns: none
** \ingroup QT
*/
void slaterdetset_add(SlaterDetSet *sdset, int index, int alphastring,
  int betastring)
{
  SlaterDet *det;
  StringSet *alphaset = sdset->alphastrings;
  StringSet *betaset = sdset->betastrings;

  if (index < sdset->size && index >= 0) {
    det = sdset->dets + index;
  }
  det->index = index;
  if (alphastring < alphaset->size && alphastring >= 0)
    det->alphastring = alphastring;
  if (betastring < betaset->size && betastring >= 0)
    det->betastring = betastring;
}

/*!
** slaterdetset_write(): Write a Slater Determinant Set to disk.
**
** \param unit      = file number to write to
** \param prefix    = prefix string to come before libpsio entry keys
** \param sdset     = pointer to SlaterDetSet to write
**
** Returns: none
** \ingroup QT
*/
void slaterdetset_write(ULI unit, const char *prefix, SlaterDetSet *sdset)
{
  int i;
  int need_to_init_psio = 0;
  int unit_opened = 1;
  char *size_key, *set_key;
  char *alphaprefix, *betaprefix;
  psio_address ptr;

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  alphaprefix = (char *) malloc( strlen(prefix) +
    strlen(SDSET_KEY_ALPHASTRINGS) + 2);
  sprintf(alphaprefix,"%s:%s",prefix,SDSET_KEY_ALPHASTRINGS);
  betaprefix = (char *) malloc( strlen(prefix) +
    strlen(SDSET_KEY_BETASTRINGS) + 2);
  sprintf(betaprefix,"%s:%s",prefix,SDSET_KEY_BETASTRINGS);

  stringset_write( unit, alphaprefix, sdset->alphastrings);
  stringset_write( unit, betaprefix, sdset->betastrings);

  free(alphaprefix);
  free(betaprefix);

  size_key = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_SIZE) + 3);
  sprintf(size_key,":%s:%s",prefix,SDSET_KEY_SIZE);
  set_key = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_DETERMINANTS)
    + 3);
  sprintf(set_key,":%s:%s",prefix,SDSET_KEY_DETERMINANTS);

  psio_write_entry( unit, size_key, (char *)&sdset->size, sizeof(int));
  psio_write_entry( unit, set_key, (char *)sdset->dets,
    sdset->size*sizeof(SlaterDet));

PSIO_CLOSE(unit)
PSIO_DONE

  free(size_key);
  free(set_key);
}

/*!
** slaterdetset_read(): Read a Slater Determinant Set
**
** \param unit      = file number of the PSIO file
** \param prefix    = prefix string to come before libpsio entry keys
** \param sdset     = pointer to SlaterDetSet to read into
**
** Returns: none
** \ingroup QT
*/
void slaterdetset_read(ULI unit, const char *prefix, SlaterDetSet **slaterdetset)
{
  int i, size;
  int need_to_init_psio = 0;
  int unit_opened = 1;
  char *size_key, *set_key;
  char *alphaprefix, *betaprefix;
  psio_address ptr;
  StringSet *alphastrings, *betastrings;
  SlaterDetSet *sdset = (SlaterDetSet *) malloc(sizeof(SlaterDetSet));

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  alphaprefix = (char *) malloc( strlen(prefix) +
    strlen(SDSET_KEY_ALPHASTRINGS) + 2);
  sprintf(alphaprefix,"%s:%s",prefix,SDSET_KEY_ALPHASTRINGS);
  betaprefix = (char *) malloc( strlen(prefix) +
    strlen(SDSET_KEY_BETASTRINGS) + 2);
  sprintf(betaprefix,"%s:%s",prefix,SDSET_KEY_BETASTRINGS);

  stringset_read( unit, alphaprefix, &alphastrings);
  stringset_read( unit, betaprefix, &betastrings);

  free(alphaprefix);
  free(betaprefix);

  size_key = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_SIZE) + 3);
  sprintf(size_key,":%s:%s",prefix,SDSET_KEY_SIZE);
  set_key = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_DETERMINANTS)
    + 3);
  sprintf(set_key,":%s:%s",prefix,SDSET_KEY_DETERMINANTS);

  psio_read_entry( unit, size_key, (char *)&size, sizeof(int));
  slaterdetset_init(sdset,size,alphastrings,betastrings);
  psio_read_entry( unit, set_key, (char *)sdset->dets,
    sdset->size*sizeof(SlaterDet));

PSIO_CLOSE(unit)
PSIO_DONE

  free(size_key);
  free(set_key);

  *slaterdetset = sdset;
}


/*!
** slaterdetvector_init(): Initialize a vector of coefficients
**   corresponding to a Slater Determinant set
**
** \param sdvector = pointer to SlaterDetVector to initialize (coeffs
**   member will be allocated)
** \param sdset    = pointer to SlaterDetSet the vector is associated with
**
** Returns: none
** \ingroup QT
*/
void slaterdetvector_init(SlaterDetVector *sdvector, SlaterDetSet *sdset)
{
  sdvector->size = sdset->size;
  sdvector->sdset = sdset;
  sdvector->coeffs = init_array(sdvector->size);
}

/*!
** slaterdetvector_delete(): De-allocate a SlaterDetVector
**
** \param sdvector = pointer to SlaterDetVector to de-allocate
**
** Note: does NOT fully free the associated SlaterDetSet.  For that, see
** function slaterdetvector_delete_full()
**
** Returns: none
** \ingroup QT
*/
void slaterdetvector_delete(SlaterDetVector *sdvector)
{
  sdvector->size = 0;
  sdvector->sdset = NULL;
  if (sdvector->coeffs) {
    free(sdvector->coeffs);
    sdvector->coeffs = NULL;
  }
}


/*!
** slaterdetvector_delete_full(): De-allocate a SlaterDetVector and its
**   associated SlaterDetSet.
**
** To keep the SlaterDetSet itself, use similar function
** slaterdetvector_delete().
**
** \param sdvector = pointer to SlaterDetVector to delete
**
** Returns: none
** \ingroup QT
*/
void slaterdetvector_delete_full(SlaterDetVector *sdvector)
{
  sdvector->size = 0;
  if (sdvector->sdset) {
    slaterdetset_delete_full(sdvector->sdset);
    sdvector->sdset = NULL;
  }
  if (sdvector->coeffs) {
    free(sdvector->coeffs);
    sdvector->coeffs = NULL;
  }
}


/*!
** slaterdetvector_add(): Add a coefficient to a SlaterDetVector
**
** \param sdvector = Pointer to SlaterDetVector to add to
** \param index    = location in vector for writing the coefficient
** \param coeff    = the coefficient to write to location index
**
** Returns: none
** \ingroup QT
*/
void slaterdetvector_add(SlaterDetVector *sdvector, int index, double coeff)
{
  if (index < sdvector->size && index >= 0) {
    sdvector->coeffs[index] = coeff;
  }
}


/*!
** slaterdetvector_set(): Set a SlaterDetVector's vector to a set of
**   coefficients supplied by array coeffs
**
** \param sdvector = pointer to SlaterDetVector for writing coefficients
** \param coeffs   = array of coefficients to write to sdvector
**
** \ingroup QT
*/
void slaterdetvector_set(SlaterDetVector *sdvector, double *coeffs)
{
  int i;
  const int size = sdvector->size;
  double *v = sdvector->coeffs;
  if (v) {
    for(i=0; i<size; i++)
      v[i] = coeffs[i];
  }
}


/*!
** slaterdetvector_write(): Write a SlaterDetVector to disk.
**
** This writes a vector in the space of Slater determinants, along with
** the set of determinants itself, to a PSIO file.
**
** Use this if we only need to write a single vector.  Otherwise, call
** slaterdetset_write(); slaterdetset_write_vect();
** to allow for multiple vectors per slaterdetset to be written to disk.
**
** \param unit      = file number of the UNINITIALIZED PSIO file
** \param prefix    = prefix string to come before libpsio entry keys
** \param vector    = SlaterDetVector to write to disk
**
** Returns: none
** \ingroup QT
*/
void slaterdetvector_write(ULI unit, const char *prefix, SlaterDetVector *vector)
{
  int need_to_init_psio = 0;
  int unit_opened = 1;

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  slaterdetset_write(unit, prefix, vector->sdset);
  slaterdetset_write_vect(unit, prefix, vector->coeffs, vector->size, 0);

PSIO_CLOSE(unit)
PSIO_DONE

}


/*!
** slaterdetset_write_vect(): Write to disk the coefficients for a single
** vector associated with a set of Slater determinants.
**
** This function already assumes we've already called slaterdetset_write()
** to write out the string and determinant information.  This is only
** going to write out the coefficients.  This has been split out because
** we might want to write several roots for a given determinant setup.
** This does not actually dpend on the presence of a SlaterDetVector object
** so it is called a SlaterDetSet function.
**
** \param unit      = file number of the UNINITIALIZED PSIO file
** \param prefix    = prefix string to come before libpsio entry keys
** \param coeffs    = array of coefficients to write
** \param size      = number of elements in coeffs array
** \param vectnum   = the vector number (to make a PSIO key).  Start
**                    numbering from zero.
**
** Returns: none
**
** CDS 8/03
** \ingroup QT
*/
void slaterdetset_write_vect(ULI unit, const char *prefix,
  double *coeffs, int size, int vectnum)
{
  int need_to_init_psio = 0;
  int unit_opened = 1;
  char *vector_key;

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  if (vectnum < 0 || vectnum > 99) {
    outfile->Printf( "(slaterdetset_write_vect): vectnum out of bounds\n");
    abort();
  }

  vector_key = (char *) malloc(strlen(prefix)+strlen(SDVECTOR_KEY_VECTOR)+5);
  sprintf(vector_key,":%s:%s%2d",prefix,SDVECTOR_KEY_VECTOR,vectnum);

  psio_write_entry(unit, vector_key, (char *)coeffs, size*sizeof(double));

PSIO_CLOSE(unit)
PSIO_DONE

  free(vector_key);
}



/*!
** slaterdetvector_read(): Read a SlaterDetVector from disk
**
** Use this if we only need to read a single vector.  Otherwise, call
** slaterdetset_read(); slaterdetset_read_vect();
** to allow for multiple vectors per slaterdetset to be read from disk.
**
** \param unit      = file number to read from
** \param prefix    = prefix string to come before libpsio entry keys
** \param sdvector  = pointer to hold pointer to SlaterDetVector allocated
**                    by this function
**
** Returns: none
** \ingroup QT
*/
void slaterdetvector_read(ULI unit, const char *prefix, SlaterDetVector **sdvector)
{
  int need_to_init_psio = 0;
  int unit_opened = 1;
  SlaterDetSet *sdset;
  SlaterDetVector *vector = (SlaterDetVector *) malloc(sizeof(SlaterDetVector));

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  slaterdetset_read(unit, prefix, &sdset);
  slaterdetvector_init(vector,sdset);
  slaterdetset_read_vect(unit, prefix, vector->coeffs, vector->size, 0);

PSIO_CLOSE(unit)
PSIO_DONE

  *sdvector = vector;
}


/*!
** slaterdetset_read_vect(): Read in the coefficients for a single vector
** associated with a SlaterDetSet.
**
** This function already assumes we've already called slaterdetset_read()
** to read in the string and determinant information.  This is only
** going to read in the coefficients.  This has been split out because
** we might want to read several roots for a given determinant setup.
** This does not actually depend on the presence of a SlaterDetVector
** object and is called a SlaterDetSet function.
**
** \param unit      = file number of the UNINITIALIZED PSIO file
** \param prefix    = prefix string to come before libpsio entry keys
** \param coeffs    = array to hold coefficients read
** \param size      = number of elements in coeffs array
** \param vectnum   = the vector number (for the PSIO key).  Start from 0.
**
** Returns: none
**
** CDS 8/03
** \ingroup QT
*/
void slaterdetset_read_vect(ULI unit, const char *prefix, double *coeffs,
  int size, int vectnum)
{
  int need_to_init_psio = 0;
  int unit_opened = 1;
  char *vector_key;

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  vector_key = (char *) malloc(strlen(prefix)+strlen(SDVECTOR_KEY_VECTOR)+5);
  sprintf(vector_key,":%s:%s%2d",prefix,SDVECTOR_KEY_VECTOR,vectnum);

  psio_read_entry(unit, vector_key, (char *)coeffs, size*sizeof(double));

PSIO_CLOSE(unit)
PSIO_DONE

  free(vector_key);

}

}
