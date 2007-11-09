
/*!
  \file slaterdset.c
  Edward Valeev, June 2002
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include "slaterdset.h"

extern "C" {
	
#define PSIO_INIT if (!psio_state()) { \
    psio_init(); psio_ipv1_config(); \
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


/*! stringset_init()
 */
void stringset_init(StringSet *sset, int size, int nelec, int nfzc,
  short int *frozen_occ)
{
  int i;

  sset->size = size;
  sset->nelec = nelec;
  sset->nfzc = nfzc;
  sset->strings = (String *) malloc(size*sizeof(String));
  memset(sset->strings,0,size*sizeof(String));
  if (nfzc > 0) {
    sset->fzc_occ = (short int *) malloc(nfzc * sizeof(short int));
    for (i=0; i<nfzc; i++) {
      sset->fzc_occ[i] = frozen_occ[i];  
    }
  }
}


/*! stringset_delete()
 */
void stringset_delete(StringSet *sset)
{
  if (sset->nfzc > 0) free(sset->fzc_occ);
  sset->size = 0;
  sset->nelec = 0;
  sset->nfzc = 0;
  if (sset->strings) free(sset->strings);
  sset->strings = NULL;
}

/*! stringset_add
 */
void stringset_add(StringSet *sset, int index, unsigned char *Occ)
{
  int i;
  int nact = sset->nelec - sset->nfzc;
  String *s;

  if (index < sset->size && index >= 0) {
    s = sset->strings + index;
  }
  s->index = index;
  s->occ = (short int*) malloc(nact*sizeof(short int));
  for(i=0;i<nact;i++)
    s->occ[i] = Occ[i];
}

/*! stringset_reindex
 */
void stringset_reindex(StringSet* sset, short int* mo_map)
{
  int s, mo, core;
  short int* occ;
  int nstrings = sset->size;
  int nact = sset->nelec - sset->nfzc;

  for (core=0; core<sset->nfzc; core++) {
    sset->fzc_occ[core] = mo_map[sset->fzc_occ[core]];
  }

  for(s=0; s<nstrings; s++) {
    occ = (sset->strings + s)->occ;
    for(mo=0; mo<nact; mo++)
      occ[mo] = mo_map[occ[mo]];
  }
}

void stringset_write(ULI unit, char *prefix, StringSet *sset)
{
  int i, size, nact;
  int need_to_init_psio = 0;
  int unit_opened = 1;
  char *size_key, *nelec_key, *nfzc_key, *strings_key, *fzc_occ_key;
  psio_address ptr;

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  size_key = (char *) malloc(strlen(prefix) + strlen(STRINGSET_KEY_SIZE) + 3);
  sprintf(size_key,":%s:%s",prefix,STRINGSET_KEY_SIZE);
  nelec_key = (char *) malloc(strlen(prefix) + strlen(STRINGSET_KEY_NELEC) + 3);
  sprintf(nelec_key,":%s:%s",prefix,STRINGSET_KEY_NELEC);
  nfzc_key = (char *) malloc(strlen(prefix) + strlen(STRINGSET_KEY_NFZC) + 3);
  sprintf(nfzc_key,":%s:%s",prefix,STRINGSET_KEY_NFZC);
  fzc_occ_key = (char *) malloc(strlen(prefix) +
    strlen(STRINGSET_KEY_FZC_OCC) + 3);
  sprintf(fzc_occ_key,":%s:%s",prefix,STRINGSET_KEY_FZC_OCC);
  strings_key = (char *) malloc(strlen(prefix) + strlen(STRINGSET_KEY_STRINGS) + 3);
  sprintf(strings_key,":%s:%s",prefix,STRINGSET_KEY_STRINGS);

  psio_write_entry( unit, size_key, (char *)&sset->size, sizeof(int));
  psio_write_entry( unit, nelec_key, (char *)&sset->nelec, sizeof(int));
  psio_write_entry( unit, nfzc_key, (char *)&sset->nfzc, sizeof(int));
  if (sset->nfzc) {
    psio_write_entry( unit, fzc_occ_key, (char *)sset->fzc_occ, 
      sset->nfzc*sizeof(short int));
  }

  ptr = PSIO_ZERO;
  size = sset->size;
  nact = sset->nelec - sset->nfzc;
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
  free(nfzc_key);
  free(strings_key);
  free(fzc_occ_key);
}


void stringset_read(ULI unit, char *prefix, StringSet **stringset)
{
  int i, size, nelec, nfzc, nact;
  int need_to_init_psio = 0;
  int unit_opened = 1;
  char *size_key, *nelec_key, *nfzc_key, *fzc_occ_key, *strings_key;
  short int *fzc_occ;
  psio_address ptr;
  StringSet *sset = (StringSet *) malloc(sizeof(StringSet));

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  size_key = (char *) malloc( strlen(prefix) + strlen(STRINGSET_KEY_SIZE) + 3);
  sprintf(size_key,":%s:%s",prefix,STRINGSET_KEY_SIZE);
  nelec_key = (char *) malloc( strlen(prefix) + strlen(STRINGSET_KEY_NELEC) + 3);
  sprintf(nelec_key,":%s:%s",prefix,STRINGSET_KEY_NELEC);
  nfzc_key = (char *) malloc( strlen(prefix) + strlen(STRINGSET_KEY_NFZC) + 3);
  sprintf(nfzc_key,":%s:%s",prefix,STRINGSET_KEY_NFZC);
  fzc_occ_key = (char *) malloc(strlen(prefix) +
    strlen(STRINGSET_KEY_FZC_OCC) + 3);
  sprintf(fzc_occ_key,":%s:%s",prefix,STRINGSET_KEY_FZC_OCC);
  strings_key = (char *) malloc( strlen(prefix) + strlen(STRINGSET_KEY_STRINGS) + 3);
  sprintf(strings_key,":%s:%s",prefix,STRINGSET_KEY_STRINGS);

  psio_read_entry( unit, size_key, (char *)&size, sizeof(int));
  psio_read_entry( unit, nelec_key, (char *)&nelec, sizeof(int));
  psio_read_entry( unit, nfzc_key, (char *)&nfzc, sizeof(int));
  if (nfzc > 0) {
    fzc_occ = (short int *) malloc(nfzc*sizeof(short int));
    psio_read_entry( unit, fzc_occ_key, (char *)fzc_occ, 
      nfzc*sizeof(short int));
  }
  else fzc_occ = NULL;

  stringset_init(sset, size, nelec, nfzc, fzc_occ);

  nact = nelec - nfzc;
  ptr = PSIO_ZERO;
  for(i=0; i<size; i++) {
    psio_read( unit, strings_key, (char *) &(sset->strings[i].index), sizeof(int), ptr, &ptr);
    sset->strings[i].occ = (short int*) malloc(nact*sizeof(short int));
    psio_read( unit, strings_key, (char *) sset->strings[i].occ, nact*sizeof(short int), ptr, &ptr);
  }

PSIO_CLOSE(unit)
PSIO_DONE

  free(size_key);
  free(nelec_key);
  free(nfzc_key);
  free(fzc_occ_key);
  free(strings_key);
  if (nfzc > 0) free(fzc_occ);
  *stringset = sset;
}


/*! slaterdetset_init()
 */
void slaterdetset_init(SlaterDetSet *sdset, int size, StringSet *alphastrings, StringSet *betastrings)
{
  sdset->size = size;
  sdset->dets = (SlaterDet *) malloc(size*sizeof(SlaterDet));
  memset(sdset->dets,0,size*sizeof(SlaterDet));
  sdset->alphastrings = alphastrings;
  sdset->betastrings = betastrings;
}

/*! slaterdetset_delete()
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

/*! slaterdetset_delete_full()
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

/*! slaterdetset_add
 */
void slaterdetset_add(SlaterDetSet *sdset, int index, int alphastring, int betastring)
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

void slaterdetset_write(ULI unit, char *prefix, SlaterDetSet *sdset)
{
  int i;
  int need_to_init_psio = 0;
  int unit_opened = 1;
  char *size_key, *set_key;
  char *alphaprefix, *betaprefix;
  psio_address ptr;

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  alphaprefix = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_ALPHASTRINGS) + 2);
  sprintf(alphaprefix,"%s:%s",prefix,SDSET_KEY_ALPHASTRINGS);
  betaprefix = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_BETASTRINGS) + 2);
  sprintf(betaprefix,"%s:%s",prefix,SDSET_KEY_BETASTRINGS);

  stringset_write( unit, alphaprefix, sdset->alphastrings);
  stringset_write( unit, betaprefix, sdset->betastrings);
  
  free(alphaprefix);
  free(betaprefix);

  size_key = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_SIZE) + 3);
  sprintf(size_key,":%s:%s",prefix,SDSET_KEY_SIZE);
  set_key = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_DETERMINANTS) + 3);
  sprintf(set_key,":%s:%s",prefix,SDSET_KEY_DETERMINANTS);

  psio_write_entry( unit, size_key, (char *)&sdset->size, sizeof(int));
  psio_write_entry( unit, set_key, (char *)sdset->dets, sdset->size*sizeof(SlaterDet));

PSIO_CLOSE(unit)
PSIO_DONE

  free(size_key);
  free(set_key);
}

void slaterdetset_read(ULI unit, char *prefix, SlaterDetSet **slaterdetset)
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

  alphaprefix = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_ALPHASTRINGS) + 2);
  sprintf(alphaprefix,"%s:%s",prefix,SDSET_KEY_ALPHASTRINGS);
  betaprefix = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_BETASTRINGS) + 2);
  sprintf(betaprefix,"%s:%s",prefix,SDSET_KEY_BETASTRINGS);

  stringset_read( unit, alphaprefix, &alphastrings);
  stringset_read( unit, betaprefix, &betastrings);
  
  free(alphaprefix);
  free(betaprefix);

  size_key = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_SIZE) + 3);
  sprintf(size_key,":%s:%s",prefix,SDSET_KEY_SIZE);
  set_key = (char *) malloc( strlen(prefix) + strlen(SDSET_KEY_DETERMINANTS) + 3);
  sprintf(set_key,":%s:%s",prefix,SDSET_KEY_DETERMINANTS);

  psio_read_entry( unit, size_key, (char *)&size, sizeof(int));
  slaterdetset_init(sdset,size,alphastrings,betastrings);
  psio_read_entry( unit, set_key, (char *)sdset->dets, sdset->size*sizeof(SlaterDet));

PSIO_CLOSE(unit)
PSIO_DONE

  free(size_key);
  free(set_key);

  *slaterdetset = sdset;
}


/*! slaterdetvector_init()
 */
void slaterdetvector_init(SlaterDetVector *sdvector, SlaterDetSet *sdset)
{
  sdvector->size = sdset->size;
  sdvector->sdset = sdset;
  sdvector->coeffs = init_array(sdvector->size);
}

/*! slaterdetvector_delete()
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


/*! slaterdetvector_delete_full()
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


/*! slaterdetvector_add
 */
void slaterdetvector_add(SlaterDetVector *sdvector, int index, double coeff)
{
  if (index < sdvector->size && index >= 0) {
    sdvector->coeffs[index] = coeff;
  }
}


/*! slaterdetvector_set
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


/*
** Use this if we only need to write a single vector.  Otherwise, call
** slaterdetset_write(); slaterdetset_write_vect();
** to allow for multiple vectors per slaterdetset to be written to disk.
*/
void slaterdetvector_write(ULI unit, char *prefix, SlaterDetVector *vector)
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


/*
** This function already assumes we've already called slaterdetset_write()
** to write out the string and determinant information.  This is only 
** going to write out the coefficients.  This has been split out because
** we might want to write several roots for a given determinant setup.
**
** CDS 8/03
*/
void slaterdetset_write_vect(ULI unit, char *prefix, 
  double *coeffs, int size, int vectnum)
{
  int need_to_init_psio = 0;
  int unit_opened = 1;
  char *vector_key;

PSIO_INIT
PSIO_OPEN(unit,PSIO_OPEN_OLD)

  if (vectnum < 0 || vectnum > 99) {
    fprintf(stderr, "(slaterdetset_write_vect): vectnum out of bounds\n");
    abort();
  }

  vector_key = (char *) malloc(strlen(prefix)+strlen(SDVECTOR_KEY_VECTOR)+5);
  sprintf(vector_key,":%s:%s%2d",prefix,SDVECTOR_KEY_VECTOR,vectnum);

  psio_write_entry(unit, vector_key, (char *)coeffs, size*sizeof(double));

PSIO_CLOSE(unit)
PSIO_DONE

  free(vector_key);
}



/*
** Use this if we only need to read a single vector.  Otherwise, call
** slaterdetset_read(); slaterdetset_read_vect();
** to allow for multiple vectors per slaterdetset to be read from disk.
*/
void slaterdetvector_read(ULI unit, char *prefix, SlaterDetVector **sdvector)
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


/*
** This function already assumes we've already called slaterdetset_read()
** to read in the string and determinant information.  This is only 
** going to read in the coefficients.  This has been split out because
** we might want to read several roots for a given determinant setup.
**
** CDS 8/03
*/
void slaterdetset_read_vect(ULI unit, char *prefix, double *coeffs, 
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

} /* extern "C" */