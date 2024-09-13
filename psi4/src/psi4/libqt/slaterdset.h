/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*!
  \file
  \brief Header file for SlaterDetSets
  \ingroup QT
  Edward Valeev, June 2002
*/

#ifndef _psi3_src_lib_libqt_slaterdset_h_
#define _psi3_src_lib_libqt_slaterdset_h_

#include "psi4/libpsio/psio.h"

namespace psi {

/*!
  String is an orbital occupation string
*/
typedef struct {
    int index;
    short int *occ; /* Orbital indices in QT order */
                    /* CDS: I'm gonna use Pitzer order actually */
} String;

/*!
  StringSet is a set of strings
*/
typedef struct {
    String *strings;
    int size;
    int nelec;
    int ndrc;
    short int *drc_occ;
} StringSet;

void stringset_init(StringSet *stringset, int size, int nelec, int ndrc, short int *frozen_occ);
void stringset_delete(StringSet *stringset);
void stringset_add(StringSet *stringset, int index, unsigned char *Occ);
void stringset_write(size_t unit, const char *prefix, StringSet *sset);
void stringset_read(size_t unit, const char *prefix, StringSet **sset);
void stringset_reindex(StringSet *stringset, short int *mo_map);

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
typedef struct _SlaterDetSet {
    int size;
    SlaterDet *dets;
    StringSet *alphastrings, *betastrings;
} SlaterDetSet;

void slaterdetset_init(SlaterDetSet *sdset, int size, StringSet *alphastrings, StringSet *betastrings);
void slaterdetset_delete(SlaterDetSet *sdset);
void slaterdetset_delete_full(SlaterDetSet *sdset);
void slaterdetset_add(SlaterDetSet *sdset, int index, int alphastring, int betastring);
void slaterdetset_write(size_t unit, const char *prefix, SlaterDetSet *sdset);
void slaterdetset_read(size_t unit, const char *prefix, SlaterDetSet **sdset);

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

void slaterdetvector_write(size_t unit, const char *prefix, SlaterDetVector *vector);
void slaterdetset_write_vect(size_t unit, const char *prefix, double *coeffs, int size, int vectnum);
void slaterdetvector_read(size_t unit, const char *prefix, SlaterDetVector **vector);
void slaterdetset_read_vect(size_t unit, const char *prefix, double *coeffs, int size, int vectnum);

#define STRINGSET_KEY_SIZE "StringSet Size"
#define STRINGSET_KEY_NELEC "StringSet Num. of Electrons"
#define STRINGSET_KEY_NDRC "StringSet Num. of Dropped DOCCs"
#define STRINGSET_KEY_DRC_OCC "StringSet Dropped Core Occs"
#define STRINGSET_KEY_STRINGS "StringSet Strings"
#define SDSET_KEY_SIZE "SlaterDetSet Size"
#define SDSET_KEY_DETERMINANTS "SlaterDetSet Determinants"
#define SDSET_KEY_ALPHASTRINGS "Alpha Strings"
#define SDSET_KEY_BETASTRINGS "Beta Strings"
#define SDVECTOR_KEY_VECTOR "Vector"

}  // namespace psi

#endif /* _psi3_src_lib_libqt_slaterdset_h_ */
