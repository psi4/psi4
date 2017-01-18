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

#include <iostream>
#include <cstdlib>

#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libciomr/libciomr.h"

#include "index.h"
#include "matrix.h"

namespace psi{

    namespace psimrcc{

using namespace std;


//////////////////////////////////////////////////////////////////////////////
//
//  Two-Index Addressing Routines
//
//    Case A                   Case B                   Case C
//    The matrix looks like    The matrix looks like    The matrix looks like
//     -                        -------------            ------
//    |*|                      |*************|          |**|  |
//    |*|                       -------------           |**|  |
//    |*|                                               |  |**|
//    |*|                                               |  |**|
//    |*|                                                ------
//     -
//////////////////////////////////////////////////////////////////////////////

/**
 * Get the element \f$ M_{pq} \f$ where the indices are absolute
 * @param p
 * @param q
 * @return
 */
double CCMatrix::get_two_address_element(short p, short q)
{
  if(left->get_nelements() == 2){
    // Case A
    return(matrix[0][left->get_tuple_rel_index(p,q)][0]);
  }else if (left->get_nelements() == 0){
    // Case B
    return(matrix[0][0][right->get_tuple_rel_index(p,q)]);
  }else if (left->get_nelements() == 1){
    // Case C
    return(matrix[left->get_tuple_irrep(p)][left->get_tuple_rel_index(p)][right->get_tuple_rel_index(q)]);
  }
  outfile->Printf("\n\n\tdouble CCMatrix::get_two_address_element(int p, int q) Critical Error!!!");

  exit(1);
  return(0.0);
}

/**
 * Set the element \f$ M_{pq} \f$ where the indices are absolute
 * @param p
 * @param q
 * @return
 */
void CCMatrix::set_two_address_element(short p, short q,double value)
{
  if(left->get_nelements() == 2){
    // Case A
    matrix[0][left->get_tuple_rel_index(p,q)][0]=value;
  }else if (left->get_nelements() == 0){
    // Case B
    matrix[0][0][right->get_tuple_rel_index(p,q)]=value;
  }else if (left->get_nelements() == 1){
    // Case C
    matrix[left->get_tuple_irrep(p)][left->get_tuple_rel_index(p)][right->get_tuple_rel_index(q)]=value;
  }
}

/**
 * Add value to the element \f$ M_{pq} \f$ where the indices p and q are absolute
 * @param p
 * @param q
 * @param value
 * @return
 */
void CCMatrix::add_two_address_element(short p, short q,double value)
{
  if(left->get_nelements() == 2){
    // Case A
    matrix[0][left->get_tuple_rel_index(p,q)][0]=value;
  }else if (left->get_nelements() == 0){
    // Case B
    matrix[0][0][right->get_tuple_rel_index(p,q)]=value;
  }else if (left->get_nelements() == 1){
    // Case C
    matrix[left->get_tuple_irrep(p)][left->get_tuple_rel_index(p)][right->get_tuple_rel_index(q)]=value;
  }
}

void CCMatrix::get_two_indices(short*& pq,int irrep, int i, int j)
{
  // Returns the two indices of tuples i and j belonging to irrep
  if(left->get_nelements() == 2)
  {
    short* left_tuple = left->get_tuple(i+left->get_first(irrep));
    pq[0]=left_tuple[0];
    pq[1]=left_tuple[1];
  }else if (left->get_nelements() == 0){
    // Case B
    short* right_tuple = right->get_tuple(j+right->get_first(irrep));
    pq[0]=right_tuple[0];
    pq[1]=right_tuple[1];
  }else if (left->get_nelements() == 1){
    // Case C
    short* left_tuple = left->get_tuple(i+left->get_first(irrep));
    short* right_tuple = right->get_tuple(j+right->get_first(irrep));
    pq[0]=left_tuple[0];
    pq[1]=right_tuple[0];
  }
}

void CCMatrix::get_two_indices_pitzer(short*& pq,int irrep, int i, int j)
{
  // Returns the two indices of tuples i and j belonging to irrep
  if(left->get_nelements() == 2)
  {
    vecvecint& left_tuple_to_pitzer = left->get_indices_to_pitzer();
    short* left_tuple = left->get_tuple(i+left->get_first(irrep));
    pq[0]=left_tuple_to_pitzer[0][left_tuple[0]];
    pq[1]=left_tuple_to_pitzer[1][left_tuple[1]];
  }else if (left->get_nelements() == 0){
    // Case B
    vecvecint& right_tuple_to_pitzer = right->get_indices_to_pitzer();
    short* right_tuple = right->get_tuple(j+right->get_first(irrep));
    pq[0]=right_tuple_to_pitzer[0][right_tuple[0]];
    pq[1]=right_tuple_to_pitzer[1][right_tuple[1]];
  }else if (left->get_nelements() == 1){
    // Case C
    vecvecint& left_tuple_to_pitzer = left->get_indices_to_pitzer();
    vecvecint& right_tuple_to_pitzer = right->get_indices_to_pitzer();
    short* left_tuple = left->get_tuple(i+left->get_first(irrep));
    short* right_tuple = right->get_tuple(j+right->get_first(irrep));
    pq[0]=left_tuple_to_pitzer[0][left_tuple[0]];
    pq[1]=right_tuple_to_pitzer[0][right_tuple[0]];
  }
}


//////////////////////////////////////////////////////////////////////////////
//
//  Two-Index Addressing Routines
//
//    Case A                   Case B                   Case C
//    The matrix looks like    The matrix looks like    The matrix looks like
//     -                        -------------            ------
//    |*|                      |*************|          |**|  |
//    |*|                       -------------           |**|  |
//    |*|                                               |  |**|
//    |*|                                               |  |**|
//    |*|                                                ------
//     -
//////////////////////////////////////////////////////////////////////////////
double CCMatrix::get_four_address_element(short p, short q, short r, short s)
{
  if(left->get_nelements() == 1)
  {
    // Case A
    return(matrix[left->get_tuple_irrep(p)][left->get_tuple_rel_index(p)][right->get_tuple_rel_index(q,r,s)]);
  }else if (left->get_nelements() == 2){
    // Case B
    return(matrix[left->get_tuple_irrep(p,q)][left->get_tuple_rel_index(p,q)][right->get_tuple_rel_index(r,s)]);
  }else if (left->get_nelements() == 3){
    // Case C
    return(matrix[right->get_tuple_irrep(s)][left->get_tuple_rel_index(p,q,r)][right->get_tuple_rel_index(s)]);
  }
  outfile->Printf("\n\n\tdouble CCMatrix::get_four_address_element(int p, int q, int r, int s) Critical Error!!!");

  exit(1);
  return(0.0);
}

void CCMatrix::set_four_address_element(short p, short q, short r, short s,double value)
{
  if(left->get_nelements() == 1)
  {
    // Case A
    matrix[left->get_tuple_irrep(p)][left->get_tuple_rel_index(p)][right->get_tuple_rel_index(q,r,s)]=value;
  }else if (left->get_nelements() == 2){
    // Case B
    matrix[left->get_tuple_irrep(p,q)][left->get_tuple_rel_index(p,q)][right->get_tuple_rel_index(r,s)]=value;
  }else if (left->get_nelements() == 3){
    // Case C
    matrix[right->get_tuple_irrep(s)][left->get_tuple_rel_index(p,q,r)][right->get_tuple_rel_index(s)]=value;
  }
}

void CCMatrix::add_four_address_element(short p, short q, short r, short s,double value)
{
  if(left->get_nelements() == 1){
    // Case A
    matrix[left->get_tuple_irrep(p)][left->get_tuple_rel_index(p)][right->get_tuple_rel_index(q,r,s)]+=value;
  }else if (left->get_nelements() == 2){
    // Case B
    matrix[left->get_tuple_irrep(p,q)][left->get_tuple_rel_index(p,q)][right->get_tuple_rel_index(r,s)]+=value;
  }else if (left->get_nelements() == 3){
    // Case C
    matrix[right->get_tuple_irrep(s)][left->get_tuple_rel_index(p,q,r)][right->get_tuple_rel_index(s)]+=value;
  }
}

void CCMatrix::get_four_indices(short*& pqrs,int irrep, int i, int j)
{
  short*  left_tuple =  left->get_tuple(i+left->get_first(irrep));
  short* right_tuple = right->get_tuple(j+right->get_first(irrep));
  // Returns the two indices of pairs i and j belonging to irrep
  if(left->get_nelements() == 1)
  {
    pqrs[0]=left_tuple[0];
    pqrs[1]=right_tuple[0];
    pqrs[2]=right_tuple[1];
    pqrs[3]=right_tuple[2];
  }else if (left->get_nelements() == 2){
    pqrs[0]=left_tuple[0];
    pqrs[1]=left_tuple[1];
    pqrs[2]=right_tuple[0];
    pqrs[3]=right_tuple[1];
  }else if (left->get_nelements() == 3){
    pqrs[0]=left_tuple[0];
    pqrs[1]=left_tuple[1];
    pqrs[2]=left_tuple[2];
    pqrs[3]=right_tuple[0];
  }
}

void CCMatrix::get_four_indices_pitzer(short*& pqrs,int irrep, int i, int j)
{
  short*  left_tuple =  left->get_tuple(i+left->get_first(irrep));
  short* right_tuple = right->get_tuple(j+right->get_first(irrep));

  vecvecint& left_tuple_to_pitzer = left->get_indices_to_pitzer();
  vecvecint& right_tuple_to_pitzer = right->get_indices_to_pitzer();

  // Returns the two indices of pairs i and j belonging to irrep
  if(left->get_nelements() == 1)
  {
    pqrs[0]=left_tuple_to_pitzer[0][left_tuple[0]];
    pqrs[1]=right_tuple_to_pitzer[0][right_tuple[0]];
    pqrs[2]=right_tuple_to_pitzer[1][right_tuple[1]];
    pqrs[3]=right_tuple_to_pitzer[2][right_tuple[2]];
  }else if (left->get_nelements() == 2){
    pqrs[0]=left_tuple_to_pitzer[0][left_tuple[0]];
    pqrs[1]=left_tuple_to_pitzer[1][left_tuple[1]];
    pqrs[2]=right_tuple_to_pitzer[0][right_tuple[0]];
    pqrs[3]=right_tuple_to_pitzer[1][right_tuple[1]];
  }else if (left->get_nelements() == 3){
    pqrs[0]=left_tuple_to_pitzer[0][left_tuple[0]];
    pqrs[1]=left_tuple_to_pitzer[1][left_tuple[1]];
    pqrs[2]=left_tuple_to_pitzer[2][left_tuple[2]];
    pqrs[3]=right_tuple_to_pitzer[0][right_tuple[0]];
  }
}


//////////////////////////////////////////////////////////////////////////////
//
//  Six-Index Addressing Routines
//
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Retrieve \f$ Z_{ijk}^{abc} \f$
 */
double CCMatrix::get_six_address_element(short i, short j, short k, short a, short b, short c)
{
  // Assume three-three splitting
  return(matrix[left->get_tuple_irrep(i,j,k)][left->get_tuple_rel_index(i,j,k)][right->get_tuple_rel_index(a,b,c)]);
}

/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ \bar{H}_{ijk}^{[abc]} += Z_{ijk}^{[abc]}  \f]
 */
void CCMatrix::add_six_address_element(short i, short j, short k, short a, short b, short c,double value)
{
  // Assume three-three splitting
  matrix[left->get_tuple_irrep(i,j,k)][left->get_tuple_rel_index(i,j,k)][right->get_tuple_rel_index(a,b,c)]+=value;
}

/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ \bar{H}_{ijk}^{[abc]} += Z_{ijk}^{[abc]}  \f]
 */
void CCMatrix::add_six_address_element_abc(short i, short j, short k, size_t abc,double value)
{
  // Assume three-three splitting
  matrix[left->get_tuple_irrep(i,j,k)][left->get_tuple_rel_index(i,j,k)][abc]+=value;
}

/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ \bar{H}_{[ijk]}^{abc} += Z_{[ijk]}^{abc}  \f]
 */
void CCMatrix::add_six_address_element_ijk(size_t ijk, short a, short b, short c,double value)
{
  // Assume three-three splitting
  matrix[right->get_tuple_irrep(a,b,c)][ijk][right->get_tuple_rel_index(a,b,c)]+=value;
}

/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ Z_{ijk}^{abc} \rightarrow \{ \bar{H}_{ijk}^{abc}, -\bar{H}_{jik}^{abc} \f]
 */
void CCMatrix::add_six_address_element_Pij(short i, short j, short k, short a, short b, short c,double value)
{
  // Assume three-three splitting
  int  irrep = left->get_tuple_irrep(i,j,k);
  size_t ijk = left->get_tuple_rel_index(i,j,k);
  size_t jik = left->get_tuple_rel_index(j,i,k);

  size_t abc = right->get_tuple_rel_index(a,b,c);
  matrix[irrep][ijk][abc]+=value;
  matrix[irrep][jik][abc]-=value;
}

/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ Z_{ijk}^{abc} \rightarrow \{ \bar{H}_{ijk}^{abc}, -\bar{H}_{jik}^{abc} \f]
 */
void CCMatrix::add_six_address_element_Pij_abc(short i, short j, short k, size_t abc,double value)
{
  // Assume three-three splitting
  int  irrep = left->get_tuple_irrep(i,j,k);
  size_t ijk = left->get_tuple_rel_index(i,j,k);
  size_t jik = left->get_tuple_rel_index(j,i,k);

  matrix[irrep][ijk][abc]+=value;
  matrix[irrep][jik][abc]-=value;
}

/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ Z_{ijk}^{abc} \rightarrow \{ \bar{H}_{ijk}^{abc}, -\bar{H}_{ijk}^{cba}, -\bar{H}_{ijk}^{acb} \f]
 */
void CCMatrix::add_six_address_element_Pab_c(size_t ijk, short a, short b, short c, double value)
{
  // Assume three-three splitting
  int  irrep = right->get_tuple_irrep(a,b,c);
  size_t abc = right->get_tuple_rel_index(a,b,c);
  size_t cba = right->get_tuple_rel_index(c,b,a);
  size_t acb = right->get_tuple_rel_index(a,c,b);
  matrix[irrep][ijk][abc]+=value;
  matrix[irrep][ijk][cba]-=value;
  matrix[irrep][ijk][acb]-=value;
}


/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ Z_{ijk}^{abc} \rightarrow \{ \bar{H}_{ijk}^{abc}, -\bar{H}_{ikj}^{abc} \f]
 */
void CCMatrix::add_six_address_element_Pjk(short i, short j, short k, short a, short b, short c,double value)
{
  // Assume three-three splitting
  int  irrep = left->get_tuple_irrep(i,j,k);
  size_t ijk = left->get_tuple_rel_index(i,j,k);
  size_t ikj = left->get_tuple_rel_index(i,k,j);

  size_t abc = right->get_tuple_rel_index(a,b,c);
  matrix[irrep][ijk][abc]+=value;
  matrix[irrep][ikj][abc]-=value;
}


/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ Z_{ijk}^{abc} \rightarrow \{ \bar{H}_{ijk}^{abc}, -\bar{H}_{jik}^{abc} \f]
 */
void CCMatrix::add_six_address_element_Pjk_abc(short i, short j, short k, size_t abc,double value)
{
  // Assume three-three splitting
  int  irrep = left->get_tuple_irrep(i,j,k);
  size_t ijk = left->get_tuple_rel_index(i,j,k);
  size_t ikj = left->get_tuple_rel_index(i,k,j);

  matrix[irrep][ijk][abc]+=value;
  matrix[irrep][ikj][abc]-=value;
}

/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ Z_{ijk}^{abc} \rightarrow \{ \bar{H}_{ijk}^{abc}, -\bar{H}_{ijk}^{bac} \f]
 */
void CCMatrix::add_six_address_element_Pab(short i, short j, short k, short a, short b, short c,double value)
{
  // Assume three-three splitting
  int  irrep = left->get_tuple_irrep(i,j,k);
  size_t ijk = left->get_tuple_rel_index(i,j,k);

  size_t abc = right->get_tuple_rel_index(a,b,c);
  matrix[irrep][ijk][abc]+=value;

  size_t bac = right->get_tuple_rel_index(b,a,c);
  matrix[irrep][ijk][bac]-=value;
}

/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ Z_{ijk}^{abc} \rightarrow \{ \bar{H}_{ijk}^{abc}, -\bar{H}_{ijk}^{bac} \f]
 */
void CCMatrix::add_six_address_element_Pab_ijk(size_t ijk, short a, short b, short c,double value)
{
  // Assume three-three splitting
  int  irrep = right->get_tuple_irrep(a,b,c);

  size_t abc = right->get_tuple_rel_index(a,b,c);
  matrix[irrep][ijk][abc]+=value;

  size_t bac = right->get_tuple_rel_index(b,a,c);
  matrix[irrep][ijk][bac]-=value;
}

/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ Z_{ijk}^{abc} \rightarrow \{ \bar{H}_{ijk}^{abc}, -\bar{H}_{ijk}^{acb} \f]
 */
void CCMatrix::add_six_address_element_Pbc(short i, short j, short k, short a, short b, short c,double value)
{
  // Assume three-three splitting
  int  irrep = left->get_tuple_irrep(i,j,k);
  size_t ijk = left->get_tuple_rel_index(i,j,k);

  size_t abc = right->get_tuple_rel_index(a,b,c);
  matrix[irrep][ijk][abc]+=value;

  size_t acb = right->get_tuple_rel_index(a,c,b);
  matrix[irrep][ijk][acb]-=value;
}

/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ Z_{ijk}^{abc} \rightarrow \{ \bar{H}_{ijk}^{abc}, -\bar{H}_{ijk}^{acb} \f]
 */
void CCMatrix::add_six_address_element_Pbc_ijk(size_t ijk, short a, short b, short c,double value)
{
  // Assume three-three splitting
  int  irrep = right->get_tuple_irrep(a,b,c);

  size_t abc = right->get_tuple_rel_index(a,b,c);
  matrix[irrep][ijk][abc]+=value;

  size_t acb = right->get_tuple_rel_index(a,c,b);
  matrix[irrep][ijk][acb]-=value;
}

/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ Z_{ijk}^{abc} \rightarrow \{ \bar{H}_{ijk}^{abc}, -\bar{H}_{kji}^{abc}, -\bar{H}_{ikj}^{abc} \f]
 */
void CCMatrix::add_six_address_element_Pij_k(short i, short j, short k, size_t abc, double value)
{
  // Assume three-three splitting
  int  irrep = left->get_tuple_irrep(i,j,k);
  size_t ijk = left->get_tuple_rel_index(i,j,k);
  size_t kji = left->get_tuple_rel_index(k,j,i);
  size_t ikj = left->get_tuple_rel_index(i,k,j);

  matrix[irrep][ijk][abc]+=value;
  matrix[irrep][kji][abc]-=value;
  matrix[irrep][ikj][abc]-=value;
}

/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ Z_{ijk}^{abc} \rightarrow \{ \bar{H}_{ijk}^{abc}, \bar{H}_{kij}^{abc}, \bar{H}_{jki}^{abc},-\bar{H}_{jik}^{abc},-\bar{H}_{kji}^{abc},-\bar{H}_{ikj}^{abc} \f]
 */
void CCMatrix::add_six_address_element_Pijk(short i, short j, short k, short a, short b, short c,double value)
{
  // Assume three-three splitting
  int  irrep = left->get_tuple_irrep(i,j,k);
  size_t abc = right->get_tuple_rel_index(a,b,c);

  size_t ijk = left->get_tuple_rel_index(i,j,k);
  matrix[irrep][ijk][abc]+=value;

  size_t kij = left->get_tuple_rel_index(k,i,j);
  matrix[irrep][kij][abc]+=value;

  size_t jki = left->get_tuple_rel_index(j,k,i);
  matrix[irrep][jki][abc]+=value;

  size_t jik = left->get_tuple_rel_index(j,i,k);
  matrix[irrep][jik][abc]-=value;

  size_t kji = left->get_tuple_rel_index(k,j,i);
  matrix[irrep][kji][abc]-=value;

  size_t ikj = left->get_tuple_rel_index(i,k,j);
  matrix[irrep][ikj][abc]-=value;
}

/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ Z_{ijk}^{abc} \rightarrow \{ \bar{H}_{ijk}^{abc}, -\bar{H}_{ijk}^{bac},-\bar{H}_{jik}^{abc},  \bar{H}_{jik}^{bac} \f]
 */
void CCMatrix::add_six_address_element_Pij_Pab(short i, short j, short k, short a, short b, short c,double value)
{
  // Assume three-three splitting
  int  irrep = left->get_tuple_irrep(i,j,k);
  size_t ijk = left->get_tuple_rel_index(i,j,k);
  size_t jik = left->get_tuple_rel_index(j,i,k);

  size_t abc = right->get_tuple_rel_index(a,b,c);
  matrix[irrep][ijk][abc]+=value;
  matrix[irrep][jik][abc]-=value;

  size_t bac = right->get_tuple_rel_index(b,a,c);
  matrix[irrep][ijk][bac]-=value;
  matrix[irrep][jik][bac]+=value;
}

/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ Z_{ijk}^{abc} \rightarrow \{ \bar{H}_{ijk}^{abc}, -\bar{H}_{ijk}^{acb},-\bar{H}_{ikj}^{abc},  \bar{H}_{ikj}^{acb} \f]
 */
void CCMatrix::add_six_address_element_Pjk_Pbc(short i, short j, short k, short a, short b, short c,double value)
{
  // Assume three-three splitting
  int  irrep = left->get_tuple_irrep(i,j,k);
  size_t ijk = left->get_tuple_rel_index(i,j,k);
  size_t ikj = left->get_tuple_rel_index(i,k,j);

  size_t abc = right->get_tuple_rel_index(a,b,c);
  matrix[irrep][ijk][abc]+=value;
  matrix[irrep][ikj][abc]-=value;

  size_t acb = right->get_tuple_rel_index(a,c,b);
  matrix[irrep][ijk][acb]-=value;
  matrix[irrep][ikj][acb]+=value;
}

/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ Z_{ijk}^{abc} \rightarrow
     \left\{
       \begin{array}{ccc}
           \bar{H}_{ijk}^{abc}, & - \bar{H}_{ijk}^{bac}, & - \bar{H}_{ijk}^{cba} \\
         - \bar{H}_{kji}^{abc}, &   \bar{H}_{kji}^{bac}, &   \bar{H}_{kji}^{cba} \\
         - \bar{H}_{ikj}^{abc}, &   \bar{H}_{ikj}^{bac}, &   \bar{H}_{ikj}^{cba}
       \end{array}
     \right.
   \f]
 */
void CCMatrix::add_six_address_element_Pij_k_Pa_bc(short i, short j, short k, short a, short b, short c,double value)
{
  // Assume three-three splitting
  int  irrep = left->get_tuple_irrep(i,j,k);
  size_t ijk = left->get_tuple_rel_index(i,j,k);
  size_t kji = left->get_tuple_rel_index(k,j,i);
  size_t ikj = left->get_tuple_rel_index(i,k,j);

  size_t abc = right->get_tuple_rel_index(a,b,c);
  matrix[irrep][ijk][abc]+=value;
  matrix[irrep][kji][abc]-=value;
  matrix[irrep][ikj][abc]-=value;

  size_t bac = right->get_tuple_rel_index(b,a,c);
  matrix[irrep][ijk][bac]-=value;
  matrix[irrep][kji][bac]+=value;
  matrix[irrep][ikj][bac]+=value;

  size_t cba = right->get_tuple_rel_index(c,b,a);
  matrix[irrep][ijk][cba]-=value;
  matrix[irrep][kji][cba]+=value;
  matrix[irrep][ikj][cba]+=value;
}

/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ Z_{ijk}^{abc} \rightarrow
     \left\{
       \begin{array}{ccc}
           \bar{H}_{ijk}^{abc}, & - \bar{H}_{ijk}^{cba}, & - \bar{H}_{ijk}^{acb} \\
         - \bar{H}_{jik}^{abc}, &   \bar{H}_{jik}^{cba}, &   \bar{H}_{jik}^{acb} \\
         - \bar{H}_{kji}^{abc}, &   \bar{H}_{kji}^{cba}, &   \bar{H}_{kji}^{acb}
       \end{array}
     \right.
   \f]
 */
void CCMatrix::add_six_address_element_Pi_jk_Pab_c(short i, short j, short k, short a, short b, short c,double value)
{
  // Assume three-three splitting
  int  irrep = left->get_tuple_irrep(i,j,k);
  size_t ijk = left->get_tuple_rel_index(i,j,k);
  size_t jik = left->get_tuple_rel_index(j,i,k);
  size_t kji = left->get_tuple_rel_index(k,j,i);

  size_t abc = right->get_tuple_rel_index(a,b,c);
  matrix[irrep][ijk][abc]+=value;
  matrix[irrep][jik][abc]-=value;
  matrix[irrep][kji][abc]-=value;

  size_t cba = right->get_tuple_rel_index(c,b,a);
  matrix[irrep][ijk][cba]-=value;
  matrix[irrep][jik][cba]+=value;
  matrix[irrep][kji][cba]+=value;

  size_t acb = right->get_tuple_rel_index(a,c,b);
  matrix[irrep][ijk][acb]-=value;
  matrix[irrep][jik][acb]+=value;
  matrix[irrep][kji][acb]+=value;
}


/**
 * @brief Given a value \f$ Z_{ijk}^{abc} \f$ do
 * \f[ Z_{ijk}^{abc} \rightarrow
     \left\{
       \begin{array}{ccc}
           \bar{H}_{ijk}^{abc}, & - \bar{H}_{ijk}^{bac}, & - \bar{H}_{ijk}^{cba} \\
         - \bar{H}_{jik}^{abc}, &   \bar{H}_{jik}^{bac}, &   \bar{H}_{jik}^{cba} \\
         - \bar{H}_{kji}^{abc}, &   \bar{H}_{kji}^{bac}, &   \bar{H}_{kji}^{cba}
       \end{array}
     \right.
   \f]
 */
void CCMatrix::add_six_address_element_Pi_jk_Pa_bc(short i, short j, short k, short a, short b, short c,double value)
{
  // Assume three-three splitting
  int  irrep = left->get_tuple_irrep(i,j,k);
  size_t ijk = left->get_tuple_rel_index(i,j,k);
  size_t jik = left->get_tuple_rel_index(j,i,k);
  size_t kji = left->get_tuple_rel_index(k,j,i);

  size_t abc = right->get_tuple_rel_index(a,b,c);
  matrix[irrep][ijk][abc]+=value;
  matrix[irrep][jik][abc]-=value;
  matrix[irrep][kji][abc]-=value;

  size_t bac = right->get_tuple_rel_index(b,a,c);
  matrix[irrep][ijk][bac]-=value;
  matrix[irrep][jik][bac]+=value;
  matrix[irrep][kji][bac]+=value;

  size_t cba = right->get_tuple_rel_index(c,b,a);
  matrix[irrep][ijk][cba]-=value;
  matrix[irrep][jik][cba]+=value;
  matrix[irrep][kji][cba]+=value;
}

/*
double CCMatrix::get_two_address_element_check(short p, short q)
{
  if(left->get_nelements() == 2){
    // Case A
    outfile->Printf("double CCMatrix::get_two_address_element_check(int p, int q) Case A non implemented");

    exit(1);
  }else if (left->get_nelements() == 0){
    // Case B
    outfile->Printf("double CCMatrix::get_two_address_element_check(int p, int q) Case B non implemented");

    exit(1);
  }else if (left->get_nelements() == 1){
    // Case C
    if((left->get_tuple_rel_index(p)!=-1) && (right->get_tuple_rel_index(q)!=-1))
      return(matrix[left->get_tuple_irrep(p)][left->get_tuple_rel_index(p)][right->get_tuple_rel_index(q)]);
    else return(0.0);
  }
  outfile->Printf("\n\n\tdouble CCMatrix::get_two_address_element_check(int p, int q) Critical Error!!!");

  exit(1);
  return(0.0);
}
*/

}} /* End Namespaces */
