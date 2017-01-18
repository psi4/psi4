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

#ifndef _jk_independent_impl_h
#define _jk_independent_impl_h

using namespace psi;

template<class JDriver, class KDriver>
void JKIndependent<JDriver, KDriver>::preiterations()
{
  /*
  if (do_J_)
    j_driver_.initialize();

  if (do_K_ && do_separately_)
    k_driver_.initialize();
*/
}


template<class JDriver, class KDriver>
void JKIndependent<JDriver, KDriver>::compute_JK()
{

  //printf("compute_JK independent.\n");

  if (do_J_)
  {
    std::cout << "Doing independent J computation\n";

    timer_on("Independent_J");
    j_driver_.Update(D_ao_);
    j_driver_.Compute();
    timer_off("Independent_J");
    //J_ao_ = j_driver_.J();
    J_ = j_driver_.J();

    //outfile->Printf( "J_ao in jk_independent\n");
    //J_ao_[0]->print(outfile);
    if (!do_separately_ && do_K_) {
      std::cout << "Doing independent K computation with J driver.\n";
      //K_ao_ = j_driver_.K();
      //outfile->Printf( "K_ao in jk_independent\n");
      //K_ao_[0]->print(outfile);
      K_ = j_driver_.K();
    }
  }

  if (do_K_ && do_separately_)
  {
    std::cout << "Doing independent K computation\n";
    timer_on("Independent_K");
    k_driver_.Update(D_ao_);
    k_driver_.Compute();
    timer_off("Independent_K");
    //K_ao_ = k_driver_.K();
    K_ = k_driver_.K();
  }

}

template<class JDriver, class KDriver>
void JKIndependent<JDriver, KDriver>::postiterations()
{
  /*
  if (do_J_)
    j_driver_.finalize();

  if (do_K_ && do_separately_)
    k_driver_.finalize();
   */
}

template<class JDriver, class KDriver>
void JKIndependent<JDriver, KDriver>::common_init()
{
  /*
  if (do_J_)
    j_driver_.common_init();

  if (do_K_ && do_separately_)
    k_driver_.common_init();
   */
}

template<class JDriver, class KDriver>
JKIndependent<JDriver, KDriver>::JKIndependent(std::shared_ptr<BasisSet> primary,
                                               bool do_separately)
:
JK(primary),
j_driver_(primary, D_ao_),
k_driver_(primary, D_ao_),
do_separately_(do_separately)
{

  j_driver_.set_do_J(do_J_);
  j_driver_.set_do_K(do_separately_ ? false : do_K_);

  k_driver_.set_do_J(false);
  k_driver_.set_do_K(do_separately_ ? do_K_ : false);

}

template<class JDriver, class KDriver>
JKIndependent<JDriver, KDriver>::~JKIndependent()
{
  // should all be taken care of automatically
}

template<class JDriver, class KDriver>
void JKIndependent<JDriver, KDriver>::print_header() const
{

  if (print_) {
    outfile->Printf( "  ==> Independent J and K computations <==\n\n");
    outfile->Printf( "Coulomb computation: \n");

    j_driver_.print_header();

    outfile->Printf( "Exchange computation: \n");

    k_driver_.print_header();

  }
}

#endif
