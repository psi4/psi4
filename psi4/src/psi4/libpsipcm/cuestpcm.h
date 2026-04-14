/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

 #ifndef CUESTPCM_H
 #define CUESTPCM_H
 
 #ifdef USING_cuEST
 #include "cuest.h"
 #include <memory>
 #include "psi4/libmints/typedefs.h"
 
 namespace psi {
     class MintsHelper;
     class Options;
 
     class cuestPCM final{
         /// cuEST PCM integral plan
         cuestPCMIntPlan_t pcm_integral_plan_ = nullptr;
         cuestWorkspace_t* pcm_integral_ws_ptr_ = nullptr;
         double *d_q1_matrix_ = nullptr;
         double *d_q2_matrix_ = nullptr;
         double *d_D_matrix_ = nullptr;
         double *d_V_matrix_ = nullptr;
         std::shared_ptr<MintsHelper> mintshelper_;
         double pcm_convergence_;
 
     public:
         cuestPCM() = delete;
         cuestPCM(const Options &options, const std::shared_ptr<MintsHelper> &mintshelper);
         ~cuestPCM();
         cuestPCM(const cuestPCM &) = delete;
         cuestPCM(cuestPCM &&) = delete;
         cuestPCM& operator=(const cuestPCM &) = delete;
         cuestPCM& operator=(cuestPCM &&) = delete;
 
         std::pair<double, SharedMatrix> compute_PCM_terms(const SharedMatrix &D);
         SharedMatrix compute_PCM_gradient(const SharedMatrix &D);
         SharedMatrix compute_V(const SharedMatrix &D);
     };
 } // namespace psi
 #endif // USING_cuEST
 
 #endif // CUESTPCM_H
 