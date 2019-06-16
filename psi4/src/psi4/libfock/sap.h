/*
  Copyright (c) 2019, Susi Lehtola
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
  * Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
  * Neither the name of the <organization> nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
  USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
  SUCH DAMAGE.
*/

#ifndef SAP_POTENTIAL
#define SAP_POTENTIAL
/*
  Routines for the implementation of the superposition of atomic
  potentials guess for electronic structure calculations, see

  S. Lehtola, "Assessment of Initial Guesses for Self-Consistent Field
  Calculations. Superposition of Atomic Potentials: Simple yet
  Efficient", J. Chem. Theory Comput. 15, 1593 (2019).
  DOI: 10.1021/acs.jctc.8b01089

  This function evaluates the repulsive part of the LDA exchange-only
  potential of a neutral atom. The potentials have been calculated for
  the ground-states of spherically symmetric atoms at the non-relativistic
  level of theory, using accurate finite-element calculations as described
  in

  S. Lehtola, "Fully numerical Hartree-Fock and density functional
  calculations. I. Atoms", Int J Quantum Chem. e25945 (2019).
  DOI: 10.1002/qua.25945
*/
double sap_effective_charge(int Z, double r);
#endif
