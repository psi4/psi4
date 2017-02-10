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

#include "sapt.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/potential.h"
#include "psi4/libmints/basisset.h"

namespace psi { namespace sapt {

SAPT::SAPT(SharedWavefunction Dimer, SharedWavefunction MonomerA,
           SharedWavefunction MonomerB, Options& options,
           std::shared_ptr<PSIO> psio)
            : Wavefunction(options)
{
  shallow_copy(Dimer);

  if ((Dimer->nirrep() != 1) || (MonomerA->nirrep() != 1) || (MonomerA->nirrep() != 1)){
    throw PSIEXCEPTION("SAPT must be run in C1 symmetry. Period.");
  }

  if ((Dimer->soccpi().sum() != 0) || (MonomerA->soccpi().sum() != 0) ||
       (MonomerA->soccpi().sum() != 0)){
    throw PSIEXCEPTION("This is a RHF SAPT constructor. Pair those electrons up cracker!");
  }


  psio_ = psio;

#ifdef USING_LAPACK_MKL
  mkl_set_dynamic(1);
#endif

#ifdef _OPENMP
  omp_set_nested(0);
#endif

  initialize(MonomerA, MonomerB);
  get_denom();
}

SAPT::~SAPT()
{
  if (evalsA_ != NULL) free(evalsA_);
  if (evalsB_ != NULL) free(evalsB_);
  if (diagAA_ != NULL) free(diagAA_);
  if (diagBB_ != NULL) free(diagBB_);
  if (CA_ != NULL) free_block(CA_);
  if (CB_ != NULL) free_block(CB_);
  if (CHFA_ != NULL) free_block(CHFA_);
  if (CHFB_ != NULL) free_block(CHFB_);
  if (sAB_ != NULL) free_block(sAB_);
  if (vABB_ != NULL) free_block(vABB_);
  if (vBAA_ != NULL) free_block(vBAA_);
  if (vAAB_ != NULL) free_block(vAAB_);
  if (vBAB_ != NULL) free_block(vBAB_);
  zero_.reset();
}

void SAPT::initialize(SharedWavefunction MonomerA, SharedWavefunction MonomerB)
{
  evalsA_ = NULL;
  evalsB_ = NULL;
  diagAA_ = NULL;
  diagBB_ = NULL;
  CA_ = NULL;
  CB_ = NULL;
  CHFA_ = NULL;
  CHFB_ = NULL;
  sAB_ = NULL;
  vABB_ = NULL;
  vBAA_ = NULL;
  vAAB_ = NULL;
  vBAB_ = NULL;

  // We inherit from the dimer basis
  ribasis_ = get_basisset("DF_BASIS_SAPT");
  elstbasis_ = get_basisset("DF_BASIS_ELST");

  // Compare pointers
  if (ribasis_ == elstbasis_){
      elst_basis_ = 0;
  } else {
      elst_basis_ = 1;
  }

  zero_ = std::shared_ptr<BasisSet>(BasisSet::zero_ao_basis_set());

  if(options_.get_str("EXCH_SCALE_ALPHA") == "FALSE") {
      exch_scale_alpha_ = 0.0;
  } else if (options_.get_str("EXCH_SCALE_ALPHA") == "TRUE") {
      exch_scale_alpha_ = 1.0;                // Default value for alpha
  } else {
      exch_scale_alpha_ = std::atof(options_.get_str("EXCH_SCALE_ALPHA").c_str());
  }
  Process::environment.globals["SAPT ALPHA"] = exch_scale_alpha_;
  print_ = options_.get_int("PRINT");
  debug_ = options_.get_int("DEBUG");
  schwarz_ = options_.get_double("INTS_TOLERANCE");
  mem_ = (long int) ((double) memory_*options_.get_double("SAPT_MEM_SAFETY"));
  mem_ /= 8L;

  std::vector<int> realsA;
  realsA.push_back(0);
  std::vector<int> ghostsA;
  ghostsA.push_back(1);
  std::shared_ptr<Molecule> monomerA = molecule_->extract_subsets(realsA,
    ghostsA);
  foccA_ = monomerA->nfrozen_core(options_.get_str("FREEZE_CORE"));

  std::vector<int> realsB;
  realsB.push_back(1);
  std::vector<int> ghostsB;
  ghostsB.push_back(0);
  std::shared_ptr<Molecule> monomerB = molecule_->extract_subsets(realsB,
    ghostsB);
  foccB_ = monomerB->nfrozen_core(options_.get_str("FREEZE_CORE"));

  natomsA_ = 0;
  natomsB_ = 0;

  for (int n=0; n<monomerA->natom(); n++)
    if (monomerA->Z(n)) natomsA_++;
  for (int n=0; n<monomerB->natom(); n++)
    if (monomerB->Z(n)) natomsB_++;

  ndf_ = ribasis_->nbf();

  double enucD, enucA, enucB;
  double eHFD, eHFA, eHFB;

  eHFD = energy_;
  enucD = molecule_->nuclear_repulsion_energy();

  // Monomer A info
  nsoA_  = MonomerA->nso();
  nmoA_  = MonomerA->nmo();
  noccA_ = MonomerA->doccpi().sum();
  nvirA_ = nmoA_ - noccA_;
  NA_    = 2 * noccA_;
  eHFA   = MonomerA->reference_energy();
  enucA  = MonomerA->molecule()->nuclear_repulsion_energy();
  aoccA_ = noccA_ - foccA_;

  // Monomer B info
  nsoB_  = MonomerB->nso();
  nmoB_  = MonomerB->nmo();
  noccB_ = MonomerB->doccpi().sum();
  nvirB_ = nmoB_ - noccB_;
  NB_    = 2 * noccB_;
  eHFB   = MonomerB->reference_energy();
  enucB  = MonomerB->molecule()->nuclear_repulsion_energy();
  aoccB_ = noccB_ - foccB_;

  enuc_ = enucD - enucA - enucB;
  eHF_ =  eHFD - eHFA - eHFB;

  evalsA_ = init_array(nmoA_);
  std::memcpy(evalsA_, MonomerA->epsilon_a()->pointer(), sizeof(double) * nmoA_);

  evalsB_ = init_array(nmoB_);
  std::memcpy(evalsB_, MonomerB->epsilon_a()->pointer(), sizeof(double) * nmoB_);

  CA_ = block_matrix(nso_,nmoA_);
  double **tempA = block_matrix(nsoA_,nmoA_);
  std::memcpy(tempA[0], MonomerA->Ca()->pointer()[0], sizeof(double) * nmoA_ * nsoA_);

  if (nsoA_ != nso_) {
    for (int n=0; n<nsoA_; n++)
      C_DCOPY(nmoA_,tempA[n],1,CA_[n],1);
  }
  else
    C_DCOPY(nso_*nmoA_,tempA[0],1,CA_[0],1);
  free_block(tempA);

  CB_ = block_matrix(nso_,nmoB_);
  double **tempB = block_matrix(nsoB_,nmoB_);
  std::memcpy(tempB[0], MonomerB->Ca()->pointer()[0], sizeof(double) * nmoB_ * nsoB_);
  if (nsoB_ != nso_) {
    for (int n=0; n<nsoB_; n++)
      C_DCOPY(nmoB_,tempB[n],1,CB_[n+nsoA_],1);
  }
  else
    C_DCOPY(nso_*nmoB_,tempB[0],1,CB_[0],1);
  free_block(tempB);

  int nbf[8];
  nbf[0] = nso_;
  std::shared_ptr<MatrixFactory> fact =
    std::shared_ptr<MatrixFactory>(new MatrixFactory);
  fact->init_with(1, nbf, nbf);

  std::shared_ptr<IntegralFactory> intfact =
    std::shared_ptr<IntegralFactory>(new IntegralFactory(basisset_,
    basisset_, basisset_, basisset_));

  std::shared_ptr<OneBodyAOInt> Sint(intfact->ao_overlap());
  SharedMatrix Smat = SharedMatrix
    (fact->create_matrix("Overlap"));
  Sint->compute(Smat);

  double **sIJ = Smat->pointer();
  double **sAJ = block_matrix(nmoA_,nso_);
  sAB_ = block_matrix(nmoA_,nmoB_);

  C_DGEMM('T','N',nmoA_,nso_,nso_,1.0,CA_[0],nmoA_,sIJ[0],nso_,
      0.0,sAJ[0],nso_);
  C_DGEMM('N','N',nmoA_,nmoB_,nso_,1.0,sAJ[0],nso_,CB_[0],nmoB_,
      0.0,sAB_[0],nmoB_);

  free_block(sAJ);

  std::shared_ptr<PotentialInt> potA(static_cast<PotentialInt*>(
    intfact->ao_potential()));
  SharedMatrix ZxyzA(new Matrix("Charges A (Z,x,y,z)", natomsA_, 4));
  for (int n=0, p=0; n<monomerA->natom(); n++) {
    if (monomerA->Z(n)) {
      double Z = (double) monomerA->Z(n);
      double x = monomerA->x(n);
      double y = monomerA->y(n);
      double z = monomerA->z(n);
      ZxyzA->set(0, p, 0, Z);
      ZxyzA->set(0, p, 1, x);
      ZxyzA->set(0, p, 2, y);
      ZxyzA->set(0, p, 3, z);
      p++;
    }
  }
  potA->set_charge_field(ZxyzA);
  SharedMatrix VAmat = SharedMatrix
    (fact->create_matrix("Nuclear Attraction (Monomer A)"));
  potA->compute(VAmat);

  std::shared_ptr<PotentialInt> potB(static_cast<PotentialInt*>(
    intfact->ao_potential()));
  SharedMatrix ZxyzB(new Matrix("Charges B (Z,x,y,z)", natomsB_, 4));
  for (int n=0, p=0; n<monomerB->natom(); n++) {
    if (monomerB->Z(n)) {
      double Z = (double) monomerB->Z(n);
      double x = monomerB->x(n);
      double y = monomerB->y(n);
      double z = monomerB->z(n);
      ZxyzB->set(0, p, 0, Z);
      ZxyzB->set(0, p, 1, x);
      ZxyzB->set(0, p, 2, y);
      ZxyzB->set(0, p, 3, z);
      p++;
    }
  }
  potB->set_charge_field(ZxyzB);
  SharedMatrix VBmat = SharedMatrix
    (fact->create_matrix("Nuclear Attraction (Monomer B)"));
  potB->compute(VBmat);

  double **vIB = block_matrix(nso_,nmoB_);
  double **vAJ = block_matrix(nmoA_,nso_);
  vAAB_ = block_matrix(nmoA_,nmoB_);
  vABB_ = block_matrix(nmoB_,nmoB_);
  vBAA_ = block_matrix(nmoA_,nmoA_);
  vBAB_ = block_matrix(nmoA_,nmoB_);

  double **vIJ = VAmat->pointer();

  C_DGEMM('N','N',nso_,nmoB_,nso_,1.0,vIJ[0],nso_,CB_[0],nmoB_,0.0,
    vIB[0],nmoB_);
  C_DGEMM('T','N',nmoA_,nmoB_,nso_,1.0,CA_[0],nmoA_,vIB[0],nmoB_,0.0,
    vAAB_[0],nmoB_);
  C_DGEMM('T','N',nmoB_,nmoB_,nso_,1.0,CB_[0],nmoB_,vIB[0],nmoB_,0.0,
    vABB_[0],nmoB_);

  vIJ = VBmat->pointer();

  C_DGEMM('T','N',nmoA_,nso_,nso_,1.0,CA_[0],nmoA_,vIJ[0],nso_,0.0,
    vAJ[0],nso_);
  C_DGEMM('N','N',nmoA_,nmoA_,nso_,1.0,vAJ[0],nso_,CA_[0],nmoA_,0.0,
    vBAA_[0],nmoA_);
  C_DGEMM('N','N',nmoA_,nmoB_,nso_,1.0,vAJ[0],nso_,CB_[0],nmoB_,0.0,
    vBAB_[0],nmoB_);

  free_block(vIB);
  free_block(vAJ);
}

void SAPT::get_denom()
{
  std::shared_ptr<Vector> evals_aoccA(new Vector(aoccA_));
  std::shared_ptr<Vector> evals_virA(new Vector(nvirA_));
  std::shared_ptr<Vector> evals_aoccB(new Vector(aoccB_));
  std::shared_ptr<Vector> evals_virB(new Vector(nvirB_));

  for (int a=0; a<aoccA_; a++)
    evals_aoccA->set(0,a,evalsA_[a+foccA_]);
  for (int r=0; r<nvirA_; r++)
    evals_virA->set(0,r,evalsA_[r+noccA_]);
  for (int b=0; b<aoccB_; b++)
    evals_aoccB->set(0,b,evalsB_[b+foccB_]);
  for (int s=0; s<nvirB_; s++)
    evals_virB->set(0,s,evalsB_[s+noccB_]);

  denom_ = SAPTDenominator::buildDenominator(
    options_.get_str("DENOMINATOR_ALGORITHM"),
    evals_aoccA, evals_virA, evals_aoccB, evals_virB,
    options_.get_double("DENOMINATOR_DELTA"), debug_);

  if (debug_ > 1)
    denom_->debug();

  SharedMatrix tauAR = denom_->denominatorA();
  SharedMatrix tauBS = denom_->denominatorB();

  dAR_ = tauAR->pointer();
  dBS_ = tauBS->pointer();

  nvec_ = denom_->nvector();
}

CPHFDIIS::CPHFDIIS(int length, int maxvec)
{
  max_diis_vecs_ = maxvec;
  vec_length_ = length;

  curr_vec_ = 0;
  num_vecs_ = 0;

  t_vecs_ = block_matrix(maxvec,length);
  err_vecs_ = block_matrix(maxvec,length);
}

CPHFDIIS::~CPHFDIIS()
{
  free_block(t_vecs_);
  free_block(err_vecs_);
}

void CPHFDIIS::store_vectors(double *t_vec, double *err_vec)
{
  C_DCOPY(vec_length_,t_vec,1,t_vecs_[curr_vec_],1);
  C_DCOPY(vec_length_,err_vec,1,err_vecs_[curr_vec_],1);

  curr_vec_ = (curr_vec_+1)%max_diis_vecs_;
  num_vecs_++;
  if (num_vecs_ > max_diis_vecs_) num_vecs_ = max_diis_vecs_;
}

void CPHFDIIS::get_new_vector(double *t_vec)
{
  int *ipiv;
  double *Cvec;
  double **Bmat;

  ipiv = init_int_array(num_vecs_+1);
  Bmat = block_matrix(num_vecs_+1,num_vecs_+1);
  Cvec = (double *) malloc((num_vecs_+1)*sizeof(double));

  for (int i=0; i<num_vecs_; i++) {
    for (int j=0; j<=i; j++) {
      Bmat[i][j] = Bmat[j][i] = C_DDOT(vec_length_,err_vecs_[i],1,
        err_vecs_[j],1);
  }}

  for (int i=0; i<num_vecs_; i++) {
    Bmat[num_vecs_][i] = -1.0;
    Bmat[i][num_vecs_] = -1.0;
    Cvec[i] = 0.0;
  }

  Bmat[num_vecs_][num_vecs_] = 0.0;
  Cvec[num_vecs_] = -1.0;

  C_DGESV(num_vecs_+1,1,&(Bmat[0][0]),num_vecs_+1,&(ipiv[0]),&(Cvec[0]),
    num_vecs_+1);

  for (int i=0; i<num_vecs_; i++) {
    C_DAXPY(vec_length_,Cvec[i],t_vecs_[i],1,t_vec,1);
  }

  free(ipiv);
  free(Cvec);
  free_block(Bmat);
}

}}
