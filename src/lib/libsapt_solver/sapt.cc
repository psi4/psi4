#include "sapt.h"

namespace psi { namespace sapt {

SAPT::SAPT(Options& options, boost::shared_ptr<PSIO> psio, 
  boost::shared_ptr<Chkpt> chkpt) : Wavefunction(options, psio, chkpt)
{
#ifdef HAVE_MKL
  mkl_set_dynamic(1);
#endif

#ifdef _OPENMP
  omp_set_nested(0);
#endif

  initialize();
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
  ribasis_.reset();
  zero_.reset();
}

void SAPT::initialize()
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

  boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
  ribasis_ = boost::shared_ptr<BasisSet>(BasisSet::construct(parser, molecule_, 
    "DF_BASIS_SAPT"));
  elst_basis_ = 0;
  if (options_.get_str("DF_BASIS_ELST") != "") {
    elstbasis_ = boost::shared_ptr<BasisSet>(BasisSet::construct(parser, 
      molecule_,"DF_BASIS_ELST"));
    elst_basis_ = 1;
  }
  zero_ = boost::shared_ptr<BasisSet>(BasisSet::zero_ao_basis_set());
  parser.reset();

  print_ = options_.get_int("PRINT");
  debug_ = options_.get_int("DEBUG");
  schwarz_ = options_.get_double("INTS_TOLERANCE");
  mem_ = (long int) ((double) memory_*options_.get_double("SAPT_MEM_SAFETY"));
  mem_ /= 8L;

  std::vector<int> realsA;
  realsA.push_back(0);
  std::vector<int> ghostsA;
  ghostsA.push_back(1);
  boost::shared_ptr<Molecule> monomerA = molecule_->extract_subsets(realsA,
    ghostsA);
  foccA_ = monomerA->nfrozen_core(options_.get_str("FREEZE_CORE"));

  std::vector<int> realsB;
  realsB.push_back(1);
  std::vector<int> ghostsB;
  ghostsB.push_back(0);
  boost::shared_ptr<Molecule> monomerB = molecule_->extract_subsets(realsB,
    ghostsB);
  foccB_ = monomerB->nfrozen_core(options_.get_str("FREEZE_CORE"));

  natomsA_ = 0;
  natomsB_ = 0;

  for (int n=0; n<monomerA->natom(); n++)
    if (monomerA->Z(n)) natomsA_++;
  for (int n=0; n<monomerB->natom(); n++)
    if (monomerB->Z(n)) natomsB_++;

  ndf_ = ribasis_->nbf();

  psio_->open(PSIF_SAPT_DIMER,PSIO_OPEN_OLD);
  psio_->open(PSIF_SAPT_MONOMERA,PSIO_OPEN_OLD);
  psio_->open(PSIF_SAPT_MONOMERB,PSIO_OPEN_OLD);

  double enucD, enucA, enucB;
  double eHFD, eHFA, eHFB;

  psio_->read_entry(PSIF_SAPT_DIMER,"Dimer NSO",(char *) &nso_,sizeof(int));
  psio_->read_entry(PSIF_SAPT_DIMER,"Dimer NMO",(char *) &nmo_,sizeof(int));
  psio_->read_entry(PSIF_SAPT_DIMER,"Dimer HF Energy",(char *) &eHFD, 
    sizeof(double));
  psio_->read_entry(PSIF_SAPT_DIMER,"Dimer Nuclear Repulsion Energy",(char *)
    &enucD, sizeof(double));

  psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer NSO",(char *) &nsoA_, 
    sizeof(int));
  psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer NMO",(char *) &nmoA_, 
    sizeof(int));
  psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer NOCC",(char *) &noccA_, 
    sizeof(int));
  psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer NVIR",(char *) &nvirA_, 
    sizeof(int));
  psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer Number of Electrons",(char *)
    &NA_, sizeof(int));
  psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer HF Energy",(char *) &eHFA, 
    sizeof(double));
  psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer Nuclear Repulsion Energy",
    (char *) &enucA, sizeof(double));

  aoccA_ = noccA_ - foccA_;

  psio_->read_entry(PSIF_SAPT_MONOMERB,"Monomer NSO",(char *) &nsoB_, 
    sizeof(int));
  psio_->read_entry(PSIF_SAPT_MONOMERB,"Monomer NMO",(char *) &nmoB_, 
    sizeof(int));
  psio_->read_entry(PSIF_SAPT_MONOMERB,"Monomer NOCC",(char *) &noccB_, 
    sizeof(int));
  psio_->read_entry(PSIF_SAPT_MONOMERB,"Monomer NVIR",(char *) &nvirB_, 
    sizeof(int));
  psio_->read_entry(PSIF_SAPT_MONOMERB,"Monomer Number of Electrons",(char *)
    &NB_, sizeof(int));
  psio_->read_entry(PSIF_SAPT_MONOMERB,"Monomer HF Energy",(char *) &eHFB, 
    sizeof(double));
  psio_->read_entry(PSIF_SAPT_MONOMERB,"Monomer Nuclear Repulsion Energy",
    (char *) &enucB, sizeof(double));

  aoccB_ = noccB_ - foccB_;

  enuc_ = enucD - enucA - enucB;
  eHF_ =  eHFD - eHFA - eHFB;

  evalsA_ = init_array(nmoA_);
  psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer HF Eigenvalues",(char *)
    &(evalsA_[0]), sizeof(double)*nmoA_);

  evalsB_ = init_array(nmoB_);
  psio_->read_entry(PSIF_SAPT_MONOMERB,"Monomer HF Eigenvalues",(char *)
    &(evalsB_[0]), sizeof(double)*nmoB_);

  CA_ = block_matrix(nso_,nmoA_);
  double **tempA = block_matrix(nsoA_,nmoA_);
  psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer HF Coefficients",(char *)
    &(tempA[0][0]), sizeof(double)*nmoA_*nsoA_);
  if (nsoA_ != nso_) {
    for (int n=0; n<nsoA_; n++)
      C_DCOPY(nmoA_,tempA[n],1,CA_[n],1);
  }
  else
    C_DCOPY(nso_*nmoA_,tempA[0],1,CA_[0],1);
  free_block(tempA);

  CB_ = block_matrix(nso_,nmoB_);
  double **tempB = block_matrix(nsoB_,nmoB_);
  psio_->read_entry(PSIF_SAPT_MONOMERB,"Monomer HF Coefficients",(char *)
    &(tempB[0][0]), sizeof(double)*nmoB_*nsoB_);
  if (nsoB_ != nso_) {
    for (int n=0; n<nsoB_; n++)
      C_DCOPY(nmoB_,tempB[n],1,CB_[n+nsoA_],1);
  }
  else
    C_DCOPY(nso_*nmoB_,tempB[0],1,CB_[0],1);
  free_block(tempB);

  psio_->close(PSIF_SAPT_DIMER,1);
  psio_->close(PSIF_SAPT_MONOMERA,1);
  psio_->close(PSIF_SAPT_MONOMERB,1);

  int nbf[8];
  nbf[0] = nso_;
  boost::shared_ptr<MatrixFactory> fact = 
    boost::shared_ptr<MatrixFactory>(new MatrixFactory);
  fact->init_with(1, nbf, nbf);

  boost::shared_ptr<IntegralFactory> intfact = 
    boost::shared_ptr<IntegralFactory>(new IntegralFactory(basisset_, 
    basisset_, basisset_, basisset_));

  boost::shared_ptr<OneBodyAOInt> Sint(intfact->ao_overlap());
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

  boost::shared_ptr<PotentialInt> potA(static_cast<PotentialInt*>(
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

  boost::shared_ptr<PotentialInt> potB(static_cast<PotentialInt*>(
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
  boost::shared_ptr<Vector> evals_aoccA(new Vector(aoccA_));
  boost::shared_ptr<Vector> evals_virA(new Vector(nvirA_));
  boost::shared_ptr<Vector> evals_aoccB(new Vector(aoccB_));
  boost::shared_ptr<Vector> evals_virB(new Vector(nvirB_));

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
