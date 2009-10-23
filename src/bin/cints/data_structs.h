#ifndef _psi_src_bin_cints_data_structs_h
#define _psi_src_bin_cints_data_structs_h

/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
//#ifndef DATA_STRUCTS_H
//#define DATA_STRUCTS_H
#include"defines.h"

namespace psi {
  namespace cints {
    /*-----------------------
      structure declarations
      -----------------------*/
    
    /*- 2-index quantity -*/
    struct oebuf {
      double val;
      int ij;
    };
    
    /*- 4-index quantity with index range of (0,65535) -*/
    struct tebuf {
      double val;
      short int i;
      short int j;
      short int k;
      short int l;
    };
    
    
    struct gaussian_function{
      double exp;		/* orbital exponent */
      double ccoeff[CINTS_MAX_AM];	/* comb. of contraction coeff and normalization */
    };
    
    struct shell_def{
      int center;		/* atom on which shell is centered */
      int am;		/* angular momentum of shell (+1?)*/
      int n_prims;		/* number of primitives in shell */
      int fprim;		/* pointer to the first primitive in shell */
      int fbf;             /* pointer to the first basis function from shell */
      int fao;             /* pointer to the first AO from shell */
      double rad_extent;   /* radial extent of the shell (the distance at
			      which the radial part of the first basis function
			      drops below BasisSet.thresh */
      int *trans_vec;      /* shell symmetry transformation vector */
    };
    
    struct coordinates{
      double x;  /*  what do you think these are? */
      double y;
      double z;
      double Z_nuc; /* nuclear charge */
    };
    
    struct shell_pair{
      int i, j;
      double ***P;
      double AB[3];
      double ***PA;
      double ***PB;
      double *a1, *a2, **gamma;
      double *inorm, *jnorm;
      double **Sovlp;
      /*--- Having open-shell and closed-shell and alpha and beta densities
        in UHF case is redundant. Thus only one pair is used at a given time, and the
	other pair is empty (ptr=NULL)!
	
	In fact, the difference open- and closed-shell densities in UHF case are
	used in hf_fock to form the Fock matrix incrementally.
	the total alpha and beta densities in UHF case are used in DFT to form the total
	XC matrix. ---*/
      double **dmat;
      double **dmato;
      double **dmatb;
      double **dmata;
      /*--- SCF Lagrangian ---*/
      double **lagr;
      double Smax;
      double Dmax;
    };
    
    struct unique_shell_pair{
      int *SOpair_npi;
      int **SOpair_so_i;
      int **SOpair_so_j;
      int **SOpair_bf_i;
      int **SOpair_bf_j;
    };
    
    /*--- Matrix element (to be used to handle sparse matrices) ---*/
    typedef struct {
      int row;
      int column;
      double value;
    } mat_elem;
    
    
    /*--- These are basically identical, left from Justin ---*/
    struct double_array {
      int n1;
      double *d;
    };
    typedef struct double_array double_array_t;
    
    struct struct_double_matrix {
      int n1;
      int n2;
      double **d /*[n1]*/ /*[n2]*/;
    };
    typedef struct struct_double_matrix double_matrix_t;
    
    
    enum scftype {rhf = 0, uhf = 1, rohf = 2, twocon = 3};
    enum frametype {canonical = 0, reference = 1};
    
    typedef struct {
      char *wfn;                         /* Wavefunction */
      char *dertype;                     /* Derivative type */
      double cutoff;                     /* Cutoff on ERIs/Fock matrix elements */
      double hf_exch;                    /* Fraction of exact HF exchange in the Fock matrix */
      int  make_dft;                     /* Use DFT? */
      int print_lvl;                     /* Print level */
      long int max_memory;               /* Maximum amount of memory to use, in double words */
      long int memory;                   /* Amount left available */
      int make_oei;                      /* Flag to compute one-electron integrals */
      int make_fock;                     /* Flag to compute Fock matrix */
      int make_eri;                      /* Flag to compute two-electron integrals */
      int make_deriv1;                   /* Flag to compute first derivatives of one- and two-electron integrals */
      int make_deriv1_mvd;               /* Flag to compute first derivatives of SCF+MVD relativisitic correction */
      int make_deriv2;                   /* Flag to compute second derivatives of one- and two-electron integrals */
      int make_oeprop;                   /* Flag tp compute one-electron property integrals */
      int make_mp2;                      /* Flag to compute MP2 energy directly */
      int make_r12ints;                  /* Compute integrals for linear R12 methods */
      int make_mkpt2_ints;               /* Flag to compute MKPT2 integrals directly and dump them to disk */
      int make_mp2r12;                   /* Flag to compute MP2-R12 energy directly */
      int make_cc_bt2;                   /* Flag to compute CC four-virtuals T2 term directly */
      int make_giao_deriv;               /* Flag to compute derivative integrals WRT B and E fields over
					    GIAO Gaussians */
      int symm_ints;                     /* This flag should be set whe individual integrals over SO need to be computed */
      int scf_only;                      /* Means that ERIs will be used only in SCF calculations
					    (may save some space) */
      int num_threads;                   /* Number of threads */
      enum scftype reftype;              /* Reference type, e.g. RHF, ROHF, UHF */
      int restart;                       /* Is this a restart? */
      int restart_task;                  /* Where to restart? */
      struct coordinates origin;         /* user-selected origin for magnetic dipole integrals */
      double fine_structure_alpha;       /* scalar to multiply fine-structure constant */
      double E[3];                       /* electric field vector */
      bool E_given;                      /* Was EFIELD given? */
      enum frametype E_frame;            /* in which frame is the field given? the default is canonical */
      int empirical_dispersion;          /* add grad to empirical dispersion terms? */
    } UserOptions_t;
    
    typedef struct {
      int num_prims;                     /* number of primitive gaussians */
      int num_shells;                    /* number of shells */
      int max_num_prims;                 /* maximum number of primitives per shell */
      int puream;                        /* pure angular momentum flag */
      int num_ao;                        /* number of AO's */
      int max_am;                        /* maximum angular momentum in the basis + 1 */
      int *am2shell;                     /* Mapping array for am ordering to shell ordering */
      int *shells_per_am;                /* Number of shells per am type */
      double thresh;                     /* Threshold used to evaluate radial extents of shells ---*/
      double **schwartz_eri;             /* the matrix num_shells by num_shells:
					    [si][sj] = max(ij|ij) i in si, j in sj  */
      struct shell_def *shells;          /* shell info */
      struct gaussian_function *cgtos;   /* cartesian gaussian information */
      struct shell_pair **shell_pairs;   /* shell pair info */
      
    } BasisSet_t;
    
    typedef struct {
      int nirreps;                       /* number of irreps */
      int max_stab_index;                /* maximum stabilizer index */
      int num_unique_atoms;              /* number of symmetry unique atoms */
      int num_unique_shells;             /* number of symmetry unique shells */
      int num_so;                        /* number of SO's */
      int *atom_positions;               /* symmetry positions/stabilizers of atoms */
      int **ict;                         /* transformation properties of nuclei under symmetry operations */
      int *ua2a;                         /* unique atom number to full atom number mapping array */
      int *us2s;                         /* unique shell number to full shell number mapping array */
      int *sopi;                         /* number of SO per irrep */
      int *sym_oper;                     /* mapping array between "canonical" and symmetry.h-defined
					    ordering of symmetry operations */
      int *so2symblk;                    /* SO number to symmetry block mapping array */
      int **dp_table;                    /* direct product multiplication table */
      int ***dcr;                        /* double coset representatives */
      int **dcr_dim;                     /* dimensions of double coset representatives */
      int **dcr_deg;
      int **GnG;
      int *cdsalcpi;                     /* Number of cartesian displacement SALCs per irrep */
      int *cdsalc_ioffset;               /* offsets for the above */
      char *symlabel;                    /* symmetry label */
      char **irr_labels;                 /* labels of irreps */
      double **cartrep;                  /* cartesian representation matrices */
      double **usotao;                   /* SO to (basis functions if puream && !make_fock, AO otherwise)
					    transformation matrix */
      double **cdsalc2cd;                /* Cartesian displacement SALCs (in columns) */
      struct unique_shell_pair **us_pairs; /* unique shell symmetry info */
    } SymmetryInfo_t;
    

    typedef struct {
      int num_atoms;                     /* number of atoms */
      double Enuc;                       /* nuclear repulsion energy */
      char *label;                       /* calculation label */
      struct coordinates *centers;       /* nuclear centers */
      double **Rref;                     /* rotation back to reference frame */
    } Molecule_t;
    
    typedef struct {
      int itap30;               /* Checkpoint file */
      int itap33;               /* SO ERI file in IWL format */
      int itapS;                /* SO Overlap IWL file */
      int itapT;                /* SO Kinetic energy IWL file */
      int itapV;                /* SO Potential energy IWL file */
      int itapOEInt_Misc;       /* File for all miscellaneous one-electron integrals */
      int itapS_AO;             /* AO Overlap IWL file */
      int itapMX_AO;            /* AO mu(x) IWL file */
      int itapMY_AO;            /* AO mu(y) IWL file */
      int itapMZ_AO;            /* AO mu(z) IWL file */
      int itapQXX_AO;           /* AO q(xx) IWL file */
      int itapQXY_AO;           /* AO q(xy) IWL file */
      int itapQXZ_AO;           /* AO q(xz) IWL file */
      int itapQYY_AO;           /* AO q(yy) IWL file */
      int itapQYZ_AO;           /* AO q(yz) IWL file */
      int itapQZZ_AO;           /* AO q(zz) IWL file */
      int itapNablaX_AO;        /* AO nabla(x) IWL file */
      int itapNablaY_AO;        /* AO nabla(y) IWL file */
      int itapNablaZ_AO;        /* AO nabla(z) IWL file */
      int itapDSCF;             /* "Interface" file between DSCF and CINTS */
      int itapD;                /* Correlated AO OPDM and Lagrangian from transqt */
      int itapG;                /* Correlated AO TPDM from transqt */
      int itapR12;              /* SO integrals of r12 operator */
      int itapT1;               /* SO integrals of [r12,T1] operator */
      int itapERI_MO;           /* MO ERI integrals */
      int itapR12_MO;           /* MO R12 integrals */
      int itapR12T2_MO;         /* MO [r12,T2] integrals */
      int itapdgdB[3];          /* AO dgd/Bi integrals over GIAO Gaussians */
      int itapD1ERI_SO;         /* SO derivative ERI integrals are stored in files itapD1ERI_SO, itapD1ERI_SO+1, ... itapD1ERI_SO+3*natoms */
    } IOUnits_t;
    
    typedef struct {
      double **bf_norm;                  /* "angular" parts of the normalization constants for cartesian GTOs of each
					    angular momentum level */
      double ***cart2pureang;            /* cartesian to pure angular momentum transformation matrices */
      double ****cc2pp;                  /* composite (CxC) cartesian to pure angular momentum transformation
					    matrices */
      mat_elem ****cc2pp_sparse;         /* sparse representation (row-compressed) of cc2pp */
      mat_elem ****pp2cc_sparse;         /* sparse representation (row-compressed) of the reverse of cc2pp */
      int ***cc2pp_rowlength;            /* this holds lengths of "compressed" rows in cc2pp_sparse */
      int ***pp2cc_rowlength;            /* see above */
    } GTOs_t;
    
    typedef struct {
      double **grid;            /* Table of "exact" Fm(T) values. Row index corresponds to
				   values of T (max_T+1 rows), column index to values
				   of m (max_m+1 columns) */
      double delT;              /* The step size for T, depends on cutoff */
      double cutoff;            /* Tolerance cutoff used in all computations of Fm(T) */
      int order_interp;         /* Order of (Taylor) interpolation */
      int max_m;                /* Maximum value of m in the table, depends on cutoff
				   and the number of terms in Taylor interpolation */
      int max_T;                /* Maximum index of T in the table, depends on cutoff
				   and m */
      double *T_crit;           /* Maximum T for each row, depends on cutoff;
				   for a given m and T_idx <= max_T_idx[m] use Taylor interpolation,
				   for a given m and T_idx > max_T_idx[m] use the asymptotic formula */
      void (*compute_Fm)(double *, double, unsigned int);  /* The function which computes a set of Fm(T), 0<=m<=l
							      for given T and l */
    } Fm_Eval_t;
    
    typedef struct {
      double Escf;              /* SCF energy */
      double Ecorr;             /* Correlation energy */
      double Eref;              /* Reference energy (if not SCF reference) */
      double *scf_evals[2];     /* SCF eigenvalues (alpha and beta spin) */
      double *scf_evals_occ[2]; /* Eigenvalues for active occupied orbitals in QTS order */
      double *scf_evals_uocc[2];/* Eigenvalues for active virtual orbitals in QTS order */
      double **scf_evec[2];     /* SCF eigenvectors in AO basis (alpha and beta spin)
				   NOTE: MOs are arranged in rows!!!!! */
      double **scf_evec_occ[2]; /* SCF eigenvectors in AO basis for all
				   doubly-occupied MOs in QTS order:
				   frozen DOCC MOs for each symmetry block come first,
				   then active DOCC MOs for each symmetry block */
      double **scf_evec_uocc[2];/* SCF eigenvectors in AO basis for all
				   vacant MOs in QTS order:
				   active UOCC MOs for each symmetry block come first,
				   then frozen UOCC MOs for each symmetry block */
      double tcscf_occ[2];      /* Squared coefficients of determinants in TCSCF wavefunction */
      double **Alpha, **Beta;   /* Alpha and Beta energy coupling coeffcients */
      int *mo2symblk;           /* Array that maps MO Pitzer index to its symblk number;
				   useful in manipulating Pitzer-indexed MOs */
      int *mo2symblk_occ[2];    /* Array that maps docc index to its symblk number;
				   useful in manipulating QTS-indexed MOs */
      int *mo2symblk_uocc[2];   /* Array that maps uocc index to its symblk number;
				   useful in manipulating QTS-indexed MOs */
      int *orbspi;              /* number of MOs per irrep */
      int *clsdpi;              /* number of doubly-occupied MOs per irrep */
      int *openpi;              /* number of singly-occupied MOs per irrep */
      int *virtpi;              /* number of vacant MOs per irrep */
      int *frozen_docc;         /* number of frozen doubly-occupied MOs per irrep */
      int *frozen_uocc;         /* number of frozen vacant MOs per irrep */
      int num_mo;               /* number of MOs */
      int ndocc;                /* number of doubly-occupied MOs */
      int nfrdocc;              /* number of "frozen" doubly occupied MOs */
      int nactdocc;             /* number of correlated doubly occupied MOs */
      int nsocc;                /* number of singly-occupied MOs */
      int nuocc;                /* number of vacant MOs */
      int nfruocc;              /* number of "frozen" vacant MOs */
      int nactuocc;             /* number of "active" vacant MOs */
      int num_moshells;         /* number of shells of MOs */
      int num_openmoshells;     /* number of shells of singly-occupied MOs */
      int alpha_occ;            /* number of alpha occupied orbitals */
      int beta_occ;             /* number of beta  occupied orbitals */
      int *occ_to_pitzer;       /* The occupied (frozen docc + docc + active) to Pitzer
                                   array for the MkPT2 routine */
      int *vir_to_pitzer;       /* The virtual (active + virtual) to Pitzer
                                   array for the MkPT2 routine */
//these were added by ACS for the direct ump2r12 routine (01/06)
    double **scf_evec_alpha;     /* alpha SCF eigenvectors in AO basis
                                 NOTE: MOs are arranged in rows!!!!! */
    double **scf_evec_beta;     /*  beta SCF eigenvectors in AO basis
                                 NOTE: MOs are arranged in rows!!!!! */
    double **scf_evec_occ_alpha; /* alpha SCF eigenvectors in AO basis for all
                                 doubly-occupied MOs in QTS order:
                                 frozen DOCC MOs for each symmetry block come first,
                                 then active DOCC MOs for each symmetry block */
    double **scf_evec_occ_beta; /* beta SCF eigenvectors in AO basis for all
                                 doubly-occupied MOs in QTS order:
                                 frozen DOCC MOs for each symmetry block come first,
                                 then active DOCC MOs for each symmetry block */
    int *virtpi_alpha;        /* number of vacant alpha MOs per irrep */
    int *virtpi_beta;         /* number of vacant beta MOs per irrep */
    int alpha_act_occ;        /* number of active alpha occupied orbitals */
    int beta_act_occ;         /* number of active beta  occupied orbitals */
    int *mo2symblk_occ_alpha; /* Array that maps alpha occ index to its symblk number;
                                 useful in manipulating QTS-indexed MOs */
    int *mo2symblk_occ_beta;  /* Array that maps beta occ index to its symblk number;
                                 useful in manipulating QTS-indexed MOs */

    } MOInfo_t;
    
    
    typedef struct {
      double **T2_s;  /* source T2's */
      double **T2_t;  /* target T2's */
      int nvirt;  /* no. active virtuals */
      int nocc;   /* no. active occupieds */
    } CCInfo_t;
    
    /* Cartesian derivative SALCs */
    typedef struct {
      int nsalcs;
      char *atom_irreps;            /* Bit-packed irreps of all derivatives wrt coordinates of a given atom.
				       Assume Abelian groups, i.e. at most 8 irreps */
#ifdef __cplusplus
      typedef struct {
	int nsalcs;
	int* salcs;
      } cd2salc_map_t;
      cd2salc_map_t* cd2salc_map;   /* Maps cartesian derivative to the list of all SALCs to which it contributes */
#else
      struct cd2salc_map_t {
	int nsalcs;
	int* salcs;
      } *cd2salc_map;   /* Maps cartesian derivative to the list of all SALCs to which it contributes */
#endif
      
      int* salc2irrep;              /* Maps SALC to its irrep */
      
    } CDSALC_t;
    
    /* -------------------------------------------------------
       
    DFT Data Structures
    
    -------------------------------------------------------*/
    /* -------------------------
       Pruned Grid parameters
       ------------------------*/
    struct param_set_s { 
      int n_ang_grids;               /* Number of different lebedev spheres 
					in atomic grid */
      double *alpha;                 /* The cutoff parameters */
      int *angpoints;                /* the different angular grids */
    };
    
    struct pruned_info_s {
      int *a2param;                  /* Elements tell which parameter set 
					is used for each atom */
      int n_diff_ang_grids;          /* If there are more than one sets 
					of angular grids, this is how many */
      int n_tot_ang_grids;         /* Total number of different angular grids */
      int n_param_sets;              /* Number of total parameter sets */
      struct param_set_s *param_set; /* The parameter sets */
    };
    
    /*------------------------------
      Primitve class types
      -----------------------------*/
    
    typedef struct{
      struct coordinates p_cart;
      double ang_weight;
    } leb_point_t;
    
    typedef struct{
      double r;
      double drdq;
      int n_ang_points;
      leb_point_t *points;
    } leb_sphere_t;
    
    typedef struct{
      int size;
      int radial_start;
      int radial_end;
      leb_sphere_t *spheres;
    } prim_leb_chunk_t;
    
    typedef struct{
      int chunk_num;
      prim_leb_chunk_t *leb_chunk;
    } prim_atomic_grid_t;
    
    /* ----------------------
       Concrete classes
       ----------------------*/

    struct close_shell_info_s{
      int num_close_aos;
      int num_close_shells;
      int *shells_close_to_chunk;
      int *aos_close_to_chunk;
      int *close_shells_per_am;
      double **close_COCC;
      double **close_COCC_a;
      double **close_COCC_b;
    };
    
    struct leb_chunk_s{
      int size;
      int radial_start;
      int radial_end;
      leb_sphere_t *spheres;
      int *shells_close_to_chunk;
      int *close_shells_per_am;
      double **close_COCC;
    };
    
    struct atomic_grid_s{
      int atom_num;
      int atom_degen;
      struct coordinates atom_center;
      double Bragg_radii;
      int chunk_num;
      struct leb_chunk_s *leb_chunk;
    };
    
    typedef struct{
      int n_rad_points;
      int pruned_flag;
      char *label;
      struct atomic_grid_s *atomic_grid;
      prim_atomic_grid_t prim_atomic_grid;
      prim_atomic_grid_t *prim_pruned_atomic_grids;
      struct pruned_info_s pruned_info;
    } grid_t;
    
    struct den_info_s{
      double den;
      double dena;
      double denb;
      double gradx;
      double grady;
      double gradz;
      double gamma;
    };
    
    struct fun_info_s{
      double eval;
      double dval;
      double dvala;
      double dvalb;
      double dpval;
      double dgval;
      double ddval;
      double ddvala;
      double ddvalb;
    };
    
    struct xc_info_s{
      struct fun_info_s exch_info;
      struct fun_info_s corr_info;
    };
    
    
    typedef struct{
      int prtflag;                /* dft printing flag */
      
      double *basis;              /* This is an array to hold the value of
				     basis functions at a given point */
      
      double *gamma_basis;
      double *gradx;
      double *grady;
      double *gradz;
      
      double *Bragg;
      double XC_energy;           /* Exchange Correlation Energy */
      double X_energy;            /* Exchange Energy */
      double C_energy;            /* Correlation Energy */
      
    /* All function pointers */
      
      struct fun_info_s (*exchange_func)(struct den_info_s);
      struct fun_info_s (*correlation_func)(struct den_info_s);
      
      struct den_info_s (*den_calc)(struct coordinates,int);
      
      struct close_shell_info_s close_shell_info;
      
      grid_t grid;
      
    } DFT_options_t;
    
    
  }
}
    
#endif
//#endif

