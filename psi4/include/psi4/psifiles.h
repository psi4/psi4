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

/*
** PSIFILES.H
**
** This header file contains the definitions of the numbers assigned
**  to various binary files in PSI.  This was created primarily to
**  help avoid conflicts in the numbering of new PSI files in developmental
**  programs but will grow to encompass some older binary files.
**
** This additional level of abstraction will aid in the maintenance of
**  code.  You are strongly encouraged to refer to files using these
**  definitions rather than the actual numbers; the numbers may change
**  in the future but the names will not.
**
** Created by C. David Sherrill on 29 April 1998
**
** Make sure that you comment any new additions to this, so they are inlined 
**  into the manual as well as the psifiles.py python file.  Note the format 
**  of the existing comment markers, and keep comments on one line.
*/

#ifndef _psi_include_psifiles_h_
#define _psi_include_psifiles_h_

#define PSI_DEFAULT_FILE_PREFIX "psi"

#define PSIF_OPTKING           1    /*-  -*/
#define PSIF_GRAD              11   /*- geometry optimization, geometry, and gradient; currently is an ASCII file like output.grad -*/
#define PSIF_INTCO             12   /*- internal coordinates file, currently is ASCII file like output.intco -*/
#define PSIF_3INDEX            16   /*-  -*/
#define PSIF_DSCF              31   /*-  -*/
#define PSIF_CHKPT             32   /*- new libpsio checkpoint file number -*/
#define PSIF_SO_TEI            33   /*-  -*/
#define PSIF_SO_PK             34   /*-  -*/
#define PSIF_OEI               35   /*-  -*/
#define PSIF_SO_ERF_TEI        36   /*-  -*/
#define PSIF_SO_ERFC_TEI       37   /*-  -*/
#define PSIF_SO_R12            38   /*-  -*/
#define PSIF_SO_R12T1          39   /*-  -*/
#define PSIF_DERINFO           40   /*-  -*/
#define PSIF_SO_PRESORT        41   /*-  -*/
#define PSIF_OLD_CHKPT         42   /*- Until we have flexible PSIF_CHKPT this will store previous calculation info -*/
#define PSIF_CIVECT            43   /*- CI vector from DETCI along with string and determinant info -*/

#define PSIF_AO_DGDBX          44   /*- B-field derivative AO integrals over GIAO Gaussians -- only bra-ket permutational symmetry holds -*/
#define PSIF_AO_DGDBY          45   /*-  -*/
#define PSIF_AO_DGDBZ          46   /*-  -*/
/* PSIMRCC files */
#define PSIF_PSIMRCC_INTEGRALS 50   /*-  -*/
#define PSIF_PSIMRCC_RESTART   51   /*-  -*/
/* MCSCF files */
#define PSIF_MCSCF             52   /*-  -*/
#define PSIF_TPDM_HALFTRANS    53   /*-  -*/
#define PSIF_DETCAS            60   /*-  -*/
// The integral files used by libtrans
#define PSIF_LIBTRANS_DPD      61   /*- libtrans: All transformed integrals in DPD format are sent here by default -*/
#define PSIF_LIBTRANS_A_HT     62   /*- libtrans: Alpha half-transformed integrals in DPD format -*/
#define PSIF_LIBTRANS_B_HT     63   /*- libtrans: Beta half-tranformed integrals in DPD format -*/

#define PSIF_LIBDIIS           64   /*- Storage for libdiis -*/
#define PSIF_DFT_GRID          65   /*- Storage for DFT/pseudospectral grid -*/
// DFCC 3-index files
#define PSIF_DF_TENSOR         66   /*-  -*/
#define PSIF_PS_TENSOR         67   /*-  -*/

#define PSIF_TPDM_PRESORT      71   /*-  -*/
#define PSIF_MO_TEI            72   /*-  -*/
#define PSIF_MO_OPDM           73   /*-  -*/
#define PSIF_MO_TPDM           74   /*-  -*/
#define PSIF_MO_LAG            75   /*-  -*/
#define PSIF_AO_OPDM           76   /*- PSIF_AO_OPDM also contains AO Lagrangian -*/
#define PSIF_AO_TPDM           77   /*-  -*/

#define PSIF_DBOC              78   /*- dboc: -*/

#define PSIF_MO_R12            79   /*-  -*/
#define PSIF_MO_R12T2          80   /*-  -*/

/*
** Additions for UHF-based transformations.
** -TDC, 6/01
*/
#define PSIF_MO_AA_TEI         81   /*-  -*/
#define PSIF_MO_BB_TEI         82   /*-  -*/
#define PSIF_MO_AB_TEI         83   /*-  -*/
#define PSIF_MO_AA_TPDM        84   /*-  -*/
#define PSIF_MO_BB_TPDM        85   /*-  -*/
#define PSIF_MO_AB_TPDM        86   /*-  -*/
#define PSIF_AA_PRESORT        87   /*- AA UHF twopdm presort file -*/
#define PSIF_BB_PRESORT        88   /*- BB UHF twopdm presort file -*/
#define PSIF_AB_PRESORT        89   /*- AB UHF twopdm presort file -*/
/*
** MO Hessian File (also contains specialized integral and Fock lists.
** See programs STABLE and CPHF for more info.
** -TDC, 7/00
*/
#define PSIF_MO_HESS           90   /*-  -*/
#define PSIF_CPHF              91   /*-  -*/

#define PSIF_SO_PKSUPER1       92   /*-  -*/
#define PSIF_SO_PKSUPER2       93   /*-  -*/

// The half-transformed integrals
#define PSIF_HALFT0            94   /*-  -*/
#define PSIF_HALFT1            95   /*-  -*/

#define PSIF_DFSCF_A           96   /*- B Matrix containing 3-index tensor in AOs for use with DF-SCF -*/
#define PSIF_DFSCF_BJ          97   /*- B Matrix containing 3-index tensor in AOs with J^-1/2 for use with DF-SCF -*/
#define PSIF_DFSCF_K           98   /*- Exchange tensor for DF-SCF -*/
#define PSIF_DFSCF_BJI         99   /*- The three-center integrals for DF-SCF -*/

/* The miscellaneous CC information file */
#define PSIF_CC_INFO        100   /*-  -*/
/* One-electron integral files */
#define PSIF_CC_OEI         101   /*-  -*/
/* Two-electron integral files */   /* pqnum  rsnum */
#define PSIF_CC_AINTS       102   /*-  -*/
#define PSIF_CC_BINTS       103   /*-  -*/
#define PSIF_CC_CINTS       104   /*-  -*/
#define PSIF_CC_DINTS       105   /*-  -*/
#define PSIF_CC_EINTS       106   /*-  -*/
#define PSIF_CC_FINTS       107   /*-  -*/

/* Two-electron amplitudes, intermediates, and densities */
#define PSIF_CC_DENOM       108   /*-  -*/
#define PSIF_CC_TAMPS       109   /*-  -*/
#define PSIF_CC_GAMMA       110   /*-  -*/
#define PSIF_CC_MISC        111   /*-  -*/
#define PSIF_CC_HBAR        112   /*-  -*/

#define PSIF_CC_OEI_NEW     113   /*-  -*/
#define PSIF_CC_GAMMA_NEW   114   /*-  -*/
#define PSIF_CC_AINTS_NEW   115   /*-  -*/
#define PSIF_CC_BINTS_NEW   116   /*-  -*/
#define PSIF_CC_CINTS_NEW   117   /*-  -*/
#define PSIF_CC_DINTS_NEW   118   /*-  -*/
#define PSIF_CC_EINTS_NEW   119   /*-  -*/
#define PSIF_CC_FINTS_NEW   120   /*-  -*/

/* ground state lambda and intermediates for excited states */
#define PSIF_CC_LAMBDA      121   /*-  -*/

/* converged eigenvectors of hbar */
#define PSIF_CC_RAMPS       122   /*-  -*/
#define PSIF_CC_LAMPS       123   /*-  -*/

#define PSIF_CC_LR          124   /*-  -*/

#define PSIF_CC_DIIS_ERR    125   /*-  -*/
#define PSIF_CC_DIIS_AMP    126   /*-  -*/

#define PSIF_CC_TMP         127   /*-  -*/
#define PSIF_CC_TMP0        128   /*-  -*/
#define PSIF_CC_TMP1        129   /*-  -*/
#define PSIF_CC_TMP2        130   /*-  -*/
#define PSIF_CC_TMP3        131   /*-  -*/
#define PSIF_CC_TMP4        132   /*-  -*/
#define PSIF_CC_TMP5        133   /*-  -*/
#define PSIF_CC_TMP6        134   /*-  -*/
#define PSIF_CC_TMP7        135   /*-  -*/
#define PSIF_CC_TMP8        135   /*-  -*/
#define PSIF_CC_TMP9        137   /*-  -*/
#define PSIF_CC_TMP10       138   /*-  -*/
#define PSIF_CC_TMP11       139   /*-  -*/
/* temporary files for CCEOM and CCLAMBDA */
#define PSIF_EOM_D          140   /*-  -*/
#define PSIF_EOM_CME        141   /*-  -*/
#define PSIF_EOM_Cme        142   /*-  -*/
#define PSIF_EOM_CMNEF      143   /*-  -*/
#define PSIF_EOM_Cmnef      144   /*-  -*/
#define PSIF_EOM_CMnEf      145   /*-  -*/
#define PSIF_EOM_SIA        146   /*-  -*/
#define PSIF_EOM_Sia        147   /*-  -*/
#define PSIF_EOM_SIJAB      148   /*-  -*/
#define PSIF_EOM_Sijab      149   /*-  -*/
#define PSIF_EOM_SIjAb      150   /*-  -*/
#define PSIF_EOM_R          151 /*- holds residual -*/
#define PSIF_CC_GLG         152 /*- left-hand psi for g.s. parts of cc-density -*/
#define PSIF_CC_GL          153 /*- left-hand psi for e.s. parts of cc-density -*/
#define PSIF_CC_GR          154 /*- right-hand eigenvector for cc-density -*/
#define PSIF_EOM_TMP1       155 /*- intermediates just for single contractions -*/
#define PSIF_EOM_TMP0       156 /*- temporary copies of density -*/
#define PSIF_EOM_TMP_XI     157 /*- intermediates for xi computation -*/
#define PSIF_EOM_XI         158 /*- xi = dE/dt amplitudes -*/
#define PSIF_EOM_TMP        159 /*- intermediates used more than once -*/
#define PSIF_CC3_HET1       160 /*- [H,e^T1] -*/
#define PSIF_CC3_HC1        161 /*- [H,C1] -*/
#define PSIF_CC3_HC1ET1     162 /*- [[H,e^T1],C1] -*/
#define PSIF_CC3_MISC       163 /*- various intermediates needed in CC3 codes -*/
#define PSIF_CC2_HET1       164 /*- [H,e^T1] -*/

/* Markers for the first and last CC file numbers */
#define PSIF_CC_MIN  PSIF_CC_INFO  /*-  -*/
#define PSIF_CC_MAX  PSIF_CC2_HET1 /*-  -*/

#define PSIF_WK_PK          165  /*- File to contain wK pre-sorted integrals for PK -*/

#define PSIF_SCF_MOS           180  /*- Save SCF orbitals for re-use later as guess, etc. -*/
#define PSIF_DFMP2_AIA         181  /*- Unfitted three-index MO ints for DFMP2 -*/
#define PSIF_DFMP2_QIA         182  /*- Fitted-three index MO ints for DFMP2 -*/

#define PSIF_ADC               183  /*- ADC -*/
#define PSIF_ADC_SEM           184  /*- ADC -*/

#define PSIF_SAPT_DIMER        190  /*- SAPT Two-Body Dimer -*/
#define PSIF_SAPT_MONOMERA     191  /*- SAPT Two-Body Mon A -*/
#define PSIF_SAPT_MONOMERB     192  /*- SAPT Two-Body Mon B -*/

#define PSIF_SAPT_AA_DF_INTS   193  /*- SAPT AA DF Ints -*/
#define PSIF_SAPT_AB_DF_INTS   194  /*- SAPT AB DF Ints -*/
#define PSIF_SAPT_BB_DF_INTS   195  /*- SAPT BB DF Ints -*/
#define PSIF_SAPT_AMPS         196  /*- SAPT Amplitudes -*/
#define PSIF_SAPT_TEMP         197  /*- SAPT Temporary worlds fastest code file -*/

#define PSIF_SAPT_LRINTS       198  /*- SAPT0 2-Body linear response LDA integrals -*/

#define PSIF_SO_D1OEI          199  /*- Derivative OEIs are stored in file 199 -*/
#define PSIF_SO_D1ERI          200  /*- Derivative ERIs are stored in files 200, 201, 202, etc. File 200 -*/

#define PSIF_3B_SAPT_TRIMER              220  /*- SAPT Three-Body Trimer -*/
#define PSIF_3B_SAPT_DIMER_AB            221  /*- SAPT Three-Body Dimer AB -*/
#define PSIF_3B_SAPT_DIMER_AC            222  /*- SAPT Three-Body Dimer AC -*/
#define PSIF_3B_SAPT_DIMER_BC            223  /*- SAPT Three-Body Dimer BC -*/
#define PSIF_3B_SAPT_MONOMER_A           224  /*- SAPT Three-Body Mon A -*/
#define PSIF_3B_SAPT_MONOMER_B           225  /*- SAPT Three-Body Mon B -*/
#define PSIF_3B_SAPT_MONOMER_C           226  /*- SAPT Three-Body Mon C -*/
#define PSIF_3B_SAPT_AA_DF_INTS          227  /*-  -*/
#define PSIF_3B_SAPT_BB_DF_INTS          228  /*-  -*/
#define PSIF_3B_SAPT_CC_DF_INTS          229  /*-  -*/
#define PSIF_3B_SAPT_AMPS                230  /*-  -*/

#define PSIF_DCC_IJAK          250  /*- CEPA/CC (ij|ak) -*/
#define PSIF_DCC_IJAK2         251  /*- CEPA/CC (ij|ak) -*/
#define PSIF_DCC_ABCI          252  /*- CEPA/CC (ia|bc) -*/
#define PSIF_DCC_ABCI2         253  /*- CEPA/CC (ia|bc) -*/
#define PSIF_DCC_ABCI3         254  /*- CEPA/CC (ia|bc) -*/
#define PSIF_DCC_ABCI4         255  /*- CEPA/CC (ia|bc) -*/
#define PSIF_DCC_ABCI5         256  /*- CEPA/CC (ia|bc) -*/
#define PSIF_DCC_ABCD1         257  /*- CEPA/CC (ab|cd)+ -*/
#define PSIF_DCC_ABCD2         258  /*- CEPA/CC (ab|cd)- -*/
#define PSIF_DCC_IJAB          259  /*- CEPA/CC (ij|ab) -*/
#define PSIF_DCC_IAJB          260  /*- CEPA/CC (ia|jb) -*/
#define PSIF_DCC_IJKL          261  /*- CEPA/CC (ij|kl) -*/
#define PSIF_DCC_OVEC          262  /*- CEPA/CC old vectors for diis -*/
#define PSIF_DCC_EVEC          263  /*- CEPA/CC error vectors for diis -*/
#define PSIF_DCC_R2            264  /*- CEPA/CC residual -*/
#define PSIF_DCC_TEMP          265  /*- CEPA/CC temporary storage -*/
#define PSIF_DCC_T2            266  /*- CEPA/CC t2 amplitudes -*/
#define PSIF_DCC_QSO           267  /*- DFCC 3-index integrals -*/
#define PSIF_DCC_SORT_START    270  /*- CEPA/CC integral sort starting file number -*/

#define PSIF_SAPT_CCD          271  /*- SAPT2+ CCD Utility File -*/

#define PSIF_HESS              272  /*- Hessian Utility File -*/

#define PSIF_OCC_DPD           273  /*- OCC DPD -*/
#define PSIF_OCC_DENSITY       274  /*- OCC Density -*/
#define PSIF_OCC_IABC          275  /*- OCC out-of-core <IA|BC> -*/
#define PSIF_DFOCC_INTS        276  /*- DFOCC Integrals -*/ 
#define PSIF_DFOCC_AMPS        277  /*- DFOCC Amplitudes -*/
#define PSIF_DFOCC_DENS        278  /*- DFOCC PDMs -*/
#define PSIF_DFOCC_IABC        279  /*- DFOCC (IA|BC) -*/ 
#define PSIF_DFOCC_TEMP        280  /*- DFOCC temporary storage -*/

#define PSIF_SAD               300  /*- A SAD file (File for SAD related quantities -*/

/* All of these one-electron quantities have been moved into PSIF_OEI
   Most integrals are real Hermitian hence only lower triangle of the matrix is written out */
/* These macros give libpsio TOC strings for easy identification.     */
#define PSIF_SO_S           "SO-basis Overlap Ints"
#define PSIF_SO_T           "SO-basis Kinetic Energy Ints"
#define PSIF_SO_V           "SO-basis Potential Energy Ints"
#define PSIF_SO_H           "SO-basis Core Hamiltonian"
#define PSIF_AO_S           "AO-basis Overlap Ints"
#define PSIF_AO_MX          "AO-basis Mu-X Ints"
#define PSIF_AO_MY          "AO-basis Mu-Y Ints"
#define PSIF_AO_MZ          "AO-basis Mu-Z Ints"
#define PSIF_MO_MX          "MO-basis Mu-X Ints"
#define PSIF_MO_MY          "MO-basis Mu-Y Ints"
#define PSIF_MO_MZ          "MO-basis Mu-Z Ints"
#define PSIF_AO_QXX         "AO-basis Q-XX Ints"    /* Electric quadrupole moment integrals */
#define PSIF_AO_QXY         "AO-basis Q-XY Ints"
#define PSIF_AO_QXZ         "AO-basis Q-XZ Ints"
#define PSIF_AO_QYY         "AO-basis Q-YY Ints"
#define PSIF_AO_QYZ         "AO-basis Q-YZ Ints"
#define PSIF_AO_QZZ         "AO-basis Q-ZZ Ints"
#define PSIF_AO_TXX         "AO-basis T-XX Ints"    /* Traceless electric quadrupole moment integrals */
#define PSIF_AO_TXY         "AO-basis T-XY Ints"    /* Traceless electric quadrupole moment integrals */
#define PSIF_AO_TXZ         "AO-basis T-XZ Ints"    /* Traceless electric quadrupole moment integrals */
#define PSIF_AO_TYY         "AO-basis T-YY Ints"    /* Traceless electric quadrupole moment integrals */
#define PSIF_AO_TYZ         "AO-basis T-YZ Ints"    /* Traceless electric quadrupole moment integrals */
#define PSIF_AO_TZZ         "AO-basis T-ZZ Ints"    /* Traceless electric quadrupole moment integrals */

/* These integrals are anti-Hermitian -- upper triangle has sign opposite of that of the lower triangle */
#define PSIF_AO_NablaX      "AO-basis Nabla-X Ints" /* integrals of nabla operator */
#define PSIF_AO_NablaY      "AO-basis Nabla-Y Ints"
#define PSIF_AO_NablaZ      "AO-basis Nabla-Z Ints"

/* These integrals are pure imaginary Hermitian. We write the full matrix of the imaginary part of these
integrals out (i.e. multiply by i=sqrt(-1) to get the integrals) */
#define PSIF_AO_LX          "AO-basis LX Ints"      /* integrals of angular momentum operator */
#define PSIF_AO_LY          "AO-basis LY Ints"
#define PSIF_AO_LZ          "AO-basis LZ Ints"
#define PSIF_AO_DSDB_X      "AO-basis dS/dBx Ints"      /* Overlap derivative integrals WRT B field */
#define PSIF_AO_DSDB_Y      "AO-basis dS/dBy Ints"
#define PSIF_AO_DSDB_Z      "AO-basis dS/dBz Ints"
#define PSIF_AO_DHDB_X      "AO-basis dh/dBx Ints"      /* One-electron derivative integrals WRT B field */
#define PSIF_AO_DHDB_Y      "AO-basis dh/dBy Ints"
#define PSIF_AO_DHDB_Z      "AO-basis dh/dBz Ints"
#define PSIF_AO_D2HDBDE_XX  "AO-basis d2h/dBxdEx Ints"  /* One-electron derivative integrals WRT E and B fields */
#define PSIF_AO_D2HDBDE_XY  "AO-basis d2h/dBxdEy Ints"
#define PSIF_AO_D2HDBDE_XZ  "AO-basis d2h/dBxdEz Ints"
#define PSIF_AO_D2HDBDE_YX  "AO-basis d2h/dBydEx Ints"
#define PSIF_AO_D2HDBDE_YY  "AO-basis d2h/dBydEy Ints"
#define PSIF_AO_D2HDBDE_YZ  "AO-basis d2h/dBydEz Ints"
#define PSIF_AO_D2HDBDE_ZX  "AO-basis d2h/dBzdEx Ints"
#define PSIF_AO_D2HDBDE_ZY  "AO-basis d2h/dBzdEy Ints"
#define PSIF_AO_D2HDBDE_ZZ  "AO-basis d2h/dBzdEz Ints"
#define PSIF_MO_DFDB_X      "AO-basis dF/dBx Ints"      /* Fock operator derivative integrals WRT B field */
#define PSIF_MO_DFDB_Y      "AO-basis dF/dBy Ints"
#define PSIF_MO_DFDB_Z      "AO-basis dF/dBz Ints"

#define PSIF_MO_OEI         "MO-basis One-electron Ints"
#define PSIF_MO_A_OEI       "MO-basis Alpha One-electron Ints"
#define PSIF_MO_B_OEI       "MO-basis Beta One-electron Ints"
#define PSIF_MO_FZC         "MO-basis Frozen-Core Operator"
#define PSIF_MO_FOCK        "MO-basis Fock Matrix"
#define PSIF_MO_A_FZC       "MO-basis Alpha Frozen-Core Oper"
#define PSIF_MO_B_FZC       "MO-basis Beta Frozen-Core Oper"
#define PSIF_MO_A_FOCK      "MO-basis Alpha Fock Matrix"
#define PSIF_MO_B_FOCK      "MO-basis Beta Fock Matrix"

/* More macros */
#define PSIF_AO_OPDM_TRIANG "AO-basis OPDM triang"
#define PSIF_AO_LAG_TRIANG  "AO-basis Lagrangian triang"
#define PSIF_AO_OPDM_SQUARE "AO-basis OPDM square"
#define PSIF_SO_OPDM        "SO-basis OPDM"
#define PSIF_SO_OPDM_TRIANG "SO-basis triang"

/* PSI return codes --- for new PSI driver           */
#define PSI_RETURN_SUCCESS      0
#define PSI_RETURN_FAILURE      1
#define PSI_RETURN_ENDLOOP      2
#define PSI_RETURN_BALK         3

//Added by ACS (01/06) for the UMP2R12 routines
#define PSIF_MO_A_MX        "MO-basis Alpha Mu-X Ints"
#define PSIF_MO_A_MY        "MO-basis Alpha Mu-Y Ints"
#define PSIF_MO_A_MZ        "MO-basis Alpha Mu-Z Ints"
#define PSIF_MO_B_MX        "MO-basis Beta Mu-X Ints"
#define PSIF_MO_B_MY        "MO-basis Beta Mu-Y Ints"
#define PSIF_MO_B_MZ        "MO-basis Beta Mu-Z Ints"
#define PSIF_MO_A_QXX       "MO-basis Alpha Q-XX Ints"
#define PSIF_MO_A_QYY       "MO-basis Alpha Q-YY Ints"
#define PSIF_MO_A_QZZ       "MO-basis Alpha Q-ZZ Ints"
#define PSIF_MO_B_QXX       "MO-basis Beta Q-XX Ints"
#define PSIF_MO_B_QYY       "MO-basis Beta Q-YY Ints"
#define PSIF_MO_B_QZZ       "MO-basis Beta Q-ZZ Ints"
#define PSIF_AO_QRR         "AO-basis Q-XX + Q-YY + Q-ZZ Ints"
#define PSIF_MO_QRR         "MO-basis Q-XX + Q-YY + Q-ZZ Ints"
#define PSIF_MO_A_QRR       "MO-basis Alpha Q-XX + Q-YY + Q-ZZ Ints"
#define PSIF_MO_B_QRR       "MO-basis Beta Q-XX + Q-YY + Q-ZZ Ints"
// end ACS additions

#endif /* header guard */
