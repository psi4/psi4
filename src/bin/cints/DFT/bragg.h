#ifndef _psi_src_bin_cints_DFT_bragg_h
#define _psi_src_bin_cints_DFT_bragg_h

/*! \file bragg.h
    \ingroup CINTS
    \brief Enter brief description of file here 

    These are the Slater-Bragg Radii in Angstoms */
namespace psi { namespace CINTS {
  double Bragg_radii[] =
    { 1.00, /* No Value*/ 
      1.00 /*0.50*/ , /* H */
      0.5882352941,                 /* He */
  /*3.0769,*/     3.0767153861,            /* Li */
      /*2.0513*/2.0511769221,                 /* Be */
      1.5385,                 /* B */
      1.2308,                 /* C */
      1.0256,                 /* N */
      0.8791/*0.87904725401*/,         /* O */ 
      0.7692,                 /* F */
      0.6838,                 /* Ne */ 
      4.0909,                 /* Na */
      3.1579,                 /* Mg */ 
      2.5714,                 /* Al */
      2.1687,                 /* Si */
      1.8750,                 /* P */ 
      1.6514,                 /* S */
      1.4754,                 /* Cl */
      1.3333,                 /* Ar */ 
      
      /* These are not right, need to fix */
      2.20, 
      1.80, 
      1.60, 
      1.40, 
      1.35, 
      1.40, 
      1.40,  /* K-Mn */
      1.40, 
      1.35, 
      1.35, 
      1.35, 
      1.35, 
      1.30, 
      1.25,  /* Fe-Ge */
      1.15, 
      1.15, 
      1.15, /* As-Br */
      2.35, 
      2.00, 
      1.80, 
      1.55, 
      1.45, 
      1.45, 
      1.35,  /* Rb-Tc */
      1.30, 
      1.35, 
      1.40, 
      1.60, 
      1.55, 
      1.55, 
      1.45,  /* Ru-Sn */
      1.45, 
      1.40,
      1.40  /* Sb-I  */
    };
};};
#endif
