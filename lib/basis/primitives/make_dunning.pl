#!/usr/bin/perl

$emslmanipulate = "./emsl_manipulate.pl";
@HORDINAL = ("d", "t", "q", "5", "6");

$excl = "exclusive";
$hhe  = "insert_hhe";
$alar = "insert_alar";


foreach $ord (@HORDINAL) {

   # access gbs component files that can be used directly from EMSL
   $basis            =              "basis-cc-pv"   . $ord . "z.gbs";
   $basis_ri         =              "basis-cc-pv"   . $ord . "z-ri.gbs";
   $basis_jk         =       "molpro-basis-cc-pv"   . $ord . "z-jkfit.gbs";
   $basis_dk         =              "basis-cc-pv"   . $ord . "z-dk.gbs";
   $basis_w_core_dk  =            "basis-cc-pwcv"   . $ord . "z-dk.gbs";

   $diffuse          =        "diffuse-aug-cc-pv"   . $ord . "z.gbs";
   $d_diffuse        =      "diffuse-d-aug-cc-pv"   . $ord . "z.gbs";
   $diffuse_ri       =        "diffuse-aug-cc-pv"   . $ord . "z-ri.gbs";
   $diffuse_jk       = "molpro-diffuse-aug-cc-pv"   . $ord . "z-jkfit.gbs";
   $diffuse_dk       =        "diffuse-aug-cc-pv"   . $ord . "z-dk.gbs";

   $diffuse_jun      =   "hold-diffuse-jun-cc-pv"   . $ord . "z.gbs";
   $diffuse_may      =   "hold-diffuse-may-cc-pv"   . $ord . "z.gbs";
   $diffuse_apr      =   "hold-diffuse-apr-cc-pv"   . $ord . "z.gbs";
   $diffuse_mar      =   "hold-diffuse-mar-cc-pv"   . $ord . "z.gbs";
   $diffuse_feb      =   "hold-diffuse-feb-cc-pv"   . $ord . "z.gbs";

   $diffuse_jun_ri   =   "hold-diffuse-jun-cc-pv"   . $ord . "z-ri.gbs";
   $diffuse_may_ri   =   "hold-diffuse-may-cc-pv"   . $ord . "z-ri.gbs";
   $diffuse_apr_ri   =   "hold-diffuse-apr-cc-pv"   . $ord . "z-ri.gbs";
   $diffuse_mar_ri   =   "hold-diffuse-mar-cc-pv"   . $ord . "z-ri.gbs";
   $diffuse_feb_ri   =   "hold-diffuse-feb-cc-pv"   . $ord . "z-ri.gbs";

   $diffuse_jun_jk   =   "hold-diffuse-jun-cc-pv"   . $ord . "z-jkfit.gbs";
   $diffuse_may_jk   =   "hold-diffuse-may-cc-pv"   . $ord . "z-jkfit.gbs";
   $diffuse_apr_jk   =   "hold-diffuse-apr-cc-pv"   . $ord . "z-jkfit.gbs";
   $diffuse_mar_jk   =   "hold-diffuse-mar-cc-pv"   . $ord . "z-jkfit.gbs";
   $diffuse_feb_jk   =   "hold-diffuse-feb-cc-pv"   . $ord . "z-jkfit.gbs";

   # append blank H-He to core component files from EMSL so that all defined elements are present
   $frag_core        =        "corevalence-cc-pcv"  . $ord . "z.gbs";
   $frag_w_core      =              "tight-cc-pwcv" . $ord . "z.gbs";
   $frag_w_core_ri   =              "tight-cc-pwcv" . $ord . "z-ri.gbs";
   $frag_hhe         = "basis-blankHHe.gbs";

   $core = "corevalence-cc-pcv" . $ord . "z_autogen.gbs";
   mergegbs($core, $frag_core, $frag_hhe, $hhe);

   $w_core = "tight-cc-pwcv" . $ord . "z_autogen.gbs";
   mergegbs($w_core, $frag_w_core, $frag_hhe, $hhe);

   $w_core_ri = "tight-cc-pwcv" . $ord . "z-ri_autogen.gbs";
   mergegbs($w_core_ri, $frag_w_core_ri, $frag_hhe, $hhe);

   # append plain basis and diffuse files to xpd, for which only Al-Ar are defined
   $frag_basis_xpd   =                   "basis-cc-pv_"  . $ord . "pd_z.gbs";
   $frag_diffuse_xpd =             "diffuse-aug-cc-pv_"  . $ord . "pd_z.gbs";
   $frag_diffuse_jun_xpd =    "hold-diffuse-jun-cc-pv_"  . $ord . "pd_z.gbs";
   $frag_diffuse_may_xpd =    "hold-diffuse-may-cc-pv_"  . $ord . "pd_z.gbs";
   $frag_diffuse_apr_xpd =    "hold-diffuse-apr-cc-pv_"  . $ord . "pd_z.gbs";
   $frag_diffuse_mar_xpd =    "hold-diffuse-mar-cc-pv_"  . $ord . "pd_z.gbs";
   $frag_diffuse_feb_xpd =    "hold-diffuse-feb-cc-pv_"  . $ord . "pd_z.gbs";

   $basis_xpd = "basis-cc-pv_" . $ord . "pd_z_autogen.gbs";
   mergegbs($basis_xpd, $basis, $frag_basis_xpd, $alar);

   $diffuse_xpd = "diffuse-aug-cc-pv_" . $ord . "pd_z_autogen.gbs";
   mergegbs($diffuse_xpd, $diffuse, $frag_diffuse_xpd, $alar);

   $diffuse_jun_xpd = "diffuse-jun-cc-pv_" . $ord . "pd_z_autogen.gbs";
   mergegbs($diffuse_jun_xpd, $diffuse_jun, $frag_diffuse_jun_xpd, $alar);

   $diffuse_may_xpd = "diffuse-may-cc-pv_" . $ord . "pd_z_autogen.gbs";
   mergegbs($diffuse_may_xpd, $diffuse_may, $frag_diffuse_may_xpd, $alar);

   $diffuse_apr_xpd = "diffuse-apr-cc-pv_" . $ord . "pd_z_autogen.gbs";
   mergegbs($diffuse_apr_xpd, $diffuse_apr, $frag_diffuse_apr_xpd, $alar);

   $diffuse_mar_xpd = "diffuse-mar-cc-pv_" . $ord . "pd_z_autogen.gbs";
   mergegbs($diffuse_mar_xpd, $diffuse_mar, $frag_diffuse_mar_xpd, $alar);

   $diffuse_feb_xpd = "diffuse-feb-cc-pv_" . $ord . "pd_z_autogen.gbs";
   mergegbs($diffuse_feb_xpd, $diffuse_feb, $frag_diffuse_feb_xpd, $alar);


   # start forming final gbs files, numbering according to chart in NOTES

   # ordinary basis sets

   $gbs01 = "cc-pv" . $ord . "z.gbs";
   copygbs($gbs01, $basis);

   $gbs02 = "cc-pv_" . $ord . "pd_z.gbs";
   copygbs($gbs02, $basis_xpd);

   $gbs03 = "cc-pcv" . $ord . "z.gbs";
   mergegbs($gbs03, $gbs01, $core, $excl);

   $gbs04 = "cc-pcv_" . $ord . "pd_z.gbs";
   mergegbs($gbs04, $gbs02, $core, $excl);

   $gbs05 = "cc-pwcv" . $ord . "z.gbs";
   mergegbs($gbs05, $gbs01, $w_core, $excl);

   $gbs06 = "cc-pwcv_" . $ord . "pd_z.gbs";
   mergegbs($gbs06, $gbs02, $w_core, $excl);

   $gbs07 = "aug-cc-pv" . $ord . "z.gbs";
   mergegbs($gbs07, $gbs01, $diffuse, $excl);

   $gbs08 = "aug-cc-pv_" . $ord . "pd_z.gbs";
   mergegbs($gbs08, $gbs02, $diffuse_xpd, $excl);

   $gbs09 = "aug-cc-pcv" . $ord . "z.gbs";
   mergegbs($gbs09, $gbs03, $diffuse, $excl);
  
   $gbs10 = "aug-cc-pcv_" . $ord . "pd_z.gbs";
   mergegbs($gbs10, $gbs04, $diffuse_xpd, $excl);
  
   $gbs11 = "aug-cc-pwcv" . $ord . "z.gbs";
   mergegbs($gbs11, $gbs05, $diffuse, $excl);

   $gbs12 = "aug-cc-pwcv_" . $ord . "pd_z.gbs";
   mergegbs($gbs12, $gbs06, $diffuse_xpd, $excl);

   $gbs73 = "heavy-aug-cc-pv" . $ord . "z.gbs";
   mergegbs($gbs73, $gbs07, $gbs01, $hhe);

   $gbs74 = "heavy-aug-cc-pv_" . $ord . "pd_z.gbs";
   mergegbs($gbs74, $gbs08, $gbs01, $hhe);

   $gbs75 = "heavy-aug-cc-pcv" . $ord . "z.gbs";
   mergegbs($gbs75, $gbs09, $gbs01, $hhe);

   $gbs76 = "heavy-aug-cc-pcv_" . $ord . "pd_z.gbs";
   mergegbs($gbs76, $gbs10, $gbs01, $hhe);

   $gbs77 = "heavy-aug-cc-pwcv" . $ord . "z.gbs";
   mergegbs($gbs77, $gbs11, $gbs01, $hhe);

   $gbs78 = "heavy-aug-cc-pwcv_" . $ord . "pd_z.gbs";
   mergegbs($gbs78, $gbs12, $gbs01, $hhe);

   $gbs121 = "jun-cc-pv" . $ord . "z.gbs";
   mergegbs($gbs121, $gbs01, $diffuse_jun, $excl);

   $gbs122 = "jun-cc-pv_" . $ord . "pd_z.gbs";
   mergegbs($gbs122, $gbs02, $diffuse_jun_xpd, $excl);

   $gbs123 = "jun-cc-pcv" . $ord . "z.gbs";
   mergegbs($gbs123, $gbs03, $diffuse_jun, $excl);
  
   $gbs124 = "jun-cc-pcv_" . $ord . "pd_z.gbs";
   mergegbs($gbs124, $gbs04, $diffuse_jun_xpd, $excl);
  
   $gbs125 = "jun-cc-pwcv" . $ord . "z.gbs";
   mergegbs($gbs125, $gbs05, $diffuse_jun, $excl);

   $gbs126 = "jun-cc-pwcv_" . $ord . "pd_z.gbs";
   mergegbs($gbs126, $gbs06, $diffuse_jun_xpd, $excl);

   $gbs127 = "may-cc-pv" . $ord . "z.gbs";
   mergegbs($gbs127, $gbs01, $diffuse_may, $excl);

   $gbs128 = "may-cc-pv_" . $ord . "pd_z.gbs";
   mergegbs($gbs128, $gbs02, $diffuse_may_xpd, $excl);

   $gbs129 = "may-cc-pcv" . $ord . "z.gbs";
   mergegbs($gbs129, $gbs03, $diffuse_may, $excl);
  
   $gbs130 = "may-cc-pcv_" . $ord . "pd_z.gbs";
   mergegbs($gbs130, $gbs04, $diffuse_may_xpd, $excl);
  
   $gbs131 = "may-cc-pwcv" . $ord . "z.gbs";
   mergegbs($gbs131, $gbs05, $diffuse_may, $excl);

   $gbs132 = "may-cc-pwcv_" . $ord . "pd_z.gbs";
   mergegbs($gbs132, $gbs06, $diffuse_may_xpd, $excl);

   $gbs133 = "apr-cc-pv" . $ord . "z.gbs";
   mergegbs($gbs133, $gbs01, $diffuse_apr, $excl);

   $gbs134 = "apr-cc-pv_" . $ord . "pd_z.gbs";
   mergegbs($gbs134, $gbs02, $diffuse_apr_xpd, $excl);

   $gbs135 = "apr-cc-pcv" . $ord . "z.gbs";
   mergegbs($gbs135, $gbs03, $diffuse_apr, $excl);
  
   $gbs136 = "apr-cc-pcv_" . $ord . "pd_z.gbs";
   mergegbs($gbs136, $gbs04, $diffuse_apr_xpd, $excl);
  
   $gbs137 = "apr-cc-pwcv" . $ord . "z.gbs";
   mergegbs($gbs137, $gbs05, $diffuse_apr, $excl);

   $gbs138 = "apr-cc-pwcv_" . $ord . "pd_z.gbs";
   mergegbs($gbs138, $gbs06, $diffuse_apr_xpd, $excl);

   $gbs139 = "mar-cc-pv" . $ord . "z.gbs";
   mergegbs($gbs139, $gbs01, $diffuse_mar, $excl);

   $gbs140 = "mar-cc-pv_" . $ord . "pd_z.gbs";
   mergegbs($gbs140, $gbs02, $diffuse_mar_xpd, $excl);

   $gbs141 = "mar-cc-pcv" . $ord . "z.gbs";
   mergegbs($gbs141, $gbs03, $diffuse_mar, $excl);
  
   $gbs142 = "mar-cc-pcv_" . $ord . "pd_z.gbs";
   mergegbs($gbs142, $gbs04, $diffuse_mar_xpd, $excl);
  
   $gbs143 = "mar-cc-pwcv" . $ord . "z.gbs";
   mergegbs($gbs143, $gbs05, $diffuse_mar, $excl);

   $gbs144 = "mar-cc-pwcv_" . $ord . "pd_z.gbs";
   mergegbs($gbs144, $gbs06, $diffuse_mar_xpd, $excl);

   $gbs145 = "feb-cc-pv" . $ord . "z.gbs";
   mergegbs($gbs145, $gbs01, $diffuse_feb, $excl);

   $gbs146 = "feb-cc-pv_" . $ord . "pd_z.gbs";
   mergegbs($gbs146, $gbs02, $diffuse_feb_xpd, $excl);

   $gbs147 = "feb-cc-pcv" . $ord . "z.gbs";
   mergegbs($gbs147, $gbs03, $diffuse_feb, $excl);
  
   $gbs148 = "feb-cc-pcv_" . $ord . "pd_z.gbs";
   mergegbs($gbs148, $gbs04, $diffuse_feb_xpd, $excl);
  
   $gbs149 = "feb-cc-pwcv" . $ord . "z.gbs";
   mergegbs($gbs149, $gbs05, $diffuse_feb, $excl);

   $gbs150 = "feb-cc-pwcv_" . $ord . "pd_z.gbs";
   mergegbs($gbs150, $gbs06, $diffuse_feb_xpd, $excl);

   $gbs13 = "d-aug-cc-pv" . $ord . "z.gbs";
   mergegbs($gbs13, $gbs07, $d_diffuse, $excl);

   $gbs15 = "d-aug-cc-pcv" . $ord . "z.gbs";
   mergegbs($gbs15, $gbs09, $d_diffuse, $excl);

   $gbs17 = "d-aug-cc-pwcv" . $ord . "z.gbs";
   mergegbs($gbs17, $gbs11, $d_diffuse, $excl);



   # mp2-fitting basis sets

   $gbs19 = "cc-pv" . $ord . "z-ri.gbs";
   copygbs($gbs19, $basis_ri);
   
   $gbs20 = "cc-pv_" . $ord . "pd_z-ri.gbs";
   copygbs($gbs20, $gbs19);

   $gbs23 = "cc-pwcv" . $ord . "z-ri.gbs";
   mergegbs($gbs23, $gbs19, $w_core_ri, $excl);

   $gbs24 = "cc-pwcv_" . $ord . "pd_z-ri.gbs";
   copygbs($gbs24, $gbs23);

   $gbs25 = "aug-cc-pv" . $ord . "z-ri.gbs";
   mergegbs($gbs25, $gbs19, $diffuse_ri, $excl);

   $gbs26 = "aug-cc-pv_" . $ord . "pd_z-ri.gbs";
   copygbs($gbs26, $gbs25);

   $gbs29 = "aug-cc-pwcv" . $ord . "z-ri.gbs";
   mergegbs($gbs29, $gbs23, $diffuse_ri, $excl);

   $gbs30 = "aug-cc-pwcv_" . $ord . "pd_z-ri.gbs";
   copygbs($gbs30, $gbs29);

   $gbs79 = "heavy-aug-cc-pv" . $ord . "z-ri.gbs";
   mergegbs($gbs79, $gbs25, $gbs19, $hhe);

   $gbs80 = "heavy-aug-cc-pv_" . $ord . "pd_z-ri.gbs";
   copygbs($gbs80, $gbs79);

   $gbs83 = "heavy-aug-cc-pwcv" . $ord . "z-ri.gbs";
   mergegbs($gbs83, $gbs29, $gbs19, $hhe);

   $gbs84 = "heavy-aug-cc-pwcv_" . $ord . "pd_z-ri.gbs";
   copygbs($gbs84, $gbs83);

   $gbs151 = "jun-cc-pv" . $ord . "z-ri.gbs";
   mergegbs($gbs151, $gbs19, $diffuse_jun_ri, $excl);

   $gbs152 = "jun-cc-pv_" . $ord . "pd_z-ri.gbs";
   copygbs($gbs152, $gbs151);

   $gbs155 = "jun-cc-pwcv" . $ord . "z-ri.gbs";
   mergegbs($gbs155, $gbs23, $diffuse_jun_ri, $excl);

   $gbs156 = "jun-cc-pwcv_" . $ord . "pd_z-ri.gbs";
   copygbs($gbs156, $gbs155);

   $gbs157 = "may-cc-pv" . $ord . "z-ri.gbs";
   mergegbs($gbs157, $gbs19, $diffuse_may_ri, $excl);

   $gbs158 = "may-cc-pv_" . $ord . "pd_z-ri.gbs";
   copygbs($gbs158, $gbs157);

   $gbs161 = "may-cc-pwcv" . $ord . "z-ri.gbs";
   mergegbs($gbs161, $gbs23, $diffuse_may_ri, $excl);

   $gbs162 = "may-cc-pwcv_" . $ord . "pd_z-ri.gbs";
   copygbs($gbs162, $gbs161);

   $gbs163 = "apr-cc-pv" . $ord . "z-ri.gbs";
   mergegbs($gbs163, $gbs19, $diffuse_apr_ri, $excl);

   $gbs164 = "apr-cc-pv_" . $ord . "pd_z-ri.gbs";
   copygbs($gbs164, $gbs163);

   $gbs167 = "apr-cc-pwcv" . $ord . "z-ri.gbs";
   mergegbs($gbs167, $gbs23, $diffuse_apr_ri, $excl);

   $gbs168 = "apr-cc-pwcv_" . $ord . "pd_z-ri.gbs";
   copygbs($gbs168, $gbs167);

   $gbs169 = "mar-cc-pv" . $ord . "z-ri.gbs";
   mergegbs($gbs169, $gbs19, $diffuse_mar_ri, $excl);

   $gbs170 = "mar-cc-pv_" . $ord . "pd_z-ri.gbs";
   copygbs($gbs170, $gbs169);

   $gbs173 = "mar-cc-pwcv" . $ord . "z-ri.gbs";
   mergegbs($gbs173, $gbs23, $diffuse_mar_ri, $excl);

   $gbs174 = "mar-cc-pwcv_" . $ord . "pd_z-ri.gbs";
   copygbs($gbs174, $gbs173);

   $gbs175 = "feb-cc-pv" . $ord . "z-ri.gbs";
   mergegbs($gbs175, $gbs19, $diffuse_feb_ri, $excl);

   $gbs176 = "feb-cc-pv_" . $ord . "pd_z-ri.gbs";
   copygbs($gbs176, $gbs175);

   $gbs179 = "feb-cc-pwcv" . $ord . "z-ri.gbs";
   mergegbs($gbs179, $gbs23, $diffuse_feb_ri, $excl);

   $gbs180 = "feb-cc-pwcv_" . $ord . "pd_z-ri.gbs";
   copygbs($gbs180, $gbs179);



   # hf-fitting basis sets

   $gbs37 = "cc-pv" . $ord . "z-jkfit.gbs";
   copygbs($gbs37, $basis_jk);

   $gbs38 = "cc-pv_" . $ord . "pd_z-jkfit.gbs";
   copygbs($gbs38, $gbs37);

   $gbs43 = "aug-cc-pv" . $ord . "z-jkfit.gbs";
   mergegbs($gbs43, $gbs37, $diffuse_jk, $excl);

   $gbs44 = "aug-cc-pv_" . $ord . "pd_z-jkfit.gbs";
   copygbs($gbs44, $gbs43);

   $gbs85 = "heavy-aug-cc-pv" . $ord . "z-jkfit.gbs";
   mergegbs($gbs85, $gbs43, $gbs37, $hhe);

   $gbs86 = "heavy-aug-cc-pv_" . $ord . "pd_z-jkfit.gbs";
   copygbs($gbs86, $gbs85);

   $gbs181 = "jun-cc-pv" . $ord . "z-jkfit.gbs";
   mergegbs($gbs181, $gbs37, $diffuse_jun_jk, $excl);

   $gbs182 = "jun-cc-pv_" . $ord . "pd_z-jkfit.gbs";
   copygbs($gbs182, $gbs181);

   $gbs187 = "may-cc-pv" . $ord . "z-jkfit.gbs";
   mergegbs($gbs187, $gbs37, $diffuse_may_jk, $excl);

   $gbs188 = "may-cc-pv_" . $ord . "pd_z-jkfit.gbs";
   copygbs($gbs188, $gbs187);

   $gbs193 = "apr-cc-pv" . $ord . "z-jkfit.gbs";
   mergegbs($gbs193, $gbs37, $diffuse_apr_jk, $excl);

   $gbs194 = "apr-cc-pv_" . $ord . "pd_z-jkfit.gbs";
   copygbs($gbs194, $gbs193);

   $gbs199 = "mar-cc-pv" . $ord . "z-jkfit.gbs";
   mergegbs($gbs199, $gbs37, $diffuse_mar_jk, $excl);

   $gbs200 = "mar-cc-pv_" . $ord . "pd_z-jkfit.gbs";
   copygbs($gbs200, $gbs199);

   $gbs205 = "feb-cc-pv" . $ord . "z-jkfit.gbs";
   mergegbs($gbs205, $gbs37, $diffuse_feb_jk, $excl);

   $gbs206 = "feb-cc-pv_" . $ord . "pd_z-jkfit.gbs";
   copygbs($gbs206, $gbs205);



   # Douglas-Kroll basis sets

   $gbs97  = "cc-pv" . $ord . "z-dk.gbs";
   copygbs($gbs97, $basis_dk);

   $gbs99  = "cc-pcv" . $ord . "z-dk.gbs";
   mergegbs($gbs99, $gbs97, $core, $excl);

   $gbs101 = "cc-pwcv" . $ord . "z-dk.gbs";
   copygbs($gbs101, $basis_w_core_dk);

   $gbs103 = "aug-cc-pv" . $ord . "z-dk.gbs";
   mergegbs($gbs103, $gbs97, $diffuse_dk, $excl);

   $gbs105 = "aug-cc-pcv" . $ord . "z-dk.gbs";
   mergegbs($gbs105, $gbs99, $diffuse_dk, $excl);

   $gbs107 = "aug-cc-pwcv" . $ord . "z-dk.gbs";
   mergegbs($gbs107, $gbs101, $diffuse_dk, $excl);

   $gbs109 = "heavy-aug-cc-pv" . $ord . "z-dk.gbs";
   mergegbs($gbs109, $gbs103, $gbs97, $hhe);

   $gbs111 = "heavy-aug-cc-pcv" . $ord . "z-dk.gbs";
   mergegbs($gbs111, $gbs105, $gbs97, $hhe);

   $gbs113 = "heavy-aug-cc-pwcv" . $ord . "z-dk.gbs";
   mergegbs($gbs113, $gbs107, $gbs97, $hhe);

}


   # Misc. basis sets that need hand editing from others

   copygbs("aug-cc-pvdzp.gbs",           "hold-aug-cc-pvdzp.gbs");
   copygbs("aug-cc-pvdzp-jkfit.gbs",     "hold-aug-cc-pvdzp-jkfit.gbs");
   copygbs("aug-cc-pvdzp-ri.gbs",        "hold-aug-cc-pvdzp-ri.gbs");

   copygbs("cc-pvtz-dual.gbs",           "hold-cc-pvtz-dual.gbs");
   copygbs("cc-pvqz-dual.gbs",           "hold-cc-pvqz-dual.gbs");
   copygbs("aug-cc-pvdz-dual.gbs",       "hold-aug-cc-pvdz-dual.gbs");
   copygbs("aug-cc-pvtz-dual.gbs",       "hold-aug-cc-pvtz-dual.gbs");
   copygbs("aug-cc-pvqz-dual.gbs",       "hold-aug-cc-pvqz-dual.gbs");
   copygbs("heavy-aug-cc-pvtz-dual.gbs", "hold-heavy-aug-cc-pvtz-dual.gbs");
   copygbs("heavy-aug-cc-pvqz-dual.gbs", "hold-heavy-aug-cc-pvqz-dual.gbs");




sub mergegbs {
   
   my($dest, $base, $append, $mode) = @_;

   if ( (-e $base) && (-e $append) ) {

      system ("$emslmanipulate $base $append $mode quiet");

      if (-e "merged.gbs") {

         system ("mv merged.gbs $dest");
         printf ("formed   %-38s from   %-38s and   %-38s by   $mode\n", $dest, $base, $append, $mode);
      }
      else {

         printf ("FAILED2  %-38s\n", $dest);
      }
   }
   else {

      printf ("FAILED   %-38s\n", $dest);
   }
}


sub copygbs {

   my($dest, $origin) = @_;

   if (-e $origin) {

      system ("cp $origin $dest");
      printf ("copied   %-38s from   %-38s\n", $dest, $origin);
   }
   else {

      printf ("FAILED   %-38s\n", $dest);
   }
}


