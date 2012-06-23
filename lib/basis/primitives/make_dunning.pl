#!/usr/bin/perl

$emslmanipulate = "./emsl_manipulate.pl";
@HORDINAL = ("D", "T", "Q", "5", "6");

$excl = "exclusive";
$hhe  = "insert_hhe";
$alar = "insert_alar";

open(FAM_OUT,">basislistdunning.py");

print FAM_OUT "\"\"\"Module (auto-generated from make_dunning.pl script)\n";
print FAM_OUT "with commands building :py:class:`~basislist.BasisFamily` objects that\n";
print FAM_OUT "encode the Dunning basis set orbital definitions in\n";
print FAM_OUT ":source:`lib/basis/NOTES` and fitting bases designed for those\n";
print FAM_OUT "orbital bases.\n\n\"\"\"\n";

print FAM_OUT "from basislist import *\n\n\n";
print FAM_OUT "def load_basfam_dunning():\n\n";

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

   # append alkali to main group basis sets since EMSL stores F12 in two separate downloads
   $frag_mgp_f12       =  "basis-maingroup-cc-pv"   . $ord . "z-f12.gbs";
   $frag_alk_f12       =     "basis-alkali-cc-pv"   . $ord . "z-f12.gbs";
   $frag_mgp_f12_optri =  "basis-maingroup-cc-pv"   . $ord . "z-f12-optri.gbs";
   $frag_alk_f12_optri =     "basis-alkali-cc-pv"   . $ord . "z-f12-optri.gbs";

   $basis_f12 = "basis-cc-pv" . $ord . "z-f12_autogen.gbs";
   mergegbs($basis_f12, $frag_mgp_f12, $frag_alk_f12, $excl, 0);

   $basis_f12_optri = "basis-cc-pv" . $ord . "z-f12-optri_autogen.gbs";
   mergegbs($basis_f12_optri, $frag_mgp_f12_optri, $frag_alk_f12_optri, $excl, 0);

   # append blank H-He to core component files from EMSL so that all defined elements are present
   $frag_core        =        "corevalence-cc-pcv"  . $ord . "z.gbs";
   $frag_w_core      =              "tight-cc-pwcv" . $ord . "z.gbs";
   $frag_w_core_ri   =              "tight-cc-pwcv" . $ord . "z-ri.gbs";
   $frag_hhe         = "basis-blankHHe.gbs";

   $core = "corevalence-cc-pcv" . $ord . "z_autogen.gbs";
   mergegbs($core, $frag_core, $frag_hhe, $hhe, 0);

   $w_core = "tight-cc-pwcv" . $ord . "z_autogen.gbs";
   mergegbs($w_core, $frag_w_core, $frag_hhe, $hhe, 0);

   $w_core_ri = "tight-cc-pwcv" . $ord . "z-ri_autogen.gbs";
   mergegbs($w_core_ri, $frag_w_core_ri, $frag_hhe, $hhe, 0);

   # append plain basis and diffuse files to xpd, for which only Al-Ar are defined
   $frag_basis_xpd   =                   "basis-cc-pv_"  . $ord . "pd_z.gbs";
   $frag_diffuse_xpd =             "diffuse-aug-cc-pv_"  . $ord . "pd_z.gbs";
   $frag_diffuse_jun_xpd =    "hold-diffuse-jun-cc-pv_"  . $ord . "pd_z.gbs";
   $frag_diffuse_may_xpd =    "hold-diffuse-may-cc-pv_"  . $ord . "pd_z.gbs";
   $frag_diffuse_apr_xpd =    "hold-diffuse-apr-cc-pv_"  . $ord . "pd_z.gbs";
   $frag_diffuse_mar_xpd =    "hold-diffuse-mar-cc-pv_"  . $ord . "pd_z.gbs";
   $frag_diffuse_feb_xpd =    "hold-diffuse-feb-cc-pv_"  . $ord . "pd_z.gbs";

   $basis_xpd = "basis-cc-pv_" . $ord . "pd_z_autogen.gbs";
   mergegbs($basis_xpd, $basis, $frag_basis_xpd, $alar, 0);

   $diffuse_xpd = "diffuse-aug-cc-pv_" . $ord . "pd_z_autogen.gbs";
   mergegbs($diffuse_xpd, $diffuse, $frag_diffuse_xpd, $alar, 0);

   $diffuse_jun_xpd = "diffuse-jun-cc-pv_" . $ord . "pd_z_autogen.gbs";
   mergegbs($diffuse_jun_xpd, $diffuse_jun, $frag_diffuse_jun_xpd, $alar, 0);

   $diffuse_may_xpd = "diffuse-may-cc-pv_" . $ord . "pd_z_autogen.gbs";
   mergegbs($diffuse_may_xpd, $diffuse_may, $frag_diffuse_may_xpd, $alar, 0);

   $diffuse_apr_xpd = "diffuse-apr-cc-pv_" . $ord . "pd_z_autogen.gbs";
   mergegbs($diffuse_apr_xpd, $diffuse_apr, $frag_diffuse_apr_xpd, $alar, 0);

   $diffuse_mar_xpd = "diffuse-mar-cc-pv_" . $ord . "pd_z_autogen.gbs";
   mergegbs($diffuse_mar_xpd, $diffuse_mar, $frag_diffuse_mar_xpd, $alar, 0);

   $diffuse_feb_xpd = "diffuse-feb-cc-pv_" . $ord . "pd_z_autogen.gbs";
   mergegbs($diffuse_feb_xpd, $diffuse_feb, $frag_diffuse_feb_xpd, $alar, 0);


   # start forming final gbs files, numbering according to chart in NOTES

   # ordinary basis sets

   $gbs01 = "cc-pV" . $ord . "Z.gbs";
   copygbs($gbs01, $basis, 1);

   $gbs02 = "cc-pV(" . $ord . "+d)Z.gbs";
   copygbs($gbs02, $basis_xpd, 1);

   $gbs03 = "cc-pCV" . $ord . "Z.gbs";
   mergegbs($gbs03, $gbs01, $core, $excl, 1);

   $gbs04 = "cc-pCV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs04, $gbs02, $core, $excl, 1);

   $gbs05 = "cc-pwCV" . $ord . "Z.gbs";
   mergegbs($gbs05, $gbs01, $w_core, $excl, 1);

   $gbs06 = "cc-pwCV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs06, $gbs02, $w_core, $excl, 1);

   $gbs07 = "aug-cc-pV" . $ord . "Z.gbs";
   mergegbs($gbs07, $gbs01, $diffuse, $excl, 1);

   $gbs08 = "aug-cc-pV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs08, $gbs02, $diffuse_xpd, $excl, 1);

   $gbs09 = "aug-cc-pCV" . $ord . "Z.gbs";
   mergegbs($gbs09, $gbs03, $diffuse, $excl, 1);
  
   $gbs10 = "aug-cc-pCV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs10, $gbs04, $diffuse_xpd, $excl, 1);
  
   $gbs11 = "aug-cc-pwCV" . $ord . "Z.gbs";
   mergegbs($gbs11, $gbs05, $diffuse, $excl, 1);

   $gbs12 = "aug-cc-pwCV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs12, $gbs06, $diffuse_xpd, $excl, 1);

   $gbs73 = "heavy-aug-cc-pV" . $ord . "Z.gbs";
   mergegbs($gbs73, $gbs07, $gbs01, $hhe, 1);

   $gbs74 = "heavy-aug-cc-pV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs74, $gbs08, $gbs01, $hhe, 1);

   $gbs75 = "heavy-aug-cc-pCV" . $ord . "Z.gbs";
   mergegbs($gbs75, $gbs09, $gbs01, $hhe, 1);

   $gbs76 = "heavy-aug-cc-pCV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs76, $gbs10, $gbs01, $hhe, 1);

   $gbs77 = "heavy-aug-cc-pwCV" . $ord . "Z.gbs";
   mergegbs($gbs77, $gbs11, $gbs01, $hhe, 1);

   $gbs78 = "heavy-aug-cc-pwCV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs78, $gbs12, $gbs01, $hhe, 1);

   $gbs121 = "jun-cc-pV" . $ord . "Z.gbs";
   mergegbs($gbs121, $gbs01, $diffuse_jun, $excl, 1);

   $gbs122 = "jun-cc-pV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs122, $gbs02, $diffuse_jun_xpd, $excl, 1);

   $gbs123 = "jun-cc-pCV" . $ord . "Z.gbs";
   mergegbs($gbs123, $gbs03, $diffuse_jun, $excl, 1);
  
   $gbs124 = "jun-cc-pCV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs124, $gbs04, $diffuse_jun_xpd, $excl, 1);
  
   $gbs125 = "jun-cc-pwCV" . $ord . "Z.gbs";
   mergegbs($gbs125, $gbs05, $diffuse_jun, $excl, 1);

   $gbs126 = "jun-cc-pwCV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs126, $gbs06, $diffuse_jun_xpd, $excl, 1);

   $gbs127 = "may-cc-pV" . $ord . "Z.gbs";
   mergegbs($gbs127, $gbs01, $diffuse_may, $excl, 1);

   $gbs128 = "may-cc-pV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs128, $gbs02, $diffuse_may_xpd, $excl, 1);

   $gbs129 = "may-cc-pCV" . $ord . "Z.gbs";
   mergegbs($gbs129, $gbs03, $diffuse_may, $excl, 1);
  
   $gbs130 = "may-cc-pCV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs130, $gbs04, $diffuse_may_xpd, $excl, 1);
  
   $gbs131 = "may-cc-pwCV" . $ord . "Z.gbs";
   mergegbs($gbs131, $gbs05, $diffuse_may, $excl, 1);

   $gbs132 = "may-cc-pwCV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs132, $gbs06, $diffuse_may_xpd, $excl, 1);

   $gbs133 = "apr-cc-pV" . $ord . "Z.gbs";
   mergegbs($gbs133, $gbs01, $diffuse_apr, $excl, 1);

   $gbs134 = "apr-cc-pV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs134, $gbs02, $diffuse_apr_xpd, $excl, 1);

   $gbs135 = "apr-cc-pCV" . $ord . "Z.gbs";
   mergegbs($gbs135, $gbs03, $diffuse_apr, $excl, 1);
  
   $gbs136 = "apr-cc-pCV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs136, $gbs04, $diffuse_apr_xpd, $excl, 1);
  
   $gbs137 = "apr-cc-pwCV" . $ord . "Z.gbs";
   mergegbs($gbs137, $gbs05, $diffuse_apr, $excl, 1);

   $gbs138 = "apr-cc-pwCV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs138, $gbs06, $diffuse_apr_xpd, $excl, 1);

   $gbs139 = "mar-cc-pV" . $ord . "Z.gbs";
   mergegbs($gbs139, $gbs01, $diffuse_mar, $excl, 1);

   $gbs140 = "mar-cc-pV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs140, $gbs02, $diffuse_mar_xpd, $excl, 1);

   $gbs141 = "mar-cc-pCV" . $ord . "Z.gbs";
   mergegbs($gbs141, $gbs03, $diffuse_mar, $excl, 1);
  
   $gbs142 = "mar-cc-pCV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs142, $gbs04, $diffuse_mar_xpd, $excl, 1);
  
   $gbs143 = "mar-cc-pwCV" . $ord . "Z.gbs";
   mergegbs($gbs143, $gbs05, $diffuse_mar, $excl, 1);

   $gbs144 = "mar-cc-pwCV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs144, $gbs06, $diffuse_mar_xpd, $excl, 1);

   $gbs145 = "feb-cc-pV" . $ord . "Z.gbs";
   mergegbs($gbs145, $gbs01, $diffuse_feb, $excl, 1);

   $gbs146 = "feb-cc-pV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs146, $gbs02, $diffuse_feb_xpd, $excl, 1);

   $gbs147 = "feb-cc-pCV" . $ord . "Z.gbs";
   mergegbs($gbs147, $gbs03, $diffuse_feb, $excl, 1);
  
   $gbs148 = "feb-cc-pCV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs148, $gbs04, $diffuse_feb_xpd, $excl, 1);
  
   $gbs149 = "feb-cc-pwCV" . $ord . "Z.gbs";
   mergegbs($gbs149, $gbs05, $diffuse_feb, $excl, 1);

   $gbs150 = "feb-cc-pwCV(" . $ord . "+d)Z.gbs";
   mergegbs($gbs150, $gbs06, $diffuse_feb_xpd, $excl, 1);

   $gbs13 = "d-aug-cc-pV" . $ord . "Z.gbs";
   mergegbs($gbs13, $gbs07, $d_diffuse, $excl, 1);

   $gbs15 = "d-aug-cc-pCV" . $ord . "Z.gbs";
   mergegbs($gbs15, $gbs09, $d_diffuse, $excl, 1);

   $gbs17 = "d-aug-cc-pwCV" . $ord . "Z.gbs";
   mergegbs($gbs17, $gbs11, $d_diffuse, $excl, 1);



   # mp2-fitting basis sets

   $gbs19 = "cc-pV" . $ord . "Z-RI.gbs";
   copygbs($gbs19, $basis_ri, 1);
   
   $gbs20 = "cc-pV(" . $ord . "+d)Z-RI.gbs";
   copygbs($gbs20, $gbs19, 2);

   $gbs23 = "cc-pwCV" . $ord . "Z-RI.gbs";
   mergegbs($gbs23, $gbs19, $w_core_ri, $excl, 1);

   $gbs24 = "cc-pwCV(" . $ord . "+d)Z-RI.gbs";
   copygbs($gbs24, $gbs23, 2);

   $gbs25 = "aug-cc-pV" . $ord . "Z-RI.gbs";
   mergegbs($gbs25, $gbs19, $diffuse_ri, $excl, 1);

   $gbs26 = "aug-cc-pV(" . $ord . "+d)Z-RI.gbs";
   copygbs($gbs26, $gbs25, 2);

   $gbs29 = "aug-cc-pwCV" . $ord . "Z-RI.gbs";
   mergegbs($gbs29, $gbs23, $diffuse_ri, $excl, 1);

   $gbs30 = "aug-cc-pwCV(" . $ord . "+d)Z-RI.gbs";
   copygbs($gbs30, $gbs29, 2);

   $gbs79 = "heavy-aug-cc-pV" . $ord . "Z-RI.gbs";
   mergegbs($gbs79, $gbs25, $gbs19, $hhe, 1);

   $gbs80 = "heavy-aug-cc-pV(" . $ord . "+d)Z-RI.gbs";
   copygbs($gbs80, $gbs79, 2);

   $gbs83 = "heavy-aug-cc-pwCV" . $ord . "Z-RI.gbs";
   mergegbs($gbs83, $gbs29, $gbs19, $hhe, 1);

   $gbs84 = "heavy-aug-cc-pwCV(" . $ord . "+d)Z-RI.gbs";
   copygbs($gbs84, $gbs83, 2);

   $gbs151 = "jun-cc-pV" . $ord . "Z-RI.gbs";
   mergegbs($gbs151, $gbs19, $diffuse_jun_ri, $excl, 1);

   $gbs152 = "jun-cc-pV(" . $ord . "+d)Z-RI.gbs";
   copygbs($gbs152, $gbs151, 2);

   $gbs155 = "jun-cc-pwCV" . $ord . "Z-RI.gbs";
   mergegbs($gbs155, $gbs23, $diffuse_jun_ri, $excl, 1);

   $gbs156 = "jun-cc-pwCV(" . $ord . "+d)Z-RI.gbs";
   copygbs($gbs156, $gbs155, 2);

   $gbs157 = "may-cc-pV" . $ord . "Z-RI.gbs";
   mergegbs($gbs157, $gbs19, $diffuse_may_ri, $excl, 1);

   $gbs158 = "may-cc-pV(" . $ord . "+d)Z-RI.gbs";
   copygbs($gbs158, $gbs157, 2);

   $gbs161 = "may-cc-pwCV" . $ord . "Z-RI.gbs";
   mergegbs($gbs161, $gbs23, $diffuse_may_ri, $excl, 1);

   $gbs162 = "may-cc-pwCV(" . $ord . "+d)Z-RI.gbs";
   copygbs($gbs162, $gbs161, 2);

   $gbs163 = "apr-cc-pV" . $ord . "Z-RI.gbs";
   mergegbs($gbs163, $gbs19, $diffuse_apr_ri, $excl, 1);

   $gbs164 = "apr-cc-pV(" . $ord . "+d)Z-RI.gbs";
   copygbs($gbs164, $gbs163, 2);

   $gbs167 = "apr-cc-pwCV" . $ord . "Z-RI.gbs";
   mergegbs($gbs167, $gbs23, $diffuse_apr_ri, $excl, 1);

   $gbs168 = "apr-cc-pwCV(" . $ord . "+d)Z-RI.gbs";
   copygbs($gbs168, $gbs167, 2);

   $gbs169 = "mar-cc-pV" . $ord . "Z-RI.gbs";
   mergegbs($gbs169, $gbs19, $diffuse_mar_ri, $excl, 1);

   $gbs170 = "mar-cc-pV(" . $ord . "+d)Z-RI.gbs";
   copygbs($gbs170, $gbs169, 2);

   $gbs173 = "mar-cc-pwCV" . $ord . "Z-RI.gbs";
   mergegbs($gbs173, $gbs23, $diffuse_mar_ri, $excl, 1);

   $gbs174 = "mar-cc-pwCV(" . $ord . "+d)Z-RI.gbs";
   copygbs($gbs174, $gbs173, 2);

   $gbs175 = "feb-cc-pV" . $ord . "Z-RI.gbs";
   mergegbs($gbs175, $gbs19, $diffuse_feb_ri, $excl, 1);

   $gbs176 = "feb-cc-pV(" . $ord . "+d)Z-RI.gbs";
   copygbs($gbs176, $gbs175, 2);

   $gbs179 = "feb-cc-pwCV" . $ord . "Z-RI.gbs";
   mergegbs($gbs179, $gbs23, $diffuse_feb_ri, $excl, 1);

   $gbs180 = "feb-cc-pwCV(" . $ord . "+d)Z-RI.gbs";
   copygbs($gbs180, $gbs179, 2);



   # hf-fitting basis sets

   $gbs37 = "cc-pV" . $ord . "Z-JKFIT.gbs";
   copygbs($gbs37, $basis_jk, 1);

   $gbs38 = "cc-pV(" . $ord . "+d)Z-JKFIT.gbs";
   copygbs($gbs38, $gbs37, 2);

   $gbs43 = "aug-cc-pV" . $ord . "Z-JKFIT.gbs";
   mergegbs($gbs43, $gbs37, $diffuse_jk, $excl, 1);

   $gbs44 = "aug-cc-pV(" . $ord . "+d)Z-JKFIT.gbs";
   copygbs($gbs44, $gbs43, 2);

   $gbs85 = "heavy-aug-cc-pV" . $ord . "Z-JKFIT.gbs";
   mergegbs($gbs85, $gbs43, $gbs37, $hhe, 1);

   $gbs86 = "heavy-aug-cc-pV(" . $ord . "+d)Z-JKFIT.gbs";
   copygbs($gbs86, $gbs85, 2);

   $gbs181 = "jun-cc-pV" . $ord . "Z-JKFIT.gbs";
   mergegbs($gbs181, $gbs37, $diffuse_jun_jk, $excl, 1);

   $gbs182 = "jun-cc-pV(" . $ord . "+d)Z-JKFIT.gbs";
   copygbs($gbs182, $gbs181, 2);

   $gbs187 = "may-cc-pV" . $ord . "Z-JKFIT.gbs";
   mergegbs($gbs187, $gbs37, $diffuse_may_jk, $excl, 1);

   $gbs188 = "may-cc-pV(" . $ord . "+d)Z-JKFIT.gbs";
   copygbs($gbs188, $gbs187, 2);

   $gbs193 = "apr-cc-pV" . $ord . "Z-JKFIT.gbs";
   mergegbs($gbs193, $gbs37, $diffuse_apr_jk, $excl, 1);

   $gbs194 = "apr-cc-pV(" . $ord . "+d)Z-JKFIT.gbs";
   copygbs($gbs194, $gbs193, 2);

   $gbs199 = "mar-cc-pV" . $ord . "Z-JKFIT.gbs";
   mergegbs($gbs199, $gbs37, $diffuse_mar_jk, $excl, 1);

   $gbs200 = "mar-cc-pV(" . $ord . "+d)Z-JKFIT.gbs";
   copygbs($gbs200, $gbs199, 2);

   $gbs205 = "feb-cc-pV" . $ord . "Z-JKFIT.gbs";
   mergegbs($gbs205, $gbs37, $diffuse_feb_jk, $excl, 1);

   $gbs206 = "feb-cc-pV(" . $ord . "+d)Z-JKFIT.gbs";
   copygbs($gbs206, $gbs205, 2);



   # Douglas-Kroll basis sets

   $gbs97  = "cc-pV" . $ord . "Z-DK.gbs";
   copygbs($gbs97, $basis_dk, 1);

   $gbs99  = "cc-pCV" . $ord . "Z-DK.gbs";
   mergegbs($gbs99, $gbs97, $core, $excl, 1);

   $gbs101 = "cc-pwCV" . $ord . "Z-DK.gbs";
   copygbs($gbs101, $basis_w_core_dk, 1);

   $gbs103 = "aug-cc-pV" . $ord . "Z-DK.gbs";
   mergegbs($gbs103, $gbs97, $diffuse_dk, $excl, 1);

   $gbs105 = "aug-cc-pCV" . $ord . "Z-DK.gbs";
   mergegbs($gbs105, $gbs99, $diffuse_dk, $excl, 1);

   $gbs107 = "aug-cc-pwCV" . $ord . "Z-DK.gbs";
   mergegbs($gbs107, $gbs101, $diffuse_dk, $excl, 1);

   $gbs109 = "heavy-aug-cc-pV" . $ord . "Z-DK.gbs";
   mergegbs($gbs109, $gbs103, $gbs97, $hhe, 1);

   $gbs111 = "heavy-aug-cc-pCV" . $ord . "Z-DK.gbs";
   mergegbs($gbs111, $gbs105, $gbs97, $hhe, 1);

   $gbs113 = "heavy-aug-cc-pwCV" . $ord . "Z-DK.gbs";
   mergegbs($gbs113, $gbs107, $gbs97, $hhe, 1);

}


   # Misc. basis sets that need hand editing from others

#   copygbs("aug-cc-pvdzp.gbs",           "hold-aug-cc-pvdzp.gbs");
#   copygbs("aug-cc-pvdzp-jkfit.gbs",     "hold-aug-cc-pvdzp-jkfit.gbs");
#   copygbs("aug-cc-pvdzp-ri.gbs",        "hold-aug-cc-pvdzp-ri.gbs");

   copygbs("cc-pvtz-dual.gbs",           "hold-cc-pvtz-dual.gbs", 1);
   copygbs("cc-pvqz-dual.gbs",           "hold-cc-pvqz-dual.gbs", 1);
   copygbs("aug-cc-pvdz-dual.gbs",       "hold-aug-cc-pvdz-dual.gbs", 1);
   copygbs("aug-cc-pvtz-dual.gbs",       "hold-aug-cc-pvtz-dual.gbs", 1);
   copygbs("aug-cc-pvqz-dual.gbs",       "hold-aug-cc-pvqz-dual.gbs", 1);
   copygbs("heavy-aug-cc-pvtz-dual.gbs", "hold-heavy-aug-cc-pvtz-dual.gbs", 1);
   copygbs("heavy-aug-cc-pvqz-dual.gbs", "hold-heavy-aug-cc-pvqz-dual.gbs", 1);

close(FAM_OUT);



sub mergegbs {
   
   my($dest, $base, $append, $mode, $pywrite) = @_;

   $plaindest = lc($dest);
   $plaindest =~ s/\(/_/g;
   $plaindest =~ s/\)/_/g;
   $plaindest =~ s/,/_/g;
   $plaindest =~ s/\+/p/g;
   $plaindest =~ s/\*/s/g;

   $plainbase = lc($base);
   $plainbase =~ s/\(/_/g;
   $plainbase =~ s/\)/_/g;
   $plainbase =~ s/,/_/g;
   $plainbase =~ s/\+/p/g;
   $plainbase =~ s/\*/s/g;

   $plainappend = lc($append);
   $plainappend =~ s/\(/_/g;
   $plainappend =~ s/\)/_/g;
   $plainappend =~ s/,/_/g;
   $plainappend =~ s/\+/p/g;
   $plainappend =~ s/\*/s/g;

   if ( (-e $plainbase) && (-e $plainappend) ) {

      system ("$emslmanipulate $plainbase $plainappend $mode quiet");

      if (-e "merged.gbs") {

         system ("mv merged.gbs $plaindest");
         printf "formed   %-38s from   %-38s and   %-38s by   $mode\n", $plaindest, $plainbase, $plainappend, $mode;
         if ($pywrite == 1) { writepython($plaindest, $dest); }
      }
      else { printf ("FAILED2  %-38s\n", $plaindest); }
   }
   else { printf ("FAILED   %-38s\n", $plaindest); }
}


sub copygbs {

   my($dest, $origin, $pywrite) = @_;

   $plaindest = $dest;
   $plaindest =~ s/\(/_/g;
   $plaindest =~ s/\)/_/g;
   $plaindest =~ s/,/_/g;
   $plaindest =~ s/\+/p/g;
   $plaindest =~ s/\*/s/g;

   $plainorigin = $origin;
   $plainorigin =~ s/\(/_/g;
   $plainorigin =~ s/\)/_/g;
   $plainorigin =~ s/,/_/g;
   $plainorigin =~ s/\+/p/g;
   $plainorigin =~ s/\*/s/g;

   if (-e $plainorigin) {

      system ("cp $plainorigin $plaindest");
      printf "copied   %-38s from   %-38s\n", $plaindest, $plainorigin;
      if ($pywrite == 1) { writepython($plaindest, $dest); }
      if ($pywrite == 2) { writepython($plaindest, $origin); }
   }
   else { printf ("FAILED   %-38s\n", $plaindest); }
}


sub writepython {

   my($dest, $origin) = @_;

   $originroot = $origin;
   $originroot =~ s/\.gbs//g;

   $destorbinstance = "basis_" . lc($dest);
   $destorbinstance =~ s/\.gbs//g;
   $destorbinstance =~ s/-//g;
   $destorbinstance =~ s/jkfit$//g;
   $destorbinstance =~ s/ri$//g;
   $destorbinstance =~ s/dual$//g;

   if ($dest =~ /jkfit\.gbs/i) {
      printf FAM_OUT "    %s.add_jkfit('%s')\n", $destorbinstance, $originroot;
   }
   elsif ($dest =~ /ri\.gbs$/i) {
      printf FAM_OUT "    %s.add_rifit('%s')\n", $destorbinstance, $originroot;
   }
   elsif ($dest =~ /dual\.gbs$/i) {
      printf FAM_OUT "    %s.add_dualfit('%s')\n", $destorbinstance, $originroot;
   }
   else {
      printf FAM_OUT "    %s = BasisFamily('%s')\n", $destorbinstance, $originroot;
      printf FAM_OUT "    basisfamily_list.append(%s)\n", $destorbinstance;
   }
}

