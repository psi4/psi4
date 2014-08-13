#
# f77funcgen.pl
#
# find all pnga_ functions and output as many F77_FUNC_ macros as possible
# output pnga_ functions that don't quite fit
# 
if ($#ARGV+1 != 1) {
    die "usage: f77funcgen.pl filename";
}
open FILE, "<$ARGV[0]" or die $!;
while (<FILE>) {
    if (/pnga_/) {
        chomp;
        s/^.*pnga_(.*)\(.*$/\1/;
        $big = uc $_;
        print "#define ga_${_}_  F77_FUNC_(ga_$_, GA_$big)\n";
        print "#define ga_c${_}_ F77_FUNC_(ga_c$_,GA_C$big)\n";
        print "#define ga_d${_}_ F77_FUNC_(ga_d$_,GA_D$big)\n";
        print "#define ga_i${_}_ F77_FUNC_(ga_i$_,GA_I$big)\n";
        print "#define ga_s${_}_ F77_FUNC_(ga_s$_,GA_S$big)\n";
        print "#define ga_z${_}_ F77_FUNC_(ga_z$_,GA_Z$big)\n";
        print "#define nga_${_}_  F77_FUNC_(nga_$_, NGA_$big)\n";
        print "#define nga_c${_}_ F77_FUNC_(nga_c$_,NGA_C$big)\n";
        print "#define nga_d${_}_ F77_FUNC_(nga_d$_,NGA_D$big)\n";
        print "#define nga_i${_}_ F77_FUNC_(nga_i$_,NGA_I$big)\n";
        print "#define nga_s${_}_ F77_FUNC_(nga_s$_,NGA_S$big)\n";
        print "#define nga_z${_}_ F77_FUNC_(nga_z$_,NGA_Z$big)\n";
    }
}
print "/* the missing functions are either complex type or strangely named */\n";
print "\n";
print "#define gai_cdot_patch_     F77_FUNC_(gai_cdot_patch,GAI_CDOT_PATCH)\n";
print "#define gai_zdot_patch_     F77_FUNC_(gai_zdot_patch,GAI_ZDOT_PATCH)\n";
print "#define ngai_cdot_patch_    F77_FUNC_(ngai_cdot_patch,NGAI_CDOT_PATCH)\n";
print "#define ngai_zdot_patch_    F77_FUNC_(ngai_zdot_patch,NGAI_ZDOT_PATCH)\n";
print "\n";
print "#define ga_cscal_patch_     F77_FUNC_(ga_cscal_patch,GA_CSCAL_PATCH)\n";
print "#define ga_dscal_patch_     F77_FUNC_(ga_dscal_patch,GA_DSCAL_PATCH)\n";
print "#define ga_iscal_patch_     F77_FUNC_(ga_iscal_patch,GA_ISCAL_PATCH)\n";
print "#define ga_sscal_patch_     F77_FUNC_(ga_sscal_patch,GA_SSCAL_PATCH)\n";
print "#define ga_zscal_patch_     F77_FUNC_(ga_zscal_patch,GA_ZSCAL_PATCH)\n";
print "\n";
print "#define ga_cscal_           F77_FUNC_(ga_cscal,GA_CSCAL)\n";
print "#define ga_dscal_           F77_FUNC_(ga_dscal,GA_DSCAL)\n";
print "#define ga_iscal_           F77_FUNC_(ga_iscal,GA_ISCAL)\n";
print "#define ga_sscal_           F77_FUNC_(ga_sscal,GA_SSCAL)\n";
print "#define ga_zscal_           F77_FUNC_(ga_zscal,GA_ZSCAL)\n";
print "\n";
print "#define gai_cdot_           F77_FUNC_(gai_cdot,GAI_CDOT)\n";
print "#define gai_zdot_           F77_FUNC_(gai_zdot,GAI_ZDOT)\n";
print "#define ngai_cdot_          F77_FUNC_(ngai_cdot,NGAI_CDOT)\n";
print "#define ngai_zdot_          F77_FUNC_(ngai_zdot,NGAI_ZDOT)\n";
print "\n";
print "#define ga_cgemm_           F77_FUNC_(ga_cgemm,GA_CGEMM)\n";
print "#define ga_dgemm_           F77_FUNC_(ga_dgemm,GA_DGEMM)\n";
print "#define ga_sgemm_           F77_FUNC_(ga_sgemm,GA_SGEMM)\n";
print "#define ga_zgemm_           F77_FUNC_(ga_zgemm,GA_ZGEMM)\n";
print "\n";
print "#define nga_periodic_get_   F77_FUNC_(nga_periodic_get,NGA_PERIODIC_GET)\n";
print "#define nga_periodic_put_   F77_FUNC_(nga_periodic_put,NGA_PERIODIC_PUT)\n";
print "#define nga_periodic_acc_   F77_FUNC_(nga_periodic_acc,NGA_PERIODIC_ACC)\n";
print "\n";
print "#define ga_access_          F77_FUNC_(ga_access,GA_ACCESS)\n";
print "#define nga_access_         F77_FUNC_(nga_access,NGA_ACCESS)\n";
print "#define nga_access_block_   F77_FUNC_(nga_access_block,NGA_ACCESS_BLOCK)\n";
print "#define nga_access_block_grid_      F77_FUNC_(nga_access_block_grid,NGA_ACCESS_BLOCK_GRID)\n";
print "#define nga_access_block_segment_   F77_FUNC_(nga_access_block_segment,NGA_ACCESS_BLOCK_SEGMENT)\n";
print "\n";
print "#define ga_pgroup_cgop_     F77_FUNC_(ga_pgroup_cgop,GA_PGROUP_CGOP)\n";
print "#define ga_pgroup_dgop_     F77_FUNC_(ga_pgroup_dgop,GA_PGROUP_DGOP)\n";
print "#define ga_pgroup_igop_     F77_FUNC_(ga_pgroup_igop,GA_PGROUP_IGOP)\n";
print "#define ga_pgroup_sgop_     F77_FUNC_(ga_pgroup_sgop,GA_PGROUP_SGOP)\n";
print "#define ga_pgroup_zgop_     F77_FUNC_(ga_pgroup_zgop,GA_PGROUP_ZGOP)\n";
print "\n";
print "#define nga_pgroup_cgop_    F77_FUNC_(nga_pgroup_cgop,NGA_PGROUP_CGOP)\n";
print "#define nga_pgroup_dgop_    F77_FUNC_(nga_pgroup_dgop,NGA_PGROUP_DGOP)\n";
print "#define nga_pgroup_igop_    F77_FUNC_(nga_pgroup_igop,NGA_PGROUP_IGOP)\n";
print "#define nga_pgroup_sgop_    F77_FUNC_(nga_pgroup_sgop,NGA_PGROUP_SGOP)\n";
print "#define nga_pgroup_zgop_    F77_FUNC_(nga_pgroup_zgop,NGA_PGROUP_ZGOP)\n";
