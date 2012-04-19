#!/usr/bin/perl

use strict;
use warnings;
use File::Path qw(remove_tree);

# This script reads the driver's options setup to provide a list of 
# keywords expected by each module. The results of this parsing is put into TeX files, which are
# inlined into the manual.


#
# Read the options for each module
#
my $DriverPath = "";
if ($#ARGV == 0) { $DriverPath = $ARGV[0] . "/"; }
my $DriverFile = $DriverPath . "../../src/bin/psi4/read_options.cc";

my $CurrentModule;
my %Keywords;
my %ModuleDescriptions;
my %ModuleSubsections;
my %Expert;
open(DRIVER, $DriverFile) or die "\nI can't read the PSI driver file\n";

my $CurrentSubsection = "";
my @OrderedSubsection = ();

# Start looping over the driver file
while(<DRIVER>){
    # Look for the name of the module for which options are being added
    if(/name\s*\=\=\s*\"(\w+)\"/){
        $CurrentModule = $1;
        $CurrentSubsection = "";
        @OrderedSubsection = ();
    }
    my $CommentString;
    my $SphinxCommentString;
    my $Type;
    my $Default;
    my $Keyword;
    my $Possibilities;
    my $Expert;
    # If we find a /*- it means the comment block has started, but we
    # don't know if it's a multi-line comment, let's find out
    if(/\/\*-\s*SUBSECTION\s+(.*)-\*\// and $CurrentModule){
        $CurrentSubsection = $1;
        #print "$CurrentSubsection\n";
        push @{$ModuleSubsections{$CurrentModule}}, $CurrentSubsection;
    }elsif(/\/\*-\s*MODULEDESCRIPTION/ and $CurrentModule){
        $ModuleDescriptions{$CurrentModule} = get_description($_);
    }elsif(/\/\*-/ and $CurrentModule){
        ($CommentString, $Expert) = determine_comment($_);
        $CommentString =~ s/_/\\_/g;
        # process @@ as math mode subscript in tex
        $CommentString =~ s/@@/_/g;
        $SphinxCommentString = $CommentString;
        $SphinxCommentString =~ s/ \$/ :math:`/g;
        $SphinxCommentString =~ s/\(\$/(\\ :math:`/g;
        $SphinxCommentString =~ s/\$ /` /g;
        $SphinxCommentString =~ s/\$\./`./g;
        $SphinxCommentString =~ s/\$,/`,/g;
        $SphinxCommentString =~ s/\$\)/`\\ \)/g;
        $SphinxCommentString =~ s/\\_/_/g;
        ($Keyword, $Type, $Default, $Possibilities) = determine_keyword_type_and_default();
        if($Expert){
            $Expert{$CurrentModule}{$CurrentSubsection}{$Keyword}{"Type"}    = $Type;
            $Expert{$CurrentModule}{$CurrentSubsection}{$Keyword}{"Default"} = $Default;
            $Expert{$CurrentModule}{$CurrentSubsection}{$Keyword}{"Comment"} = $CommentString;
            $Expert{$CurrentModule}{$CurrentSubsection}{$Keyword}{"SphComment"} = $SphinxCommentString;
            $Expert{$CurrentModule}{$CurrentSubsection}{$Keyword}{"Possibilities"} = $Possibilities;
        }else{
            $Keywords{$CurrentModule}{$CurrentSubsection}{$Keyword}{"Type"}    = $Type;
            $Keywords{$CurrentModule}{$CurrentSubsection}{$Keyword}{"Default"} = $Default;
            $Keywords{$CurrentModule}{$CurrentSubsection}{$Keyword}{"Comment"} = $CommentString;
            $Keywords{$CurrentModule}{$CurrentSubsection}{$Keyword}{"SphComment"} = $SphinxCommentString;
            $Keywords{$CurrentModule}{$CurrentSubsection}{$Keyword}{"Possibilities"} = $Possibilities;
        }
    }
}
close DRIVER;


my @temp = ();
print_hash(\%Keywords, "keywords.tex", 1, "kw");
my @PSIMODULES = @temp;
print_hash(\%Expert, "expert_keywords.tex", 0, "ekw");

sub print_hash
{
 my %hash     = %{$_[0]};
 my $filename = $_[1];
 my $print_description = $_[2];
 my $label_string = $_[3];
 open(OUT,">$filename") or die "\nI can't write to $filename\n";
 print OUT "{\n \\footnotesize\n";
 my $ttsout = "source/autodoc_abbr_options_c.rst";
 my $tsout = "source/autodoc_glossary_options_c.rst";
 my $sout = "source/autodoc_options_c_bymodule.rst";
 if ($print_description) { 
    open(TTSOUT,">$ttsout") or die "\nI can't write to $ttsout\n";
    open(TSOUT,">$tsout") or die "\nI can't write to $tsout\n";
    print TSOUT "\n.. include:: autodoc_abbr_options_c.rst\n\n";
    print TSOUT "\n.. _`apdx:options_c_alpha`:\n\n";
    print TSOUT "Keywords by Alpha\n=================\n\n";
    print TSOUT ".. glossary::\n   :sorted:\n\n";
    open(SOUT,">$sout") or die "\nI can't write to $sout\n";
    print SOUT "\n.. _`apdx:options_c_module`:\n\n";
    print SOUT "Keywords by Module\n==================\n\n.. toctree::\n   :maxdepth: 1\n\n";
 }
 else { 
    open(TTSOUT,">>$ttsout") or die "\nI can't write to $ttsout\n";
    open(TSOUT,">>$tsout") or die "\nI can't write to $tsout\n";
    open(SOUT,">/dev/null");
 }
 my @RearrModules = sort {$a cmp $b} keys %hash;
 @RearrModules = grep { $_ ne "GLOBALS"} @RearrModules;
 unshift(@RearrModules, "GLOBALS");
 foreach my $Module (@RearrModules){
     $Module =~ s/_/-/g; # Things like plugin_module_name will screw things up...
     push(@temp, $Module);
     my $Moddivider = "=" x length($Module);
     printf OUT "\n\\subsection{%s}\\label{%s-%s}\n",$Module,$label_string, $Module;
     print SOUT "   autodir_options_c/module__" . lc($Module) . "\n";
     my $ssout = "source/autodir_options_c/module__" . lc($Module) . ".rst";
     if ($print_description) { 
        printf "Auto-documenting options in module %s\n", lc($Module);
        open(SSOUT,">$ssout") or die "\nI can't write to $ssout\n";
        printf SSOUT ".. _`apdx:%s`:\n\n", lc($Module);
        printf SSOUT "\n%s\n%s\n\n", uc($Module), $Moddivider;
        if (exists $ModuleDescriptions{$Module}) { printf SSOUT "$ModuleDescriptions{$Module}\n\n"; }
        printf SSOUT ".. toctree::\n   :hidden:\n   :glob:\n\n   %s__*\n\n", lc($Module);
     }
     else { 
        open(SSOUT,">>$ssout") or die "\nI can't write to $ssout\n";
     }
     if (exists $ModuleDescriptions{$Module} and $print_description){
         printf OUT "\n{\\normalsize $ModuleDescriptions{$Module}}\\\\\n";
         # Insert an empty table entry as a spacer
         printf OUT '\\begin{tabular*}{\\textwidth}[tb]{c}';
         printf OUT "\n\t  \\\\ \n";
         print OUT "\\end{tabular*}\n";
     }
     my @OrderedSubsection = ("");
     if (exists $ModuleSubsections{$Module}) {
          @OrderedSubsection = @{$ModuleSubsections{$Module}};
     }
     foreach my $Subsection (@OrderedSubsection) {
       if (defined(%{$hash{$Module}{$Subsection}})) { 
         if($Subsection){
             if ($print_description) { 
                 my $Secdivider = "_" x (length($Subsection)-1);
                 print OUT "\\subsubsection{$Subsection}\n";
                 print SSOUT "\n$Subsection\n$Secdivider\n\n";
             }
             else {
                 my $Secdivider = "_" x (8+length($Subsection));
                 print OUT "\\subsubsection{$Subsection}\n";
                 print SSOUT "\n*Expert* $Subsection\n$Secdivider\n\n";
             }
         }
         else {
            if ($print_description) { print SSOUT "\nGeneral\n_______\n\n"; }
            else                    { print SSOUT "\n*Expert*\n________\n\n"; }
         }    
         my %SectionHash = %{$hash{$Module}{$Subsection}};
         foreach my $Keyword (sort {$a cmp $b} keys %SectionHash){
             my %KeyHash = %{$SectionHash{$Keyword}};
             my $DashedKeyword = $Keyword;
             $DashedKeyword =~ s/\\_/-/g;
             my $UnderscoredKeyword = $Keyword;
             $UnderscoredKeyword =~ s/\\_/_/g;
             my $Keydivider = "\"" x (14+length($Module)+2*length($UnderscoredKeyword));
             my $keywordfilename = lc($Module) . "__" . lc($UnderscoredKeyword);
             my $fullkeywordfilename = "source/autodir_options_c/" . $keywordfilename . ".rst";
             print SSOUT ".. include:: $keywordfilename.rst\n";
             open(SSSOUT,">$fullkeywordfilename") or die "\nI can't write to $fullkeywordfilename\n";
             printf SSSOUT ":term:`%s <%s (%s)>`\n%s\n\n", uc($UnderscoredKeyword), uc($UnderscoredKeyword), uc($Module), $Keydivider;
             my $SphCommentSubstit = $KeyHash{"SphComment"};
             $SphCommentSubstit = substitute_comment($SphCommentSubstit);
             printf SSSOUT "      %s\n\n", $SphCommentSubstit;
             printf TTSOUT ".. |%s__%s| replace:: :term:`%s <%s (%s)>`\n", lc($Module), lc($UnderscoredKeyword), uc($UnderscoredKeyword), uc($UnderscoredKeyword), uc($Module);
             if ($print_description) {
                printf TSOUT "   %s (%s)\n      :ref:`apdx:%s` |w---w| %s\n\n", uc($UnderscoredKeyword), uc($Module), uc($Module), $KeyHash{"SphComment"};
             } else {
                printf TSOUT "   %s (%s)\n      :ref:`apdx:%s` **(Expert)** |w---w| %s\n\n", uc($UnderscoredKeyword), uc($Module), uc($Module), $KeyHash{"SphComment"};
             }
             if(($KeyHash{"Type"} eq "bool") || ($KeyHash{"Type"} eq "boolean")) {
                printf SSSOUT "      * **Type**: :ref:`boolean <op_c_boolean>`\n";
                printf TSOUT  "      * **Type**: :ref:`boolean <op_c_boolean>`\n";
             }
             elsif (($KeyHash{"Type"} eq "double") && ((lc($Keyword) =~ /conv/) || (lc($Keyword) =~ /tol/))) {
                printf SSSOUT "      * **Type**: :ref:`conv double <op_c_conv>`\n";
                printf TSOUT  "      * **Type**: :ref:`conv double <op_c_conv>`\n";
             }
             elsif (($KeyHash{"Type"} eq "string") && ((lc($Keyword) eq "basis") || (index(lc($Keyword), "df\\_basis") == 0))) {
                printf SSSOUT "      * **Type**: %s\n", $KeyHash{"Type"};
                printf TSOUT  "      * **Type**: %s\n", $KeyHash{"Type"};
                printf SSSOUT "      * **Possible Values**: :ref:`basis string <apdx:basisElement>`\n";
                printf TSOUT  "      * **Possible Values**: :ref:`basis string <apdx:basisElement>`\n";
             }
             else {
                printf SSSOUT "      * **Type**: %s\n", $KeyHash{"Type"};
                printf TSOUT  "      * **Type**: %s\n", $KeyHash{"Type"};
             }
             printf OUT "\\paragraph{%s}\\label{op-%s-%s} \n", $Keyword, $Module, $DashedKeyword;
             printf OUT '\\begin{tabular*}{\\textwidth}[tb]{p{0.05\\textwidth}p{0.95\\textwidth}}';
             printf OUT "\n\t & %s \\\\ \n", $KeyHash{"Comment"};
             if($KeyHash{"Possibilities"}){
                  my @Options = split(/ +/, $KeyHash{"Possibilities"});
                  printf OUT "\n\t  & {\\bf Possible Values:} %s \\\\ \n", join(", ", @Options);
                  printf SSSOUT "      * **Possible Values**: %s\n", join(", ", @Options);
                  printf TSOUT "      * **Possible Values**: %s\n", join(", ", @Options);
             }
             print OUT "\\end{tabular*}\n";
             printf OUT '\\begin{tabular*}{\\textwidth}[tb]{p{0.3\\textwidth}p{0.35\\textwidth}p{0.35\\textwidth}}';
             printf OUT "\n\t   & {\\bf Type:} %s &  {\\bf Default:} %s\\\\\n\t & & \\\\\n", $KeyHash{"Type"}, $KeyHash{"Default"};
             printf SSSOUT "      * **Default**: %s\n\n", $KeyHash{"Default"};
             printf TSOUT "      * **Default**: %s\n\n", $KeyHash{"Default"};
             print OUT "\\end{tabular*}\n";
             close SSSOUT;
         }  # keyword
       }
     }  # subsection
     print SSOUT "\n";
     close SSOUT;
 }  # module
 print SOUT "\n";
 close SOUT;
 print TTSOUT "\n";
 close TTSOUT;
 print TSOUT "\n";
 close TSOUT;
 print OUT "}\n";
 close OUT;
}


sub get_description
{
 my $Line = shift;
 chomp $Line;
 my $String;
 #Process the inputted line
 if($Line =~ /(?:\/\*\-\s*MODULEDESCRIPTION)\s+(.*)/){
     $String = $1;
 }else{
     die "\nI don't know what to do with $Line\n";
 }
 if(!($String =~ s/\s*-\*\///g)){
     # 'twas more than a one-liner, let's keep searching
     while(<DRIVER>){
         # Add on the current line, sans the newline
         chomp;
         $String .= $_;
         # Attempt to nuke any -*/ patterns, if successful we're done
         last if $String =~ s/\s*-\*\///g;
     }
 }
 # Search and replace multiple spaces with a single space, not that LaTeX cares
 $String =~ s/ +/ /g;
 return $String;
}


# including the options abbr substitutions file in every SSSOUT option file slows
#   compilation by a factor of ten. so, back-translate |%s__%s| into :term:`%s`
sub substitute_comment
{
    my $Line = shift;
    while (1) {
        (my $before, my $pattern_module, my $pattern_keyword, my $after) = 
          $Line =~ m/^(.*?)[\s\(]\|(\w+)__(\w+)\|[\s\).,](.*?)$/;
        if ($pattern_module) {
            $Line = $before . " :term:`" . uc($pattern_keyword) . " <" . 
              uc($pattern_keyword) . " (" . uc($pattern_module) . ")>` " . $after;
        } else { last; }
    }
    return $Line
}


sub determine_comment
{
 my $Line = shift;
 chomp $Line;
 my $String;
 #Process the inputted line
 if($Line =~ /(?:\/\*\-)\s+(.*)/){
     $String = $1;
 }else{
     die "\nI don't know what to do with $Line\n";
 }
 if(!($String =~ s/\s*-\*\///g)){
     # 'twas more than a one-liner, let's keep searching
     while(<DRIVER>){
         # Add on the current line, sans the newline
         chomp;
         $String .= $_;
         # Attempt to nuke any -*/ patterns, if successful we're done
         last if $String =~ s/\s*-\*\///g;
     }
 }
 # Search and replace multiple spaces with a single space, not that LaTeX cares
 $String =~ s/ +/ /g;
 # Look for the expert flag and nuke it, if found
 my $Expert = ($String =~ s/!expert//g);
 ($String, $Expert);
}


sub determine_keyword_type_and_default
{
 my $Type    = "";
 my $Default = "No Default";
 my $Keyword = "";
 my $Possibilities = "";
 while(<DRIVER>){
     # Ignore blank lines
     next unless /\w+/;
     if(/add_str\(\s*\"(.*)\"\s*\,\s*\"(.*)\"\s*\,\s*\"(.*)\"\s*\)/){
         # This is a string, with default and options
         $Type = "str";
         $Keyword = $1;
         $Default = $2;
         $Possibilities = $3;
         $Default = "No Default" unless $Default =~ /\w+/;
     }elsif(/add_str_i\(\s*\"(.*)\"\s*\,\s*\"(.*)\"\s*\,\s*\"(.*)\"\s*\)/){
         # This is a string, with default and options
         $Type = "str";
         $Keyword = $1;
         $Default = $2;
         $Possibilities = $3;
         $Default = "No Default" unless $Default =~ /\w+/;
     }elsif(/add_str\(\s*\"(.*)\"\s*,\s*\"(.*)\"\s*\)/){
         # This is a string, with default, but without options
         $Type = "str";
         $Keyword = $1;
         $Default = $2;
         $Default = "No Default" unless $Default =~ /\w+/;
     }elsif(/add_str_i\(\s*\"(.*)\"\s*,\s*\"(.*)\"\s*\)/){
         # This is a string, with default, but without options
         $Type = "str";
         $Keyword = $1;
         $Default = $2;
         $Default = "No Default" unless $Default =~ /\w+/;
     }elsif(/add_bool\(\s*\"(.*)\"\s*,\s*(?:\")?([-\w]+)(?:\")?/){
         # This is a boolean with a default
         $Type = "bool";
         $Keyword = $1;
         $Default = lc ($2);
         if ($Default eq "1") { $Default = "true"; }
         if ($Default eq "0") { $Default = "false"; }
     }elsif(/add_double\(\s*\"(.*)\"\s*,\s*(?:\")?([-\/.\w]+)(?:\")?/){
         # This is a double with a default
         $Type = "double";
         $Keyword = $1;
         $Default = lc ($2);
     }elsif(/add_(\w+)\(\s*\"(\w+)\"\s*\,\s*(?:\")?([-\w]+)(?:\")?/){
         # This is a keyword with a default
         $Type = $1;
         $Keyword = $2;
         $Default = $3;
     }elsif(/add\(\s*\"(\w*)\"\,\s*new\s+(\w+)\(\)/){
         # This is a custom DataType thingy
         $Keyword = $1;
         if($2 eq "ArrayType"){
             $Type = "array";
         }elsif($2 eq "MapType"){
             $Type = "map";
         }elsif($2 eq "PythonDataType"){
             $Type = "python";
         }else{
             print $_;
             die "\nUnrecognized type: $2\n";
         }
     }
     last if /\;/;
 }
 # Make the underscores LaTeX-friendly
 $Keyword =~ s/_/\\_/g;
 $Default =~ s/_/\\_/g;
 $Possibilities =~ s/_/\\_/g;
 if($Type eq "str"){
     $Type = "string";
 }elsif($Type eq "int"){
     $Type = "integer";
 }elsif($Type eq "bool"){
     $Type = "boolean";
 }elsif($Type eq "double"){
 }elsif($Type eq "array"){
 }elsif($Type eq "map"){
 }elsif($Type eq "python"){
 }else{
     print $_;
     die "\nUnrecognized type: $Type\n";
 }
 ($Keyword, $Type, $Default, $Possibilities);
}

