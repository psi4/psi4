#!/usr/bin/perl

use strict;
use warnings;
use File::Path qw(remove_tree);

# This script does two things: a) reads the driver's options setup to provide a list of 
# keywords expected by each module; b) looks for comments in each test case to identify
# the type of calculation.  The results of this parsing is put into TeX files, which are
# inlined into the manual.
#
# Then, we look at the test cases, using the specially formatted comments to add a listing
# of samples in the users' manual; the tests are copied (sans validation stuff) to the appropriate
# location in the samples folder.
#
# Then, for each module detected in read_options.cc, we collect the process environment
# variables set in src/bin/module and src/lib/libmodule_solver and place them in a TeX file,
# which is inlined into the manual.

#
# First, read the options for each module
#
my $DriverPath = "";
if ($#ARGV == 0) { $DriverPath = $ARGV[0] . "/"; }
my $DriverFile = $DriverPath . "../../src/bin/psi4/read_options.cc";

my $CurrentModule;
my %Keywords;
my %ModuleDescriptions;
my %Expert;
open(DRIVER, $DriverFile) or die "\nI can't read the PSI driver file\n";

my $CurrentSubsection = "";

# Start looping over the driver file
while(<DRIVER>){
    # Look for the name of the module for which options are being added
    if(/name\s*\=\=\s*\"(\w+)\"/){
        $CurrentModule = $1;
        $CurrentSubsection = "";
    }
    my $CommentString;
    my $Type;
    my $Default;
    my $Keyword;
    my $Possibilities;
    my $Expert;
    # If we find a /*- it means the comment block has started, but we
    # don't know if it's a multi-line comment, let's find out
    if(/\/\*-\s*SUBSECTION\s+(.*)-\*\// and $CurrentModule){
        $CurrentSubsection = $1;
    }elsif(/\/\*-\s*MODULEDESCRIPTION/ and $CurrentModule){
        $ModuleDescriptions{$CurrentModule} = get_description($_);
    }elsif(/\/\*-/ and $CurrentModule){
        ($CommentString, $Expert) = determine_comment($_);
        $CommentString =~ s/_/\\_/g;
        # process @@ as math mode subscript in tex
        $CommentString =~ s/@@/_/g;
        ($Keyword, $Type, $Default, $Possibilities) = determine_keyword_type_and_default();
        if($Expert){
            $Expert{$CurrentModule}{$CurrentSubsection}{$Keyword}{"Type"}    = $Type;
            $Expert{$CurrentModule}{$CurrentSubsection}{$Keyword}{"Default"} = $Default;
            $Expert{$CurrentModule}{$CurrentSubsection}{$Keyword}{"Comment"} = $CommentString;
            $Expert{$CurrentModule}{$CurrentSubsection}{$Keyword}{"Possibilities"} = $Possibilities;
        }else{
            $Keywords{$CurrentModule}{$CurrentSubsection}{$Keyword}{"Type"}    = $Type;
            $Keywords{$CurrentModule}{$CurrentSubsection}{$Keyword}{"Default"} = $Default;
            $Keywords{$CurrentModule}{$CurrentSubsection}{$Keyword}{"Comment"} = $CommentString;
            $Keywords{$CurrentModule}{$CurrentSubsection}{$Keyword}{"Possibilities"} = $Possibilities;
        }
    }
}
close DRIVER;

my @temp = ();
print_hash(\%Keywords, "keywords.tex", 1);
my @PSIMODULES = @temp;
print_hash(\%Expert, "expert_keywords.tex", 0);

sub print_hash
{
 my %hash     = %{$_[0]};
 my $filename = $_[1];
 my $print_description = $_[2];
 open(OUT,">$filename") or die "\nI can't write to $filename\n";
 print OUT "{\n \\footnotesize\n";
 foreach my $Module (sort {$a cmp $b} keys %hash){
     push(@temp, $Module);
     printf OUT "\n\\subsection{%s}\n",$Module;
     if (exists $ModuleDescriptions{$Module} and $print_description){
         printf OUT "\n{\\normalsize $ModuleDescriptions{$Module}}\\\\\n";
         # Insert an empty table entry as a spacer
         printf OUT '\\begin{tabular*}{\\textwidth}[tb]{c}';
         printf OUT "\n\t  \\\\ \n";
         print OUT "\\end{tabular*}\n";
     }
     foreach my $Subsection (sort {$a gt $b} keys %{$hash{$Module}}){
         if($Subsection){
             print OUT "\\subsubsection{$Subsection}\n";
         }
         my %SectionHash = %{$hash{$Module}{$Subsection}};
         foreach my $Keyword (sort {$a gt $b} keys %SectionHash){
             my %KeyHash = %{$SectionHash{$Keyword}};
             printf OUT '\\begin{tabular*}{\\textwidth}[tb]{p{0.3\\textwidth}p{0.7\\textwidth}}';
             printf OUT "\n\t %s & %s \\\\ \n", $Keyword, exists $KeyHash{"Comment"} ? $KeyHash{"Comment"} : "";
             if(exists $KeyHash{"Possibilities"}){
                  my @Options = split(/ +/, $KeyHash{"Possibilities"});
                  printf OUT "\n\t  & {\\bf Possible Values:} %s \\\\ \n", join(", ", @Options);
             }
             print OUT "\\end{tabular*}\n";
             printf OUT '\\begin{tabular*}{\\textwidth}[tb]{p{0.3\\textwidth}p{0.35\\textwidth}p{0.35\\textwidth}}';
             printf OUT "\n\t   & {\\bf Type:} %s &  {\\bf Default:} %s\\\\\n\t & & \\\\\n", $KeyHash{"Type"}, $KeyHash{"Default"};
             print OUT "\\end{tabular*}\n";
         }
     }
 }
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
         $Type = "bool";
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


#
# Secondly we raid the test cases looking for tags
#

my $SamplesFolder = $DriverPath . "../../samples";
opendir(SAMPLES, $SamplesFolder) or die "I can't read $SamplesFolder\n";
foreach my $File(readdir SAMPLES){
    next unless $File =~ /\w+/; # Make sure we don't nuke .. or . !
    next if $File =~ /example_psi4rc_file/; # Keep the example psi4rc file
    remove_tree("$SamplesFolder/$File");
}

my $TestsFolder = $DriverPath . "../../tests";
opendir(TESTS, $TestsFolder) or die "I can't read $TestsFolder\n";
my $TexSummary = "tests_descriptions.tex";
# Create a plain-text summary in the samples directory
my $Summary = $SamplesFolder."/SUMMARY";
open(SUMMARY,">$Summary") or die "I can't write to $Summary\n";
# Make a LaTeX version for the manual, too
open(TEXSUMMARY,">$TexSummary") or die "I can't write to $TexSummary\n";
foreach my $Dir(readdir TESTS){
    my $Input = $TestsFolder."/".$Dir."/input.dat";
    # Look for an input file in each subdirectory, or move on
    open(INPUT, "<$Input") or next;
    #
    # Remember, we want to copy the test suite over to psi4/samples, omitting the test functions
    #
    my $SamplesDirectory = $SamplesFolder."/".$Dir;
    unless(-d $SamplesDirectory){
        # This directory doesn't exist in psi4/samples, make it now
        mkdir $SamplesDirectory or die "\nI can't create $SamplesDirectory\n";
    }
    my $SampleInput = $SamplesDirectory."/input.dat";
    open(SAMPLE, ">$SampleInput") or die "\nI can't write to $SampleInput\n";
    my $Description;
    while(<INPUT>){
        # If this line isn't associated with testing, put it in the sample
        print SAMPLE unless /\#TEST/;
        # Now we only want to grab the comments.  Move on if this is not a comment.
        next unless s/\#\!//;
        $Description .= $_;
        chomp $Description;
    }
    close INPUT;
    close SAMPLE;

    # Process the comment that we grabbed from the input.
    if($Description){
        # make directory name tex safe
        my $Dir_tex = $Dir;
        $Dir_tex =~ s/_/\\_/g;
        my $Description_tex = $Description;
        $Description_tex =~ s/_/\\_/g;
        print TEXSUMMARY "\\begin{tabular*}{\\textwidth}[tb]{p{0.2\\textwidth}p{0.8\\textwidth}}\n";
        print TEXSUMMARY "{\\bf $Dir_tex} & $Description_tex \\\\\n\\\\\n";
        print TEXSUMMARY "\\end{tabular*}\n";
        printf SUMMARY "%-12s %s\n\n\n", $Dir.":", $Description;
    }else{
        warn "Warning!!! Undocumented input: $Input\n";
    }
}
close TEXSUMMARY;
close SUMMARY;
closedir TESTS;


#
# Thirdly, scan the source for Process::Environment variables
#

my $SrcFolder = $DriverPath . "../../src";
$TexSummary = "variables_list.tex";
open(TEXOUT,">$TexSummary") or die "I can't write to $TexSummary\n";
print TEXOUT "\\section{Environment Variables Set by Each Module}\\label{variableslist}\n";
print TEXOUT "{\n \\footnotesize\n";

# Grab psi modules and ordering from options parsing above
foreach my $Module (@PSIMODULES) {
    # Set path for each module of bin/module and lib/libmodule_solver
    #     Assign stray variables as for OEPROP below
    my @RelevantDirs = ($SrcFolder . "/bin/" . lc($Module), $SrcFolder . "/lib/lib" . lc($Module) . "_solver");
    if ($Module eq "OEPROP") { push(@RelevantDirs, $SrcFolder . "/lib/libmints"); }
    my @EnvVariables = ();
    printf TEXOUT "\n\\subsection{%s}\n",$Module;
    # Search each line in each file in each module-relevant directory for environment variables
    foreach my $Dir (@RelevantDirs) {
        if (opendir(SRC, $Dir)) {
            while (my $file = readdir(SRC)) {
                if (open(CODE, "<$Dir/$file")) {
                    my @text = <CODE>;
                    foreach my $line (@text) {
                        if ($line =~ /\QProcess::environment.globals\E/) {
                            my @ltemp = split( /"/, $line);
                            if ($ltemp[0] =~ /\QProcess::environment.globals\E/) {
                                push(@EnvVariables, $ltemp[1]);
                            }
                        }
                    }
                    close(CODE);
                }
            }
            closedir SRC;
        }
    }
    # Remove duplicate env variables, sort into alphabetical order, and print to tex file
    my %EnvHash;
    foreach my $EnvVar (@EnvVariables){
         $EnvHash{$EnvVar} = 1 if $EnvVar;
    }
    foreach my $Var (sort keys %EnvHash) {
        printf TEXOUT '\\begin{tabular*}{\\textwidth}[tb]{p{1.0\\textwidth}}';
        printf TEXOUT "\n\t %s \\\\ \n", $Var;
        print TEXOUT "\\end{tabular*}\n";
    }
}
print TEXOUT "}\n";
close TEXOUT;

#
# Now, grab the physical constants
#
my $PhysconstFile = $DriverPath . "../../include/physconst.h";
my $PyPhysconstFile = $DriverPath . "../../lib/python/physconst.py";
open(PHYSCONST, "<$PhysconstFile") or die "I can't open $PhysconstFile\n";
open(TEXOUT, ">physconst.tex") or die "I can't write to phyconst.tex\n";
open(PYOUT, ">$PyPhysconstFile") or die "I can't write to $PyPhysconstFile\n";
print PYOUT "# Do not modify this file! It is auto-generated by the document_options_and_tests\n".
            "# script, from psi4topdir/include/physconst.h\n";
while(<PHYSCONST>){
    next unless /\s*#define\s+(\w+)\s+([-Ee0-9.]+)\s+\/\*-(.*)-\*\//;
    my $Var     = $1;
    my $Val     = $2;
    my $Comment = $3;
    printf PYOUT "%-25s = %-20s #%-40s\n", $Var, $Val, $Comment;
    $Var =~ s/_/\\_/g; # Make things TeX-friendly
    printf TEXOUT "%-25s & %-20s & %-40s\\\\\n", $Var, $Val, $Comment;
}
close PHYSCONST;
close PYOUT;
close TEXOUT;
