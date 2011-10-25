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
# First, read the options for each module
#
my $DriverPath = "";
if ($#ARGV == 0) { $DriverPath = $ARGV[0] . "/"; }
my $DriverFile = $DriverPath . "../../src/bin/psi4/read_options.cc";

my $CurrentModule;
my %Keywords;
my %Expert;
open(DRIVER, $DriverFile) or die "\nI can't read the PSI driver file\n";

# Start looping over the driver file
while(<DRIVER>){
    # Look for the name of the module for which options are being added
    $CurrentModule = $1 if(/name\s*\=\=\s*\"(\w+)\"/);
    my $CommentString;
    my $Type;
    my $Default;
    my $Keyword;
    my $Possibilities;
    my $Expert;
    # If we find a /*- it means the comment block has started, but we
    # don't know if it's a multi-line comment, let's find out
    if(/\/\*-/ and $CurrentModule){
        ($CommentString, $Expert) = determine_comment($_);
        $CommentString =~ s/_/\\_/g;
        ($Keyword, $Type, $Default, $Possibilities) = determine_keyword_type_and_default();
        if($Expert){
            $Expert{$CurrentModule}{$Keyword}{"Type"}    = $Type;
            $Expert{$CurrentModule}{$Keyword}{"Default"} = $Default;
            $Expert{$CurrentModule}{$Keyword}{"Comment"} = $CommentString;
            $Expert{$CurrentModule}{$Keyword}{"Possibilities"} = $Possibilities;
        }else{
            $Keywords{$CurrentModule}{$Keyword}{"Type"}    = $Type;
            $Keywords{$CurrentModule}{$Keyword}{"Default"} = $Default;
            $Keywords{$CurrentModule}{$Keyword}{"Comment"} = $CommentString;
            $Keywords{$CurrentModule}{$Keyword}{"Possibilities"} = $Possibilities;
        }
    }
}
close DRIVER;

print_hash(\%Keywords, "keywords.tex", "Keywords Recognized by Each Module", "keywords");
print_hash(\%Expert, "expert_keywords.tex", "Expert Keywords Recognized by Each Module, for Advanced Users", "expertkeywords");



sub print_hash
{
 my %hash     = %{$_[0]};
 my $filename = $_[1];
 my $title    = $_[2];
 my $label    = $_[3];
 open(OUT,">$filename") or die "\nI can't write to $filename\n";
 print OUT "\\section{$title}\\label{$label}\n";
 print OUT "{\n \\footnotesize\n";
 foreach my $Module (sort {$a gt $b} keys %hash){
     printf OUT "\n\\subsection{%s}\n",$Module;
     foreach my $Keyword (sort {$a gt $b} keys %{$hash{$Module}}){
         printf OUT '\\begin{tabular*}{\\textwidth}[tb]{p{0.3\\textwidth}p{0.7\\textwidth}}';
         printf OUT "\n\t %s & %s \\\\ \n", $Keyword, $hash{$Module}{$Keyword}{"Comment"};
         my $Options = $hash{$Module}{$Keyword}{"Possibilities"};
         if($Options =~ /\w+/){
              my @Options = split(/ +/,$Options);
              printf OUT "\n\t  & {\\bf Possible Values:} %s \\\\ \n", join(", ", @Options);
         }
         print OUT "\\end{tabular*}\n";
         printf OUT '\\begin{tabular*}{\\textwidth}[tb]{p{0.3\\textwidth}p{0.35\\textwidth}p{0.35\\textwidth}}';
         printf OUT "\n\t   & {\\bf Type:} %s &  {\\bf Default:} %s\\\\\n\t & & \\\\\n", 
                     $hash{$Module}{$Keyword}{"Type"}, $hash{$Module}{$Keyword}{"Default"}; 
         print OUT "\\end{tabular*}\n";
     }
 }
 print OUT "}\n";
 close OUT;
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
 }else{
     print $_;
     die "\nUnrecognized type: $Type\n";
 }
 ($Keyword, $Type, $Default, $Possibilities);
}


#
# Now we raid the test cases looking for tags
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


