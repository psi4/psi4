#!/usr/bin/perl

use strict;
use warnings;

# This script does two things: a) reads the driver's options setup to provide a list of 
# keywords expected by each module; b) looks for comments in each test case to identify
# the type of calculation.  The results of this parsing is put into TeX files, which are
# inlined into the manual.


#
# First, read the options for each module
#
my $DriverFile = "../../src/bin/psi4/read_options.cc";

my $CurrentModule;
my %Data;
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
    # If we find a /*- it means the comment block has started, but we
    # don't know if it's a multi-line comment, let's find out
    if(/\/\*-/){
        $CommentString = determine_comment($_);
        ($Keyword, $Type, $Default, $Possibilities) = determine_keyword_type_and_default();
        $Data{$CurrentModule}{$Keyword}{"Type"}    = $Type;
        $Data{$CurrentModule}{$Keyword}{"Default"} = $Default;
        $Data{$CurrentModule}{$Keyword}{"Comment"} = $CommentString;
        $Data{$CurrentModule}{$Keyword}{"Possibilities"} = $Possibilities;
    }
}
close DRIVER;

open(OUT,">keywords.tex") or die "\nI can't write to keywords.tex\n";
print OUT "\\section{Keywords Recognized by Each Module}\n";
print OUT "{\n \\footnotesize\n";
foreach my $Module (sort {$a gt $b} keys %Data){
    printf OUT "\n\\subsection{%s}\n",$Module;
    foreach my $Keyword (sort {$a gt $b} keys %{$Data{$Module}}){
        printf OUT '\\begin{tabular*}{\\textwidth}[tb]{p{0.3\\textwidth}p{0.7\\textwidth}}';
        printf OUT "\n\t %s & %s \\\\ \n", $Keyword, $Data{$Module}{$Keyword}{"Comment"};
        my $Options = $Data{$Module}{$Keyword}{"Possibilities"};
        if($Options =~ /\w+/){
             my @Options = split(/ +/,$Options);
             printf OUT "\n\t  & {\\bf Possible Values:} %s \\\\ \n", join(", ", @Options);
        }
        print OUT "\\end{tabular*}\n";
        printf OUT '\\begin{tabular*}{\\textwidth}[tb]{p{0.3\\textwidth}p{0.35\\textwidth}p{0.35\\textwidth}}';
        printf OUT "\n\t   & {\\bf Type:} %s &  {\\bf Default:} %s\\\\\n\t & & \\\\\n", 
                    $Data{$Module}{$Keyword}{"Type"}, $Data{$Module}{$Keyword}{"Default"}; 
        print OUT "\\end{tabular*}\n";
    }
}
print OUT "}\n";
close OUT;

sub determine_comment
{
 my $Line = shift;
 my $String;
 #Process the inputted line
 if($Line =~ /(?:\/\*\-)\s+(.*)/){
     $String = $1;
     # If we find -*/ we're done
     return $String if $String =~ s/ -\*\///g;
 }else{
     die "\nI don't know what to do with $Line\n";
 }
 if($Line !~ /-\*\//){
     # 'twas more than a one-liner, let's keep searching
     while(<DRIVER>){
         # Add on the current line, sans the newline
         chomp;
         $String .= $_;
         # Attempt to nuke any -*/ patterns, if successful we're done
         last if $String =~ s/ -\*\///g;
     }
 }
 # Search and replace multiple spaces with a single space, not that LaTeX cares
 $String =~ s/ +/ /g;
 $String;
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
     }elsif(/add_(\w+)\(\s*\"(\w+)\"\s*\,\s*(?:\")?(\w+)(?:\")?/){
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

my $TestsFolder = "../../tests";
opendir(TESTS, $TestsFolder) or die "I can't read $TestsFolder\n";
my $Out = "tests_descriptions.tex";
open(OUT,">$Out") or die "I can't write to $Out\n";
printf OUT "\\begin{tabular*}{\\textwidth}[tb]{p{0.15\\textwidth}p{0.75\\textwidth}}\n";
foreach my $Dir(readdir TESTS){
    my $Input = join("/",$TestsFolder,$Dir,"input.dat");
    # Look for an input file in each subdirectory, or move on
    open(INPUT, "<$Input") or next;
    my $Description;
    while(<INPUT>){
        next unless s/\%\!//;
        $Description .= $_;
        chomp $Description;
    }
    close INPUT;
    if($Description){
        print OUT "{\\bf $Dir} & $Description\\\\\n";
    }else{
        warn "Warning!!! Undocumented input: $Input\n";
    }
}
print OUT "\\end{tabular*}";
close OUT;
closedir TESTS;


