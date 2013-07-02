#!/usr/bin/perl

system ("ls -1 *cc-*.gbs | grep -v 'autogen' | grep -v 'tight' | grep -v 'polarization' | grep -v 'molpro' | grep -v 'diffuse' | grep -v 'basis' | grep -v 'corevalence' | grep -v 'hold' > realdunnings.txt");

system ("echo 'changed dunning basis sets\n\n' > basisdunningdiffer.txt");
open IN, "realdunnings.txt";
@text = <IN>;
close(IN);

foreach $file (@text) {

   chomp($file);
   print "$file\n";

   open BASFILE, $file;
   @bastext = <BASFILE>;
   close(BASFILE);

   if (!($bastext[0] =~ /spherical/)) {
      system ("./prepend_puream_onefile.sh $file");
   }

   #system("diff -bw -I '^! merged' $file ../../$file");  # less restrictive- shows changes in comment lines that arise when basis sets change but affected atoms not present
   system("diff -bw -I '^!' $file ../../$file");
   #system("diff -bw -I '^!' $file ../$file");
   system("diff -qbw -I '^!' $file ../../$file >> basisdunningdiffer.txt");
   #system("diff -qbw -I '^!' $file ../$file >> basisdunningdiffer.txt");
   print "\n\n\n";
}


