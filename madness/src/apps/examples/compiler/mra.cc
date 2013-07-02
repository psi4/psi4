#include <iostream>
#include <fstream>
#include "mra-driver.hh"

std::map<std::string, bool> dectab;        // keeps track of CXX declared symbols

int
main (int argc, char *argv[])
{
  mra_driver driver;
  for (++argv; argv[0]; ++argv)
    if (*argv == std::string ("-p"))
      driver.trace_parsing = true;
    else if (*argv == std::string ("-s"))
      driver.trace_scanning = true;
    else if (!driver.parse (*argv))
        ;
  
  std::ofstream treefile("prog.tree",std::ios_base::trunc);
  std::ostream treeout(treefile.rdbuf());
  driver.print_tree(treeout);
  treefile.close();

  std::ofstream regfile("prog.reg",std::ios_base::trunc);
  std::ostream regout(regfile.rdbuf());
  driver.regenerate(regout);
  regfile.close();

  std::ofstream texfile("prog.tex",std::ios_base::trunc);
  std::ostream texout(texfile.rdbuf());
  driver.generate_tex(texout);
  texfile.close();

  std::ofstream ccfile("prog.cc",std::ios_base::trunc);
  std::ostream ccout(ccfile.rdbuf());
  driver.generate_cxx(ccout);
  ccfile.close();

  return 0;
}
