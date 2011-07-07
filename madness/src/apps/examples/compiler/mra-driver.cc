#include "mra-driver.hh"
#include "mra-parser.hh"

mra_driver::mra_driver ()
    : scopedepth(0)
    , tmpvarcnt(0)
    , use_k_default(true)
    , use_eps_default(true)
    , trace_scanning (false)
    , trace_parsing (false)
{
    insert_sym("built in", new Exp("pi",Exp::REAL));
}

mra_driver::~mra_driver ()
{
}

int
mra_driver::parse (const std::string &f)
{
  file = f;
  scan_begin ();
  yy::mra_parser parser (*this);
  parser.set_debug_level (trace_parsing);
  int res = parser.parse ();
  scan_end ();
  return res;
}

void
mra_driver::error (const yy::location& l, const std::string& m)
{
  std::cerr << l << ": " << m << std::endl;
}

void
mra_driver::error (const std::string& m)
{
  std::cerr << m << std::endl;
}
