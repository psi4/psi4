#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <vector>

#include "Comm.h"
#include "FuncReg.h"

using namespace std;
using namespace tascel;
using namespace tascel::comm;


TslFuncRegTbl::TslFuncRegTbl()
  : ftbl()
{}

TslFuncRegTbl::TslFuncRegTbl(const TslFuncRegTbl &that)
  : ftbl(that.ftbl)
{}

TslFuncRegTbl& TslFuncRegTbl::operator = (const TslFuncRegTbl &that)
{
  if (this == &that) {
      return *this;
  }
  ftbl = that.ftbl;
  return *this;
}

TslFuncRegTbl::~TslFuncRegTbl() {}

TslFunc
TslFuncRegTbl::add(TslFunc_t f) {
  barrier();
  ftbl.push_back(f);
  return ftbl.size() - 1;
}

TslFunc_t
TslFuncRegTbl::get(TslFunc fn) const {
  return ftbl.at(fn);
}

