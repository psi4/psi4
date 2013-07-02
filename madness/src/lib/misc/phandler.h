/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


*/

#ifndef MADNESS_MISC_PHANDLER_H__INCLUDED
#define MADNESS_MISC_PHANDLER_H__INCLUDED

/// \file misc/phandler.h
/// \brief Interface for the muParser library for turning user-defined functions into bytecode.

/* Example:

   #include <misc/phandler.h>

    typedef FunctionFactory<double,3> factoryT;
    typedef std::shared_ptr< FunctionFunctorInterface<double, 3> > functorT;
    Function<double,3> pFunc = factoryT(world).functor(functorT(
                      new ParserHandler<double,3>("exp(-abs(r))")));

    pfunc.compress();

    ...
*/

#include <muParser/muParser.h>
#include <string>
#include <mra/mra.h>

//using namespace madness;
//using namespace mu;

// T can be double or complex, but muParser results will always be real
// NDIM can be 1 to 6
// Varables allowed in strings are x,y,z,u,v,w, and r
// "r" will always mean magnitude of given vector

template <typename T, int NDIM>
class ParserHandler : public madness::FunctionFunctorInterface<T, NDIM> {

  private:
    mu::Parser parser;                   // the string-to-function parser
    static const int MAX_DIM = 6;
    mutable double vars[MAX_DIM];              // variables used in expression
    mutable double r;                    // distance to origin
    typedef madness::Vector<double, NDIM> coordT;

  public:
    ParserHandler(std::string expr)   {
        if (NDIM > MAX_DIM) MADNESS_EXCEPTION("too many dim for parser!",0);
        try {
          if (NDIM >= 1) parser.DefineVar("x", &vars[0]);
          if (NDIM >= 2) parser.DefineVar("y", &vars[1]);
          if (NDIM >= 3) parser.DefineVar("z", &vars[2]);
          if (NDIM >= 4) parser.DefineVar("u", &vars[3]);
          if (NDIM >= 5) parser.DefineVar("v", &vars[4]);
          if (NDIM >= 6) parser.DefineVar("w", &vars[5]);
          parser.DefineVar("r", &r);
          parser.SetExpr(expr);
        } catch (mu::Parser::exception_type &e) {
          std::cout << "muParser: " << e.GetMsg() << std::endl;
        }
    }

    virtual T operator() (const coordT &vals_in) const {
      r = 0;
      for (int i = 0; i < NDIM; ++i) {
        vars[i] = vals_in[i];
        r += vals_in[i]*vals_in[i];
      }
      r = sqrt(r);

      try {
        return parser.Eval();
      } catch (mu::Parser::exception_type &e) {
        std::cout << "muParser: " << e.GetMsg() << std::endl;
        return 0.0;
      }
    } // end operator()

    void changeExpr(const std::string newExpr) {
      try {
        parser.SetExpr(newExpr);
      } catch (mu::Parser::exception_type &e) {
        std::cout << "muParser: " << e.GetMsg() << std::endl;
      }
    }
}; // end class parserhandler

#endif // MADNESS_MISC_PHANDLER_H__INCLUDED
