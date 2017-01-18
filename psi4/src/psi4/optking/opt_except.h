/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include <vector>

namespace opt {

/**
  \ingroup optking
  \brief   Exception class for problems with internal coordinates.
  */

/* If they relate to values of the coordinates and derivatives, then try new
   coordinates ; if it looks like user error in definition, then
   quit right away.  */
class INTCO_EXCEPT {
  private:
   const char * message;
   bool really_quit;

  public:

    static int dynamic_level; // defined in optking.cc
    //static bool override_fragment_mode;
    static std::vector<int> linear_angles;


    INTCO_EXCEPT(const char * m) {
      message = m;
      really_quit = false;
      //try_other_intcos = false;
    }

    INTCO_EXCEPT(const char * m, bool t) {
      message = m;
      really_quit = t;
      //try_other_intcos = t;
    }

    INTCO_EXCEPT(const char * m, std::vector<int> l) {
      message = m;
      linear_angles = l;
      really_quit = false;
      //try_other_intcos = t;
    }

    void increment_dynamic_level(void) {
      if (dynamic_level == 0) // turning 'on' dynamic
        dynamic_level = 1;

      dynamic_level += 1;
    }

    ~INTCO_EXCEPT() {};

    //bool try_again() { return try_other_intcos; }
    const char *g_message(void) { return message; }
    bool g_really_quit(void) { return really_quit; }
};

class BAD_STEP_EXCEPT {
  private:
   const char * message;

  public:
    BAD_STEP_EXCEPT(const char * m) { message = m; }

    ~BAD_STEP_EXCEPT() {};

    const char *g_message(void) { return message; }
};

class BROKEN_SYMMETRY_EXCEPT {
  private:
    const char * message;

  public:
    BROKEN_SYMMETRY_EXCEPT(const char * m) { message = m; }

    ~BROKEN_SYMMETRY_EXCEPT() {};

    const char *g_message(void) { return message; }
};

}