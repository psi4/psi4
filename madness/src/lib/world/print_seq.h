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

  $Id: print_seq.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/
#ifndef MADNESS_WORLD_PRINT_SEQ_H__INCLUDED
#define MADNESS_WORLD_PRINT_SEQ_H__INCLUDED

/// \file print_seq.h
/// \brief Implements print_seq ... included by world.h

namespace madness {
    /// Sequentially ordered printing of (serializable) data from every process ... collective no fence
    template <typename A, typename B, typename C, typename D>
    void print_seq(World& world, const A& a, const B& b, const C& c, const D& d) {
        if (world.rank() == 0) {
            printf("%6d : ",0);
            print(a, b, c, d);
            for (int p=1; p<world.size(); ++p) {
                A aa;
                B bb;
                C cc;
                D dd;
                archive::MPIOutputArchive(world,p) & 1;
                archive::MPIInputArchive(world, p) & aa & bb & cc & dd;
                printf("%6d : ",p);
                print(aa, bb, cc, dd);
            }
        }
        else {
            int i;
            archive::MPIInputArchive(world, 0) & i;
            archive::MPIOutputArchive(world, 0) & a & b & c & d;
        }
    }

    /// Sequentially ordered printing of (serializable) data from every process ... collective no fence
    template <typename A, typename B, typename C>
    void print_seq(World& world, const A& a, const B& b, const C& c) {
        if (world.rank() == 0) {
            printf("%6d : ",0);
            print(a, b, c);
            for (int p=1; p<world.size(); ++p) {
                A aa;
                B bb;
                C cc;
                archive::MPIOutputArchive(world,p) & 1;
                archive::MPIInputArchive(world, p) & aa & bb & cc;
                printf("%6d : ",p);
                print(aa, bb, cc);
            }
        }
        else {
            int i;
            archive::MPIInputArchive(world, 0) & i;
            archive::MPIOutputArchive(world, 0) & a & b & c;
        }
    }

    /// Sequentially ordered printing of (serializable) data from every process ... collective no fence
    template <typename A, typename B>
    void print_seq(World& world, const A& a, const B& b) {
        if (world.rank() == 0) {
            printf("%6d : ",0);
            print(a, b);
            for (int p=1; p<world.size(); ++p) {
                A aa;
                B bb;
                archive::MPIOutputArchive(world,p) & 1;
                archive::MPIInputArchive(world, p) & aa & bb;
                printf("%6d : ",p);
                print(aa, bb);
            }
        }
        else {
            int i;
            archive::MPIInputArchive(world, 0) & i;
            archive::MPIOutputArchive(world, 0) & a & b;
        }
    }

    /// Sequentially ordered printing of (serializable) data from every process ... collective no fence
    template <typename A>
    void print_seq(World& world, const A& a) {
        if (world.rank() == 0) {
            printf("%6d : ",0);
            print(a);
            for (int p=1; p<world.size(); ++p) {
                A aa;
                archive::MPIOutputArchive(world,p) & 1;
                archive::MPIInputArchive(world, p) & aa;
                printf("%6d : ",p);
                print(aa);
            }
        }
        else {
            int i;
            archive::MPIInputArchive(world, 0) & i;
            archive::MPIOutputArchive(world, 0) & a;
        }
    }
}

#endif // MADNESS_WORLD_PRINT_SEQ_H__INCLUDED
