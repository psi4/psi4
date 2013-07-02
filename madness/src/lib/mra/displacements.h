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

  $Id: displacements.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/
#ifndef MADNESS_MRA_DISPLACEMENTS_H__INCLUDED
#define MADNESS_MRA_DISPLACEMENTS_H__INCLUDED

namespace madness {
    /// Holds displacements for applying operators to avoid replicating for all operators
    template <std::size_t NDIM>
    class Displacements {

        static std::vector< Key<NDIM> > disp;
        static std::vector< Key<NDIM> > disp_periodicsum[64];

        static int bmax_default() {
            int bmax;
            if      (NDIM == 1) bmax = 7;
            else if (NDIM == 2) bmax = 5;
            else if (NDIM == 3) bmax = 3;
            else                bmax = 2;
            return bmax;
        }

        static bool cmp_keys(const Key<NDIM>& a, const Key<NDIM>& b) {
            return a.distsq() < b.distsq();
        }

        static bool cmp_keys_periodicsum(const Key<NDIM>& a, const Key<NDIM>& b) {
            Translation twonm1 = (Translation(1)<<a.level())>>1;

            uint64_t suma=0, sumb=0;
            for (std::size_t d=0; d<NDIM; ++d) {
                Translation la = a.translation()[d];
                if (la > twonm1) la -= twonm1*2;
                if (la <-twonm1) la += twonm1*2;
                suma += la*la;

                Translation lb = b.translation()[d];
                if (lb > twonm1) lb -= twonm1*2;
                if (lb <-twonm1) lb += twonm1*2;
                sumb += lb*lb;
            }
            return suma < sumb;
        }

        static void make_disp(int bmax) {
            // Note newer loop structure in make_disp_periodic_sum
            Vector<Translation,NDIM> d(0);

            int num = 1;
            for (std::size_t i=0; i<NDIM; ++i) num *= (2*bmax + 1);
            disp.resize(num,Key<NDIM>(0));

            num = 0;
            if (NDIM == 1) {
                for (d[0]=-bmax; d[0]<=bmax; ++d[0])
                    disp[num++] = Key<NDIM>(0,d);
            }
            else if (NDIM == 2) {
                for (d[0]=-bmax; d[0]<=bmax; ++d[0])
                    for (d[1]=-bmax; d[1]<=bmax; ++d[1])
                        disp[num++] = Key<NDIM>(0,d);
            }
            else if (NDIM == 3) {
                for (d[0]=-bmax; d[0]<=bmax; ++d[0])
                    for (d[1]=-bmax; d[1]<=bmax; ++d[1])
                        for (d[2]=-bmax; d[2]<=bmax; ++d[2])
                            disp[num++] = Key<NDIM>(0,d);
            }
            else if (NDIM == 4) {
                for (d[0]=-bmax; d[0]<=bmax; ++d[0])
                    for (d[1]=-bmax; d[1]<=bmax; ++d[1])
                        for (d[2]=-bmax; d[2]<=bmax; ++d[2])
                            for (d[3]=-bmax; d[3]<=bmax; ++d[3])
                                disp[num++] = Key<NDIM>(0,d);
            }
            else if (NDIM == 5) {
                for (d[0]=-bmax; d[0]<=bmax; ++d[0])
                    for (d[1]=-bmax; d[1]<=bmax; ++d[1])
                        for (d[2]=-bmax; d[2]<=bmax; ++d[2])
                            for (d[3]=-bmax; d[3]<=bmax; ++d[3])
                                for (d[4]=-bmax; d[4]<=bmax; ++d[4])

                                    disp[num++] = Key<NDIM>(0,d);
            }
            else if (NDIM == 6) {
                for (d[0]=-bmax; d[0]<=bmax; ++d[0])
                    for (d[1]=-bmax; d[1]<=bmax; ++d[1])
                        for (d[2]=-bmax; d[2]<=bmax; ++d[2])
                            for (d[3]=-bmax; d[3]<=bmax; ++d[3])
                                for (d[4]=-bmax; d[4]<=bmax; ++d[4])
                                    for (d[5]=-bmax; d[5]<=bmax; ++d[5])
                                        disp[num++] = Key<NDIM>(0,d);
            }
            else {
                MADNESS_EXCEPTION("_make_disp: hard dimension loop",NDIM);
            }

            std::sort(disp.begin(), disp.end(), cmp_keys);
        }

        static void make_disp_periodicsum(int bmax, Level n) {
            Translation twon = Translation(1)<<n;

            if (bmax > (twon-1)) bmax=twon-1;

            // Make permissible 1D translations
            Translation b[4*bmax+1];
            int i=0;
            for (Translation lx=-bmax; lx<=bmax; ++lx) {
                b[i++] = lx;
                if ((lx < 0) && (lx+twon > bmax)) b[i++] = lx + twon;
                if ((lx > 0) && (lx-twon <-bmax)) b[i++] = lx - twon;
            }
            MADNESS_ASSERT(i <= 4*bmax+1);
            int numb = i;

            disp_periodicsum[n] = std::vector< Key<NDIM> >();
            Vector<long,NDIM> lim(numb);
            for (IndexIterator index(lim); index; ++index) {
                Vector<Translation,NDIM> d;
                for (std::size_t i=0; i<NDIM; ++i) {
                    d[i] = b[index[i]];
                }
                disp_periodicsum[n].push_back(Key<NDIM>(n,d));
            }

            std::sort(disp_periodicsum[n].begin(), disp_periodicsum[n].end(), cmp_keys_periodicsum);
//             print("KEYS AT LEVEL", n);
//             print(disp_periodicsum[n]);
        }


    public:
        Displacements() {
            if (disp.size() == 0) {
                make_disp(bmax_default());

                if (NDIM <= 3) {
                    Level nmax = 8*sizeof(Translation) - 2;
                    for (Level n=0; n<nmax; ++n) make_disp_periodicsum(bmax_default(), n);
                }
            }
        }

        const std::vector< Key<NDIM> >& get_disp(Level n, bool isperiodicsum) {
            if (isperiodicsum) {
                MADNESS_ASSERT(NDIM <= 3);
                return disp_periodicsum[n];
            }
            else {
                return disp;
            }
        }

    };
}
#endif // MADNESS_MRA_DISPLACEMENTS_H__INCLUDED
