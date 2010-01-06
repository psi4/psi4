/*! \file
\ingroup INPUT
\brief Enter brief description of file here 
*/
#define EXTERN
#include <cstdio>
#include <libciomr/libciomr.h>
#include <cstdlib>
#include <cmath>
#include <libqt/qt.h>
#include "input.h"
#include "global.h"
#include "defines.h"
#include <physconst.h>
#include <masses.h>

namespace {
    void if_to_invert_axis(double* v1, int* must_invert, int* should_invert, double* maxproj);
}

namespace psi { namespace input {

    void reorient()
    {
        int i,j;
        int degen;
        int deg_IM1, deg_IM2;
        int nspher_set;
        int prox_i, prox_j;
        int axis, axis1, axis2, xyz;
        int nshould, nmust, must_invert[3], should_invert[3];
        double Xcm = 0.0;
        double Ycm = 0.0;
        double Zcm = 0.0;
        double mass = 0.0;
        double tmp, abs, rel;
        double Zmax, r, min_ij;
        double **IT, *IM, **ITAxes;
        double *v1, *v2, *v3, **R;
        double vabs, maxproj[3];
        double **sset_geom, *sset_dist;
        double origin[] = {0.0, 0.0, 0.0};

        double cos_theta, sin_theta, theta;
        double cos_phix, cos_phiy, phix, sin_phix;
        double v2norm;
        double **Rz, **Rx, **Rzt;

        const double im2rotconst = 0.25/(M_PI*_c_au*_bohr2cm);

        v1 = init_array(3);
        v2 = init_array(3);
        v3 = init_array(3);
        IT = block_matrix(3,3);
        IM = init_array(3);
        ITAxes = block_matrix(3,3);

        for(i=0;i<num_atoms;i++) {
            tmp = an2masses[(int)nuclear_charges[i]];
            Xcm += tmp*geometry[i][0];
            Ycm += tmp*geometry[i][1];
            Zcm += tmp*geometry[i][2];
            mass += tmp;
        }
        Xcm /= mass; Ycm /= mass; Zcm /= mass;

        if (!no_comshift) {
            for(i=0;i<num_allatoms;i++) {
                full_geom[i][0] -= Xcm;
                full_geom[i][1] -= Ycm;
                full_geom[i][2] -= Zcm;
            }
        }

        if (num_atoms > 1) {
            for(i=0;i<num_atoms;i++) {
                tmp = an2masses[(int)nuclear_charges[i]]/_au2amu;
                IT[0][0] += tmp*(geometry[i][1]*geometry[i][1] + geometry[i][2]*geometry[i][2]);
                IT[1][1] += tmp*(geometry[i][0]*geometry[i][0] + geometry[i][2]*geometry[i][2]);
                IT[2][2] += tmp*(geometry[i][0]*geometry[i][0] + geometry[i][1]*geometry[i][1]);
                IT[0][1] -= tmp*geometry[i][0]*geometry[i][1];
                IT[0][2] -= tmp*geometry[i][0]*geometry[i][2];
                IT[1][2] -= tmp*geometry[i][1]*geometry[i][2];
            }
            IT[1][0] = IT[0][1];
            IT[2][0] = IT[0][2];
            IT[2][1] = IT[1][2];
            sq_rsp(3,3,IT,IM,1,ITAxes,1.0E-14);
            IM[0] = fabs(IM[0]); /*Fixing numerical errors in the linear case*/
            fprintf(outfile,"\n  -Rotational constants (cm-1) :\n");
            if (IM[0] < ZERO_MOMENT_INERTIA) /* Linear molecule */
                fprintf(outfile,"    A = **********  ");
            else   /* Regular molecule */
                fprintf(outfile,"    A = %10.5lf  ",im2rotconst/IM[0]);
            if (IM[1] < ZERO_MOMENT_INERTIA)  /* Atom */
                fprintf(outfile,"B = **********  C = **********\n");
            else /* molecule */
                fprintf(outfile,"B = %10.5lf  C = %10.5lf\n",im2rotconst/IM[1],im2rotconst/IM[2]);

            degen = 0;
            for(i=0;i<2;i++)
            for(j=i+1;j<3 && degen<2;j++) {
                abs = fabs(IM[i] - IM[j]);
                tmp = (IM[i] > IM[j]) ? IM[i] : IM[j];
                if (abs > 1.0E-14)
                    rel = abs/tmp;
                else
                    rel = 0.0;
                if (rel < ZERO_MOMENT_INERTIA) {
                    degen++;
                    deg_IM1 = i;
                    deg_IM2 = j;
                }
            }

            v1[0] = ITAxes[0][1];
            v1[1] = ITAxes[1][1];
            v1[2] = ITAxes[2][1];
            v2[0] = ITAxes[0][2];
            v2[1] = ITAxes[1][2];
            v2[2] = ITAxes[2][2];
            cross_prod(v1,v2,v3);
            ITAxes[0][0] = v3[0];
            ITAxes[1][0] = v3[1];
            ITAxes[2][0] = v3[2];

            nmust = 0;
            nshould = 0;
            for(axis=0; axis<3; axis++) {

                v1[0] = ITAxes[0][axis];
                v1[1] = ITAxes[1][axis];
                v1[2] = ITAxes[2][axis];

                if_to_invert_axis(v1,&(must_invert[axis]),&(should_invert[axis]),&(maxproj[axis]));
                nmust += must_invert[axis];
                nshould += should_invert[axis];

            }

            R = block_matrix(3,3);
            if (nmust == 2) {
                for(axis=0; axis<3; axis++) {
                    if (must_invert[axis])
                        R[axis][axis] = -1.0;
                    else
                        R[axis][axis] = 1.0;
                }
            }
            else if (nmust == 1 && nshould > 0) {
                if (nshould == 2) {
                    for(axis=0; axis<3; axis++)
                    if (should_invert[axis]) {
                        axis1 = axis;
                        axis++;
                        break;
                    }
                    for(; axis<3; axis++)
                    if (should_invert[axis]) {
                        axis2 = axis;
                        break;
                    }
                    if (fabs(maxproj[axis1]) > fabs(maxproj[axis2])) {
                        nshould = 1;
                        should_invert[axis2] = 0;
                    }
                    else {
                        nshould = 1;
                        should_invert[axis1] = 0;
                    }
                }

                for(axis=0; axis<3; axis++) {
                    if (must_invert[axis])
                        R[axis][axis] = -1.0;
                    else if (should_invert[axis])
                        R[axis][axis] = -1.0;
                    else
                        R[axis][axis] = 1.0;
                }
            }
            else if (nmust == 3) {
                R[0][0] = -1.0;
                R[1][1] = -1.0;
                R[2][2] = 1.0;
            }
            else if (nmust == 0 && nshould > 1) {
                if (nshould == 3) {
                    tmp = fabs(maxproj[0]);
                    i = 0;
                    for(axis=1; axis<3; axis++) {
                        if (fabs(maxproj[axis]) < fabs(tmp)) {
                            i = axis;
                            tmp = fabs(maxproj[axis]);
                        }
                    }
                    should_invert[i] = 0;
                    nshould = 2;
                }
                for(axis=0; axis<3; axis++) {
                    if (should_invert[axis])
                        R[axis][axis] = -1.0;
                    else
                        R[axis][axis] = 1.0;
                }
            }
            else {
                R[0][0] = 1.0;
                R[1][1] = 1.0;
                R[2][2] = 1.0;
            }

            if (degen == 0 && !no_reorient) {
                rotate_full_geom(ITAxes);
                rotate_full_geom(R);
            }
            free_block(R);

            if (degen == 1 && !no_reorient) {

                int must_invert, should_invert, unique_axis;
                double maxproj, invert_pfac;

                if (deg_IM1 + deg_IM2 == 3)
                    unique_axis = 0;
                else
                    unique_axis = 2;

                v1[0] = ITAxes[0][unique_axis];
                v1[1] = ITAxes[1][unique_axis];
                v1[2] = ITAxes[2][unique_axis];

                if_to_invert_axis(v1, &must_invert, &should_invert, &maxproj);
                if (must_invert || should_invert)
                    invert_pfac = -1.0;
                else
                    invert_pfac = 1.0;

                v1[0] *= invert_pfac;
                v1[1] *= invert_pfac;
                v1[2] *= invert_pfac;

                cos_theta = v1[2];
                if ( (1.0 - fabs(cos_theta)) > ZERO_MOMENT_INERTIA) {
                    theta = acos(cos_theta);
                    sin_theta = sin(theta);

                    v3[0] = 0.0;
                    v3[1] = 0.0;
                    v3[2] = 1.0;
                    cross_prod(v1,v3,v2);
                    v2norm = sqrt(v2[0] * v2[0] +
                        v2[1] * v2[1] +
                        v2[2] * v2[2]);
                    v2[0] /= v2norm;
                    v2[1] /= v2norm;
                    v2[2] /= v2norm;

                    cos_phix = v2[0];
                    cos_phiy = v2[1];
                    phix = acos(cos_phix);

                    if (cos_phiy > 0.0) {
                        phix *= -1.0; 
                    }
                    sin_phix = sin(phix);

                    Rz = block_matrix(3,3);
                    Rz[2][2] = 1.0;
                    Rz[0][0] = cos_phix;
                    Rz[1][1] = cos_phix;
                    Rz[0][1] = sin_phix;
                    Rz[1][0] = -sin_phix;
                    rotate_full_geom(Rz);
                    free_block(Rz);

                    Rx = block_matrix(3,3);
                    Rx[0][0] = 1.0;
                    Rx[1][1] = cos_theta;
                    Rx[2][2] = cos_theta;
                    Rx[1][2] = sin_theta;
                    Rx[2][1] = -sin_theta;
                    rotate_full_geom(Rx);
                    free_block(Rx);

                    Rzt = block_matrix(3,3);
                    Rzt[2][2] = 1.0;
                    Rzt[0][0] = cos_phix;
                    Rzt[1][1] = cos_phix;
                    Rzt[0][1] = -sin_phix;
                    Rzt[1][0] = sin_phix;
                    rotate_full_geom(Rzt);
                    free_block(Rzt);
                }
            }

            if (degen == 0) {
                fprintf(outfile,"    It is an asymmetric top.\n");
                rotor = asymmtop;
            }
            else if (degen == 1)
            switch (deg_IM1 + deg_IM2) {
                case 3:
                if (IM[0] < ZERO_MOMENT_INERTIA) {
                    fprintf(outfile,"    It is a linear molecule.\n");
                    rotor = linear;
                }
                else {
                    fprintf(outfile,"    It is a prolate symmetric top.\n");
                    rotor = symmtop;
                }
                break;
                case 1:
                fprintf(outfile,"    It is an oblate symmetric top.\n");
                rotor = symmtop;
                break;
            }
            else if (degen == 2) {
                fprintf(outfile,"    It is a spherical top.\n");
                rotor = sphtop;
                Zmax = 0;
                for(i=0;i<num_atoms;i++) 
                    if (Zmax < (int)nuclear_charges[i] && sqrt(dot_prod(geometry[i],geometry[i])) > ZERO)
                    Zmax = (int)nuclear_charges[i];
                r = 0.0;
                for(i=0;i<num_atoms && r < ZERO;i++)
                    r = (nuclear_charges[i] == Zmax) ? sqrt(dot_prod(geometry[i],geometry[i])) : 0.0;
                sset_geom = init_matrix(num_atoms,3);
                sset_geom[0][0] = geometry[i-1][0];
                sset_geom[0][1] = geometry[i-1][1];
                sset_geom[0][2] = geometry[i-1][2];
                nspher_set = 1;
                for(j=i;j<num_atoms;j++)
                if (nuclear_charges[j] == Zmax) {
                    tmp = sqrt(dot_prod(geometry[i],geometry[i]));
                    if ( fabs(tmp - r) < ZERO) {
                        sset_geom[nspher_set][0] = geometry[j][0];
                        sset_geom[nspher_set][1] = geometry[j][1];
                        sset_geom[nspher_set][2] = geometry[j][2];
                        nspher_set++;
                    }
                }
                if (!no_reorient)
                switch(nspher_set) {
                    case 4:
                    median_vec(sset_geom[0], sset_geom[1], v1);
                    median_vec(sset_geom[1], sset_geom[2], v2);
                    cross_prod(v1, v2, v3);
                    vectors_to_matrix(v1, v2, v3, ITAxes);
                    rotate_full_geom(ITAxes);
                    break;

                    case 6:
                    if (!inv_related(sset_geom[0], sset_geom[1]))
                        median_vec(sset_geom[0], sset_geom[1], v1);
                    else 
                        median_vec(sset_geom[0], sset_geom[2], v1);
                    unit_vec(sset_geom[0],origin,v2);
                    cross_prod(v1, v2, v3);
                    unit_vec(v3,origin,v3);
                    cross_prod(v3, v1, v2);
                    vectors_to_matrix(v1, v2, v3, ITAxes);
                    rotate_full_geom(ITAxes);
                    break;

                    case 8:
                    sset_dist = init_array(nspher_set*(nspher_set+1)/2);
                    calc_distance(sset_geom,sset_dist,nspher_set);
                    min_ij = sset_dist[ioff[1]+0];
                    prox_i = 1;
                    for(i=2;i<nspher_set-2;i++)
                    if (min_ij > (tmp = sset_dist[ioff[i]+0]) && !inv_related(sset_geom[0],sset_geom[i])) {
                        min_ij = tmp;
                        prox_i = i;
                        break;
                    }
                    for(j=prox_i+1;j<nspher_set;j++)
                    if (fabs(min_ij - sset_dist[ioff[j]+0]) < ZERO) {
                        prox_j = j;
                        break;
                    }
                    unit_vec(sset_geom[0], sset_geom[prox_i], v1);
                    unit_vec(sset_geom[0], sset_geom[prox_j], v2);
                    cross_prod(v1, v2, v3);
                    vectors_to_matrix(v1, v2, v3, ITAxes);
                    rotate_full_geom(ITAxes);
                    free(sset_dist);
                    break;

                    default:
                    break;
                }
                free_matrix(sset_geom,num_atoms);
            }
        }
        else if (num_atoms == 1) {
            fprintf(outfile,"    It is a spherical top.\n");
            rotor = atom;
        }
        else if (num_atoms <= 0) {
            punt("Fewer than 1 atom");
        }

        free(v1);
        free(v2);
        free(v3);
        free(IM);
        free_block(IT);
        free_block(ITAxes);
    }

}} // namespace psi::input

namespace {
    void if_to_invert_axis(double* v1, int* must_invert, int* should_invert, double* maxproj)
    {
        using namespace psi::input;
        int xyz, nzero;
        double vabs;

        *maxproj = 0.0;
        *must_invert = 0;
        *should_invert = 0;

        nzero = 0;

        for(xyz=0; xyz<3; xyz++) {

            vabs = fabs(v1[xyz]);

            if (vabs < ZERO)
                nzero++;

            if (vabs > fabs(*maxproj)) {
                *maxproj = v1[xyz];
            }

        }

        if (nzero == 2) {
            if (*maxproj < 0.0)
                *must_invert = 1;
            else
                *must_invert = 0;
        }
        else if (nzero < 2) {
            if (*maxproj < 0.0)
                *should_invert = 1;
            else
                *should_invert = 0;
        }
    }
} // namespace
