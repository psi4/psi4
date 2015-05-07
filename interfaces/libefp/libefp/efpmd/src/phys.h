#ifndef EFPMD_PHYS_H
#define EFPMD_PHYS_H

/* Bohr radius in angstroms */
#define BOHR_RADIUS 0.52917721092

/* Boltzmann constant in [Hartree / K] */
#define BOLTZMANN 3.166811429e-6

/* Femtoseconds to atomic units of time conversion */
#define FS_TO_AU (1.0 / 2.41888432650212e-2)

/* AMU to atomic units of mass conversion */
#define AMU_TO_AU (1.0 / 5.485799094622e-4)

/* Hertree energy in Joules */
#define HARTREE 4.35974434e-18

/* Bar to atomic units of pressure */
#define BAR_TO_AU (1.0e-25 * BOHR_RADIUS * BOHR_RADIUS * BOHR_RADIUS / HARTREE)

/* Fine structure constant */
#define FINE_CONST 7.297352569824e-3

#endif /* EFPMD_PHYS_H */
