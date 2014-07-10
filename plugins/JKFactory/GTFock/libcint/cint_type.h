#ifndef __CINT_TYPE_H__
#define __CINT_TYPE_H__


struct OED
{
    int nalpha;
    int ncoeff;
    int ncgto1;
    int ncgto2;
    int npgto1;
    int npgto2;
    int shell1;
    int shell2;
    int natoms;
    int ncsum;
    int spheric;
    int screen;  
    double x1;
    double y1;
    double z1;
    double x2;
    double y2;
    double z2;
    double *xn;
    double *yn;
    double *zn;
    double *ncharge;
    double *cc;
    double *alpha;
    int cc_beg[2];
    int cc_end[2];
    int imax;
    int zmax;
    double *zcore;
    double *zcore2;
    int *icore;

    int fp_memory_opt;
    int int_memory_opt;
    int *coef_offset;
    int *exp_offset;
};


struct ERD
{
    int shell1;
    int shell2;
    int shell3;
    int shell4;
    int ncoef;
    int nalpha;
    int ncsum;
    int ncgto1;
    int ncgto2;
    int ncgto3;
    int ncgto4;
    int npgto1;
    int npgto2;
    int npgto3;
    int npgto4;
    int imax;
    int zmax;
    int spheric;
    int screen;
    int cc_beg[4];
    int cc_end[4];
    double *cc;
    double *alpha;
    double x1;
    double y1;
    double z1;
    double x2;
    double y2;
    double z2;
    double x3;
    double y3;
    double z3;
    double x4;
    double y4;
    double z4;
    int *icore;
    double *zcore;

    int fp_memory_opt;
    int int_memory_opt;
    int *coef_offset;
    int *exp_offset;
};


struct BasisSet
{
    //Number of atoms
    int natoms;
    //Atomic numbers (indexed from 0 to natoms-1), actual Zs, not Z-1, i.e. if atom i is H, then eid[i]=1
    int *eid;
    //These are the x,y,z coordinates of the atoms (indexed from 0 to natoms-1) and in bohr
    double *xn;
    double *yn;
    double *zn;
    //These are the nuclear charges of the atoms (indexed from 0 to natoms-1), really just a double cast of eid
    double *ncharge;
    //The total number of electrons, obtained by summing eid
    int nelectrons;

    //*************************These are details for the entire basis set, not just that of the molecule****************//

    //This is set to a MACRO value (either CARTESIAN=0 or SPHERICAL=1)
    int basistype;
    //This is the total number of atoms this basis set is defined for (e.g. the number of atoms 6-31G* is defined for, not this object)
    int lenatom0;
    /* This is an ELEN (ELEN is the number of elements that have been hard coded in GTFock) long array
     * that maps Z to one of the 0 to lenatom0-1 atoms this basis is defined for.   eptr[Z]=-1 if element Z
     * is not supported in this basis set.  NOte the indexes are Z and Z-1
     */
    int *eptr;
    ///This is an lenatom0+1 long array that keeps track of where shells start for a given atom
    ///it is indexed from 0 to lenatom0-1, the returned shell is indexed from 0 to lenshell0-1
    ///The last element in the array is the total number of shells in the basis, lenshell0
    int *atom_start;
    //Total number of shells this basis set is defined for (e.g. the number of shells 6-31G* is defined for, not this object)
    int lenshell0;
    //Total number of primitives this basis set is defined for (e.g. the number of primitives 6-31G* is defined for, not this object)
    int totnexp;
    //An lenshell0+1 long array that says where a shell's primitives start (both shells and primitives are indexed from 0)
    //The last element is the total number of primitives in the array
    int *ptrshell;
    //An totnexp long array that contains the number of primitives in the current shell
    int *nexp0;
    //An totnexp long array of the primitive exponents (indexed 0 to totnexp-1)
    double *exp0;
    //An totnexp long array of the primitive expansion coefficients (indexed from 0 to totnexp-1)
    double *cc0;
    //An lenshell0 long array of each shell's momentum (shells are indexed 0...lenshell0-1), momenta are s=0,p=1,etc.
    int *momentum0;
    
    //************************************These are details for the AO basis set of the molecule***********************//

    //Total number of shells for the molecule
    int nshells;    
    //Total number of basis functions for the molecule
    int nfunctions;
    //An nshells long array, where shell i's (indexed from 0 to nshells-1) basis functions start (indexed from 0 to nfunctions-1)
    int *f_start_id;
    //An nshells long array of where shell i's (indexed from 0 to nshells-1) basis functions end (indexed form 0 to nfunctions-1)
    int *f_end_id;
    //An natoms+1 long array, where element i is where atom i's (indexed from 0 to natoms-1) shells start (indexed from 0 to nshells-1)
    //element natoms is the total number of shells for the molecule
    int *s_start_id;
    //An nshells long array of the number of primitives in a shell (shells are indexed 0 to nshells-1)
    int *nexp;
    //An nshells long array where exp[i] is the starting address of exp0 for shell i(indexed from 0)
    double **exp;
    //An nshells long array where cc[i] is the starting address of cc0 for shell i (indexed from 0)
    double **cc;
    //An nshells long array where momentum[i] is the momentum of shell i (indexed from 0 to nshells-1)
    int *momentum;
    //nshells long arrays where q[i] is the shell's center's location in the q-th spatial dimension (shells are indexed starting at 0)
    double *x;
    //nshells long arrays where q[i] is the shell's center's location in the q-th spatial dimension (shells are indexed starting at 0)
    double *y;
    //nshells long arrays where q[i] is the shell's center's location in the q-th spatial dimension (shells are indexed starting at 0)
    double *z;

    int maxdim; // max number of functions of a shell 
    int max_momentum;//max momentum of the molecule
    int max_nexp;//maximum number of primitives in a shell
    int max_nexp_id;//shell number with that many primitives (indexed from 0)
    
    char str_buf[512];  
};


#endif /* __CINT_TYPE_H__ */
