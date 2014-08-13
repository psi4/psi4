#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>
#include <string.h>
#include <ctype.h>

#include "config.h"
#include "basisset.h"


#define ELEN         50
#define SLEN         5
#define MAXNS        3
#define MAXATOMNAME  2
#define A2BOHR       1.889726
#define CARTESIAN    0
#define SPHERICAL    1


static char etable[ELEN][MAXATOMNAME + 1] =
{
  "H",  "He", "Li", "Be", "B",
  "C",  "N",  "O",  "F",  "Ne",
  "Na", "Mg", "Al", "Si", "P",
  "S",  "Cl", "Ar", "K",  "Ca",
  "Sc", "Ti", "V",  "Cr", "Mn",
  "Fe", "Co", "Ni", "Cu", "Zn",
  "Ga", "Ge", "As", "Se", "Br",
  "Kr", "Rb", "Sr", "Y",  "Zr",
  "Nb", "Mo", "Tc", "Ru", "Rh",
  "Pd", "Ag", "Cd", "In", "Sn"
};

static char mtable[SLEN] =
{
  'S',  'P', 'D', 'F', 'G' 
};


static void normalization (BasisSet_t basis)
{
    double sum;
    double temp;
    double temp2;
    double temp3;
    double xnorm;
    double a1;
    double a2;
    int i;
    int j;
    int k;
    int offset;

    for (i = 0; i < basis->lenshell0; i++)
    {
        sum = 0.0;
        offset = basis->ptrshell[i];
        for (j = 0; j < basis->nexp0[i]; j++)
        {
            for (k = 0; k <= j; k++)
            {
                a1 = basis->exp0[offset + j];
                a2 = basis->exp0[offset + k];
                temp = basis->cc0[offset + j] * basis->cc0[offset + k]; 
                temp2 = basis->momentum0[i] + 1.5;
                temp3 = 2.0 * sqrt (a1 * a2) / (a1 + a2);
                temp3 = pow (temp3, temp2);
                temp = temp * temp3;
                sum = sum + temp;
                if (j != k)
                {
                    sum = sum + temp;
                }
            }
        }
        xnorm = 1.0 / sqrt (sum);
        for (j = 0; j < basis->nexp0[i]; j++)
        {
            basis->cc0[offset + j] *= xnorm;
        }
    }
}


void _maxMomentum (BasisSet_t basis, int *max_momentum)
{
    *max_momentum = basis->max_momentum;
}


void _maxPrimid (BasisSet_t basis, int *max_primid)
{
    *max_primid = basis->max_nexp_id;
}


void _maxnumExp (BasisSet_t basis, int *max_nexp)
{
    *max_nexp = basis->max_nexp;   
}


CIntStatus_t CInt_createBasisSet (BasisSet_t *_basis)
{
    BasisSet_t basis;
    basis = (BasisSet_t )malloc (sizeof(struct BasisSet));
    if (NULL == basis)
    {
        CINT_PRINTF (1, "memory allocation failed\n");
        return CINT_STATUS_ALLOC_FAILED;
    }
    memset (basis, 0, sizeof(struct BasisSet));

    *_basis = basis;
    return CINT_STATUS_SUCCESS;
}


static CIntStatus_t import_molecule (char *file, BasisSet_t basis)
{
    FILE *fp;
    char line[1024];
    char str[1024];
    int natoms;
    int nelectrons;
    int nsc;
    int i;

    fp = fopen (file, "r");
    if (fp == NULL)
    {
        CINT_PRINTF (1, "failed to open molecule file %s\n", file);
        return CINT_STATUS_FILEIO_FAILED;
    }
    if (fgets (line, 1024, fp) == NULL)
    {
        CINT_PRINTF (1, "file %s has a wrong format\n", file);
        return CINT_STATUS_FILEIO_FAILED;
    }

    // number of atoms
    sscanf (line, "%d", &(basis->natoms));
    basis->xn = (double *)malloc (sizeof(double) * basis->natoms);
    basis->yn = (double *)malloc (sizeof(double) * basis->natoms);
    basis->zn = (double *)malloc (sizeof(double) * basis->natoms);
    basis->ncharge = (double *)malloc (sizeof(double) * basis->natoms); 
    basis->eid = (int *)malloc (sizeof(int) * basis->natoms);
    if (NULL == basis->xn ||
        NULL == basis->yn ||
        NULL == basis->zn ||
        NULL == basis->ncharge ||
        NULL == basis->eid)
    {
        CINT_PRINTF (1, "memory allocation failed\n");
        return CINT_STATUS_ALLOC_FAILED;
    }

    // read x, y and z
    natoms = 0;
    nelectrons = 0;
    while (fgets (line, 1024, fp) != NULL)
    {
        nsc = sscanf (line, "%s %lf %lf %lf",
                      str, &(basis->xn[natoms]), 
                      &(basis->yn[natoms]), &(basis->zn[natoms]));
        basis->xn[natoms] = basis->xn[natoms] * A2BOHR;
        basis->yn[natoms] = basis->yn[natoms] * A2BOHR;
        basis->zn[natoms] = basis->zn[natoms] * A2BOHR;   
        if (strlen(str) > MAXATOMNAME || nsc == EOF)
        {
            CINT_PRINTF (1, "atom %s in %s is not supported\n", str, file);
            return CINT_STATUS_INVALID_VALUE;
        }
        for (i = 0; i < ELEN; i++)
        {
            if (strcmp (str, etable[i]) == 0)
            {
                basis->eid[natoms] = i + 1;
                break;
            }
        }
        if (i == ELEN)
        {
            CINT_PRINTF (1, "atom %s is not supported\n", str);
            return CINT_STATUS_INVALID_VALUE;
        }
        basis->ncharge[natoms] = (double)(basis->eid[natoms]);
        nelectrons += basis->eid[natoms];
        natoms++;
    }
    basis->nelectrons = nelectrons;
    if (natoms != basis->natoms)
    {
        CINT_PRINTF (1, "file %s natoms %d does not match the header\n",
            file, natoms);
        return CINT_STATUS_FILEIO_FAILED;
    }
    
    fclose (fp);
    
    return CINT_STATUS_SUCCESS;
}


static CIntStatus_t import_basis (char *file, BasisSet_t basis)
{
    FILE *fp;
    char line[1024];
    char str[1024];
    int natoms;
    int nshells;
    int i;
    int j;
    int nexp;
    int ns;
    double beta;
    double cc[MAXNS];
    long int mark;
    int totnexp;

    fp = fopen (file, "r");
    if (fp == NULL)
    {
        CINT_PRINTF (1, "failed to open molecule file %s\n", file);
        return CINT_STATUS_FILEIO_FAILED;
    }

    // read the basis type
    if (fgets (line, 1024, fp) == NULL)
    {
        CINT_PRINTF (1, "file %s has a wrong format\n", file);
        return CINT_STATUS_FILEIO_FAILED;    
    }
    sscanf (line, "%s", str);
    if (strcmp (str, "cartesian") == 0)
    {
        basis->basistype = CARTESIAN;
    }
    else if (strcmp (str, "spherical") == 0)
    {
        basis->basistype = SPHERICAL;
    }

    // get number of atoms
    natoms = 0;
    nshells = 0;
    totnexp = 0;
    while (fgets (line, 1024, fp) != NULL)
    {
        if (isalpha (line[0]))
        {
            // a new atom
            natoms++;
            while (fgets (line, 1024, fp) != NULL)
            {
                if (isalpha (line[0]))
                {
                    sscanf (line, "%s %d %lf",
                        str, &nexp, &beta);
                    ns = strlen (str);               
                    nshells += ns;
                    totnexp += ns * nexp;
                }
                if (line[0] == '*')
                {
                    break;
                }
            }
         }
    }
    //Total number of atoms this basis set is defined for
    basis->lenatom0 = natoms;
    //Total number of shells the basis set is defined for
    basis->lenshell0 = nshells;
    //Total number of primitive gaussians in the basis set
    basis->totnexp = totnexp;
    //RMR added the +1 because otherwise you are not allocating enough space since the index==Z and not Z-1
    basis->eptr = (int *)malloc (sizeof(int) * (ELEN+1));
    //This is an array of where the shells for this atom starts atom_start[i] is where the i'th element starts (i from 0)
    //presumably this would be called something like atom_start[eptr[i]] to get where the shell for atom i starts
    //The last element is the the number of shells
    basis->atom_start = (int *)malloc (sizeof(int) * (natoms + 1));
    //nexp0[i] is the number of primitives in shell i, (i from 0)
    basis->nexp0 = (int *)malloc (sizeof(int) * nshells);
    //These are the coefficients of the primitives
    basis->cc0 = (double *)malloc (sizeof(double) * totnexp);
    //These are the arguments of the exponentials
    basis->exp0 = (double *)malloc (sizeof(double) * totnexp);
    //This is the momentum of the current shell
    basis->momentum0 = (int *)malloc (sizeof(int) * nshells);
    //ptrshell[i] is where primitives in shell i start (i from 0), the last element is the total number of primitives
    basis->ptrshell = (int *)malloc (sizeof(int) * (nshells + 1)); 
    if (NULL == basis->atom_start ||
        NULL == basis->eptr ||
        NULL == basis->nexp0 ||
        NULL == basis->cc0 ||
        NULL == basis->exp0 ||
        NULL == basis->momentum0 ||
        NULL == basis->ptrshell)
    {
        CINT_PRINTF (1, "memory allocation failed\n");
        return CINT_STATUS_ALLOC_FAILED;    
    }
    for (i = 0; i < ELEN; i++)
    {
        basis->eptr[i] = -1;
    }
    
    // get nshells
    rewind (fp);
    fgets (line, 1024, fp);
    natoms = 0;
    nshells = 0;
    totnexp = 0;
    while (fgets (line, 1024, fp) != NULL)
    {
        if (isalpha (line[0]))
        {
            // a new atom
            sscanf (line, "%s", str);
            for (i = 0; i < ELEN; i++)
            {
                if (strcmp (str, etable[i]) == 0)
                {
                    basis->eptr[i + 1] = natoms;
                    break;
                }
            }
            if (i == ELEN)
            {
                CINT_PRINTF (1, "atom %s in %s is not supported\n", str, file);
                return CINT_STATUS_INVALID_VALUE;
            }
            basis->atom_start[natoms] = nshells;           
            natoms++;
            // read shells
            while (fgets (line, 1024, fp) != NULL)
            {
                if (isalpha (line[0]))
                {
                    sscanf (line, "%s %d %lf",
                        str, &nexp, &beta);
                    ns = strlen (str);
                    if (nexp <= 0 || ns <= 0 || ns > MAXNS)
                    {
                        CINT_PRINTF (1, "file %s contains invalid values\n", file);
                        return CINT_STATUS_INVALID_VALUE;                        
                    }
                    mark = ftell (fp);
                    for (i = 0; i < ns; i++)//If str was something like SP loop over S and then P
                    {
                        basis->nexp0[nshells] = nexp;
                        basis->ptrshell[nshells] = totnexp;
                        for (j = 0; j < SLEN; j++)//loop over recognized angular momentum numbers
                        {
                            if (str[i] == mtable[j])
                            {
                                basis->momentum0[nshells] = j;
                                break;
                            }
                        }//End loop over momenta
                        if (j == SLEN)
                        {
                            CINT_PRINTF (1, "shell %s in file %s is not supported\n",
                                str, file);
                            return CINT_STATUS_INVALID_VALUE;  
                        }
                        fseek (fp, mark, SEEK_SET);
                        for (j = 0; j < basis->nexp0[nshells]; j++)//Loop over primitives in the shell
                        {
                            if (fgets (line, 1024, fp) == NULL ||
                                line[0] == '*' ||
                                isalpha (line[0]))
                            {
                                CINT_PRINTF (1, "file %s has a wrong format\n", file);
                                return CINT_STATUS_FILEIO_FAILED;
                            }
                            sscanf (line, "%lf %lf %lf %lf",
                                    &(basis->exp0[totnexp + j]),
                                    &(cc[0]), &(cc[1]), &(cc[2]));
                            basis->cc0[totnexp + j] = cc[i];
                        }//End loop over primitives
                        totnexp += basis->nexp0[nshells];
                        nshells++;
                    }//End loop over shells in SP, for example
                }//End current shell
                if (line[0] == '*')
                {
                    break;
                }
            }
         }
    }
    basis->atom_start[natoms] = nshells;
    basis->ptrshell[nshells] = totnexp;
    if (natoms != basis->lenatom0 ||
        nshells != basis->lenshell0 ||
        basis->totnexp != totnexp)
    {
        CINT_PRINTF (1, "file %s has a wrong format\n", file);
        return CINT_STATUS_FILEIO_FAILED;    
    }

    fclose (fp);

    normalization (basis);
    
    return CINT_STATUS_SUCCESS;
}


static CIntStatus_t parse_molecule (BasisSet_t basis)
{
    int natoms;
    int nshells;   
    int nfunctions;
    int maxdim;
    int max_momentum;
    int max_nexp;
    int max_nexp_id;
    int i;
    int j;
    int eid;
    int atom_start;
    int atom_end;
    int offset;

    // get lengths
    natoms = basis->natoms;
    nshells = 0;
    for (i = 0; i < natoms; i++)
    {
        eid = basis->eid[i];//Atomic number with no offset (i.e. H==1, He==2,...)
        atom_start = basis->atom_start[basis->eptr[eid]];
        atom_end = basis->atom_start[basis->eptr[eid] + 1];
        nshells += atom_end - atom_start;
    }
    //Array of where an atom's shells start, last element total number of shells
    basis->s_start_id = (int *)malloc (sizeof(int) * (natoms + 1));
    //Array of where each shell starts, in terms of basis functions not primitives
    basis->f_start_id = (int *)malloc (sizeof(int) * nshells);
    //Same except where it ends
    basis->f_end_id = (int *)malloc (sizeof(int) * nshells);
    //cart coordinates of the center for each shell
    basis->x = (double *)malloc (sizeof(double) * nshells);
    basis->y = (double *)malloc (sizeof(double) * nshells);
    basis->z = (double *)malloc (sizeof(double) * nshells);
    //Number of primitives in the shell
    basis->nexp = (int *)malloc (sizeof(int) * nshells);
    //Expansion coefficients of the shell
    basis->cc = (double **)malloc (sizeof(double *) * nshells);
    //Exponent scale factor for each shell
    basis->exp = (double **)malloc (sizeof(double *) * nshells);
    //Momentum of each shell
    basis->momentum = (int *)malloc (sizeof(int) * nshells);   
    if (NULL == basis->f_start_id ||
        NULL == basis->f_end_id ||
        NULL == basis->s_start_id ||
        NULL == basis->x ||
        NULL == basis->y ||
        NULL == basis->z ||
        NULL == basis->nexp ||
        NULL == basis->momentum ||
        NULL == basis->cc ||
        NULL == basis->exp)
    {
        CINT_PRINTF (1, "memory allocation failed\n");
        return CINT_STATUS_ALLOC_FAILED;    
    }
    basis->nshells = nshells;

    // parse molecules
    nshells = 0;
    nfunctions = 0;
    maxdim = 0;
    max_momentum = 0;
    max_nexp = 0;
    max_nexp_id = 0;
    for (i = 0; i < natoms; i++)
    {
        eid = basis->eid[i]; 
        atom_start = basis->atom_start[basis->eptr[eid]];
        atom_end = basis->atom_start[basis->eptr[eid] + 1];
        basis->s_start_id[i] = nshells;
        //Loop over shells in atom i
        for (j = atom_start; j < atom_end; j++)
        {
            basis->f_start_id[nshells + j - atom_start] = nfunctions;
            basis->nexp[nshells + j - atom_start] = basis->nexp0[j];
            basis->x[nshells + j - atom_start] = basis->xn[i];
            basis->y[nshells + j - atom_start] = basis->yn[i];
            basis->z[nshells + j - atom_start] = basis->zn[i];
            basis->momentum[nshells + j - atom_start] = basis->momentum0[j];
            //Determine the highest angular momentum of the basis
            max_momentum = (max_momentum > basis->momentum0[j] ?
                max_momentum : basis->momentum0[j]);
            //Determine which shell has the most primitives
            if (max_nexp < basis->nexp0[j])
            {
                max_nexp  = basis->nexp0[j];
                max_nexp_id = nshells + j - atom_start;
            }
            offset = basis->ptrshell[j];
            basis->cc[nshells + j - atom_start] = &(basis->cc0[offset]);
            basis->exp[nshells + j - atom_start] = &(basis->exp0[offset]);
            if (basis->basistype == SPHERICAL)
            {
                nfunctions += 2 * basis->momentum0[j] + 1;
                maxdim = (2 * basis->momentum0[j] + 1) > maxdim ?
                    (2 * basis->momentum0[j] + 1) : maxdim;
            }
            else if (basis->basistype == CARTESIAN)
            {
                nfunctions += (basis->momentum0[j] + 1)*(basis->momentum0[j] + 2)/2;
                maxdim = ((basis->momentum0[j] + 1)*(basis->momentum0[j] + 2)/2) > maxdim ?
                    ((basis->momentum0[j] + 1)*(basis->momentum0[j] + 2)/2) : maxdim;
            }
            basis->f_end_id[nshells + j - atom_start] = nfunctions - 1;
        }
        nshells += atom_end - atom_start;
    }
    basis->s_start_id[natoms] = nshells;
    //This is the maximum number of orbitals any shell has
    basis->maxdim = maxdim;
    //Total number of orbitals
    basis->nfunctions = nfunctions;
    basis->max_momentum = max_momentum;
    basis->max_nexp = max_nexp;
    basis->max_nexp_id = max_nexp_id;
    
    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_loadBasisSet (BasisSet_t basis, char *bsfile, char *molfile)
{
    CIntStatus_t ret;

    // read xyz file
    if ((ret = import_molecule (molfile, basis)) != CINT_STATUS_SUCCESS)
    {
        return ret;
    }
    // read basis set
    if ((ret = import_basis (bsfile, basis)) != CINT_STATUS_SUCCESS)
    {
        return ret;
    }
    //parse xyz
    if ((ret = parse_molecule (basis)) != CINT_STATUS_SUCCESS)
    {
        return ret;
    }

    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_destroyBasisSet (BasisSet_t basis)
{
    free (basis->f_start_id);
    free (basis->f_end_id);
    free (basis->xn);
    free (basis->yn);
    free (basis->zn);
    free (basis->ncharge);
    free (basis->x);
    free (basis->y);
    free (basis->z);
    
    free (basis->eid);
    free (basis->eptr);
    free (basis->atom_start);
    free (basis->ptrshell);
    free (basis->nexp0);
    free (basis->exp0);
    free (basis->cc0);
    free (basis->momentum0);
    
    free (basis->s_start_id);
    free (basis->nexp);
    free (basis->exp);
    free (basis->cc);
    free (basis->momentum);

    free (basis);

    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_packBasisSet (BasisSet_t basis,
                                void **buf,
                                int *bufsize)
{
    int _bufsize;
    char *_buf;
    int offset;
    
    _bufsize = 6 * sizeof(int) + 4 * basis->natoms * sizeof(double) +                
               (3 * basis->lenshell0 + ELEN +
                basis->lenatom0 + basis->natoms + 2) * sizeof(int) +
                basis->totnexp * 2 * sizeof(double);
    _buf = (char *)malloc (_bufsize);
    assert (_buf != NULL);
    offset = 0;    
    memcpy (&(_buf[offset]), &(basis->natoms), sizeof(int));
    offset += sizeof(int);
    memcpy (&(_buf[offset]), &(basis->nelectrons), sizeof(int));
    offset += sizeof(int);
    memcpy (&(_buf[offset]), &(basis->basistype), sizeof(int));
    offset += sizeof(int);
    memcpy (&(_buf[offset]), &(basis->lenshell0), sizeof(int));
    offset += sizeof(int);
    memcpy (&(_buf[offset]), &(basis->totnexp), sizeof(int));
    offset += sizeof(int);
    memcpy (&(_buf[offset]), &(basis->lenatom0), sizeof(int));
    offset += sizeof(int);

    memcpy (&(_buf[offset]), basis->xn, sizeof(double) * basis->natoms);
    offset += sizeof(double) * basis->natoms;
    memcpy (&(_buf[offset]), basis->yn, sizeof(double) * basis->natoms);
    offset += sizeof(double) * basis->natoms;
    memcpy (&(_buf[offset]), basis->zn, sizeof(double) * basis->natoms);
    offset += sizeof(double) * basis->natoms;
    memcpy (&(_buf[offset]), basis->ncharge, sizeof(double) * basis->natoms);
    offset += sizeof(double) * basis->natoms;
    memcpy (&(_buf[offset]), basis->eid, sizeof(int) * basis->natoms);
    offset += sizeof(int) * basis->natoms;

    memcpy (&(_buf[offset]), basis->exp0, sizeof(double) * basis->totnexp);
    offset += sizeof(double) * basis->totnexp;
    memcpy (&(_buf[offset]), basis->cc0, sizeof(double) * basis->totnexp);
    offset += sizeof(double) * basis->totnexp;

    memcpy (&(_buf[offset]), basis->eptr, sizeof(int) * ELEN);
    offset += sizeof(int) * ELEN;
    memcpy (&(_buf[offset]), basis->atom_start, sizeof(int) * (basis->lenatom0 + 1));
    offset += sizeof(int) * (basis->lenatom0 + 1);
    
    memcpy (&(_buf[offset]), basis->momentum0, sizeof(int) * basis->lenshell0);
    offset += sizeof(int) * basis->lenshell0;
    memcpy (&(_buf[offset]), basis->nexp0, sizeof(int) * basis->lenshell0);
    offset += sizeof(int) * basis->lenshell0;
    memcpy (&(_buf[offset]), basis->ptrshell, sizeof(int) * (basis->lenshell0 + 1));
    offset += sizeof(int) * (basis->lenshell0 + 1);
    
    assert (offset == _bufsize);

    *bufsize = _bufsize;
    *buf = (char *)_buf;

    return CINT_STATUS_SUCCESS; 
}


CIntStatus_t CInt_unpackBasisSet (BasisSet_t basis,
                                  void *buf)
{
    int offset;
    CIntStatus_t ret;
    char *_buf;

    _buf = (char *)buf;
    offset = 0;    
    memcpy (&(basis->natoms), &(_buf[offset]), sizeof(int));
    offset += sizeof(int);
    memcpy (&(basis->nelectrons), &(_buf[offset]), sizeof(int));
    offset += sizeof(int);
    memcpy (&(basis->basistype), &(_buf[offset]), sizeof(int));
    offset += sizeof(int);
    memcpy (&(basis->lenshell0), &(_buf[offset]), sizeof(int));
    offset += sizeof(int);
    memcpy (&(basis->totnexp), &(_buf[offset]), sizeof(int));
    offset += sizeof(int);
    memcpy (&(basis->lenatom0), &(_buf[offset]), sizeof(int));
    offset += sizeof(int);
    
    basis->xn = (double *)malloc (sizeof(double) * basis->natoms);
    basis->yn = (double *)malloc (sizeof(double) * basis->natoms);
    basis->zn = (double *)malloc (sizeof(double) * basis->natoms);
    basis->ncharge = (double *)malloc (sizeof(double) * basis->natoms); 
    basis->eid = (int *)malloc (sizeof(int) * basis->natoms);
    if (NULL == basis->xn ||
        NULL == basis->yn ||
        NULL == basis->zn ||
        NULL == basis->ncharge ||
        NULL == basis->eid)
    {
        CINT_PRINTF (1, "memory allocation failed\n");
        return CINT_STATUS_ALLOC_FAILED;
    }

    memcpy (basis->xn, &(_buf[offset]), sizeof(double) * basis->natoms);
    offset += sizeof(double) * basis->natoms;
    memcpy (basis->yn, &(_buf[offset]), sizeof(double) * basis->natoms);
    offset += sizeof(double) * basis->natoms;
    memcpy (basis->zn, &(_buf[offset]), sizeof(double) * basis->natoms);
    offset += sizeof(double) * basis->natoms;
    memcpy (basis->ncharge, &(_buf[offset]), sizeof(double) * basis->natoms);
    offset += sizeof(double) * basis->natoms;
    memcpy (basis->eid, &(_buf[offset]), sizeof(int) * basis->natoms);
    offset += sizeof(int) * basis->natoms;

    basis->cc0 = (double *)malloc (sizeof(double) * basis->totnexp);
    basis->exp0 = (double *)malloc (sizeof(double) * basis->totnexp);
    basis->eptr = (int *)malloc (sizeof(int) * ELEN);
    basis->atom_start = (int *)malloc (sizeof(int) * (basis->lenatom0 + 1)); 
    basis->momentum0 = (int *)malloc (sizeof(int) * basis->lenshell0);
    basis->nexp0 = (int *)malloc (sizeof(int) * basis->lenshell0);
    basis->ptrshell = (int *)malloc (sizeof(int) * (basis->lenshell0 + 1));
    if (NULL == basis->atom_start ||
        NULL == basis->eptr ||
        NULL == basis->nexp0 ||
        NULL == basis->cc0 ||
        NULL == basis->exp0 ||
        NULL == basis->momentum0 ||
        NULL == basis->ptrshell)
    {
        CINT_PRINTF (1, "memory allocation failed\n");
        return CINT_STATUS_ALLOC_FAILED;    
    }

    memcpy (basis->exp0, &(_buf[offset]), sizeof(double) * basis->totnexp);
    offset += sizeof(double) * basis->totnexp;
    memcpy (basis->cc0, &(_buf[offset]), sizeof(double) * basis->totnexp);
    offset += sizeof(double) * basis->totnexp;

    memcpy (basis->eptr, &(_buf[offset]), sizeof(int) * ELEN);
    offset += sizeof(int) * ELEN;
    memcpy (basis->atom_start, &(_buf[offset]), sizeof(int) * (basis->lenatom0 + 1));
    offset += sizeof(int) * (basis->lenatom0 + 1);
    
    memcpy (basis->momentum0, &(_buf[offset]), sizeof(int) * basis->lenshell0);
    offset += sizeof(int) * basis->lenshell0;
    memcpy (basis->nexp0, &(_buf[offset]), sizeof(int) * basis->lenshell0);
    offset += sizeof(int) * basis->lenshell0;
    memcpy (basis->ptrshell, &(_buf[offset]), sizeof(int) * (basis->lenshell0 + 1));
    offset += sizeof(int) * (basis->lenshell0 + 1);

    if ((ret = parse_molecule (basis)) != CINT_STATUS_SUCCESS)
    {
        return ret;
    }

    return CINT_STATUS_SUCCESS;    
}


int CInt_getNumShells (BasisSet_t basis)
{
    return (basis->nshells);
}


int CInt_getNumFuncs (BasisSet_t basis)
{
    return (basis->nfunctions);
}


int CInt_getNumAtoms (BasisSet_t basis)
{
    return (basis->natoms);
}


int CInt_getShellDim (BasisSet_t basis, int shellid)
{
    return (basis->f_end_id[shellid] -
        basis->f_start_id[shellid] + 1);
}


int CInt_getMaxShellDim (BasisSet_t basis)
{
    return (basis->maxdim);
}


int CInt_getNumOccOrb (BasisSet_t basis)
{
    return (basis->nelectrons/2);
}


int CInt_getFuncStartInd (BasisSet_t basis, int shellid)
{
    return (basis->f_start_id[shellid]);
}


int CInt_getFuncEndInd (BasisSet_t basis, int shellid)
{
    return (basis->f_end_id[shellid]);
}


int CInt_getAtomStartInd (BasisSet_t basis, int atomid)
{
    return (basis->s_start_id[atomid]);
}
