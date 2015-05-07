#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>
#include <string.h>
#include <ctype.h>
#include <libgen.h>

#include "config.h"
#include "basisset.h"

#ifdef HAS_MALLOC_H
#include <malloc.h>
#endif

#define ELEN         50
#define SLEN         5
#define MAXNS        3
#define MAXATOMNAME  2
//#define A2BOHR       1.889726
#define A2BOHR       (1.0/0.52917720859)
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
    /*double sum;
    double temp;
    double temp2;
    double temp3;
    double xnorm;
    double a1;
    double a2;
    int i;
    int j;
    int k;
    double power;
    int shell;

    for (i = 0; i < basis->bs_nshells; i++)
    {
        sum = 0.0;
        for (j = 0; j < basis->bs_nexp[i]; j++)
        {
            for (k = 0; k <= j; k++)
            {
                a1 = basis->bs_exp[i][j];
                a2 = basis->bs_exp[i][k];
                temp = basis->bs_cc[i][j] * basis->bs_cc[i][k]; 
                temp2 = basis->bs_momentum[i] + 1.5;
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
        shell = basis->bs_momentum[i];
        power = (double) shell *0.5 + 0.75;
        for (j = 0; j < basis->bs_nexp[i]; j++)
        {
            basis->bs_cc[i][j] *= xnorm;
            if (basis->bs_exp[i][j] == 0.0)
            {
                basis->bs_norm[i][j] = 1.0;
            }
            else
            {
                basis->bs_norm[i][j] = pow (basis->bs_exp[i][j], power);
            }
        }              
    }    */
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


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

static inline double vector_min(size_t length, const double vector[restrict static length]) {
    double result = vector[0];
    for (size_t i = 1; i < length; i++) {
        const double element = vector[i];
        result = (element < result) ? element : result;
    }
    return result;
}


CIntStatus_t CInt_createBasisSet (BasisSet_t *_basis)
{
    BasisSet_t basis;
    basis = (BasisSet_t )malloc (sizeof(struct BasisSet));
    CINT_ASSERT(basis != NULL);
    memset (basis, 0, sizeof(struct BasisSet));

    *_basis = basis;
    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_destroyBasisSet (BasisSet_t basis)
{
    free (basis->f_start_id);
    free (basis->f_end_id);

    for (int i = 0; i < basis->bs_nshells; i++)
    {
        ALIGNED_FREE (basis->bs_cc[i]);
        ALIGNED_FREE (basis->bs_exp[i]);
        ALIGNED_FREE (basis->bs_norm[i]);
    }
    free (basis->bs_cc);
    free (basis->bs_exp);
    free (basis->bs_norm);

    if (basis->guess != NULL)
    {
        for (int i = 0; i < basis->bs_natoms; i++)
        {
            free (basis->guess[i]);
        }
        free (basis->guess);
    }    
    free (basis->eid);
    free (basis->xn);
    free (basis->yn);
    free (basis->zn);   
    free (basis->charge);
    free (basis->bs_eptr);
    free (basis->bs_atom_start);
    free (basis->bs_nexp);
    free (basis->bs_momentum);
    free (basis->cc);
    free (basis->exp);
    free (basis->minexp);
    free (basis->norm);

    free (basis);

    return CINT_STATUS_SUCCESS;
}


CIntStatus_t parse_molecule (BasisSet_t basis)
{
    int natoms;
    int nshells;   
    int nfunctions;
    int maxdim;
    int max_momentum;
    int max_nexp;
    int max_nexp_id;
    int eid;
    int atom_start;
    int atom_end;

    // get lengths
    natoms = basis->natoms;
    nshells = 0;
    for (uint32_t i = 0; i < natoms; i++) {
        eid = basis->eid[i];
        atom_start = basis->bs_atom_start[basis->bs_eptr[eid - 1]];
        atom_end = basis->bs_atom_start[basis->bs_eptr[eid - 1] + 1];
        nshells += atom_end - atom_start;
    }
    
    basis->s_start_id = (uint32_t *)malloc (sizeof(uint32_t) * (natoms + 1));
    basis->f_start_id = (uint32_t *)malloc (sizeof(uint32_t) * nshells);
    basis->f_end_id = (uint32_t *)malloc (sizeof(uint32_t) * nshells);
#ifdef HAS_MALLOC_H
    basis->xyz0 = (double *)memalign(32, sizeof(double) * nshells * 4);
#else
    basis->xyz0 = (double *)ALIGNED_MALLOC(sizeof(double) * nshells * 4);
#endif
    basis->nexp = (uint32_t *)malloc(sizeof(uint32_t) * nshells);
    basis->cc = (double **)malloc(sizeof(double *) * nshells);
    basis->exp = (double **)malloc(sizeof(double *) * nshells);
    basis->minexp = (double*)malloc(sizeof(double) * nshells);
    basis->norm = (double **)malloc(sizeof(double *) * nshells);
    basis->momentum = (uint32_t *)malloc(sizeof(uint32_t) * nshells);   
    CINT_ASSERT(basis->s_start_id != NULL);
    CINT_ASSERT(basis->f_start_id != NULL);
    CINT_ASSERT(basis->f_end_id != NULL);
    CINT_ASSERT(basis->xyz0 != NULL);
    CINT_ASSERT(basis->nexp != NULL);
    CINT_ASSERT(basis->cc != NULL);
    CINT_ASSERT(basis->minexp != NULL);
    CINT_ASSERT(basis->norm != NULL);
    CINT_ASSERT(basis->momentum != NULL);
    basis->nshells = nshells;

    // parse molecules
    nshells = 0;
    nfunctions = 0;
    maxdim = 0;
    max_momentum = 0;
    max_nexp = 0;
    max_nexp_id = 0;
    for (uint32_t i = 0; i < natoms; i++) {
        eid = basis->eid[i];    
        atom_start = basis->bs_atom_start[basis->bs_eptr[eid - 1]];
        atom_end = basis->bs_atom_start[basis->bs_eptr[eid - 1] + 1];
        /* Atom not supported */
        CINT_ASSERT(basis->bs_eptr[eid - 1] != -1);
        basis->s_start_id[i] = nshells;
        for (uint32_t j = atom_start; j < atom_end; j++) {
            basis->f_start_id[nshells + j - atom_start] = nfunctions;
            basis->nexp[nshells + j - atom_start] = basis->bs_nexp[j];
            basis->xyz0[(nshells + j - atom_start) * 4 + 0] = basis->xn[i];
            basis->xyz0[(nshells + j - atom_start) * 4 + 1] = basis->yn[i];
            basis->xyz0[(nshells + j - atom_start) * 4 + 2] = basis->zn[i];
            basis->momentum[nshells + j - atom_start] = basis->bs_momentum[j];
            max_momentum = (max_momentum > basis->bs_momentum[j] ?
                max_momentum : basis->bs_momentum[j]);
            if (max_nexp < basis->bs_nexp[j]) {
                max_nexp = basis->bs_nexp[j];
                max_nexp_id = nshells + j - atom_start;
            }
            basis->cc[nshells + j - atom_start] = basis->bs_cc[j];
            basis->exp[nshells + j - atom_start] = basis->bs_exp[j];
            basis->minexp[nshells + j - atom_start] = vector_min(basis->nexp[nshells + j - atom_start], basis->exp[nshells + j - atom_start]);
            basis->norm[nshells + j - atom_start] = basis->bs_norm[j];
            if (basis->basistype == SPHERICAL) {
                nfunctions += 2 * basis->bs_momentum[j] + 1;
                maxdim = (2 * basis->bs_momentum[j] + 1) > maxdim ?
                    (2 * basis->bs_momentum[j] + 1) : maxdim;
            }
            else if (basis->basistype == CARTESIAN) {
                nfunctions += (basis->bs_momentum[j] + 1)*(basis->bs_momentum[j] + 2)/2;
                maxdim = ((basis->bs_momentum[j] + 1)*(basis->bs_momentum[j] + 2)/2) > maxdim ?
                    ((basis->bs_momentum[j] + 1)*(basis->bs_momentum[j] + 2)/2) : maxdim;
            }
            basis->f_end_id[nshells + j - atom_start] = nfunctions - 1;
        }
        nshells += atom_end - atom_start;
    }
    basis->s_start_id[natoms] = nshells;
    basis->maxdim = maxdim;
    basis->nfunctions = nfunctions;
    basis->max_momentum = max_momentum;
    basis->max_nexp = max_nexp;
    basis->max_nexp_id = max_nexp_id;
    
    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_unpackBasisSet (BasisSet_t basis, void *buf)
{
    CIntStatus_t ret;
    char *_buf = (char *)buf;
    int offset = 0;
    memcpy(&(basis->natoms), &(_buf[offset]), sizeof(int));
    offset += sizeof(int);
    memcpy(&(basis->nelectrons), &(_buf[offset]), sizeof(int));
    offset += sizeof(int);
    memcpy(&(basis->bs_natoms), &(_buf[offset]), sizeof(int));
    offset += sizeof(int);
    memcpy(&(basis->basistype), &(_buf[offset]), sizeof(int));
    offset += sizeof(int);
    memcpy(&(basis->bs_nshells), &(_buf[offset]), sizeof(int));
    offset += sizeof(int);
    memcpy(&(basis->bs_totnexp), &(_buf[offset]), sizeof(int));
    offset += sizeof(int);
    memcpy(&(basis->bs_nelements), &(_buf[offset]), sizeof(int));
    offset += sizeof(int);
    basis->xn = (double *)malloc(sizeof(double) * basis->natoms);
    basis->yn = (double *)malloc(sizeof(double) * basis->natoms);
    basis->zn = (double *)malloc(sizeof(double) * basis->natoms);
    basis->charge = (double *)malloc(sizeof(double) * basis->natoms); 
    basis->eid = (int *)malloc(sizeof(int) * basis->natoms);
    CINT_ASSERT(basis->xn != NULL);
    CINT_ASSERT(basis->yn != NULL);
    CINT_ASSERT(basis->zn != NULL);
    CINT_ASSERT(basis->charge != NULL);
    CINT_ASSERT(basis->eid != NULL);
    memcpy(basis->xn, &(_buf[offset]), sizeof(double) * basis->natoms);
    offset += sizeof(double) * basis->natoms;
    memcpy(basis->yn, &(_buf[offset]), sizeof(double) * basis->natoms);
    offset += sizeof(double) * basis->natoms;
    memcpy(basis->zn, &(_buf[offset]), sizeof(double) * basis->natoms);
    offset += sizeof(double) * basis->natoms;
    memcpy(basis->charge, &(_buf[offset]), sizeof(double) * basis->natoms);
    offset += sizeof(double) * basis->natoms;
    memcpy(basis->eid, &(_buf[offset]), sizeof(int) * basis->natoms);
    offset += sizeof(int) * basis->natoms;
    basis->bs_cc = (double **)malloc(sizeof(double *) * basis->bs_nshells);
    basis->bs_exp = (double **)malloc(sizeof(double *) * basis->bs_nshells);
    basis->bs_norm = (double **)malloc(sizeof(double *) * basis->bs_nshells);
    basis->bs_eptr = (int *)malloc(sizeof(int) * basis->bs_nelements);
    basis->bs_atom_start = (int *)malloc(sizeof(int) * (basis->bs_natoms + 1)); 
    basis->bs_momentum = (int *)malloc(sizeof(int) * basis->bs_nshells);
    basis->bs_nexp = (int *)malloc(sizeof(int) * basis->bs_nshells);
    CINT_ASSERT(basis->bs_cc != NULL);
    CINT_ASSERT(basis->bs_exp != NULL);
    CINT_ASSERT(basis->bs_norm != NULL);
    CINT_ASSERT(basis->bs_eptr != NULL);
    CINT_ASSERT(basis->bs_atom_start != NULL);
    CINT_ASSERT(basis->bs_momentum != NULL);
    CINT_ASSERT(basis->bs_nexp != NULL);

    memcpy(basis->bs_nexp, &(_buf[offset]), sizeof(int) * basis->bs_nshells);
    offset += sizeof(int) * basis->bs_nshells;
    for (int i = 0; i < basis->bs_nshells; i++) {
        int nexp = basis->bs_nexp[i];
        basis->bs_cc[i] = (double *)ALIGNED_MALLOC(sizeof(double) * nexp);
        basis->bs_exp[i] = (double *)ALIGNED_MALLOC(sizeof(double) * nexp);
        basis->bs_norm[i] = (double *)ALIGNED_MALLOC(sizeof(double) * nexp);
        CINT_ASSERT(basis->bs_cc[i] != NULL);
        CINT_ASSERT(basis->bs_exp[i] != NULL);
        CINT_ASSERT(basis->bs_norm[i] != NULL);
        memcpy(basis->bs_exp[i], &(_buf[offset]), sizeof(double) * nexp);
        offset += sizeof(double) * nexp;
        memcpy(basis->bs_cc[i], &(_buf[offset]), sizeof(double) * nexp);
        offset += sizeof(double) * nexp;
        memcpy(basis->bs_norm[i], &(_buf[offset]), sizeof(double) * nexp);
        offset += sizeof(double) * nexp;
    }
    memcpy(basis->bs_eptr, &(_buf[offset]), sizeof(int) * basis->bs_nelements);
    offset += sizeof(int) * basis->bs_nelements;
    memcpy(basis->bs_atom_start, &(_buf[offset]),
           sizeof(int) * (basis->bs_natoms + 1));
    offset += sizeof(int) * (basis->bs_natoms + 1);
    memcpy(basis->bs_momentum, &(_buf[offset]), sizeof(int) * basis->bs_nshells);
    offset += sizeof(int) * basis->bs_nshells;

    if ((ret = parse_molecule(basis)) != CINT_STATUS_SUCCESS) {
        return ret;
    }

    return CINT_STATUS_SUCCESS;    
}


#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif


CIntStatus_t import_molecule (char *file, BasisSet_t basis)
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

    // number of atoms    
    if (fgets (line, 1024, fp) == NULL)
    {
        CINT_PRINTF (1, "file %s has a wrong format\n", file);
        return CINT_STATUS_FILEIO_FAILED; 
    }
    sscanf (line, "%d", &(basis->natoms));    
    if (basis->natoms <= 0)
    {
        CINT_PRINTF (1, "file %s has a wrong format\n", file);
        return CINT_STATUS_FILEIO_FAILED;
    }
        
    // skip comment line
    if (fgets (line, 1024, fp) == NULL)
    {
        CINT_PRINTF (1, "file %s has a wrong format\n", file);
        return CINT_STATUS_FILEIO_FAILED; 
    }
    
    basis->xn = (double *)malloc (sizeof(double) * basis->natoms);
    basis->yn = (double *)malloc (sizeof(double) * basis->natoms);
    basis->zn = (double *)malloc (sizeof(double) * basis->natoms);
    basis->charge = (double *)malloc (sizeof(double) * basis->natoms); 
    basis->eid = (int *)malloc (sizeof(int) * basis->natoms);
    if (NULL == basis->xn ||
        NULL == basis->yn ||
        NULL == basis->zn ||
        NULL == basis->charge ||
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
        if (isalpha(str[0]))
        {
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
            basis->charge[natoms] = (double)(basis->eid[natoms]);
            nelectrons += basis->eid[natoms];
            natoms++;
        }
    }
    basis->nelectrons = nelectrons;
    if (natoms != basis->natoms)
    {
        CINT_PRINTF (1, "file %s natoms %d does not match the header\n",
            file, natoms);
        return CINT_STATUS_FILEIO_FAILED;
    }

    // compute nuc energy
    double ene = 0.0;
    for (int A = 0; A < natoms; A++)
    {
        for (int B = A + 1; B < natoms; B++)
        {
            double dx = basis->xn[A] - basis->xn[B];
            double dy = basis->yn[A] - basis->yn[B];
            double dz = basis->zn[A] - basis->zn[B];
            double R = sqrt(dx * dx + dy * dy + dz * dz);
            ene += basis->charge[A] * basis->charge[B] / R;
        }
    }
    basis->ene_nuc = ene;
    
    fclose (fp);
    
    return CINT_STATUS_SUCCESS;
}


CIntStatus_t import_basis (char *file, BasisSet_t basis)
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
    int bs_totnexp;

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
    bs_totnexp = 0;
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
                    bs_totnexp += ns * nexp;
                }
                if (line[0] == '*')
                {
                    break;
                }
            }
         }
    }
    basis->bs_natoms = natoms;
    basis->bs_nshells = nshells;
    basis->bs_totnexp = bs_totnexp;
    basis->bs_nelements = ELEN;
    basis->bs_eptr = (int *)malloc (sizeof(int) * basis->bs_nelements);
    basis->bs_atom_start = (int *)malloc (sizeof(int) * (natoms + 1));
    basis->bs_nexp = (int *)malloc (sizeof(int) * nshells);
    basis->bs_cc = (double **)malloc (sizeof(double *) * nshells);
    basis->bs_norm = (double **)malloc (sizeof(double *) * nshells);
    basis->bs_exp = (double **)malloc (sizeof(double *) * nshells);
    basis->bs_momentum = (int *)malloc (sizeof(int) * nshells);
    basis->bs_eid = (int *)malloc (sizeof(int) * natoms);
    CINT_ASSERT(basis->bs_eptr != NULL);
    CINT_ASSERT(basis->bs_atom_start != NULL);
    CINT_ASSERT(basis->bs_nexp != NULL);
    CINT_ASSERT(basis->bs_cc != NULL);
    CINT_ASSERT(basis->bs_norm != NULL);
    CINT_ASSERT(basis->bs_exp != NULL);
    CINT_ASSERT(basis->bs_momentum != NULL);
    CINT_ASSERT(basis->bs_eid != NULL);
    for (i = 0; i < basis->bs_nelements; i++) {
        basis->bs_eptr[i] = -1;
    }
    
    // get nshells
    rewind (fp);
    fgets (line, 1024, fp);
    natoms = 0;
    nshells = 0;
    bs_totnexp = 0;
    while (fgets (line, 1024, fp) != NULL)
    {
        if (isalpha (line[0]))
        {
            // a new atom
            sscanf (line, "%s", str);
            for (i = 0; i < basis->bs_nelements; i++)
            {
                if (strcmp (str, etable[i]) == 0)
                {
                    basis->bs_eptr[i] = natoms;
                    basis->bs_eid[natoms] = i;
                    break;
                }
            }
            if (i == basis->bs_nelements)
            {
                CINT_PRINTF (1, "atom %s in %s is not supported\n", str, file);
                return CINT_STATUS_INVALID_VALUE;
            }
            basis->bs_atom_start[natoms] = nshells;           
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
                    for (i = 0; i < ns; i++)
                    {
                        basis->bs_nexp[nshells] = nexp;
                        basis->bs_cc[nshells] = (double *)ALIGNED_MALLOC (sizeof(double) * nexp);
                        basis->bs_exp[nshells] = (double *)ALIGNED_MALLOC (sizeof(double) * nexp);
                        basis->bs_norm[nshells] = (double *)ALIGNED_MALLOC (sizeof(double) * nexp);
                        assert(basis->bs_cc[nshells] != NULL);
                        assert(basis->bs_exp[nshells] != NULL);
                        assert(basis->bs_norm[nshells] != NULL);
                        for (j = 0; j < SLEN; j++)
                        {
                            if (str[i] == mtable[j])
                            {
                                basis->bs_momentum[nshells] = j;
                                break;
                            }
                        }
                        if (j == SLEN)
                        {
                            CINT_PRINTF (1, "shell %s in file %s is not supported\n",
                                str, file);
                            return CINT_STATUS_INVALID_VALUE;  
                        }
                        fseek (fp, mark, SEEK_SET);
                        for (j = 0; j < basis->bs_nexp[nshells]; j++)
                        {
                            if (fgets (line, 1024, fp) == NULL ||
                                line[0] == '*' ||
                                isalpha (line[0]))
                            {
                                CINT_PRINTF (1, "file %s has a wrong format\n", file);
                                return CINT_STATUS_FILEIO_FAILED;
                            }
                            sscanf (line, "%lf %lf %lf %lf",
                                    &(basis->bs_exp[nshells][j]),
                                    &(cc[0]), &(cc[1]), &(cc[2]));
                            basis->bs_cc[nshells][j] = cc[i];
                        }
                        bs_totnexp += basis->bs_nexp[nshells];
                        nshells++;
                    }
                }
                if (line[0] == '*')
                {
                    break;
                }
            }
         }
    }
    basis->bs_atom_start[natoms] = nshells;
    if (nshells != basis->bs_nshells || basis->bs_totnexp != bs_totnexp)
    {
        CINT_PRINTF (1, "file %s has a wrong format\n", file);
        return CINT_STATUS_FILEIO_FAILED;    
    }

    fclose (fp);

    normalization (basis);
    
    return CINT_STATUS_SUCCESS;
}


CIntStatus_t import_guess (char *file, BasisSet_t basis)
{
    char *dir;
    if (file != NULL) {
        dir = strdup(file);
        dir = dirname(dir);
    } else {
        dir = "/";
    }
    
    char fname[1024];
    char line[1024];    
    basis->guess = (double **)malloc (sizeof(double *) * basis->bs_natoms); 
    if (basis->guess == NULL)
    {
        return CINT_STATUS_ALLOC_FAILED;
    }
     
    for (int i = 0; i < basis->bs_natoms; i++)
    {
        const int atom_start = basis->bs_atom_start[i];
        const int atom_end = basis->bs_atom_start[i + 1];
        int eid = basis->bs_eid[i];
        int nfunctions = 0;
        for (int j = atom_start; j < atom_end; j++)
        {
            if (basis->basistype == SPHERICAL)
            {
                nfunctions += 2 * basis->bs_momentum[j] + 1;
            }
            else if (basis->basistype == CARTESIAN)
            {
                nfunctions += (basis->bs_momentum[j] + 1)*(basis->bs_momentum[j] + 2)/2;
            }
        }
        basis->guess[i] = (double *)malloc (sizeof(double) * nfunctions * nfunctions);

        // read guess
        eid = (eid >= ELEN ? 0 : eid);
        sprintf(fname, "%s/%s.dat", dir, etable[eid]);
        FILE *fp = fopen(fname, "r");
        int flag = 0;
        if (fp != NULL)
        {
            for (int j = 0; j < nfunctions * nfunctions; j++)
            {
                if (fgets (line, 1024, fp) == NULL)
                {
                    flag = 1;
                    goto end;
                }
                sscanf (line, "%le", &(basis->guess[i][j]));
                basis->guess[i][j] /= 2.0;
            }
            // test symmetry
            for (int j = 0; j < nfunctions; j++)
            {
                for (int k = 0; k < nfunctions; k++)
                {
                    if (basis->guess[i][j * nfunctions + k] !=
                        basis->guess[i][k * nfunctions + j])
                    {
                        flag = 1;
                        goto end;
                    }
                }
            }
        }
        else
        {
            flag = 1;
        }
end:        
        if (flag == 1)
        {
            memset (basis->guess[i], 0, sizeof(double) * nfunctions * nfunctions);          
        }
    }
    
    if (file != NULL) {
        free(dir);
    }
    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_packBasisSet (BasisSet_t basis,
                                void **buf,
                                int *bufsize)
{
    int _bufsize;
    char *_buf;
    int offset;
    int i;
    int nexp;
    
    _bufsize = 6 * sizeof(int) + 4 * basis->natoms * sizeof(double) +                
               (2 * basis->bs_nshells + basis->bs_nelements + basis->natoms
                + basis->bs_natoms + 2) * sizeof(int) +
                basis->bs_totnexp * 3 * sizeof(double);
    _buf = (char *)malloc(_bufsize);
    assert(_buf != NULL);
    offset = 0;    
    memcpy(&(_buf[offset]), &(basis->natoms), sizeof(int));
    offset += sizeof(int);
    memcpy(&(_buf[offset]), &(basis->nelectrons), sizeof(int));
    offset += sizeof(int);
    memcpy(&(_buf[offset]), &(basis->bs_natoms), sizeof(int));
    offset += sizeof(int);
    memcpy(&(_buf[offset]), &(basis->basistype), sizeof(int));
    offset += sizeof(int);
    memcpy(&(_buf[offset]), &(basis->bs_nshells), sizeof(int));
    offset += sizeof(int);
    memcpy(&(_buf[offset]), &(basis->bs_totnexp), sizeof(int));
    offset += sizeof(int);
    memcpy(&(_buf[offset]), &(basis->bs_nelements), sizeof(int));
    offset += sizeof(int);
    memcpy(&(_buf[offset]), basis->xn, sizeof(double) * basis->natoms);
    offset += sizeof(double) * basis->natoms;
    memcpy(&(_buf[offset]), basis->yn, sizeof(double) * basis->natoms);
    offset += sizeof(double) * basis->natoms;
    memcpy(&(_buf[offset]), basis->zn, sizeof(double) * basis->natoms);
    offset += sizeof(double) * basis->natoms;
    memcpy(&(_buf[offset]), basis->charge, sizeof(double) * basis->natoms);
    offset += sizeof(double) * basis->natoms;
    memcpy(&(_buf[offset]), basis->eid, sizeof(int) * basis->natoms);
    offset += sizeof(int) * basis->natoms;
    memcpy(&(_buf[offset]), basis->bs_nexp, sizeof(int) * basis->bs_nshells);
    offset += sizeof(int) * basis->bs_nshells;
    for (i = 0; i < basis->bs_nshells; i++)
    {
        nexp = basis->bs_nexp[i];
        memcpy(&(_buf[offset]), basis->bs_exp[i], sizeof(double) * nexp);
        offset += sizeof(double) * nexp;
        memcpy(&(_buf[offset]), basis->bs_cc[i], sizeof(double) * nexp);
        offset += sizeof(double) * nexp;
        memcpy(&(_buf[offset]), basis->bs_norm[i], sizeof(double) * nexp);
        offset += sizeof(double) * nexp;
        assert (offset < _bufsize);
    }
    memcpy(&(_buf[offset]), basis->bs_eptr, sizeof(int) * basis->bs_nelements);
    offset += sizeof(int) * basis->bs_nelements;
    memcpy(&(_buf[offset]), basis->bs_atom_start,
           sizeof(int) * (basis->bs_natoms + 1));
    offset += sizeof(int) * (basis->bs_natoms + 1);
    memcpy(&(_buf[offset]), basis->bs_momentum, sizeof(int) * basis->bs_nshells);
    offset += sizeof(int) * basis->bs_nshells;
    
    assert (offset == _bufsize);

    *bufsize = _bufsize;
    *buf = (void *)_buf;

    return CINT_STATUS_SUCCESS; 
}


CIntStatus_t CInt_importBasisSet(BasisSet_t basis,
                                 int natoms, int *Zs,
                                 double *X, double *Y, double *Z,
                                 int nprims, int nshells, int pure,
                                 int *shells_p_atom,
                                 int *prims_p_shell,
                                 int *L, double *cc, double *alpha)
{
    // import basis
    basis->bs_natoms = natoms;
    basis->bs_nshells = nshells;
    basis->bs_totnexp = nprims;
    if (pure == 0) {
        basis->basistype = CARTESIAN;
    } else if (pure == 1) {
        basis->basistype = SPHERICAL;
    } else {
        CINT_PRINTF(1, "pure must be 0 or 1\n");
        return CINT_STATUS_INVALID_VALUE;
    }
    basis->bs_nelements = basis->bs_natoms;
    CINT_ASSERT(natoms > 0);
    CINT_ASSERT(nshells > 0);
    CINT_ASSERT(nprims > 0);
    
    basis->bs_eptr = (int *)malloc(sizeof(int) * basis->bs_nelements);
    basis->bs_atom_start = (int *)malloc(sizeof(int) * (natoms + 1));
    basis->bs_nexp = (int *)malloc(sizeof(int) * nshells);
    basis->bs_cc = (double **)malloc(sizeof(double *) * nshells);
    basis->bs_norm = (double **)malloc(sizeof(double *) * nshells);
    basis->bs_exp = (double **)malloc(sizeof(double *) * nshells);
    basis->bs_momentum = (int *)malloc(sizeof(int) * nshells);
    basis->bs_eid = (int *)malloc(sizeof(int) * natoms);
    CINT_ASSERT(basis->bs_eptr != NULL);
    CINT_ASSERT(basis->bs_atom_start != NULL);
    CINT_ASSERT(basis->bs_nexp != NULL);
    CINT_ASSERT(basis->bs_cc != NULL);
    CINT_ASSERT(basis->bs_norm != NULL);
    CINT_ASSERT(basis->bs_exp != NULL);
    CINT_ASSERT(basis->bs_momentum != NULL);
    CINT_ASSERT(basis->bs_eid != NULL);

    basis->bs_atom_start[0] = 0;
    int start_prim = 0;
    int num_prims = 0;
    for (int i = 0; i < basis->bs_natoms; i++) {
        basis->bs_eptr[i] = i;
        basis->bs_eid[i] = i;
        basis->bs_atom_start[i + 1] =
            basis->bs_atom_start[i] + shells_p_atom[i];
        int start_shell = basis->bs_atom_start[i];
        int end_shell = basis->bs_atom_start[i + 1];    
        for (int j = start_shell; j < end_shell; j++) {
            CINT_ASSERT(j < nshells);
            int nexp = prims_p_shell[j];
            basis->bs_nexp[j] = nexp;
            basis->bs_momentum[j] = L[j];
            // copy cc, alpha and norm
            basis->bs_cc[j] =
                (double *)ALIGNED_MALLOC(sizeof(double) * nexp);
            basis->bs_exp[j] =
                (double *)ALIGNED_MALLOC(sizeof(double) * nexp);
            basis->bs_norm[j] =
                (double *)ALIGNED_MALLOC(sizeof(double) * nexp);
            CINT_ASSERT(basis->bs_cc[j] != NULL);
            CINT_ASSERT(basis->bs_exp[j] != NULL);
            CINT_ASSERT(basis->bs_norm[j] != NULL);
            memcpy(basis->bs_cc[j], &cc[start_prim], sizeof(double) * nexp);
            memcpy(basis->bs_exp[j], &alpha[start_prim], sizeof(double) * nexp);
            start_prim += nexp;
            num_prims += nexp;           
        }
    }
    CINT_ASSERT(num_prims == basis->bs_totnexp);
    normalization(basis);
    
    // import molecule
    basis->natoms = natoms;
    basis->xn = (double *)malloc(sizeof(double) * basis->natoms);
    basis->yn = (double *)malloc(sizeof(double) * basis->natoms);
    basis->zn = (double *)malloc(sizeof(double) * basis->natoms);
    basis->charge = (double *)malloc(sizeof(double) * basis->natoms); 
    basis->eid = (int *)malloc(sizeof(int) * basis->natoms);
    if (NULL == basis->xn ||
        NULL == basis->yn ||
        NULL == basis->zn ||
        NULL == basis->charge ||
        NULL == basis->eid) {
        CINT_PRINTF(1, "memory allocation failed\n");
        return CINT_STATUS_ALLOC_FAILED;
    }
    memcpy(basis->xn, X, sizeof(double) * basis->natoms);
    memcpy(basis->yn, Y, sizeof(double) * basis->natoms);
    memcpy(basis->zn, Z, sizeof(double) * basis->natoms);
    basis->nelectrons = 0;
    for (int i = 0; i < natoms; i++) {
        basis->eid[i] = i + 1;
        basis->charge[i] = (double)(Zs[i]);
        basis->nelectrons += Zs[i];
    }
    // compute nuc energy
    double ene = 0.0;
    for (int A = 0; A < natoms; A++) {
        for (int B = A + 1; B < natoms; B++) {
            double dx = basis->xn[A] - basis->xn[B];
            double dy = basis->yn[A] - basis->yn[B];
            double dz = basis->zn[A] - basis->zn[B];
            double R = sqrt(dx * dx + dy * dy + dz * dz);
            ene += basis->charge[A] * basis->charge[B] / R;
        }
    }
    basis->ene_nuc = ene;
    // parse molecule
    CIntStatus_t status;
    // parse xyz
    if ((status = parse_molecule(basis)) != CINT_STATUS_SUCCESS)
    {
        return status;
    }

    // import guess
    if ((status = import_guess(NULL, basis)) != CINT_STATUS_SUCCESS)
    {
        return status;
    }


    return CINT_STATUS_SUCCESS;
}


CIntStatus_t CInt_loadBasisSet (BasisSet_t basis, char *bsfile, char *molfile)
{
    CIntStatus_t status;

    // read xyz file
    if ((status = import_molecule (molfile, basis)) != CINT_STATUS_SUCCESS)
    {
        return status;
    }
    // read basis set
    if ((status = import_basis (bsfile, basis)) != CINT_STATUS_SUCCESS)
    {
        return status;
    }
    
    // parse xyz
    if ((status = parse_molecule (basis)) != CINT_STATUS_SUCCESS)
    {
        return status;
    }

    // import guess
    if ((status = import_guess (bsfile, basis)) != CINT_STATUS_SUCCESS)
    {
        return status;
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


void CInt_getShellxyz ( BasisSet_t basis,
                        int shellid,
                        double *x,
                        double *y,
                        double *z )
{
    *x = basis->xyz0[shellid*4 + 0];
    *y = basis->xyz0[shellid*4 + 1];
    *z = basis->xyz0[shellid*4 + 2];
}


int CInt_getShellDim (BasisSet_t basis, int shellid)
{
    return (basis->f_end_id[shellid] -
        basis->f_start_id[shellid] + 1);
}


int CInt_getFuncStartInd (BasisSet_t basis, int shellid)
{
    return (basis->f_start_id[shellid]);
}


int CInt_getMaxShellDim (BasisSet_t basis)
{
    return (basis->maxdim);
}


int CInt_getNumOccOrb (BasisSet_t basis)
{
    return (basis->nelectrons/2);
}


int CInt_getFuncEndInd (BasisSet_t basis, int shellid)
{
    return (basis->f_end_id[shellid]);
}


int CInt_getAtomStartInd (BasisSet_t basis, int atomid)
{
    return (basis->s_start_id[atomid]);
}


void CInt_getInitialGuess (BasisSet_t basis, int atomid, double **guess,
                           int *spos, int *epos)
{
    const int eid = basis->eid[atomid] - 1;
    *guess = basis->guess[basis->bs_eptr[eid]];
    const int start_shell = basis->s_start_id[atomid];
    const int end_shell = basis->s_start_id[atomid + 1];
    *spos = basis->f_start_id[start_shell];
    *epos = basis->f_end_id[end_shell - 1];
}


double CInt_getNucEnergy (BasisSet_t basis)
{
    return basis->ene_nuc;
}
