#include "psi4-dec.h"
#include <libqt/qt.h>
#include <cmath>
#include "adc.h"

namespace psi{ namespace adc{

// This source file is based on the sorting function in ccenergy, cceom and cis.    

struct onestack{
    double value;
    int i;
    int a;
};
    
//void onestack_insert(struct onestack *stack, double value, int i, int a, int level, int stacklen);
    
void 
ADC::amps_write(dpdfile2 *B, int length, FILE *outfile)
{
    struct onestack *t1stack;
    int Gia = B->my_irrep;
    
    t1stack = (struct onestack*)malloc(length*sizeof(struct onestack));
    for(int m = 0; m < length; m++) { t1stack[m].value = 0; t1stack[m].i = 0; t1stack[m].a = 0; }
    
    dpd_file2_mat_init(B);
    dpd_file2_mat_rd(B);
    
    int numt1 = 0;
    for(int h = 0;h < nirrep_;h++){
        numt1 += B->params->rowtot[h] * B->params->coltot[h^Gia];
        for(int i = 0;i < B->params->rowtot[h];i++){
            int I = B->params->roworb[h][i];
            for(int a = 0;a < B->params->coltot[h^Gia];a++){
                int A = B->params->colorb[h^Gia][a];
                double value = B->matrix[h][i][a];
                for(int m = 0;m < length;m++){
                    if((fabs(value)-fabs(t1stack[m].value)) > 1e-12){
                        onestack_insert(t1stack, value, I, A, m, length);
                        break;
                    }
                }
            }
        }
    }
    dpd_file2_mat_close(B);
    
    for(int m = 0;m < ((numt1 < length) ? numt1 : length);m++){
        if(fabs(t1stack[m].value) > 1e-6){
            fprintf(outfile, "\t        %3d %3d %20.10f\n", t1stack[m].i, t1stack[m].a, t1stack[m].value);
        }
    }
    free(t1stack);
}
    
void 
ADC::onestack_insert(struct onestack *stack, double value, int i, int a, int level, int stacklen)
{
    struct onestack temp;
    
    temp = stack[level];
    
    stack[level].value = value;
    stack[level].i = i;
    stack[level].a = a;
    
    value = temp.value;
    i = temp.i;
    a = temp.a;
    
    for(int l = level;l < stacklen-1;l++){
        temp = stack[l+1];
        
        stack[l+1].value = value;
        stack[l+1].i = i;
        stack[l+1].a = a;
        
        value = temp.value;
        i = temp.i;
        a = temp.a;
    }
}
    
}} // End Namespaces
