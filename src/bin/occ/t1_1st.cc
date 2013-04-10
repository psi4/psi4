#include <libtrans/integraltransform.h>
#include "defines.h"
#include "occwave.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace occwave{
  
void OCCWave::t1_1st_sc()
{   

    // Alpha spin case
    t1A->zero();
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int a = 0 ; a < avirtpiA[h]; ++a){
                double value = FockA->get(h, i + frzcpi_[h], i + frzcpi_[h]) - FockA->get(h, a + occpiA[h], a + occpiA[h]);
                t1A->set(h, i, a, FockA->get(h, i + frzcpi_[h], a + occpiA[h]) / value);
            }
        }
    }
    if (print_ > 1) t1A->print();

    // Beta spin case
    t1B->zero();
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < aoccpiB[h]; ++i){
            for(int a = 0 ; a < avirtpiB[h]; ++a){
                double value = FockB->get(h, i + frzcpi_[h], i + frzcpi_[h]) - FockB->get(h, a + occpiB[h], a + occpiB[h]);
                t1B->set(h, i, a, FockB->get(h, i + frzcpi_[h], a + occpiB[h]) / value);
            }
        }
    }
    if (print_ > 1) t1B->print();

} // end t1_1st_sc


//===========================================================================================
//========================= t1_1st_gen ======================================================
//===========================================================================================
void OCCWave::t1_1st_gen()
{   
     //fprintf(outfile,"\n t1_1st_gen is starting... \n"); fflush(outfile);
     // For this section the frozen-core approximation is NOT fully implemented!

    // Alpha spin case
    t1newA->zero();

    // t_I^A = F_IA
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int a = 0 ; a < avirtpiA[h]; ++a){
                t1newA->set(h, i, a, FockA->get(h, i + frzcpi_[h], a + occpiA[h]));
            }
        }
    }

    // t_I^A += \sum_{E != A} t_I^E F_AE
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int a = 0 ; a < avirtpiA[h]; ++a){
                double summ = 0.0;
                for(int e = 0 ; e < avirtpiA[h]; ++e){
                    if (e != a) summ += t1A->get(h, i, e) * FockA->get(h, a + occpiA[h], e + occpiA[h]);
                }
                t1newA->add(h, i, a, summ);
            }
        }
    }

    // t_I^A -= \sum_{M != I} t_M^A F_MI
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int a = 0 ; a < avirtpiA[h]; ++a){
                double summ = 0.0;
                for(int m = 0 ; m < aoccpiA[h]; ++m){
                    if (m != i) summ -= t1A->get(h, m, a) * FockA->get(h, m, i);
                }
                t1newA->add(h, i, a, summ);
            }
        }
    }

    // t_I^A /= D_IA
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int a = 0 ; a < avirtpiA[h]; ++a){
                double value = FockA->get(h, i, i) - FockA->get(h, a + occpiA[h], a + occpiA[h]);
                t1newA->set(h, i, a, t1newA->get(h, i, a) / value);
            }
        }
    }

    // Beta spin case
    // Off-diagonal Fock contribution
    t1newB->zero();

    // t_i^a = F_ia
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < aoccpiB[h]; ++i){
            for(int a = 0 ; a < avirtpiB[h]; ++a){
                t1newB->set(h, i, a, FockB->get(h, i, a + occpiB[h]));
            }
        }
    }

    // t_i^a += \sum_{e != a} t_i^e F_ae
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < aoccpiB[h]; ++i){
            for(int a = 0 ; a < avirtpiB[h]; ++a){
                double summ = 0.0;
                for(int e = 0 ; e < avirtpiB[h]; ++e){
                    if (e != a) summ += t1B->get(h, i, e) * FockB->get(h, a + occpiB[h], e + occpiB[h]);
                }
                t1newB->add(h, i, a, summ);
            }
        }
    }

    // t_i^a -= \sum_{m != i} t_m^a F_mi
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < aoccpiB[h]; ++i){
            for(int a = 0 ; a < avirtpiB[h]; ++a){
                double summ = 0.0;
                for(int m = 0 ; m < aoccpiB[h]; ++m){
                    if (m != i) summ -= t1B->get(h, m, a) * FockB->get(h, m, i);
                }
                t1newB->add(h, i, a, summ);
            }
        }
    }

    // t_i^a /= D_ia
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < aoccpiB[h]; ++i){
            for(int a = 0 ; a < avirtpiB[h]; ++a){
                double value = FockB->get(h, i, i) - FockB->get(h, a + occpiB[h], a + occpiB[h]);
                t1newB->set(h, i, a, t1newB->get(h, i, a) / value);
            }
        }
    }

    // RMS
    // alpha
    rms_t1A = 0.0;
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int a = 0 ; a < avirtpiA[h]; ++a){
                rms_t1A += pow(t1newA->get(h, i, a) - t1A->get(h, i, a), 2);
            }
        }
    }
    rms_t1A = sqrt(rms_t1A)/nidpA;

    // beta
    rms_t1B = 0.0;
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < aoccpiB[h]; ++i){
            for(int a = 0 ; a < avirtpiB[h]; ++a){
                rms_t1B += pow(t1newB->get(h, i, a) - t1B->get(h, i, a), 2);
            }
        }
    }
    rms_t1B = sqrt(rms_t1B)/nidpB;

    // reset
    t1A->zero();
    t1A->copy(t1newA);
    t1B->zero();
    t1B->copy(t1newB);

    if (print_ > 1) t1A->print();
    if (print_ > 1) t1B->print();

  //fprintf(outfile,"\n t1_1st_gen done. \n"); fflush(outfile);
} // end t1_1st_gen
}} // End Namespaces


