#include "scf.h"

#include <libmoinfo/libmoinfo.h>

namespace psi{ namespace MCSCF{

void SCF::construct_F()
{
  if(reference == rhf){
    Fc  = H;
    for(int batch = 0; batch < nbatch; ++batch){
      read_Raffanetti("PK",PK,batch);
      construct_G(Dc,G,PK,batch);
      Fc += G;
    }
  }else if(reference == rohf){
    Fc  = H;
    Fo  = H;
    Fo.scale(0.5);
    for(int batch = 0; batch < nbatch; ++batch){
      read_Raffanetti("PK",PK,batch);
      // Dc * PK Contributions
      construct_G(Dc,G,PK,batch);
      Fc += G;
      G.scale(0.5);
      Fo += G;

      // Do * PK Contributions
      construct_G(Do,G,PK,batch,0.5);
      Fc += G;
      G.scale(0.5);
      Fo += G;

      read_Raffanetti("K",K,batch);
      // Do * K Contributions
      construct_G(Do,G,K,batch,0.25);
      Fo += G;
    }
  }else if(reference == tcscf){
    Fc    = H;
    Favg  = H;
    for(int I = 0 ; I < nci; ++I){
      Dsum[I]  = Dc;
      Dsum[I] += Dtc[I];
      Ftc[I] = H;
      Ftc[I].scale(ci[I] * ci[I]);
      H_tcscf[I][I] = 2.0 * dot(Dsum[I],H) + moinfo_scf->get_nuclear_energy();
      for(int J = I + 1; J < nci; ++J)
        H_tcscf[I][J] = H_tcscf[J][I] = 0.0;
    }

    for(int batch = 0; batch < nbatch; ++batch){
      read_Raffanetti("PK",PK,batch);
      // Dc * PK Contributions to the Fock matrices
      construct_G(Dc,G,PK,batch);
      Fc += G;
      for(int I = 0 ; I < nci; ++I){
        T = G;
        T.scale(ci[I] * ci[I]);
        Ftc[I] += T;
      }

      // Dtc * PK Contributions to the Fock matrices
      for(int I = 0 ; I < nci; ++I){
        construct_G(Dtc[I],G,PK,batch,ci[I] * ci[I]);
        Fc += G;
        G.scale(0.5);
        Ftc[I] += G;
      }
       
      // Dsum * PK Contributions to the Hamiltonian
      for(int I = 0 ; I < nci; ++I){
        construct_G(Dsum[I],G,PK,batch);
        H_tcscf[I][I] += dot(Dsum[I],G);

        G.scale(ci[I] * ci[I]);
        Favg += G;
      }

      read_Raffanetti("K",K,batch);
      // Dtc * K Contributions
      for(int I = 0 ; I < nci; ++I){
        construct_G(Dtc[I],G,K,batch);
        T = G;
        T.scale(-0.5 * ci[I] * ci[I]);
        Ftc[I] += T;
        for(int J = 0 ; J < nci; ++J){
          if(I != J){
            T = G;
            T.scale(- ci[I] * ci[J]);
            Ftc[J] += T;
            // Compute off-diagonal elements of H
            H_tcscf[I][J] -= dot(Dtc[J],G);
          }
        }
      }

    }
  }
}

void SCF::construct_Favg()
{
  if(reference == tcscf){
    Favg  = H;
    for(int I = 0 ; I < nci; ++I){
      Dsum[I]  = Dc;
      Dsum[I] += Dtc[I];
    }

    for(int batch = 0; batch < nbatch; ++batch){
      read_Raffanetti("PK",PK,batch);
      // Dsum * PK Contributions to the Hamiltonian
      for(int I = 0 ; I < nci; ++I){
        construct_G(Dsum[I],G,PK,batch);
        G.scale(ci[I] * ci[I]);
        Favg += G;
      }
    }
  }
}

}} /* End Namespaces */
