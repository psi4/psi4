/*! \file
    \ingroup OPTKING
    \brief READ_CONSTRAINTS reads constraints from  "FIXED_INTCO" section of input
   or intco files, the constrained internals must already be present in the
   list of all of the simple internal coordinates
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"

#include <libipv1/ip_lib.h>

namespace psi { //namespace optking {

int *read_constraints(const simples_class & simples) {
  int num_type, i,j,a,b,c,d,cnt,*sign,id, sub_index, sub_index2;
  int iconstraints, *constraints;
  Intco_type intco_type;

  optinfo.constraints_present = 0;
  optinfo.nconstraints = 0;
  iconstraints = 0;

  fprintf(outfile,"\nSearching for geometrical constraints...");

  if ((ip_exist(":FIXED_INTCO",0)) || (optinfo.fix_interfragment) || (optinfo.fix_intrafragment)) {
    ip_cwk_clear(); /* search only fixed_intco */
    ip_cwk_add(":FIXED_INTCO");

    /* 1st time through we are just checking and counting */
    /* second time through allocate memory and add internals */
    for (cnt=1; cnt>=0; --cnt) {

      if (!cnt) {
        constraints = new int[optinfo.nconstraints];
        fprintf(outfile,"Coordinates to be constrained:\n");
      }

      if (ip_exist("STRE",0)) {
        num_type=0;
        ip_count("STRE",&num_type,0);
        for(i=0;i<num_type;++i) {
          ip_count("STRE",&j,1,i);
          if (j != 2) {
            fprintf(outfile,"FIXED_INTCO: Stretch %d is of wrong dimension.\n",i+1);
            exit(2);
          }
          if (!cnt) {
            ip_data("STRE","%d",&(a),2,i,0);
            ip_data("STRE","%d",&(b),2,i,1);
            swap(&a,&b);
            id = simples.get_id_from_atoms_stre(a-1,b-1);
            constraints[iconstraints++] = simples.id_to_index(id);
            fprintf(outfile,"Stretch %d: %d %d\n", id, a, b);
          }
          else {
            optinfo.nconstraints++;
          }
        }
      }

      if (ip_exist("BEND",0)) {
        num_type = 0;
        ip_count("BEND",&num_type,0);
        for(i=0;i<num_type;++i) {
          ip_count("BEND",&j,1,i);
          if (j != 3) {
            fprintf(outfile,"FIXED_INTCO: Bend %d is of wrong dimension.\n",i+1);
            exit(2);
          }
          ip_data("BEND","%d",&(a),2,i,0);
          ip_data("BEND","%d",&(b),2,i,1);
          ip_data("BEND","%d",&(c),2,i,2);
          if (!cnt) {
            swap(&a,&c);
            id = simples.get_id_from_atoms_bend(a-1,b-1,c-1);
            constraints[iconstraints++] = simples.id_to_index(id);
            fprintf(outfile,"Angle %d: %d %d %d\n", id, a, b, c);
          }
          else {
            optinfo.nconstraints++;
          }
        }
      }

      if (ip_exist("LIN1",0)) {
        num_type = 0;
        ip_count("LIN1",&num_type,0);
        for(i=0;i<num_type;++i) {
          ip_count("LIN1",&j,1,i);
          if (j != 3) {
            fprintf(outfile,"FIXED_INTCO: Linear bend(1) %d is of \
                wrong dimension.\n",i+1);
            exit(2);
          }
          ip_data("LIN1","%d",&(a),2,i,0);
          ip_data("LIN1","%d",&(b),2,i,1);
          ip_data("LIN1","%d",&(c),2,i,2);
          if (!cnt) {
            swap(&a,&c);
            id = simples.get_id_from_atoms_linb(a-1,b-1,c-1,1);
            constraints[iconstraints++] = simples.id_to_index(id);
            fprintf(outfile,"Lin1 %d: %d %d %d\n", id, a, b, c);
          }
          else {
            optinfo.nconstraints++;
          }
        }
      }

      if (ip_exist("LIN2",0)) {
        num_type = 0;
        ip_count("LIN2",&num_type,0);
        for(i=0;i<num_type;++i) {
          ip_count("LIN2",&j,1,i);
          if (j != 3) {
            fprintf(outfile,"FIXED_INTCO: Linear bend(2) %d is of \
                wrong dimension.\n",i+1);
            exit(2);
          }
          ip_data("LIN2","%d",&(a),2,i,0);
          ip_data("LIN2","%d",&(b),2,i,1);
          ip_data("LIN2","%d",&(c),2,i,2); 
          if (!cnt) {
            swap(&a,&c);
            id = simples.get_id_from_atoms_linb(a-1,b-1,c-1,2);
            constraints[iconstraints++] = simples.id_to_index(id);
            fprintf(outfile,"Lin2 %d: %d %d %d\n", id, a, b, c);
          }
          else {
            optinfo.nconstraints++;
          }
        }
      }

      if (ip_exist("TORS",0)) {
        num_type = 0;
        ip_count("TORS",&num_type,0);
        for(i=0;i<num_type;++i) {
          ip_count("TORS",&j,1,i);
          if (j != 4) {
            fprintf(outfile,"FIXED_INTCO: Torsion %d is of \
                wrong dimension.\n",i+1);
            exit(2);
          }
          ip_data("TORS","%d",&(a),2,i,0);
          ip_data("TORS","%d",&(b),2,i,1);
          ip_data("TORS","%d",&(c),2,i,2);
          ip_data("TORS","%d",&(d),2,i,3);
          if (!cnt) {
            swap_tors(&a,&b,&c,&d);
            id = simples.get_id_from_atoms_tors(a-1,b-1,c-1,d-1);
            constraints[iconstraints++] = simples.id_to_index(id);
            fprintf(outfile,"Torsion %d: %d %d %d %d\n", id, a, b, c, d);
          }
          else {
            optinfo.nconstraints++;
          }
        }
      }

      if (ip_exist("OUT",0)) {
        num_type = 0;
        ip_count("OUT",&num_type,0);
        for(i=0; i<num_type; ++i) {
          ip_count("OUT",&j,1,i);
          if (j != 4) {
            fprintf(outfile,"FIXED_INTCO: Out-of-plane %d is of \
                wrong dimension.\n",i+1);
            exit(2);
          }
          ip_data("OUT","%d",&(a),2,i,0);
          ip_data("OUT","%d",&(b),2,i,1);
          ip_data("OUT","%d",&(c),2,i,2);
          ip_data("OUT","%d",&(d),2,i,3);
          if (!cnt) {
            id = simples.get_id_from_atoms_out(a-1,b-1,c-1,d-1,sign);
            constraints[iconstraints++] = simples.id_to_index(id);
            fprintf(outfile,"Out-of-plane %d: %d %d %d %d\n", id, a, b, c, d);
          }
          else {
            optinfo.nconstraints++;
          }
        }
      }

      // freeze user-specified interfragment coordinates
      if (ip_exist("FRAG",0)) {
        num_type = 0;
        ip_count("FRAG",&num_type,0);
        for(i=0; i<num_type; ++i) {
          ip_data("FRAG","%d",&(id),2,i,0);
          if (!cnt) {
            constraints[iconstraints++] = simples.id_to_index(id);
            fprintf(outfile,"Fragment coordinate %d\n", id);
          }
          else {
            optinfo.nconstraints++;
          } 
        }
      }

      if (optinfo.fix_interfragment) {
        for (i=0; i<simples.get_num(FRAG,1); ++i) { // loop over interfragment sets
          id = simples.get_id(FRAG, i, 0); // executes once for each interfragment set
          for (int I=0; I<6; ++I) {
            if (simples.frag_get_coord_on(i,I)) {
              if (!cnt) {
                constraints[iconstraints++] = simples.id_to_index(id) + I;
                fprintf(outfile,"Fragment coordinate %d \n", id+I);
              }
              else {
                optinfo.nconstraints++;
              }
            }
          }
        }
      }

      if (optinfo.fix_intrafragment) { /* fix everything EXCEPT interfragment coordinates */
        for (i=0; i<simples.get_num(); ++i) {
          id = simples.index_to_id(i);
          simples.locate_id(id, &intco_type, &sub_index, &sub_index2);
          if (intco_type != FRAG) {
            if (!cnt) {
              constraints[iconstraints++] = simples.id_to_index(id);
              if (iconstraints == 1) fprintf(outfile,"Simple coordinates:\n");
              fprintf(outfile," %d", id);
            }
            else { 
              optinfo.nconstraints++;
            } 
          }
        }
      }

      if (cnt) fprintf(outfile,"%d found.\n",optinfo.nconstraints);
    } /* end cnt loop */

    ip_cwk_add(":DEFAULT");
    ip_cwk_add(":PSI");
    ip_cwk_add(":INTCO");
    ip_cwk_add(":OPTKING");
  } /* end if "FIXED_INTCO" */


  if (optinfo.nconstraints > 0)
    optinfo.constraints_present = 1;

  if (optinfo.constraints_present == 0)
    fprintf(outfile,"none found.\n");

  /* fixed?
  if ( (optinfo.constraints_present) && (optinfo.delocalize == 1) ) {
    fprintf(outfile,"Constraints may only be imposed on simple internal coordinates,");
    fprintf(outfile,"\nso coordinates cannot be delocalized.\n");
    exit(2);
  }
  */

  return constraints;
}

}//} /* namespace psi::optking */

