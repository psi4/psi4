/*! \defgroup INPUT input: Set up a PSI computation based on user input */

/*! 
** \file
  ** \ingroup INPUT
  ** \brief Set up a PSI computation based on user input
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cmath>

#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <chkpt_params.h>
#include <psifiles.h>
#include <liboptions/liboptions.h>
#include <physconst.h>

#include <Molecular_system.h>
#include "input.h"
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

  void print_intro();
  void print_options();
  void print_geometry(double);
  void print_full_geometry(double);
  void print_unique_geometry(double);
  void print_symm_trans();
  void print_basis_info();
  void cleanup();
  extern void build_cartdisp_salcs();
  void am_i_to_char(int am, char *am_label);

  extern void read_cart(Molecular_system &);

  using namespace psi;

/*
  void register_input_options()
  {
    // Reference wavefunction
    options.add("WFN", "SCF", "");
    // Keep the checkpoint file
    options.add("KEEP_CHKPT", (bool)false);
    // Read MOs from checkpoint and project onto new basis
    options.add("CHKPT_MOS", (bool)false);

    // These may be obseleted by psi4
    // Read geometry from checkpoin file (in findif calculations)
    options.add("CHKPT_GEOM", (bool)false);
    // Don't project MOs but simply keep them
    options.add("NOPROJECT", (bool)false);
    // Read geometry from geom.dat file (in findif calculations)
    options.add("GEOMDAT", (bool)false);
    // No center of mass shift
    options.add("NO_COM_SHIFT", (bool)false);
    // Read MOs from checkpoint file and save to a separate file
    options.add("SAVE_MOS", (bool)false);

    // Don't overwrite the output file
    // This should be removed. It is no longer input's responsiblity
    options.add("KEEP_OUTPUT", (bool)false);
    
    options.add("NO_REORIENT", (bool)false);
    options.add("LABEL","Default PSI3 Label");
    options.add("SHOWNORM", (bool)false);
    options.add("NORMALIZE", (bool)true);
    options.add("PUREAM", (bool)false);
    options.add("EXPERT", (bool)false);
    options.add("PRINT", 1);
    options.add("SUBGROUP", "", "C1 C2 CS CI C2V C2H D2 D2H");
    options.add("UNIQUE_AXIS", "", "X Y Z");
    options.add("NFRAGMENTS", 1);
    options.add("KEEP_REF_FRAME", (bool)false);
    options.add("FRAGMENT_DISTANCE_INVERSE", (bool)false);
    
    options.add("FREEZE_CORE", "FALSE", "FALSE NO TRUE YES SMALL LARGE");
    options.add("FREEZE_VIRT",0);
  }
*/
  
  int input() { return 0; }
  
  int input(Options & options, char **atom_basis, Molecular_system & molecules)
  {
   /*variables and arrays*/
    int argc = 0;
    char **argv = NULL;
    int i,j;
    int errcod;
    int is_dummy;
    int natomtri = 0;                    /*number of elements in the lower triangle of a (num_atoms)x(num_atoms) matrix */
    double *Distance=NULL;			/*Lower triangle of the internuclear distance matrix*/
    double repulsion;			/*nuclear Repulsion Energy*/
    FILE *pbasis = NULL;
    FILE *user_basis = NULL;
    FILE *local_basis = NULL;
    char *pbasis_dirname, *pbasis_filename;
    char *user_basis_filename, *user_basis_file;
    double **Rmat;
    double ***unitvecs;
    double ***bondangles;

     /*-------------------------------------
    Initialize files and parsing library
      -------------------------------------*/
      start_io(argc, argv);

    init_globals();
    wfn = options.get_str("WFN");
    parsing();
    print_intro();
     //if (options.get_int_option("PRINT"));
     //print_options();

     /* To find default basis set file first check the environment, then its location after installation */
    pbasis_dirname = getenv("PSIDATADIR");
    if (pbasis_dirname != NULL) {
      char* tmpstr = (char *) malloc(sizeof(char)*(strlen(pbasis_dirname)+12));
      sprintf(tmpstr,"%s/pbasis.dat",pbasis_dirname);
      pbasis_filename = tmpstr;
    }
    else
      pbasis_filename = strdup(INSTALLEDPSIDATADIR "/pbasis.dat");
    pbasis = fopen( pbasis_filename, "r");

    errcod = ip_string(const_cast<char*>("BASIS_FILE"),&user_basis_file,0);
    if (errcod == IPE_OK && strlen(user_basis_file) != 0) {
       /* Checking if only path to the file has been specified */
      if (!strncmp(user_basis_file+strlen(user_basis_file)-1,"/",1)) { /* Append default name - basis.dat */
        user_basis_filename = (char*)malloc(sizeof(char)*(strlen(user_basis_file)+10));
        strcpy(user_basis_filename,user_basis_file);
        user_basis_filename = strcat(user_basis_filename,"basis.dat");
      }
      else
        user_basis_filename = user_basis_file;
      user_basis = fopen(user_basis_filename, "r");
      if (user_basis != NULL) {
        ip_append(user_basis, outfile);
        fclose(user_basis);
        if (print_lvl > 0) fprintf(outfile,"\n  Parsed basis sets from %s\n",user_basis_filename);
      }
      free(user_basis_filename);
    }

    if (pbasis != NULL) {
      ip_append(pbasis, outfile);
      fclose(pbasis);
      if (print_lvl > 0) fprintf(outfile,"\n  Parsed basis sets from %s\n",pbasis_filename);
    }
    free(pbasis_filename);

    local_basis = fopen("./basis.dat", "r");
    if (local_basis != NULL) {
      ip_append(local_basis, outfile);
      fclose(local_basis);
      if (print_lvl > 0) fprintf(outfile,"\n  Parsed basis sets from basis.dat\n");
    }

    ip_cwk_add(const_cast<char*>(":BASIS"));

     // geometry is passed in

    read_cart(molecules);

    nuclear_charges = molecules.get_Z();

    // orient_fragments();

    freeze_core();
     //freeze_virt();

    repulsion = 0.0;
    if (num_atoms > 1) {
      natomtri = num_atoms*(num_atoms+1)/2;
      Distance = init_array(natomtri);
      calc_distance(geometry,Distance,num_atoms);
      Nuc_repulsion(Distance, &repulsion);
    }

    fprintf(outfile,"\n  -Geometry before Center-of-Mass shift (a.u.):\n");
    print_geometry(1.0);

    reorient();
    fprintf(outfile,"\n  -Geometry after Center-of-Mass shift and reorientation (a.u.):\n");
    print_geometry(1.0);

     /*--------------
    Get symmetric
      --------------*/
      find_symmetry();
    count_uniques();
    print_symm_trans();

     /*---------------------------------------
    Canonical frame is the reference frame
      unless read geometry from chkpt file
      ---------------------------------------*/
      if (keep_ref_frame == 0)
      canon_eq_ref_frame();

     /*---------------------
    Parse basis set data
      ---------------------*/

      read_basis(atom_basis);

     /*----------------------------------------------
    Form symmetry information arrays of all kinds
      ----------------------------------------------*/

      build_transmat();
    if (puream)
      build_cart2pureang();
    build_so_classes();
    build_usotao();
    build_cartdisp_salcs();

     /*---------------------
    Read old calculation
      ---------------------*/
    if (read_chkpt) {
      init_oldcalc();
    }

     /*-------------------------------------------------
    Write the information out to the checkpoint file
      -------------------------------------------------*/
      write_to_chkpt(repulsion);

     /*-------------------------------------------------
    Project old MOs onto new basis and write 'em out
      -------------------------------------------------*/
    if (chkpt_mos) {
      oldcalc_projection();
      write_scf_to_chkpt();
    }

     /*-------------------------------------------------
    Write old MOs and old/new basis overlap info out
      -------------------------------------------------*/
    if (save_oldcalc) {
      store_oldcalc();
    }

     /*------------------------------
    Done with the old calculation
      ------------------------------*/
    if (read_chkpt) {
      cleanup_oldcalc();
    }

     /*----------------------
    Print geometries etc.
      ----------------------*/

      print_basis_info();
    fprintf(outfile,"\n  -Unique atoms in the canonical coordinate system (a.u.):\n");
    print_unique_geometry(1.0);
    fprintf(outfile,"\n  -Geometry in the canonical coordinate system (a.u.):\n");
    print_geometry(1.0);
    fprintf(outfile,"\n  -Geometry in the canonical coordinate system (Angstrom):\n");
    print_geometry(_bohr2angstroms);
     /* Print geometry including dummy atoms, if necessary */
    is_dummy = 0;
    for(i=0;i<num_allatoms;++i) 
      if(!strncmp(full_element[i],"X\0",2) )
      is_dummy = 1;
    if(is_dummy) {
      fprintf(outfile,"\n  -Full geometry in the canonical coordinate system (a.u.):\n");
      print_full_geometry(1.0);
    }
    rotate_full_geom(Rref);
    fprintf(outfile,"\n  -Geometry in the reference coordinate system (a.u.):\n");
    print_geometry(1.0);
    fprintf(outfile,"\n  --------------------------------------------------------------------------\n");

    if (num_atoms > 1) {
       /* Print the Nuclear Repulsion Energy */
      fprintf(outfile,"\n    Nuclear Repulsion Energy (a.u.) = %20.12f\n\n",repulsion);

       /* Print the interatomic Distance Matrix in angstroms */
      for(i=0;i<natomtri;i++)
        Distance[i] *= _bohr2angstroms;
      fprintf(outfile,"  -The Interatomic Distances in angstroms:\n");
      print_array(Distance,num_atoms,outfile);
      if(print_lvl < 3) {
        fprintf(outfile,"\n    Note: To print *all* bond angles, out-of-plane\n");
        fprintf(outfile,"          angles, and torsion angles set print = 3\n");
      }
      if(num_atoms > 2 && print_lvl >= 3) { 
        Rmat = init_matrix(num_atoms,num_atoms);
        tri_to_sq(Distance,Rmat,num_atoms);
   /*print_mat(Rmat,num_atoms,num_atoms,outfile);*/
        unitvecs = unit_vectors(full_geom, Rmat);
        bondangles = calc_bond_angles(unitvecs, Rmat);
        if(num_atoms > 3 && print_lvl >= 3) {
          calc_oop_angles(unitvecs, bondangles, Rmat);
          calc_tors_angles(unitvecs, bondangles, Rmat);
        }

        free_matrix(Rmat, num_atoms);
         /* Free Memory for 3D unitvecs */
        for(i=0; i<num_atoms; i++) {
          for(j=0; j<num_atoms; j++) {
            free(unitvecs[i][j]);
          }
        }
        for(i=0; i<num_atoms; i++) {
          free(unitvecs[i]);
        }
        free(unitvecs);
         /* Free Memory for 3D bondangles */
        for(i=0; i<num_atoms; i++) {
          for(j=0; j<num_atoms; j++) {
            free(bondangles[i][j]);
          }
        }
        for(i=0; i<num_atoms; i++) {
          free(bondangles[i]);
        }
        free(bondangles);
      }

      free(Distance);
       //fprintf(outfile,"\n\n");
    }

    cleanup();
    stop_io();
    return(PSI_RETURN_SUCCESS);
  }

  void print_intro()
  {
    if (print_lvl) {
      fprintf(outfile,"                                --------------\n");
      fprintf(outfile,"                                    INPUT\n");
      fprintf(outfile,"                                --------------\n\n");
      fflush(outfile);
    }

    return;
  }


  void print_options()
  {
    if (print_lvl) {
      if (label.size())
        fprintf(outfile,"  LABEL       =\t%s\n", label.c_str());
      fprintf(outfile,"  SHOWNORM    =\t%d\n",shownorm);
      fprintf(outfile,"  PUREAM      =\t%d\n",puream);
      if (subgroup.size() > 0) {
        fprintf(outfile,"  SUBGROUP    =\t%s\n",subgroup.c_str());
        if (unique_axis.size() > 0)
          fprintf(outfile,"  UNIQUE_AXIS =\t%s\n",unique_axis.c_str());
      }
      fprintf(outfile,"  PRINT_LVL   =\t%d\n",print_lvl);
      if (chkpt_mos) {
        fprintf(outfile,"  Will read old MO vector from the checkpoint file\n");
        fprintf(outfile,"  and project onto the new basis.\n");
      }
      fflush(outfile);
    }

    return;
  }


  void print_geometry(double conv_factor)
  {
    int i,j;

    fprintf(outfile,"       Center              X                  Y                   Z\n");
    fprintf(outfile,"    ------------   -----------------  -----------------  -----------------\n");

    for(i=0;i<num_atoms;i++){
      fprintf(outfile,"    %12s ",element[i]); fflush(outfile);
      for(j=0;j<3;j++)
        fprintf(outfile,"  %17.12f",geometry[i][j]*conv_factor);
      fprintf(outfile,"\n");
    }
    fprintf(outfile,"\n");
    fflush(outfile);

    return;
  }


  void print_full_geometry(double conv_factor)
  {
    int i,j;

    fprintf(outfile,"       Center              X                  Y                   Z\n");
    fprintf(outfile,"    ------------   -----------------  -----------------  -----------------\n");

    for(i=0;i<num_allatoms;i++){
      fprintf(outfile,"  %12s ",full_element[i]); fflush(outfile);
      for(j=0;j<3;j++)
        fprintf(outfile,"  %17.12f",full_geom[i][j]*conv_factor);
      fprintf(outfile,"\n");
    }
    fflush(outfile);

    return;
  }


  void print_unique_geometry(double conv_factor)
  {
    int i,j;

    fprintf(outfile,"       Center              X                  Y                   Z\n");
    fprintf(outfile,"    ------------   -----------------  -----------------  -----------------\n");
    for(i=0;i<num_uniques;i++){
      fprintf(outfile,"    %12s ",element[u2a[i]]);
      for(j=0;j<3;j++)
        fprintf(outfile,"  %17.12f",geometry[u2a[i]][j]*conv_factor);
      fprintf(outfile,"\n");
    }
    fprintf(outfile,"\n");
    fflush(outfile);

    return;
  }


  void print_symm_trans()
  {
    int i,j;

    fprintf(outfile,"\n  -SYMMETRY INFORMATION:\n");
    fprintf(outfile,"    Computational point group is %s\n",symmetry);
    fprintf(outfile,"    Number of irr. rep.      = %d\n",nirreps);
    fprintf(outfile,"    Number of atoms          = %d\n",num_atoms);
    fprintf(outfile,"    Number of unique atoms   = %d\n\n",num_uniques);
    if (print_lvl >= PRINTSYMINFO) {
      for(i=0;i<num_uniques;i++)
        fprintf(outfile,"    Unique atom %d generates %d atoms\n",i+1,unique_degen[i]);
      fprintf(outfile,"\n    Orbits of unique atoms:\n");
      for(i=0;i<num_uniques;i++) {
        fprintf(outfile,"     %d    ",u2a[i]+1);
        for(j=0;j<nirreps;j++)
          fprintf(outfile,"%d  ",atom_orbit[u2a[i]][j]+1);
        fprintf(outfile,"\n");
      }
      fprintf(outfile,"\n");
      fprintf(outfile,"    Number of classes        = %d\n",num_classes);
      fprintf(outfile,"    Number of unique classes = %d\n\n",num_unique_classes);
      fprintf(outfile,"    Atom   Class   Symmetry position\n");
      fprintf(outfile,"    ----   -----   -----------------\n");
      for(i=0;i<num_atoms;i++)
        fprintf(outfile,"    %3d    %3d     %10d\n",i+1,atom_class[i]+1,atom_position[i]);
      fprintf(outfile,"\n");
      fprintf(outfile,"    Orbits of classes:\n");
      for(i=0;i<num_classes;i++) {
        fprintf(outfile,"     %d    ",i+1);
        for(j=0;j<nirreps;j++)
          fprintf(outfile,"%d  ",class_orbit[i][j]+1);
        fprintf(outfile,"\n");
      }
      fprintf(outfile,"\n");
    }
    fflush(outfile);

    return;
  }


  void print_basis_info()
  {
    int i,jshell,kshell,jprim,kprim,am,first,last,cnt_bf,redundant;
    double jexp, kexp;
    char am_label[2];

    fprintf(outfile,"  -BASIS SET INFORMATION:\n");
    fprintf(outfile,"    Total number of shells = %d\n",num_shells);
    fprintf(outfile,"    Number of primitives   = %d\n",num_prims);
    fprintf(outfile,"    Number of AO           = %d\n",num_ao);
    fprintf(outfile,"    Number of SO           = %d\n\n",num_so);
    fprintf(outfile,"    Irrep    Number of SO\n");
    fprintf(outfile,"    -----    ------------\n");
    for(i=0;i<nirreps;i++)
      fprintf(outfile,"    %3d      %7d\n",i+1,num_so_per_irrep[i]);

    fprintf(outfile,"\n  -Contraction Scheme:\n");
    fprintf(outfile,"    Atom     All Primitives // Unique Primitives // Shells\n");
    fprintf(outfile,"    ----     ---------------------------------------------\n");
    for(i=0;i<num_atoms;i++) {
      fprintf(outfile,"    %4d     ", i+1);  
      first = first_shell_on_atom[i];
      last = first + nshells_per_atom[i];
    /* print out # and am of primitives per atom */
      for(am=0; am<max_angmom+1; ++am) {
        cnt_bf = 0;
        for(jshell=first;jshell<last;jshell++) {
          if (am == shell_ang_mom[jshell])
            cnt_bf += nprim_in_shell[jshell];
        }
        if(cnt_bf > 0) {
          am_i_to_char(am, am_label);
          fprintf(outfile," %d%s", cnt_bf, am_label);
        }
      }
      fprintf(outfile," // ");
    /* print out # and am of unique primitives per atom */
      for(am=0; am<max_angmom+1; ++am) {
        cnt_bf = 0;
        for(jshell=first;jshell<last;jshell++) {
          if (am == shell_ang_mom[jshell]) {
            for(jprim=0; jprim<nprim_in_shell[jshell]; ++jprim) {
              jexp = exponents[first_prim_shell[jshell]+jprim];
              redundant = 0;
              for(kshell=first; kshell<jshell; kshell++) {
                for(kprim=0; kprim<nprim_in_shell[kshell]; ++kprim) {
                  kexp = exponents[first_prim_shell[kshell]+kprim];
                  if (jexp == kexp)
                    redundant = 1;
                }
              }
              if (!redundant) ++cnt_bf;
            }
          }
        }
        if(cnt_bf > 0) {
          am_i_to_char(am, am_label);
          fprintf(outfile," %d%s", cnt_bf, am_label);
        }
      }
      fprintf(outfile," // ");
    /* print out # and am of shells per atom */
      for(am=0; am<max_angmom+1; ++am) {
        cnt_bf = 0;
        for(jshell=first;jshell<last;jshell++) {
          if (am == shell_ang_mom[jshell]) ++cnt_bf;
        }
        if(cnt_bf > 0) {
          am_i_to_char(am, am_label);
          fprintf(outfile," %d%s", cnt_bf, am_label);
        }
      }
      fprintf(outfile,"\n");
    }

    fprintf(outfile,"\n");
    if (print_lvl >= DEBUGPRINT) {
      fprintf(outfile,"    Prim#     Exponent     Norm. contr. coeff.\n");
      fprintf(outfile,"    -----  --------------  -------------------\n");
      for(i=0;i<num_prims;i++)
        fprintf(outfile,"    %3d    %14.7f    %15.10f\n",i+1,exponents[i],contr_coeff[i]);
      fprintf(outfile,"\n");
      fprintf(outfile,"    Shell#  Nuc#  L  SPRIM  SLOC  SNUMG\n");
      fprintf(outfile,"    ------  ----  -  -----  ----  -----\n");
      for(i=0;i<num_shells;i++)
        fprintf(outfile,"    %4d    %3d  %2d  %3d    %3d   %3d\n",i+1,
        shell_nucleus[i]+1,shell_ang_mom[i],first_prim_shell[i]+1,first_basisfn_shell[i]+1,nprim_in_shell[i]);
      fprintf(outfile,"\n");
      fprintf(outfile,"    AMOrd. Shell#  Canon. Shell#\n");
      fprintf(outfile,"    -------------  -------------\n");
      for(i=0;i<num_shells;i++)
        fprintf(outfile,"    %7d        %7d\n",i+1,am2canon_shell_order[i]+1);
      fprintf(outfile,"\n\n");
    }
    fflush(outfile);

    return;
  }


  void cleanup()
  {
    int i,l,uc,class_curr,class_first,class_last;

  /*----------------------
    Free up global arrays
      ----------------------*/

      free(sym_oper);
    free_block(full_geom);
    free(geometry);
    free_block(Rref);
    free(nuclear_charges);
    free_char_matrix(element,num_atoms);
    free_char_matrix(full_element,num_allatoms);
    free_int_matrix(atom_orbit);
    free_int_matrix(class_orbit);
    free_int_matrix(red_unique_orbit);
    for(uc=0;uc<num_unique_classes;uc++) {
      class_curr = uc2c[uc];
      class_first = class_curr;
      class_last = class_curr + unique_class_degen[uc];
      for(i=class_first;i<class_last;i++)
      for(l=0;l<=max_angmom_class[i];l++) {
        free_matrix(class_so_coeff[i][l],ioff[l+1]*unique_class_degen[uc]);
      }
    }
    for(i=0;i<num_classes;i++) {
      free(class_so_coeff[i]);
      free_int_matrix(num_cart_so_in_class[i]);
      if (puream)
        free_int_matrix(num_pureang_so_in_class[i]);
    }
    free(class_so_coeff);
    free(num_cart_so_in_class);
    if (puream)
      free(num_pureang_so_in_class);
    free(max_angmom_class);
    free(u2a);
    free(uc2c);
    free(unique_degen);
    free(unique_class_degen);
    free(num_cart_so_per_irrep);
    free_int_matrix(num_cart_so);
    if (puream) {
      free(num_so_per_irrep);
      free_int_matrix(num_pureang_so);
      free_int_matrix(num_redun_so);
    }
    free(nshells_per_atom);
    free(first_shell_on_atom);
    free(exponents);
    free(contr_coeff);
    free(first_prim_shell);
    free(shell_nucleus);
    free(shell_ang_mom);
    free(nprim_in_shell);
    free(first_ao_shell);
    free(shells_per_am);
    free(am2canon_shell_order);
    for(i=0;i<MAX(max_angmom+1,MAXANGMOM);i++)
      free_matrix(ao_type_transmat[i],nirreps);
    free(ao_type_transmat);
    if (puream) {
      free(pureang_so_m);
      for(i=0;i<=max_angmom;i++)
        free_matrix(cart2pureang[i],2*i+1);
      free(cart2pureang);
    }
    free_matrix(usotao,num_so);
    free(atom_class);
    free(atom_position);
    free_int_matrix(ao_type_irr);
    free(ioff);
    free(df);
    free(z_geom);
    free(frag_num_atoms);
    free(frag_num_allatoms);
    free(frag_atom);
    free(frag_allatom);
    free(nref_per_fragment);
    for (i=0; i<(nfragments-1); ++i) {
      if (nref_per_fragment[i] > 0)
        free_block(ref_pts_lc[i]);
    }
    free(ref_pts_lc);

    return;
  }

  void am_i_to_char(int am, char *am_label)
  {
    if(am == 0) am_label[0] = 's';
    else if(am == 1) am_label[0] = 'p';
    else if(am == 2) am_label[0] = 'd';
    else if(am == 3) am_label[0] = 'f';
    else if(am == 4) am_label[0] = 'g';
    else if(am == 5) am_label[0] = 'h';
    else if(am == 6) am_label[0] = 'i';
    else if(am == 7) am_label[0] = 'k';
    else if(am == 8) am_label[0] = 'l';
    else if(am == 9) am_label[0] = 'm';
    else if(am == 10) am_label[0] = 'n';
    else am_label[0] = 'x';
    am_label[1] = '\0';
    return;
  }

  }} // namespace psi::input
