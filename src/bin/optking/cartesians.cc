/*! \file
    \ingroup OPTKING
    \brief cartesians.cc : contains member functions for cartesian class
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"

#include <libipv1/ip_lib.h>
#include <libchkpt/chkpt.h>
#include <masses.h>

namespace psi { namespace optking {

/*** FORCES returns forces in cartesian coordinates in aJ/Ang ***/ 
double *cartesians:: get_forces() const {
  int i;
  double *f;
  f = init_array(3*natom);
  for (i=0;i<3*natom;++i)
    f[i] = -1.0 * grad[i] * _hartree2J * 1.0E18 / _bohr2angstroms;
  return f;
}

double *cartesians:: get_fforces() const {
  int i;
  double *f;
  f = init_array(3*nallatom);
  for (i=0; i<3*nallatom; ++i) {
    f[i] = -1.0 * fgrad[i] * _hartree2J * 1.0E18 / _bohr2angstroms;
  }
  return f;
}


/*** CARTESIANS this is the class constructor for cartesian ***/ 
cartesians::cartesians() {
  int i, j, a, count, isotopes_given = 0, masses_given = 0;
  char label[MAX_LINELENGTH], line1[MAX_LINELENGTH];
  FILE *fp_11;
  double an,x,y,z,tval,tval1,tval2,tval3,tval4;
  double **geom, **fgeom, *zvals;

  natom = optinfo.natom;
  nallatom = optinfo.nallatom;

  /* All data (but maybe masses) now read from chkpt */
  chkpt_init(PSIO_OPEN_OLD);
  geom = chkpt_rd_geom();
  fgeom = chkpt_rd_fgeom();

  if ((optinfo.mode == MODE_OPT_STEP) || (optinfo.mode == MODE_GRAD_SAVE))
    grad = chkpt_rd_grad();
  else grad = init_array(3*natom);

  if ((optinfo.mode == MODE_OPT_STEP) || (optinfo.mode == MODE_ENERGY_SAVE)
    ||(optinfo.mode == MODE_DISP_NOSYMM) || (optinfo.mode == MODE_DISP_IRREP)
//    ||(optinfo.mode == MODE_DISP_LOAD) || (optinfo.mode == MODE_DISP_USER)
    || (optinfo.mode == MODE_DISP_USER)
    ||(optinfo.mode == MODE_DISP_FREQ_GRAD_CART) || (optinfo.mode == MODE_FREQ_GRAD_CART)
    ||(optinfo.mode == MODE_DISP_FREQ_ENERGY_CART) || (optinfo.mode == MODE_FREQ_ENERGY_CART)
    )

  //ACS (11/06) Allow external programs to be used
  if(optinfo.external_energies){
    energy = 0.0;
  }else{
    energy = chkpt_rd_etot();
  }

  zvals = chkpt_rd_zvals();
  chkpt_close();

  coord = new double[3*natom];
  atomic_num = new double[natom];
  mass = new double[3*natom];
  fcoord = new double[3*nallatom];
  fgrad = new double[3*nallatom];
  fatomic_num = new double[nallatom];
  fmass = new double[3*nallatom];

  count = -1;
  for(i=0; i < natom; i++) {
    atomic_num[i] = zvals[i];
    for(j=0; j < 3; j++) {
      coord[++count] = geom[i][j];
    }
  }
  free_matrix(geom);

  zero_array(fatomic_num,nallatom);

  for(i=0; i<natom; i++)
    fatomic_num[optinfo.to_dummy[i]] = zvals[i];

  count = -1;
  for(i=0; i<nallatom; i++)
    for(j=0; j<3; j++)
      fcoord[++count] = fgeom[i][j];

  free_array(zvals);
  free_matrix(fgeom);

  zero_array(fgrad, 3*nallatom);
  count = -1;
  for(i=0; i < natom; i++)
    for(j=0; j < 3; j++)
      fgrad[3*optinfo.to_dummy[i]+j]  = grad[3*i+j];

  /* read masses from input.dat or use default masses */
  count = -1;

  if (ip_exist("ISOTOPES",0)) {
    isotopes_given = 1;
    a = 0;
    ip_count("ISOTOPES", &a, 0);
    if (a != natom) {
      fprintf(outfile,"ISOTOPES array has wrong dimension.\n");
      exit(2);
    }
    for (i=0;i<natom;++i) {
      ip_data("ISOTOPES","%s", line1,1,i);
      for (j=0;j<138;j++) {
	if (strcmp(line1, mass_labels[j]) == 0) {
	  mass[++count] = atomic_masses[j];
	  mass[++count] = atomic_masses[j];
	  mass[++count] = atomic_masses[j];
	  break;
	}
      }
      if (j == 138) {
	fprintf(outfile,
		"Isotope label %s is unidentifiable.\n",mass_labels[j]);
	exit(2);
      }
    }
  }
  if (ip_exist("MASSES",0)) {
    masses_given = 1;
    if (isotopes_given)
      fprintf(outfile,"Ignoring ISOTOPES keyword, using given MASSES.\n");
    a = 0;
    ip_count("MASSES",&a,0);
    if (a != natom) {
      fprintf(outfile,"MASSES array has wrong dimension\n");
      exit(2);
    }
    else {
      for(i=0;i<natom;++i) {
        ip_data("MASSES","%lf",&tval,1,i);
        mass[++count] = tval;
        mass[++count] = tval;
        mass[++count] = tval;
      }
    }
  }
  if ((isotopes_given == 0) && (masses_given == 0)) {
    for(i=0;i<natom;++i) {
      a = (int) get_Z(i); // casting to an int for index
      mass[++count] = an2masses[a];
      mass[++count] = an2masses[a];
      mass[++count] = an2masses[a];
    }
  }

  zero_array(fmass, 3*nallatom);
  for (i=0; i<natom; ++i) {
    j = optinfo.to_dummy[i];
    fmass[3*j+0] = fmass[3*j+1] = fmass[3*j+2] = mass[3*i];
  }

/*
  fprintf(outfile,"masses without dummy atoms\n");
  for (i=0; i<nallatom; ++i)
    fprintf(outfile,"%15.10lf%15.10lf%15.10lf\n",
        fmass[3*i],fmass[3*i+1],fmass[3*i+2]);
     */

  return;
}

/*** CARTESIANS::PRINT
 * print_flag                                         to
 *    0        fcoord                                fp_out
 *    1        fatomic_num, fcoord                   fp_out
 *    2        fatomic_num, fmass, fcoord            fp_out
 *    3        fatomic_num, fmass, fcoord, fgrad     fp_out
 *    4        coord                                 fp_out
 *    5        atomic_num, coord                     fp_out
 *    6        atomic_num, mass, coord               fp_out
 *    7        atomic_num, mass, coord, grad         fp_out
 *   10        coord                                 geom.dat
 *   11        coord, grad                           file11.dat
 *   32        fcoord (implicitly coord)             chkpt
 *   12        fatomic_symb, fcoord                  fp_out
 *   13        fatomic_symb, fcoord ANGSTROM         fp_out
 *   23        fatomic_num, fmass, fcoord, fgrad, cart_rms_force, cart_max_force        fp_out
 * disp_label is only used for geom.dat writing ***/

void cartesians :: print(int print_flag, FILE *fp_out, int new_geom_file,
			 char *disp_label, int disp_num) const {
  int i,j;
  double x,y,z;
  double crms, cmax;
  int cnt = -1;
  const char *sym;

  if (print_flag == 0) {
    for (i = 0; i < natom; ++i) {
      x = coord[++cnt]; y = coord[++cnt]; z = coord[++cnt];
      fprintf(fp_out,"%20.10f%20.10f%20.10f\n",x,y,z);
    }
  }
  else if (print_flag == 1) {
    for (i = 0; i < natom; ++i) {
      x = coord[++cnt]; y = coord[++cnt]; z = coord[++cnt];
      fprintf(fp_out,"%5.1lf%15.10f%15.10f%15.10f\n",atomic_num[i],x,y,z);
    }
  }
  else if (print_flag == 2) {
    cnt = -1;
    for (i = 0; i < natom; ++i) {
      x = coord[++cnt]; y = coord[++cnt]; z = coord[++cnt];
      fprintf(fp_out,
          "%5.1lf%15.8lf%15.10f%15.10f%15.10f\n",atomic_num[i],mass[3*i],x,y,z);
    }
  }
  else if (print_flag == 3) {
    cnt = -1;
    for (i = 0; i < natom; ++i) {
      x = coord[++cnt]; y = coord[++cnt]; z = coord[++cnt];
      fprintf(fp_out,
          "%5.1lf%15.8lf%15.10f%15.10f%15.10f\n",atomic_num[i],mass[3*i],x,y,z);
    }
    cnt = -1;
    for (i = 0; i < natom; ++i) {
      x = grad[++cnt]; y = grad[++cnt]; z = grad[++cnt];
      fprintf(fp_out,
	      "%35.10lf%15.10f%15.10f\n",x,y,z);
    }
  }
  if (print_flag == 4) {
    for (i = 0; i < natom; ++i) {
      x = coord[++cnt]; y = coord[++cnt]; z = coord[++cnt];
      fprintf(fp_out,"%20.10f%20.10f%20.10f\n",x,y,z);
    }
  } 
  else if (print_flag == 5) {
    for (i = 0; i < natom; ++i) {
      x = coord[++cnt]; y = coord[++cnt]; z = coord[++cnt];
      fprintf(fp_out,"%5.1lf%15.10f%15.10f%15.10f\n",atomic_num[i],x,y,z);
    }
  }
  else if (print_flag == 6) {
    cnt = -1;
    for (i = 0; i < natom; ++i) {
      x = coord[++cnt]; y = coord[++cnt]; z = coord[++cnt];
      fprintf(fp_out,
          "%5.1lf%15.8lf%15.10f%15.10f%15.10f\n",atomic_num[i],mass[3*i],x,y,z);
    }
  }
  else if (print_flag == 7) {
    cnt = -1;
    for (i = 0; i < natom; ++i) {
      x = coord[++cnt]; y = coord[++cnt]; z = coord[++cnt];
      fprintf(fp_out,
          "%5.1lf%15.8lf%15.10f%15.10f%15.10f\n",atomic_num[i],mass[3*i],x,y,z);
    }
    cnt = -1;
    for (i = 0; i < natom; ++i) {
      x = grad[++cnt]; y = grad[++cnt]; z = grad[++cnt];
      fprintf(fp_out,
              "%35.10lf%15.10f%15.10f\n",x,y,z);
    }
  }
  else if (print_flag == 10) {
    //if (new_geom_file) fprintf(fp_out,"%%%%\n");
    fprintf(fp_out,"%% %s\n",disp_label);
    fprintf(fp_out,"geometry%1d = (\n", disp_num);
    cnt = -1;
    for (i=0; i<natom; ++i) {
      x = coord[++cnt]; y = coord[++cnt]; z = coord[++cnt];
      fprintf(fp_out,"(%20.10f%20.10f%20.10f)\n",x,y,z);
    }
    fprintf(fp_out," )\n");
  }
  else if (print_flag == 11) {
    fprintf(fp_out,"%s\n",disp_label);
    fprintf(fp_out,"%5d%20.10lf\n",natom,energy);
    cnt = -1;
    for (i = 0; i < natom; ++i) {
      x = coord[++cnt]; y = coord[++cnt]; z = coord[++cnt];
      fprintf(fp_out,
	      "%20.10lf%20.10lf%20.10lf%20.10lf\n",atomic_num[i],x,y,z);
    }
    cnt = -1;
    for (i = 0; i < natom; ++i) {
      x = grad[++cnt]; y = grad[++cnt]; z = grad[++cnt];
      fprintf(fp_out,"%40.10lf%20.10lf%20.10lf\n",x,y,z);
    }
  }
  else if (print_flag == 32) {

    fprintf(outfile,"\nGeometry written to chkpt\n");

    // send copy so print function can be 'const'
    double *fcoord_tmp = get_fcoord();
    chkpt_init(PSIO_OPEN_OLD);
    chkpt_wt_fgeom(&fcoord_tmp);
    chkpt_close();
    free_array(fcoord_tmp);

  }
  else if (print_flag == 12) { 
    for (i = 0; i < natom; ++i) {
      x = coord[++cnt]; y = coord[++cnt]; z = coord[++cnt];
      sym = atomic_labels[(int) (atomic_num[i])];
      fprintf(fp_out,"  (%3s %15.10f %15.10f %15.10f )\n",sym,x,y,z);
    }
  }
  else if (print_flag == 13) { 
    for (i = 0; i < natom; ++i) {
      x = coord[++cnt] * _bohr2angstroms; 
      y = coord[++cnt] * _bohr2angstroms; 
      z = coord[++cnt] * _bohr2angstroms;
      sym = atomic_labels[(int) (atomic_num[i])];
      fprintf(fp_out,"  (%3s %15.10f %15.10f %15.10f )\n",sym,x,y,z);
    }
  }
  else if (print_flag == 23) {
    cnt = -1;
    for (i = 0; i < natom; ++i) {
      x = coord[++cnt]; y = coord[++cnt]; z = coord[++cnt];
      fprintf(fp_out,
          "%5.1lf%15.8lf%15.10f%15.10f%15.10f\n",atomic_num[i],mass[3*i],x,y,z);
    }
    cmax = 0.0;
    crms = 0.0;
    cnt = -1;
    for (i = 0; i < natom; ++i) {
      x = grad[++cnt]; y = grad[++cnt]; z = grad[++cnt];
      fprintf(fp_out,
              "%35.10lf%15.10f%15.10f\n",x,y,z);
      crms += x*x + y*y + z*z;
      if (fabs(x) > cmax) cmax = fabs(x);
      if (fabs(y) > cmax) cmax = fabs(y);
      if (fabs(z) > cmax) cmax = fabs(z);
    }
    crms = sqrt(crms/((double) 3.0*natom));
    fprintf(outfile,"\n MAX cartesian force: %15.10lf   RMS force: %15.10lf\n", cmax, crms);
  }
  return;
}

}} /* namespace psi::optking */
