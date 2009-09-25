/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <strings.h>
#include <string.h>
#include "displacements.h"
#include "params.h"

#define EXTERN
#include "globals.h"
#undef EXTERN

#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <psifiles.h>
#include <physconst.h>
#include <masses.h>

using namespace psi::intder;

namespace psi { namespace intder {
extern Params gParams;
}}

bool Displacements::loadFromCheckPoint()
{
  int natom, nallatoms;
  int index, atomcount;
  double *zvals = NULL;
  double **coord;
  double **movedcoord;
  double *movedzvals;

  natom      = chkpt_rd_natom();          // Number of atoms
  nallatoms  = chkpt_rd_nallatom();       // Number of atoms plus dummies
  zvals      = chkpt_rd_zvals();          // Nuclear charge of atoms (does not include dummies)
  movedzvals = new double[nallatoms];

  coord      = chkpt_rd_fgeom();          // Full geom includes dummies
  movedcoord = chkpt_rd_fgeom();

  // move the geometry to the moved_dummy ordering
  for (index = 0; index < nallatoms; index++)
  {
    memcpy(movedcoord[gParams.moved_dummy[index]], coord[index], 3*sizeof(double));
    movedzvals[index] = 0.0;
  }

  // Since zvals does not include dummy atoms, must do this
  atomcount = 0;
  for (index = 0; index < nallatoms; index++)
  {
    if (!gParams.atom_dummy[index]) {
      movedzvals[gParams.moved_dummy[index]] = zvals[atomcount];
      atomcount++;
    }
  }

  Molecule ref;
  // modifed to handle dummy atoms.
  for (index = 0; index < natom; index++)
  {
    ref.addAtom((int)movedzvals[index], an2masses[(int)movedzvals[index]],
                movedcoord[index][0]*_bohr2angstroms,
                movedcoord[index][1]*_bohr2angstroms,
                movedcoord[index][2]*_bohr2angstroms);
  }
  for (index = natom; index < nallatoms; index++)
  {
    ref.addAtom((int)0, 0.00,
                movedcoord[index][0]*_bohr2angstroms,
                movedcoord[index][1]*_bohr2angstroms,
                movedcoord[index][2]*_bohr2angstroms);
  }
/*  for (index = 0; index < nallatoms; index++)
  {
    if (!gParams.atom_dummy[index]) {
      ref.addAtom((int)zvals[gParams.moved_dummy[atomcount]], an2masses[(int)zvals[gParams.moved_dummy[atomcount]]],
                  movedcoord[index][0]*_bohr2angstroms,
                  movedcoord[index][1]*_bohr2angstroms,
                  movedcoord[index][2]*_bohr2angstroms);
      atomcount++;
    }
    else {
      ref.addAtom((int)0, 0.00,
                  movedcoord[index][0]*_bohr2angstroms,
                  movedcoord[index][1]*_bohr2angstroms,
                  movedcoord[index][2]*_bohr2angstroms);
    }
  }
*/
  
  // This is the kind of C++ crap made by Justin that is nonsense. Why is this here? What does it do?
  addDisplacement(ref);
  
  free_block(coord);
  free_block(movedcoord);
  free(movedzvals);
  free(zvals);
}

// Possibly needs to be modified to handle dummy atoms
bool Displacements::loadFromOptKing()
{
  int natom = 0;
  int ndisp = 0;
  int curdisp = 0;
  int index = 0;
  double *zvals = NULL;
  double *coord = NULL;
  double **coords = NULL;

  natom = chkpt_rd_natom();
  zvals = chkpt_rd_zvals();

  psio_open(PSIF_OPTKING, PSIO_OPEN_OLD);
  psio_read_entry(PSIF_OPTKING, "OPT: Total num. of disp.",
                  (char*)&(ndisp), sizeof(int));
  printf("\n%d displacements found in OPTKING file.\n", ndisp);

  coord = new double[3*natom];
  psio_read_entry(PSIF_OPTKING, "OPT: Reference geometry",
                  (char*)&(coord[0]), 3*natom*sizeof(double));

  Molecule ref;
  for (index = 0; index < natom; index++)
  {
    ref.addAtom((int)zvals[index], an2masses[(int)zvals[index]],
                coord[3*index]*_bohr2angstroms,
                coord[3*index+1]*_bohr2angstroms,
                coord[3*index+2]*_bohr2angstroms);
  }
  addDisplacement(ref);

  coords = block_matrix(ndisp, 3*natom);
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced geometries",
                  (char*)&(coords[0][0]), 3*natom*ndisp*sizeof(double));

  for (curdisp = 0; curdisp < ndisp; curdisp++)
  {
    Molecule dispMolecule;
    for (index = 0; index < natom; index++)
    {
      dispMolecule.addAtom((int)zvals[index], an2masses[(int)zvals[index]],
                           coords[curdisp][3*index]*_bohr2angstroms,
                           coords[curdisp][3*index+1]*_bohr2angstroms,
                           coords[curdisp][3*index+2]*_bohr2angstroms);
    }
    addDisplacement(dispMolecule);
  }
  psio_close(PSIF_OPTKING, 1);

  delete[] zvals;
  delete[] coord;
  free_block(coords);

  return true; 
}

//This still doesn't work :<
bool Displacements::loadFromInput(void)
{
  FILE *geomdat;
  int natom = 0;
  int ndisp = 0;
  int curdisp = 0;
  int index = 0;
  
  int i, j, errcod;
  double Z = 0.0;
  double tmp1 = 0.0;
  double tmp2 = 0.0;
  double tmp3 = 0.0;

  char entry_name[20];
  char *atom_label;
  int entry_length = 0;
  int simple_geom = 1;
  int atomcounter = 0;
  
  ffile_noexit(&geomdat, "geom.dat", 2);

  if ( ip_exist(":GEOMETRY",0) ) {
    sprintf(entry_name, "GEOMETRY", 0);
    ip_count(entry_name, &atomcounter, 0);
    for(i=0; i < atomcounter; i++) {
      ip_count("GEOMETRY", &entry_length, 1, i);
      if(entry_length > 1) simple_geom = 0;
    }
  }
  else if (geomdat != NULL) {
    ip_append(geomdat, outfile);
    sprintf(entry_name,"GEOMETRY%d", ndisp);
    ip_count(entry_name,&atomcounter,0);
    simple_geom = 0;
    if (atomcounter == 0) {
      printf("\nThe entry in geom.dat is empty or missing!\n");
      return false; 
    }
    fprintf(outfile, "\nsimple geom? %i  entry length %i  atomcounter %i\n", simple_geom, entry_length, atomcounter);
    fclose(geomdat);
  } 
  else { 
    fprintf(outfile, "\nWhere's your geometry, bub?\n");
    return false;
  }
  
  if(simple_geom)
    natom = atomcounter / 4;
  else
    natom = entry_length;
  fflush(outfile);
  fprintf(outfile, "There are %i entries and %i atoms in %s geom is: %i (0 = array of arrays, 1 = simple array)\n", atomcounter, natom, entry_name, simple_geom);
  
  Molecule ref;
  
  for(index = 0; index < natom; index++) {   
    if(simple_geom)
      errcod = ip_string(entry_name,&atom_label,1,index * 4);
    else
      errcod = ip_string(entry_name,&atom_label,2,index,0);
    
    atom_num(atom_label, &Z);
    if (errcod != IPE_OK) {
      printf("Problem with the atom label geometry entry.");
	return false;
    }
    
    if(simple_geom) {
      errcod = ip_data(entry_name,"%lf", &tmp1,1,4*index+1);
      errcod = ip_data(entry_name,"%lf", &tmp2,1,4*index+2);
      errcod = ip_data(entry_name,"%lf", &tmp3,1,4*index+3);
    }
    else {
      errcod = ip_data(entry_name,"%lf", &tmp1,2,index,1);
      errcod = ip_data(entry_name,"%lf", &tmp2,2,index,2);
      errcod = ip_data(entry_name,"%lf", &tmp3,2,index,3);
    }
    if (errcod != IPE_OK) {
      printf("Problem with the geometry entry.");
      return false; 
    }
    //fprintf(outfile, "%i has X = %lf, Y = %lf, Z = %lf\n", (int)Z, tmp1,tmp2,tmp3);
    
    ref.addAtom((int)Z, an2masses[(int)Z],
		tmp1 * _bohr2angstroms,
		tmp2 * _bohr2angstroms,
		tmp3 * _bohr2angstroms);
  }
  printGeometries();
  fflush(outfile);

  return true;
}

//Directly lifted from input
void Displacements::atom_num(char *A, double *C)
{
  if(!strcmp(A,"GHOST") || !strcmp(A,"G")){
    *C = 0.00;
  }
  else if(!strcmp(A,"H") || !strcmp(A,"HYDROGEN")){
    *C = 1.00;
  }
  else if(!strcmp(A,"HE") || !strcmp(A,"HELIUM")){
    *C = 2.00;
  }
  else if(!strcmp(A,"LI")|| !strcmp(A,"LITHIUM")){
    *C = 3.00;
  }
  else if(!strcmp(A,"BE")|| !strcmp(A,"BERYLLIUM")){
    *C = 4.00;
  }
  else if(!strcmp(A,"B")|| !strcmp(A,"BORON")){
    *C = 5.00;
  }
  else if(!strcmp(A,"C") || !strcmp(A,"CARBON")){
    *C = 6.00;
  }
  else if(!strcmp(A,"N") || !strcmp(A,"NITROGEN")){
    *C = 7.00;
  }
  else if(!strcmp(A,"O") || !strcmp(A,"OXYGEN")){
    *C = 8.00;
  }
  else if(!strcmp(A,"F") || !strcmp(A,"FLUORINE")){
    *C = 9.00;
  }
  else if(!strcmp(A,"NE") || !strcmp(A,"NEON")){
    *C = 10.00;
  }
  else if(!strcmp(A,"NA") || !strcmp(A,"SODIUM")){
    *C = 11.00;
  }
  else if(!strcmp(A,"MG") || !strcmp(A,"MAGNESIUM")){
    *C = 12.00;
  }
  else if(!strcmp(A,"AL") || !strcmp(A,"ALUMINUM")){
    *C = 13.00;
  }
  else if(!strcmp(A,"SI") || !strcmp(A,"SILICON")){
    *C = 14.00;
  }
  else if(!strcmp(A,"P") || !strcmp(A,"PHOSPHORUS")){
  *C = 15.00;
  }
  else if(!strcmp(A,"S") || !strcmp(A,"SULPHUR") || !strcmp(A,"SULFUR")){
    *C = 16.00;
  }
  else if(!strcmp(A,"CL") || !strcmp(A,"CHLORINE")){
    *C = 17.00;
  }
  else if(!strcmp(A,"AR") || !strcmp(A,"ARGON")){
    *C = 18.00;
  }
  else if(!strcmp(A,"K") || !strcmp(A,"POTASSIUM")){
    *C = 19.00;
  }
  else if(!strcmp(A,"CA") || !strcmp(A,"CALCIUM")){
    *C = 20.00;
  }
  else if(!strcmp(A,"SC") || !strcmp(A,"SCANDIUM")){
    *C = 21.00;
  }
  else if(!strcmp(A,"TI") || !strcmp(A,"TITANIUM")){
    *C = 22.00;
  }
  else if(!strcmp(A,"V") || !strcmp(A,"VANADIUM")){
    *C = 23.00;
  }
  else if(!strcmp(A,"CR") || !strcmp(A,"CHROMIUM")){
    *C = 24.00;
  }
  else if(!strcmp(A,"MN") || !strcmp(A,"MANGANESE")){
    *C = 25.00;
  }
  else if(!strcmp(A,"FE") || !strcmp(A,"IRON")){
    *C = 26.00;
  }
  else if(!strcmp(A,"CO") || !strcmp(A,"COBALT")){
    *C = 27.00;
  }
  else if(!strcmp(A,"NI") || !strcmp(A,"NICKEL")){
    *C = 28.00;
  }
  else if(!strcmp(A,"CU") || !strcmp(A,"COPPER")){
    *C = 29.00;
  }
  else if(!strcmp(A,"ZN") || !strcmp(A,"ZINC")){
    *C = 30.00;
  }
  else if(!strcmp(A,"GA") || !strcmp(A,"GALLIUM")){
    *C = 31.00;
  }
  else if(!strcmp(A,"GE") || !strcmp(A,"GERMANIUM")){
    *C = 32.00;
  }
  else if(!strcmp(A,"AS") || !strcmp(A,"ARSENIC")){
    *C = 33.00;
  }
  else if(!strcmp(A,"SE") || !strcmp(A,"SELENIUM")){
    *C = 34.00;
  }
  else if(!strcmp(A,"BR") || !strcmp(A,"BROMINE")){
    *C = 35.00;
  }
  else if(!strcmp(A,"KR") || !strcmp(A,"KRYPTON")){
    *C = 36.00;
  }
  else if(!strcmp(A,"RB") || !strcmp(A,"RUBIDIUM")){
    *C = 37.00;
  }
  else if(!strcmp(A,"SR") || !strcmp(A,"STRONTIUM")){
    *C = 38.00;
  }
  else if(!strcmp(A,"Y") || !strcmp(A,"YTTRIUM")){
    *C = 39.00;
  }
  else if(!strcmp(A,"ZR") || !strcmp(A,"ZIRCONIUM")){
    *C = 40.00;
  }
  else if(!strcmp(A,"NB") || !strcmp(A,"NIOBIUM")){
    *C = 41.00;
  }
  else if(!strcmp(A,"MO") || !strcmp(A,"MOLYBDENUM")){
    *C = 42.00;
  }
  else if(!strcmp(A,"TC") || !strcmp(A,"TECHNETIUM")){
    *C = 43.00;
  }
  else if(!strcmp(A,"RU") || !strcmp(A,"RUTHENIUM")){
    *C = 44.00;
  }
  else if(!strcmp(A,"RH") || !strcmp(A,"RHODIUM")){
    *C = 45.00;
  }
  else if(!strcmp(A,"PD") || !strcmp(A,"PALLADIUM")){
    *C = 46.00;
  }
  else if(!strcmp(A,"AG") || !strcmp(A,"SILVER")){
    *C = 47.00;
  }
  else if(!strcmp(A,"CD") || !strcmp(A,"CADMIUM")){
    *C = 48.00;
  }
  else if(!strcmp(A,"IN") || !strcmp(A,"INDIUM")){
    *C = 49.00;
  }
  else if(!strcmp(A,"SN") || !strcmp(A,"TIN")){
    *C = 50.00;
  }
  else if(!strcmp(A,"SB") || !strcmp(A,"ANTIMONY")){
    *C = 51.00;
  }
  else if(!strcmp(A,"TE") || !strcmp(A,"TELLURIUM")){
    *C = 52.00;
  }
  else if(!strcmp(A,"I") || !strcmp(A,"IODINE")){
    *C = 53.00;
  }
  else if(!strcmp(A,"XE") || !strcmp(A,"XENON")){
    *C = 54.00;
  }
  else if(!strcmp(A,"CS") || !strcmp(A,"CESIUM")){
    *C = 55.00;
  }
  else if(!strcmp(A,"BA") || !strcmp(A,"BARIUM")){
    *C = 56.00;
  }
  else if(!strcmp(A,"LA") || !strcmp(A,"LANTHANUM")){
    *C = 57.00;
  }
  else if(!strcmp(A,"CE") || !strcmp(A,"CERIUM")){
    *C = 58.00;
  }
  else if(!strcmp(A,"PR") || !strcmp(A,"PRASEODYMIUM")){
    *C = 59.00;
  }
  else if(!strcmp(A,"ND") || !strcmp(A,"NEODYMIUM")){
    *C = 60.00;
  }
  else if(!strcmp(A,"PM") || !strcmp(A,"PROMETHIUM")){
    *C = 61.00;
  }
  else if(!strcmp(A,"SM") || !strcmp(A,"SAMARIUM")){
    *C = 62.00;
  }
  else if(!strcmp(A,"EU") || !strcmp(A,"EUROPIUM")){
    *C = 63.00;
  }
  else if(!strcmp(A,"GD") || !strcmp(A,"GADOLINIUM")){
  *C = 64.00;
  }
  else if(!strcmp(A,"TB") || !strcmp(A,"TERBIUM")){
    *C = 65.00;
  }
  else if(!strcmp(A,"DY") || !strcmp(A,"DYSPROSIUM")){
    *C = 66.00;
  }
  else if(!strcmp(A,"HO") || !strcmp(A,"HOLMIUM")){
    *C = 67.00;
  }
  else if(!strcmp(A,"ER") || !strcmp(A,"ERBIUM")){
    *C = 68.00;
  }
  else if(!strcmp(A,"TM") || !strcmp(A,"THULIUM")){
    *C = 69.00;
  }
  else if(!strcmp(A,"YB") || !strcmp(A,"YTTERBIUM")){
    *C = 70.00;
  }
  else if(!strcmp(A,"LU") || !strcmp(A,"LUTETIUM")){
  *C = 71.00;
  }
  else if(!strcmp(A,"HF") || !strcmp(A,"HAFNIUM")){
    *C = 72.00;
  }
  else if(!strcmp(A,"TA") || !strcmp(A,"TANTALUM")){
    *C = 73.00;
  }
  else if(!strcmp(A,"W") || !strcmp(A,"TUNGSTEN")){
    *C = 74.00;
  }
  else if(!strcmp(A,"RE") || !strcmp(A,"RHENIUM")){
    *C = 75.00;
  }
  else if(!strcmp(A,"OS") || !strcmp(A,"OSMIUM")){
    *C = 76.00;
  }
  else if(!strcmp(A,"IR") || !strcmp(A,"IRIDIUM")){
    *C = 77.00;
  }
  else if(!strcmp(A,"PT") || !strcmp(A,"PLATINUM")){
    *C = 78.00;
  }
  else if(!strcmp(A,"AU") || !strcmp(A,"GOLD")){
    *C = 79.00;
  }
  else if(!strcmp(A,"HG") || !strcmp(A,"MERCURY")){
    *C = 80.00;
  }
  else if(!strcmp(A,"TL") || !strcmp(A,"THALLIUM")){
    *C = 81.00;
  }
  else if(!strcmp(A,"PB") || !strcmp(A,"LEAD")){
    *C = 82.00;
 }
  else if(!strcmp(A,"BI") || !strcmp(A,"BISMUTH")){
    *C = 83.00;
  }
  else if(!strcmp(A,"PO") || !strcmp(A,"POLONIUM")){
    *C = 84.00;
  }
  else if(!strcmp(A,"AT") || !strcmp(A,"ASTATINE")){
    *C = 85.00;
 }
  else if(!strcmp(A,"RN") || !strcmp(A,"RADON")){
    *C = 86.00;
  }
  else if(!strcmp(A,"FR") || !strcmp(A,"FRANCIUM")){
    *C = 87.00;
  }
  else if(!strcmp(A,"RA") || !strcmp(A,"RADIUM")){
    *C = 88.00;
  }
  else if(!strcmp(A,"AC") || !strcmp(A,"ACTINIUM")){
    *C = 89.00;
  }
  else if(!strcmp(A,"TH") || !strcmp(A,"THORIUM")){
    *C = 90.00;
  }
  else if(!strcmp(A,"PA") || !strcmp(A,"PROTACTINIUM")){
    *C = 91.00;
  }
  else if(!strcmp(A,"U") || !strcmp(A,"URANIUM")){
    *C = 92.00;
  }
  else if(!strcmp(A,"NP") || !strcmp(A,"NEPTUNIUM")){
    *C = 93.00;
  }
  else if(!strcmp(A,"PU") || !strcmp(A,"PLUTONIUM")){
    *C = 94.00;
  }
  else if(!strcmp(A,"AM") || !strcmp(A,"AMERICIUM")){
  *C = 95.00;
  }
  else if(!strcmp(A,"CM") || !strcmp(A,"CURIUM")){
    *C = 96.00;
  }
  else if(!strcmp(A,"BK") || !strcmp(A,"BERKELIUM")){
    *C = 97.00;
  }
  else if(!strcmp(A,"CF") || !strcmp(A,"CALELSE IFORNIUM")){
    *C = 98.00;
  }
  else if(!strcmp(A,"ES") || !strcmp(A,"EINSTEINIUM")){
    *C = 99.00;
  }
  else if(!strcmp(A,"FM") || !strcmp(A,"FERMIUM")){
    *C = 100.00;
  }
  else if(!strcmp(A,"MD") || !strcmp(A,"MENDELEVIUM")){
    *C = 101.00;
  }
  else if(!strcmp(A,"NO") || !strcmp(A,"NOBELIUM")){
    *C = 102.00;
  }
  else if(!strcmp(A,"LR") || !strcmp(A,"LAWRENCIUM")){
    *C = 103.00;
  }
  else if(!strcmp(A,"UNQ")){
    *C = 104.00;
  }
  else if(!strcmp(A,"UNP")){
    *C = 105.00;
  }
  else if(!strcmp(A,"UNH")){
    *C = 106.00;
  }
  else if(!strcmp(A,"UNS")){
    *C = 107.00;
  }
  else {
    fprintf(outfile,"Unrecognized atom symbol/name!");
    exit(PSI_RETURN_FAILURE);
  }
}
