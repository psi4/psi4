/*###########################################################################*/
/*! 
  \file
  \ingroup EXTREMA
  \brief Member functions for z-matrix derived class
  
  Contains the zmat constructor, method driver functions, interfaces
  to lower level classes, and several small functions needed to
  implement z-matrix coordinates */
/*						Joseph P. Kenny 11/29/01
  ###########################################################################*/

#define EXTERN
#include"extrema.h"

using namespace psi::extrema;

/*---------------------------------------------------------------------------*/
/*! \fn zmat::zmat() 
  \brief zmat constructor

  Constructor for the top-level z-matrix derived class.  Performs setup
  of z-matrix coordinates which includes determining optimized coordinates
  and positive/negative torsion pairs.  The internals dummy constructor is
  called initially and the actual constructor named constructor_internals
  is called after the number of optimized coordinates is determined. */
/*---------------------------------------------------------------------------*/

zmat :: zmat() : internals()
   {

  int i, j, pos, dummy;        

  zmat::parse_input();

  /*read z_mat and cartesians from chkpt*/
  num_entries = chkpt_rd_nallatom();
  z_geom = chkpt_rd_zmat();
  char **temp_felement;
  felement = (char **) chkpt_rd_felement();

  dummy=0;
  for(i=0;i<num_entries;++i) 
      if(!strncmp(felement[i],"X\0",2) )
	  dummy=1;
  if(dummy) {

      free(carts);
      double** cart_temp;
      carts = init_array(3*num_entries);
      cart_temp = chkpt_rd_fgeom();

      for(i=0;i<(3*num_entries);++i) 
	  carts[i] = cart_temp[0][i];
      free_matrix(cart_temp,1);

      double* ctemp;
      ctemp = init_array(3*num_atoms);
      for(i=0;i<3*num_atoms;++i)
	      ctemp[i] = c_grads[i];
      free(c_grads);
      c_grads = init_array(3*num_entries);
      pos=0;
      for(i=0;i<num_entries;++i) {
	  if(strncmp(felement[i],"X\0",2) ) {
	      c_grads[3*i] = ctemp[3*pos];
	      c_grads[3*i+1] = ctemp[3*pos+1];
	      c_grads[3*i+2] = ctemp[3*pos+2];
	      ++pos;
	  }
      }	 
      free(ctemp);
  }

  switch(num_entries) {
    case 2: fnum_coords = 1; break;
    case 3: fnum_coords = 3; break;
    default: fnum_coords =  (num_entries*3-6); break;
  } 

  /*write z_mat to the array of simple_internal*/
  simples = (simple*) malloc(fnum_coords*sizeof(simple));  
  for(i=1;i<num_entries;++i) {
      if( i==1 ) {
	  simples[0].set_simple(0,z_geom[1].bond_val,
				2,z_geom[1].bond_atom,-1,-1,
				z_geom[1].bond_opt);
      }
      else if( i==2 ) {
	  simples[1].set_simple(0,z_geom[2].bond_val,
				3,z_geom[2].bond_atom,-1,-1,
				z_geom[2].bond_opt);
	  simples[2].set_simple(1,z_geom[2].angle_val *_pi/180.0,
				3,z_geom[2].bond_atom,z_geom[2].angle_atom,-1,
                                z_geom[2].angle_opt);
	  j=3;
	}
      else if( i>2 ) {
	  simples[j].set_simple(0,z_geom[i].bond_val,i+1,z_geom[i].bond_atom,
				-1,-1,z_geom[i].bond_opt);
	  simples[j+1].set_simple(1,z_geom[i].angle_val*_pi/180.0,
				  i+1,z_geom[i].bond_atom,
				  z_geom[i].angle_atom,-1,z_geom[i].angle_opt);
	  simples[j+2].set_simple(2,z_geom[i].tors_val*_pi/180.0,
				  i+1,z_geom[i].bond_atom,
				  z_geom[i].angle_atom,z_geom[i].tors_atom,
				  z_geom[i].tors_opt);
	  j+=3;
	}
  }

  int p=0;
  for(i=1;i<num_entries;++i) {
      if(i==1) {
	  simples[p].set_label(z_geom[1].bond_label);
	  ++p;
      }
      else if( i==2) {
          simples[p].set_label(z_geom[2].bond_label);
          simples[p+1].set_label(z_geom[2].angle_label );
	  p+=2;
      }
      else if (i>2) {
          simples[p].set_label(z_geom[i].bond_label);
          simples[p+1].set_label(z_geom[i].angle_label);
          simples[p+2].set_label(z_geom[i].tors_label);
	  p+=3;
      }
  }
  
  /*find first instance of each unique coordinate*/
  int there;
  first_unique = (int *) malloc(fnum_coords*sizeof(int));
  for(i=0;i<fnum_coords;++i) {
      there=0;
      for(j=0;j<i;++j) 
	  if( !strcmp( simples[i].get_label(), simples[j].get_label() )) {
	      first_unique[i]=0;
	      ++there;
	  }
      if( (there==0) || !strcmp(simples[i].get_label(),"") )
	  first_unique[i]=1;      
  }

  num_coords = 0;
  for(i=0;i<fnum_coords;++i) {
      if(simples[i].get_opt() && first_unique[i]) {
	  ++num_coords;
      }
  }
  if(!num_coords)
      punt("No coordinates to optimize");

  /*allocate memory now that optimized coordinate number is known*/
  internals::mem_alloc();

  for(i=0;i<fnum_coords;++i)
      fcoords[i] = fcoords_old[i] = simples[i].get_val();

  p=0;
  for(i=0;i<fnum_coords;++i) 
      if(simples[i].get_opt() && first_unique[i]) {  
	  coords[p] = simples[i].get_val();
	  ++p;
      }
	  
  /*find positive/negative torsion pairs*/
  int *is_torsion;
  is_torsion = init_int_array(num_coords);
  p=0;
  for(i=0;i<fnum_coords;++i) {
      if(simples[i].get_opt() && first_unique[i]) {
	  if(simples[i].get_type()==TORS_TYPE) 
	      is_torsion[p] = 1;
	  else is_torsion[p] = 0;
	  ++p;
      }
  }

  pos_neg_pairs = (int *) malloc(num_coords*sizeof(int));
  for(i=0;i<num_coords;++i)
      pos_neg_pairs[i] = 0;
  p=1;
  for(i=0;i<num_coords;++i) {
      if((pos_neg_pairs[i]==0) && is_torsion[i]) 
	  for(j=(i+1);j<num_coords;++j) 
	      if((pos_neg_pairs[j]==0) && is_torsion[j]) 
		  if( fabs(coords[i]+coords[j]) < POS_NEG_TORS ) {
		      pos_neg_pairs[i] = pos_neg_pairs[j] = p;
		      ++p;
		  }
  }
  
  free(is_torsion);
	
  /* test angles for extreme values, abort if hopeless */
  if(angle_abort)
      for(i=0;i<fnum_coords;++i) {
	  if( simples[i].get_type() == ANGLE_TYPE ) {
	      if( fabs(simples[i].get_val()) > NEAR_180*_pi/180.0 ) 
		  punt("Simple valence angle near 180 degrees");
	  }
	  else if (simples[i].get_type() == TORS_TYPE) 
	      if( (fabs(simples[i].get_val()) > NEAR_180*_pi/180.0 ) &&
		  (fabs(simples[i].get_val()) < NOT_180*_pi/180.0) )
		  punt("Simple torsion near 180 degrees");
      }
  
  return;
}



/*---------------------------------------------------------------------------*/
/*! \fn zmat::optimize()
  \brief z-matrix optimization driver */
/*---------------------------------------------------------------------------*/

void zmat :: optimize() {

    read_opt();
    print_carts(_bohr2angstroms);
    print_c_grads();
    grad_test();
    print_internals();
    compute_B();
    if(print_lvl >= RIDICULOUS_PRINT)
      print_B();
    grad_trans();
    switch(iteration) {
        case 1: initial_Hi(); break;
	default: update_Hi(); break; 
    }
    H_test();
    newton_step();
    back_transform();
    print_carts(_bohr2angstroms);
    print_internals();
    write_opt();
    write_chkpt();
    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn zmat::parse_input()
  \brief Performs input parsing for z-matrix coordinates. */
/*---------------------------------------------------------------------------*/
  
void zmat::parse_input() {
    
    bond_lim = BOND_LIM/_bohr2angstroms;
    errcod = ip_data("BOND_LIMIT","%lf",&bond_lim,0);
    angle_lim = ANGLE_LIM*_pi/180.0;
    if(ip_exist("ANGLE_LIMIT",0)) {
	errcod = ip_data("ANGLE_LIMIT","%lf",&angle_lim,0);
	angle_lim *= _pi/180.0;
    }
    fprintf(outfile,"\n  BOND_LIMIT:    %lf angstroms", 
	    bond_lim*_bohr2angstroms);
    fprintf(outfile,"\n  ANGLE_LIMIT:   %lf degrees", angle_lim*180.0/_pi);

    return;
}
	

/*---------------------------------------------------------------------------*/
/*! \fn zmat::initial_Hi()
  \brief Froms initial guess for inverse hessian. */
/*---------------------------------------------------------------------------*/

void zmat :: initial_Hi() {
    
    int i,j, pos=0;
    
    for(i=0;i<fnum_coords;++i) {
	if(first_unique[i] && simples[i].get_opt()) {
	    switch( simples[i].get_type() ) {
	    case 0: Hi[pos][pos] = 1.0; break;
	    case 1: Hi[pos][pos] = 4.0; break;
	    case 2: Hi[pos][pos] = 10.0; break;
	    }
	    ++pos;
	}
    }

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn zmat::print_internals()
  \brief Prints z-matrix. */
/*---------------------------------------------------------------------------*/

void zmat :: print_internals() {

  int i;

  /*print z-mat to output*/
  fprintf(outfile,"\n  Z-matrix (angstroms and degrees):\n\n");
  for(i=-1;i<(num_entries-1);++i) {
    switch(i) {
      case -1:
	  fprintf(outfile,"  %12s\n",felement[0]);
	  break;
      case 0: 
	  fprintf(outfile,"  %12s %4d %15.10lf\n",
		  felement[i+1], simples[i].get_bond(), 
		  simples[i].get_val()*_bohr2angstroms);
	  break;
      case 1: 
	  fprintf(outfile,"  %12s %4d %15.10lf %4d %15.10lf\n",
		  felement[i+1], simples[i].get_bond(), 
		  simples[i].get_val()*_bohr2angstroms,
		  simples[i+1].get_angle(), simples[i+1].get_val()*180.0/_pi); 
	  break;
      default:
	  fprintf(outfile,"  %12s %4d %15.10lf %4d %15.10lf %4d %15.10lf\n",
		  felement[i+1], simples[(i-1)*3].get_bond(), 
		  simples[(i-1)*3].get_val()*_bohr2angstroms,
		  simples[(i-1)*3+1].get_angle(), 
		  simples[(i-1)*3+1].get_val()*180.0/_pi,
		  simples[(i-1)*3+2].get_tors(), 
		  simples[(i-1)*3+2].get_val()*180.0/_pi );
	  break;
      }
  }
}



/*--------------------------------------------------------------------------*/
/*! \fn zmat::write_chkpt()
  \brief Writes geometries to chkpt. */
/*--------------------------------------------------------------------------*/

void zmat :: write_chkpt() {
    
  double **cart_matrix, **fcart_matrix;
  int i,j;

  cart_matrix = block_matrix(num_atoms,3);
  fcart_matrix = block_matrix(num_entries,3);

  int pos = 0;
  for(i=0;i<num_entries;++i) 
      for(j=0;j<3;++j) {
	  fcart_matrix[i][j] = carts[pos];
	  ++pos;
      }

  chkpt_wt_fgeom(fcart_matrix);
  
  int row=0;
  pos = 0;
  for(i=0;i<num_entries;++i) { 
      if(strncmp(felement[i],"X\0",2)) {
	  for(j=0;j<3;++j) {
	      cart_matrix[row][j] = carts[pos];
	      ++pos;
	  }
	  ++row;
      }
      else
	  pos += 3;
  }

  chkpt_wt_geom(cart_matrix);

  pos = -1;
  for(i=1;i<num_entries;++i) {
      
      if(i==1) 
	  z_geom[i].bond_val = simples[++pos].get_val(); 
      
      if(i==2) {
	  z_geom[i].bond_val = simples[++pos].get_val();
	  z_geom[i].angle_val = simples[++pos].get_val() * 180.0 / _pi;
      }
      
      if(i>2) {
	  z_geom[i].bond_val = simples[++pos].get_val();
	  z_geom[i].angle_val = simples[++pos].get_val() * 180.0 / _pi;
	  z_geom[i].tors_val = simples[++pos].get_val() * 180.0 / _pi;
      }
  }
  
  chkpt_wt_zmat(z_geom);

  return;
}



/*--------------------------------------------------------------------------*/
/*! \fn zmat::newton_step()
  \brief Computes newton-raphson z-matrix optimization step.
  
  Interface for <b>math_tools::newton_step()</b>. */
/*-------------------------------------------------------------------------*/

void zmat :: newton_step() {
    
    int i, j, k, p;
    double *s, num, con;

    fprintf(outfile,"\n\n  --------------------------------------");
    fprintf(outfile,"--------------------------------------\n");
    fprintf(outfile,"  Computing newton-raphson optimization step\n");
    fprintf(outfile,"  --------------------------------------");
    fprintf(outfile,"--------------------------------------\n");
    s = math_tools::newton_step(num_coords,Hi,grads);

    /* ensure pos/neg torsion angle pairs match */
    for(i=0;i<num_coords;++i) 
	if(pos_neg_pairs[i]) 
	    for(j=(i+1);j<num_coords;++j) 
		if(pos_neg_pairs[i]==pos_neg_pairs[j]) {
		    
		    num = (fabs(s[i]) + fabs(s[j]))/2;
		    if(s[i] < 0)
			s[i] = -1.0*num;
		    else s[i] = num;
		    if(s[j] < 0)
			s[j] = -1.0*num;
		    else s[j] = num;
		    
		    if( (s[i]+s[j]) > POS_NEG_TORS ) {
			fprintf(outfile,"\n  WARNING: positive/negative pair");
		        fprintf(outfile," %d/%d:\n",i,j);
			fprintf(outfile,"  displacements differ by more than");
			fprintf(outfile,"%lf radians\n", POS_NEG_TORS);
		    }
		}
    
    /* print before limit enforcement */
    if(print_lvl>1) {
	p=0;
	fprintf(outfile,"\n  Displacements before limit enforcement");
	fprintf(outfile," (angstroms and degrees):\n");
	for(i=0;i<fnum_coords;++i)
	    if(simples[i].get_opt() && first_unique[i]) {
		if(simples[i].get_type()==0)
		    con = _bohr2angstroms; 
		else con = 180.0/_pi;  
		fprintf(outfile,"  %i %8s: %lf\n",
			i+1,simples[i].get_label(),s[p]*con);
		++p;
	    }
    }
    
    /* enforce bond and angle limits */
    p=0;
    for(i=0;i<fnum_coords;++i) {
        if((simples[i].get_type()==0) && first_unique[i] && simples[i].get_opt() ) {
            if(s[p] > bond_lim) 
                s[p] = bond_lim; 
            else if (s[p] < (-1.0*bond_lim) ) 
                s[p] = (-1.0*bond_lim); 
	    ++p;
        }
        else if(first_unique[i] && simples[i].get_opt() ) {
            if(s[p] > angle_lim) 
                s[p] = angle_lim; 
            else if (s[p] < (-1.0*angle_lim) ) 
                s[p] = (-1.0*angle_lim);
	    ++p;
        }
    }

    /*save old coordinates to write to opt.dat and for back transform*/
    fcoord_old = init_array(fnum_coords);
    for(i=0;i<num_coords;++i) 
	coord_write[i] = coords[i];
    for(i=0;i<fnum_coords;++i)
	fcoord_old[i] = simples[i].get_val();


    /* take the step */
    for(i=0;i<num_coords;++i) 
	coords[i] += s[i];

    /*update simples*/
    p=0;
    for(i=0;i<fnum_coords;++i) {
	if(first_unique[i] && simples[i].get_opt()) {
	    simples[i].set_val(coords[p]);
	    ++p;
	    for(j=(i+1);j<fnum_coords;++j) 
		if(!strcmp(simples[i].get_label(),simples[j].get_label()) &&
		    strcmp(simples[i].get_label(),""))
		    simples[j].set_val(simples[i].get_val());
	}
    }

    /* test for extreme angle cases */
    if(angle_abort)
        for(i=0;i<fnum_coords;++i) {
            if( simples[i].get_type() == ANGLE_TYPE ) {
                if( fabs(simples[i].get_val()) > NEAR_180*_pi/180.0 )
                    punt("Simple valence angle near 180 degrees");
            }
            else if (simples[i].get_type() == TORS_TYPE)
                if( (fabs(simples[i].get_val()) > NEAR_180*_pi/180.0 ) &&
                    (fabs(simples[i].get_val()) < NOT_180*_pi/180.0) )
                    punt("Simple torsion near 180 degrees");
        }
    
    int entry;
    fprintf(outfile,"\n  Optimization Step (angstroms and degrees):\n");
    fprintf(outfile,"\n  label    initial value   gradient (a.u.)");      
    fprintf(outfile," displacement    new value");
    fprintf(outfile,"\n  -------- --------------- ---------------"); 
    fprintf(outfile," --------------- ---------------");
    entry=0;
    p=0;
    for(i=0;i<fnum_coords;++i) {
	con = _bohr2angstroms;
	if( i==2 )
	    con = 180.0/_pi;
	if( i>2 ) {
	    if(entry==0)
		++entry;
	    else if(entry==1) {
		con = 180.0/_pi;
		++entry;
	    }
	    else if(entry==2) {
		con = 180.0/_pi;
		entry = 0;
	    }
	}
	if(simples[i].get_opt() && first_unique[i]) {
	    fprintf(outfile,"\n  %-8s %15.10lf %15.10lf %15.10lf %15.10lf",
		    simples[i].get_label(),coord_write[p]*con, grads[p], 
		    s[p]*con, coords[p]*con);
	    ++p;
	}
    }

    free(s);
    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn zmat::print_carts(double)
  \brief Prints cartesians.

  \param conv conversion factor;
  either <b>1.0</b> for bohr or <b>_bohr2angstroms</b> for angstroms */
/*---------------------------------------------------------------------------*/

void zmat :: print_carts(double conv) {
 
    int i, j;
    double **temp;

    temp = init_matrix(num_entries,3);
    for(i=0;i<num_entries;++i) 
	for(j=0;j<3;++j) 
	    temp[i][j] = carts[3*i+j]*conv;
    if(conv==1.0)
	fprintf(outfile,"\n  Cartesian Coordinates (bohr):\n");
    else
	fprintf(outfile,"\n  Cartesian Coordinates (angstroms):\n");
     fprintf(outfile,"                       x              y         ");
    fprintf(outfile,"       z\n");
    fprintf(outfile,"                --------------- --------------- ");
    fprintf(outfile,"---------------\n");
    for(i=0;i<num_entries;++i)
	fprintf(outfile,"  %12s  %15.10lf %15.10lf %15.10lf\n",
		felement[i], temp[i][0], temp[i][1], temp[i][2]);
    free_matrix(temp,num_entries);

    return;
}



/*--------------------------------------------------------------------------*/
/*! \fn zmat::print_c_grads()
  \brief Prints cartesian gradients (Hartree/Bohr). */
/*--------------------------------------------------------------------------*/

void zmat :: print_c_grads() {
	
    int i, j;

    double **cgtemp;
    cgtemp = init_matrix(num_entries,3);
	
    for(i=0;i<num_entries;++i) 
	for(j=0;j<3;++j) 
	    cgtemp[i][j] = c_grads[3*i+j];
    fprintf(outfile,"\n  Cartesian Gradients (a.u):\n");
    fprintf(outfile,"                       x              y         ");
    fprintf(outfile,"       z\n");
    fprintf(outfile,"                --------------- --------------- ");
    fprintf(outfile,"---------------\n");
    for(i=0;i<num_entries;++i)
	fprintf(outfile,"  %12s  %15.10lf %15.10lf %15.10lf\n",
		felement[i], cgtemp[i][0], cgtemp[i][1], cgtemp[i][2]);
    free_matrix(cgtemp,num_entries);
    
    return;
}








