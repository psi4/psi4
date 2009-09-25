/*###########################################################################*/
/*! 
  \file
  \ingroup EXTREMA
  \brief Member functions for deloc (delocalized internals) derived class
  
  Contains the deloc constructor, method driver functions, interfaces
  to lower level classes, and several small functions needed to
  implement delocalized internal coordinates */
/*						Joseph P. Kenny 12/06/01
  ###########################################################################*/

#define EXTERN
#include"extrema.h"
#include"bondl.h"
#include"inline.h"

using namespace psi::extrema;

/*---------------------------------------------------------------------------*/
/*! \fn deloc::deloc() 

  Constructor for the top-level delocalized internals derived class. 
/*---------------------------------------------------------------------------*/

deloc :: deloc() : internals() {

  int i,j,pos;

  deloc::parse_input();  

  felement = (char **) chkpt_rd_felement();
  point_group = chkpt_rd_sym_label();

  FILE *opt_ptr;
  ffile_noexit(&opt_ptr,"opt.dat",2);
  int form_deloc=0;
  if( opt_ptr == NULL ) {
      form_deloc=1;
      if( do_opt )
	  do_deriv=0;
      iteration = 1;
  }
  
  else {
      ip_done();
      ip_set_uppercase(1);
      ip_initialize(opt_ptr,outfile);
      ip_cwk_add(":OPT_INFO");
      if(ip_exist("DO_DERIV",0)) {
	  do_deriv=1;
	  punt("Can't do derivatives yet");
      }
      else
	  do_deriv=0;
      iteration = 2; /* proper number set in read_opt */
      read_bonds();
  }
  
  deloc::init_simples(); 
  if(print_lvl > NORMAL_PRINT)
      deloc::print_simples();

  /* form B matrix */
  deloc::compute_B_simple();

  if(print_lvl >= RIDICULOUS_PRINT) {
      fprintf(outfile,"\n  B matrix in simples representation:\n");
      fnum_coords = num_simples;
      print_B();
  }

  if(form_deloc) {

      fprintf(outfile,"\n\n  Forming delocalized internal coordinates\n"); 
   
      iteration = 1;
      
      /* form and diagonalize G */
      fnum_coords = num_simples; /* for now */
      G = init_matrix(num_simples,num_simples);
      internals::compute_G();
            
      double **temp_mat, *evals, **evects;
      evals = init_array(num_simples);
      evects = init_matrix(num_simples,num_simples);
      temp_mat = init_matrix(num_simples,num_simples);
      sq_rsp(num_simples,num_simples, G, evals, 1, evects, 1.0e-14);
      mmult(G,0,evects,0,temp_mat,0,num_simples,num_simples,num_simples,0);
      for (j=0;j<num_simples;++j) {
	  error = 0;
	  for (i=0;i<num_simples;++i) {
	      if ( fabs(temp_mat[i][j] - evals[j]*evects[i][j]) > 1.0e-13)  
		  error=1;
	      if (error == 1) { 
		  fprintf(outfile,
			  "\n  WARNING error detected in evect %d\n",j); 
		  error = 0;
	      }
	  }
      }      
      free_matrix(temp_mat,num_simples);
      
      if(print_lvl>NORMAL_PRINT) {
	  pos=0;
	  fprintf(outfile,"\n  Eigenvalues of G:\n");
	  for(i=(num_simples-1);i>-1;--i) {
	      if(pos==4) {
		  fprintf(outfile,"\n");
		  pos=0;
	      }
	      fprintf(outfile,"   % 15.10lf",evals[i]);
	      ++pos;
	  }
      }

      num_nonzero = 0;
      for (i=0;i<num_simples;++i)
	  if( evals[i] > ev_tol ) ++num_nonzero;
      fprintf(outfile,
	      "\n  %d nonzero eigenvalues of G (non-redundant coordinates)\n",
	      num_nonzero);

      /* check for proper number of non-redundant coordinates */
      int rotor_type;
      rotor_type = chkpt_rd_rottype();
      switch (rotor_type) {
      case 3:
	  if(num_atoms==2)
	      degrees_of_freedom = 3*num_atoms -5;
	  else punt("Use z-matrix for linear molecules");
	  break;
      case 6:
	  degrees_of_freedom = 0;
	  break;
      default:
	  degrees_of_freedom = 3 * num_atoms - 6;
	  break;
      }
      
      if(num_nonzero == degrees_of_freedom) {
	  fprintf(outfile,
		  "  Non-zero eigenvalues of G equal degrees of freedom");
	  fprintf(outfile,"(%d)\n",num_nonzero);
      }
       else if(num_nonzero > degrees_of_freedom) {
	  fprintf(outfile,"  Only %d degrees of freedom",degrees_of_freedom);
	  fprintf(outfile,"\n  Symmetry projection and orthogonalization");
	  fprintf(outfile," will be performed\n");
      }
      if (num_nonzero < degrees_of_freedom) {
	  fprintf(outfile,"\n  Not enough non-redundant coordinates,");
	  fprintf(outfile,
		  "\n  try reducing 'ev_tol' and check atom connectivity\n");
	  punt("Not enough non-redundant coordinates");
      }
      
      /* define delocalized internal coordinates */
      deloc_define = init_matrix(num_nonzero, num_simples);
      pos=0;
      for(i=0;i<num_simples;++i)
	  if(evals[i] > ev_tol) {
	      for(j=0;j<num_simples;++j) 
		  deloc_define[pos][j] = evects[j][i];
	      ++pos;
	  }
      free_matrix(evects,num_simples);
      
      if(print_lvl>=RIDICULOUS_PRINT) {
      fprintf(outfile,"\n  Initial coordinates:\n");
      print_mat(deloc_define,num_nonzero,num_simples,outfile);
      }

      /* clean up coordinate symmetry */
      deloc::ir_project();

      if(num_coords!=degrees_of_freedom) 
	  punt("number of coordinates != degrees of freedom");
  
      for(i=0;i<num_coords;++i)
	  for(j=0;j<num_simples;++j)
	      deloc_define[i][j] *= fabs(deloc_define[i][j]);

      if(!do_deriv) {
	  
	  /* coordinates look good, allocate memory */

	  /* save old B_mat first */
	  free_matrix(B,num_simples);

	  fnum_coords = num_coords;

	  free_matrix(u,3*num_entries);

	  /* keep only symmetric coordinates for optimization */
	  int *irrep_temp;
	  irrep_temp = init_int_array(num_nonzero);
	  for(i=0;i<num_coords;++i) {
	      irrep_temp[i] = deloc_irrep[i];
	  }	      
	  free(deloc_irrep);
	  deloc_irrep = init_int_array(num_coords);

	  double **define_temp;
	  int temp_num_coords=0;
	  define_temp = init_matrix(num_coords,num_simples);
	  for(i=0;i<num_coords;++i)
	      for(j=0;j<num_simples;++j)
		  define_temp[i][j] = deloc_define[i][j];
	  free_matrix(deloc_define,num_coords);
	  deloc_define = init_matrix(num_coords,num_simples);
	  int p=0;
	  for(i=0;i<num_coords;++i) 
	      if(irrep_temp[i]==0) {
		  ++temp_num_coords;
		  for(j=0;j<num_simples;++j)
		      deloc_define[p][j] = define_temp[i][j];
		  deloc_irrep[p] = irrep_temp[i];
		  ++p;
	      }

	  num_coords = fnum_coords = temp_num_coords;
	
	  internals::mem_alloc();
	  
	  p = -1;
	  /* compute coordinate values */
	  double sum;
	  for(i=0;i<num_coords;++i) 
	      if(deloc_irrep[i]==0) {
		  sum=0;
		  for(j=0;j<num_simples;++j)
		      sum += deloc_define[i][j] * simples[j].get_val();
		  coords[++p] = sum;
	      }  	


	  
	  fprintf(outfile,"\n  Delocalized internal coordinate values:\n");
	  for(i=0;i<num_coords;++i) 
	      fprintf(outfile,"  %d: %lf\n",i+1,coords[i]);
	  
	  fflush(outfile);
	  
      }    
  }
  
  /* reform u */
  for(i=0;i<num_atoms;++i) {
      u[3*i][3*i] = 1/masses[i];
      u[3*i+1][3*i+1] = 1/masses[i];
      u[3*i+2][3*i+2] = 1/masses[i];
  }
  
  return;
}



/*---------------------------------------------------------------------------*/
/*! \fn deloc::optimize()
  \brief Delocalized internals optimization driver */
/*---------------------------------------------------------------------------*/

void deloc :: optimize() {

    read_opt();
    print_carts(_bohr2angstroms);
    print_c_grads();
    grad_test();
    compute_B();
    grad_trans();
    grads = fgrads;
    switch(iteration) {
      case 1: initial_Hi(); break;
      default: update_Hi(); break; 
    }    
    if(print_lvl > NORMAL_PRINT)
	print_H();
    deloc::newton_step();
    deloc::back_transform();
    print_carts(_bohr2angstroms);
    // print_internals();
    write_opt();
    write_chkpt();
    return;

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn deloc::back_transform()
  \brief Computes cartesians corresponding to the current delocalized coords.
 
  Interface for <b>internals::back_transform()</b>. */
/*---------------------------------------------------------------------------*/

void deloc::back_transform() {

    int i;
    double *coord_temp;
    coord_temp = init_array(num_coords);
    for(i=0;i<num_coords;++i)
	coord_temp[i] = coord_write[i];
    
    internals::back_transform(coords, coord_temp);

    free(coord_temp);

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn deloc::parse_input()
/*---------------------------------------------------------------------------*/
  
void deloc::parse_input() {

    ev_tol = DELOC_EV_TOL;
    int power;
    if(ip_exist("EV_TOL",0)) {
	errcod = ip_data("EV_TOL","%d",&power,0);
	ev_tol = 1.0 / pow(10.0,(double)power);
	fprintf(outfile,"\n  EV_TOL:        %d", power);
    }
    
    return;
}


	         
/*---------------------------------------------------------------------------*/
/*! \fn deloc::init_simples()
  \brief Initializes array of simple internal objects.
/*---------------------------------------------------------------------------*/

void deloc::init_simples() {
    
    int i,j,k,l, pos, *ioff, num;
    double internal_val, *distance;

    if(iteration==1) {
	/* compute interatomic distances */
	distance = init_array( ((num_atoms+1)*num_atoms)/2 );
	pos = -1;
	for(i=0;i<num_atoms;++i) 
	    for(j=0;j<=i;++j) 
		distance[++pos] = sqrt( (carts[3*i+0]-carts[3*j+0])
					*(carts[3*i+0]-carts[3*j+0])
					+(carts[3*i+1]-carts[3*j+1])
					*(carts[3*i+1]-carts[3*j+1])
					+(carts[3*i+2]-carts[3*j+2])
					*(carts[3*i+2]-carts[3*j+2]) );
	
	/* determine bonds */				    
	ioff = (int *) malloc (32641 * sizeof(int));
	ioff[0]=0;
	for(i=1;i<32641;++i)
	    ioff[i] = ioff[i-1] + i;
	
	bonds = init_int_matrix(num_atoms,num_atoms);
	int max, min, index;
	int warn=0;
	for(i=0;i<num_atoms;++i) 
	    for(j=0;j<=i;++j) {
		if( ((int) atomic_nums[i]) > ((int) atomic_nums[j]) ) { 
		    max = (int) atomic_nums[i];
		    min = (int) atomic_nums[j];
		}
		else {
		    max = (int) atomic_nums[j];
		    min = (int) atomic_nums[i];
		}
		index = ioff[max-1]+ (min-1);
		if(bondl[index] != 0.0) {
		    if( distance[ioff[i]+j] < (1.2*bondl[index]) )
			bonds[i][j] = bonds[j][i] = 1;
		    else if ((bondl[index] != 0.0) && (warn=0)) {
			fprintf(outfile,
				"\n\n  WARNING Bond lengths not known ");
			fprintf(outfile,"for all atoms.\n");
			fprintf(outfile,
				"  BONDS keyword can fix this.\n");
			++warn;
		    }
		}
	    }
	
	/* check input for user specified bonds or nobonds */
	int a,b;
	if (ip_exist("BONDS",0)) {
	    ip_count("BONDS",&num,0);
	    for(i=0;i<num;++i) {
		ip_data("BONDS","%d",&a,2,i,0);
		ip_data("BONDS","%d",&b,2,i,1);
		fprintf(outfile,"\n  User specified bond %d %d",a,b);
		--a;  --b;
		bonds[a][b] = 1;
		bonds[b][a] = 1;
	    }
	}
    
	if (ip_exist("NOBONDS",0)) {
	    ip_count("NOBONDS",&num,0);
	    for(i=0;i<num;++i) {
		ip_data("NOBONDS","%d",&a,2,i,0);
		ip_data("NOBONDS","%d",&b,2,i,1);
		a -= 1;  b -= 1;
		bonds[a][b] = 0;
		bonds[b][a] = 0;
	    }
	}

	free(ioff);
	free(distance);
    }

    /* count number of bonds */
    num=0;
    for(i=0;i<num_atoms;++i) 
	for(j=i+1;j<num_atoms;++j) 
	    if(bonds[i][j] == 1) 
		++num;
    count_array[0]=num;
    if( (num_atoms>1) && (num<1) )
	punt("Can't find any bonds");
    
    /* count number of bends*/
    num=0;
    for(i=0;i<num_atoms;++i) 
	for(j=0;j<num_atoms;++j) 
	    if(i!=j)
		for(k=i+1;k<num_atoms;++k) 
		    if(j!=k)
			if (bonds[i][j] && bonds[j][k]) 
			    ++num;
    count_array[1]=num;

    /* count number of torsions */
    num=0;
    for(i=0;i<num_atoms;++i) 
	for(j=0;j<num_atoms;++j) 
	    if(i!=j)
		for(k=0;k<num_atoms;++k) 
		    if(i != k && j != k) 
			for(l=i+1;l<num_atoms;++l) 
			    if( (l != j && l != k) && bonds[i][j] 
				&& bonds[j][k] && bonds[k][l]) 
				++num;
    count_array[2]=num;
    count_array[3]=0;

    if(print_lvl > NORMAL_PRINT) {
	fprintf(outfile,
		"\n\n  %d bonds, %d valence angles and %d torsion angles"
		,count_array[0],count_array[1],count_array[2]);
	fprintf(outfile," have been found");
    } 

    num_simples = (count_array[0] + count_array[1] 
		   + count_array[2] + count_array[3]);

    simples = (simple *) malloc(num_simples*sizeof(simple));

    /* set bonds */
    pos=-1;
    for(i=0;i<num_atoms;++i)
	for(j=i+1;j<num_atoms;++j)
	    if(bonds[i][j]) {
		internal_val = compute_bond(carts,i,j);
		simples[++pos].set_simple(0,internal_val,i,j,-1,-1,-1);
	    }

    /* set bends */
    for(i=0;i<num_atoms;++i) 
	for(j=0;j<num_atoms;++j) 
	    if(i!=j)
		for(k=i+1;k<num_atoms;++k) 
		    if(j!=k)
			if (bonds[i][j] && bonds[j][k]) {
			    internal_val = compute_angle(carts,i,j,k);
			    simples[++pos].set_simple(1,internal_val,
						      i,j,k,-1,-1);
			}

    /* set torsions */
    for(i=0;i<num_atoms;++i) 
	for(j=0;j<num_atoms;++j) 
	    if(i!=j)
		for(k=0;k<num_atoms;++k) 
		    if(i != k && j != k) 
			for(l=i+1;l<num_atoms;++l) 
			    if( (l != j && l != k) && bonds[i][j] 
				&& bonds[j][k] && bonds[k][l]) { 
				internal_val = compute_torsion(carts,i,j,k,l);
				simples[++pos].set_simple(2,internal_val,
							  i,j,k,l,-1);
			    }

    /* test for near 180.0 angles, abort if present */
    if(angle_abort)
	for(i=0;i<num_simples;++i) {
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
/*! \fn deloc::set_simples()
  \brief Sets values of simple internal objects.
/*---------------------------------------------------------------------------*/

void deloc::set_simples() {
    
    int i,j,k,l,pos;
    double internal_val;

    /* set values */
    for(i=0;i<num_simples;++i) {
	if(simples[i].get_type() == BOND_TYPE) {
	    internal_val = compute_bond(carts,simples[i].get_atom(),
					simples[i].get_bond());
	    simples[i].set_val(internal_val);
	    }

	else if (simples[i].get_type() == ANGLE_TYPE) {
	    internal_val = compute_angle(carts,simples[i].get_atom(),
					 simples[i].get_bond(),
					 simples[i].get_angle());
	    simples[i].set_val(internal_val);
	}

	else if( simples[i].get_type() == TORS_TYPE) { 
	    internal_val = compute_torsion(carts,simples[i].get_atom(),
					   simples[i].get_bond(),
					   simples[i].get_angle(),
					   simples[i].get_tors());
	    simples[i].set_val(internal_val);
	}

    }
    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn deloc::print_simples()
  \brief Prints simple internal coordinates.
/*---------------------------------------------------------------------------*/

void deloc::print_simples() {
    
    int i;
    
    fprintf(outfile,"\n  Generated Simple internal coordinates:");
    
    fprintf(outfile,"\n    Bonds:");
    for(i=0;i<num_simples;++i) 
	if(simples[i].get_type() == BOND_TYPE)
	    fprintf(outfile,"\n    %d  %d  %lf",
		    simples[i].get_atom()+1,simples[i].get_bond()+1,
		    simples[i].get_val());

    fprintf(outfile,"\n    Angles:");
    for(i=0;i<num_simples;++i) 
	if(simples[i].get_type() == ANGLE_TYPE)
	    fprintf(outfile,"\n    %d  %d  %d  %lf",
		    simples[i].get_atom()+1,simples[i].get_bond()+1,
		    simples[i].get_angle()+1,
		    simples[i].get_val()*180.0/_pi);
    
    fprintf(outfile,"\n    Torsions:");
    for(i=0;i<num_simples;++i) 
	if(simples[i].get_type() == TORS_TYPE)
	    fprintf(outfile,"\n    %d  %d  %d  %d  %lf",
		    simples[i].get_atom()+1,simples[i].get_bond()+1,
		    simples[i].get_angle()+1,simples[i].get_tors()+1,
		    simples[i].get_val()*180.0/_pi);

    fprintf(outfile,"\n");
    return;
}


/*--------------------------------------------------------------------------*/
/*! \fn deloc::newton_step()
  \brief Computes newton-raphson optimization step for delocalized internals.
  
  Interface for <b>math_tools::newton_step()</b>. */
/*-------------------------------------------------------------------------*/

void deloc :: newton_step() {
    
    int i, j, k, p;
    double *s, num, con;

    for(i=0;i<num_coords;++i)
	grads[i] = fgrads[i];

    fprintf(outfile,"\n\n  --------------------------------------");
    fprintf(outfile,"--------------------------------------\n");
    fprintf(outfile,"  Computing newton-raphson optimization step\n");
    fprintf(outfile,"  --------------------------------------");
    fprintf(outfile,"--------------------------------------\n");
    
    s = math_tools::newton_step(num_coords,Hi,grads);    

    /* enforce bond and angle limits */

    /*save old coordinates to write to opt.dat and for back transform*/
    for(i=0;i<num_coords;++i) 
	coord_write[i] = coords[i];

    /* limit step size */
    for(i=0;i<num_coords;++i) {
	if(s[i]>0.1)
	    s[i] = 0.1;
	if(s[i]<-0.1)
	    s[i] = -0.1;
    }

    /* ugly fix for linear fragments 
    for(i=0;i<num_coords;++i)
	if(fabs(grads[i])<1.0e-6)
	    s[i] = 0.0; */

    /* take the step */
    for(i=0;i<num_coords;++i) 
	coords[i] += s[i];

    fprintf(outfile,"\n  Optimization Step:\n");
    fprintf(outfile,"\n  label    initial value   gradient       ");      
    fprintf(outfile," displacement    new value");
    fprintf(outfile,"\n  -------- --------------- ---------------"); 
    fprintf(outfile," --------------- ---------------");
    for(i=0;i<num_coords;++i) 
	fprintf(outfile,"\n  %-4d %15.10lf %15.10lf %15.10lf %15.10lf",
		i+1,coord_write[i], grads[i], 
		s[i], coords[i]);    

    fflush(outfile);

    free(s);
    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn void compute_B_simple()
  \brief Computes B matrix in simple internal coordinates.
/*---------------------------------------------------------------------------*/

void deloc::compute_B_simple() {
    
    int i;

    B = (double**) malloc(num_simples*sizeof(double*));

    for(i=0;i<num_simples;++i) {
	if(simples[i].get_type() == 0)
	    B[i] = B_row_bond(carts,simples[i].get_atom(),
			     simples[i].get_bond());
	else if(simples[i].get_type() == 1)
	    B[i] = B_row_angle(carts,simples[i].get_atom(),
			      simples[i].get_bond(),simples[i].get_angle());
	else if(simples[i].get_type() == 2)
	    B[i] = B_row_tors(carts,simples[i].get_atom(),
			     simples[i].get_bond(),simples[i].get_angle(),
			     simples[i].get_tors());
    }

    /* throw in the u matrix for fun */
    u = init_matrix(3*num_entries,3*num_entries);
    for(i=0;i<num_atoms;++i) {
	u[3*i][3*i] = 1/masses[i];
	u[3*i+1][3*i+1] = 1/masses[i];
	u[3*i+2][3*i+2] = 1/masses[i];
    }

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn void compute_B()
  \brief Computes B matrix in delocalized internal coordinate representation.*/
/*---------------------------------------------------------------------------*/

void deloc::compute_B() {

    int i,j;
    double **B_temp;

    free(B);
    compute_B_simple();

    /* put simple B in B_temp */
    B_temp = init_matrix(num_simples,3*num_atoms);
    for(i=0;i<num_simples;++i) 
	for(j=0;j<3*num_atoms;++j)
	    B_temp[i][j] = B[i][j];

    /* compute B in internal representation */
    free_matrix(B,num_simples);
    B = init_matrix(num_coords,3*num_atoms);
    mmult(deloc_define,0,B_temp,0,B,0,num_coords,num_simples,3*num_atoms,0);

    if(print_lvl >= RIDICULOUS_PRINT) {
	fprintf(outfile,"\n B matrix: \n");
	print_mat(B,num_coords,3*num_atoms,outfile);
    }

    free_matrix(B_temp,num_simples);

    fflush(outfile);

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn void compute_block_G()
  \brief Computes block diagonal G matrix.

  A block diagonal G matrix is formed by multiplying only rows of the 
  B matrix corresponding to like simple types*/
/*---------------------------------------------------------------------------*/

void deloc::compute_block_G() {

    int i, j, first_row[4];
    first_row[0] = 0;
    first_row[1] = count_array[0];
    first_row[2] = first_row[1]+count_array[1];
    first_row[3] = first_row[2]+count_array[2];

    double **ptr;
    ptr = (double **) malloc(num_simples*sizeof(double *));
    for (i=0;i<4;++i) {
      for (j=0;j<count_array[i];++j) {
        ptr[j] = G[first_row[i]+j] + first_row[i];
      }
      if (count_array[i] != 0)
        mmult( &(B[first_row[i]]),0,
	       &(B[first_row[i]]),1,
	       ptr,0,count_array[i],num_atoms*3,count_array[i],0);
    }
    free(ptr);

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn void initial_Hi()
  \brief Computes initial Hi.*/
/*---------------------------------------------------------------------------*/

void deloc::initial_Hi() {

    int i;
    for(i=0;i<num_coords;++i)
	Hi[i][i] = 1.0;

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn  deloc::cart_to_internal(double* deloc_array)
  \brief  Computes delocalized internals from cartesian coordinates. 
  \param deloc_array values of delocalized internals */
/*--------------------------------------------------------------------------*/

void deloc :: cart_to_internal(double** deloc_array ) {

    int i,j;
    double value;

    /* recompute simples */
    deloc::set_simples();

    /* recompute delocalized internals */
    for(i=0; i<num_coords; ++i) {
	value = 0.0;
	for(j=0; j<num_simples; ++j )
	    value += deloc_define[i][j] * simples[j].get_val();
	(*deloc_array)[i] = value;
    }

    if(print_lvl >= RIDICULOUS_PRINT) {
	fprintf(outfile,"\n  intermediate internals:\n");
	for(i=0; i<num_coords; ++i)
	    fprintf(outfile,"  %d: %lf\n",i+1,(*deloc_array)[i]);
    }

    return;
}								     



/*---------------------------------------------------------------------------*/
/*! \fn deloc::read_opt()
  \brief Reads from opt.dat()
  
  Calls coord_base::read_opt() then reads deloc definitions. */
/*---------------------------------------------------------------------------*/

void deloc :: read_opt() {

    if(iteration!=1) {
	
	FILE *opt_ptr;
	
	int i, j, error=0, dim2;
	
	ffile(&opt_ptr,"opt.dat",2);
	if( opt_ptr != NULL ) {      
	    
	    ip_done();
	    ip_set_uppercase(1);
	    ip_initialize(opt_ptr,outfile);
	    ip_cwk_add(":DELOC_INFO");
	    
	    error += (!ip_exist("DELOC_DEFINE",0));
	    ip_count("DELOC_DEFINE",&num_coords,0);
	    fnum_coords = num_coords;

	    /* don't know number of coordinates until now */
	    internals::mem_alloc();
	    deloc_define = init_matrix(num_coords,num_simples);
	    deloc_irrep = init_int_array(num_coords);
	    
	    for (i=0;i<num_coords;++i) {
		//ip_count("DELOC_DEFINE",&dim2,1,i);
		//if(dim2 != num_simples)
		//    punt("Definition length != number of simples");
		error += ip_data("DELOC_DEFINE","%d",&deloc_irrep[i],2,i,0);
		for (j=0;j<num_simples;++j) {
		    error += ip_data("DELOC_DEFINE","%lf", 
				     &deloc_define[i][j], 2, i, j+1);
		}
	    }
	    
	    fclose(opt_ptr);
	}

	/* need to compute current coordinate value */
	for(i=0;i<num_coords;++i)
	    for(j=0;j<num_simples;++j)
		coords[i] += deloc_define[i][j] * simples[j].get_val(); 
	
	if( (opt_ptr==NULL) || (error>0) )
	    punt("Trouble reading delocalized coordinate definitions");
    }	

    coord_base::read_opt();

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn deloc::read_bonds()
  \brief Reads bonds matrix from opt.dat()
/*---------------------------------------------------------------------------*/

void deloc :: read_bonds() {
	
	FILE *opt_ptr;
	
	int i, j, dim;
	
	ffile(&opt_ptr,"opt.dat",2);
	if( opt_ptr != NULL ) {      
	    
	    ip_done();
	    ip_set_uppercase(1);
	    ip_initialize(opt_ptr,outfile);
	    ip_cwk_add(":DELOC_INFO");
	    
	    error += (!ip_exist("BONDS",0));
	    ip_count("BONDS",&dim,0);

	    bonds = init_int_matrix(dim,dim);
	    
	    for (i=0;i<dim;++i) 
		for (j=0;j<dim;++j) 
		    error += ip_data("BONDS","%d", 
				     &bonds[i][j], 2, i, j);
	    
	    fclose(opt_ptr);
	}

    return;
}



/*---------------------------------------------------------------------------*/
/*! \fn deloc::write_opt()
  \brief Writes to opt.dat()

  Calls coord_base::write_opt then adds delocalized definitions */
/*---------------------------------------------------------------------------*/

void deloc :: write_opt() {

    coord_base::write_opt();

    int place, i, j;
    
    FILE *opt_ptr;

    ip_done();
    ffile(&opt_ptr,"opt.dat",1);
    ip_set_uppercase(1);
    ip_initialize(opt_ptr,outfile);
    
    fprintf(opt_ptr,"\ndeloc_info: (\n\n");

    /*write coordinate vector*/
    fprintf(opt_ptr,"  deloc_define = (\n\n");
    for(i=0;i<num_coords;++i) {
	fprintf(opt_ptr,"   ( %d ", deloc_irrep[i] );
	place=0;
	for(j=0;j<num_simples;++j) {
	    if( place==8 ) {
		place = 0;
		fprintf(opt_ptr,"\n        ");
	    }
	    fprintf(opt_ptr,"%lf  ",deloc_define[i][j]);
	    ++place;
	}
    fprintf(opt_ptr,")\n\n");
    }
    fprintf(opt_ptr,"          )\n");

    /*write bonds*/
    fprintf(opt_ptr,"  bonds = (\n\n");
    for(i=0;i<num_atoms;++i) {
	fprintf(opt_ptr,"   (  " );
	place=0;
	for(j=0;j<num_atoms;++j) {
	    if( place==8 ) {
		place = 0;
		fprintf(opt_ptr,"\n        ");
	    }
	    fprintf(opt_ptr,"%d  ",bonds[i][j]);
	    ++place;
	}
    fprintf(opt_ptr,")\n\n");
    }
    fprintf(opt_ptr,"          )\n");

    fprintf(opt_ptr,"        )\n");
    
    fclose(opt_ptr);
    return;
}

    
