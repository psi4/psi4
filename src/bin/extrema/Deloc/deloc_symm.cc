/*###########################################################################*/
/*! 
  \file
  \ingroup EXTREMA
  \brief Symmetry related functions. */
/*						Joseph P. Kenny 12/06/01
  ###########################################################################*/

#define EXTERN
#include"extrema.h"

using namespace psi::extrema;

/*---------------------------------------------------------------------------*/
/*! \fn void ir_project()
  \brief Cleans up symmetry of coordinates.*/
/*---------------------------------------------------------------------------*/

void deloc::ir_project() {

    int i,j,k,l, coord, p, p2,
	nspan, nproj, nproj_coord, ncoords_irrep,
	**nuclear_ict, **simples_ict,
	*project, *proj_irreps, 
	*n_indep;

    double **reps,
	*result, **results,
	**initial_irrep,
	**proj_temp,
	**proj_coords,
	**to_ortho, ***ortho;
   
    char_table ct(point_group);


    /* -----------------------------------
       determine reducible representations 
       for each coordinate 
       -------------------------------- */


    simples_ict = init_int_matrix(ct.num_ops,num_simples);
    result = init_array(num_simples);
    reps = init_matrix(num_nonzero,ct.num_irreps);

    nuclear_ict = chkpt_rd_ict();

    /* form simples transformation matrix */
    int atdex,bodex,andex,todex;
    for(i=0;i<ct.num_ops;++i) {
	for(j=0;j<num_simples;++j) {
	    if(simples[j].get_type()==BOND_TYPE) {
		atdex = nuclear_ict[i][simples[j].get_atom()]-1;
		bodex = nuclear_ict[i][simples[j].get_bond()]-1;
		simples_ict[i][j] = get_bond_index(atdex,bodex)+1;
	    }
	    if(simples[j].get_type()==ANGLE_TYPE) {
		atdex = nuclear_ict[i][simples[j].get_atom()]-1;
		bodex = nuclear_ict[i][simples[j].get_bond()]-1;
		andex = nuclear_ict[i][simples[j].get_angle()]-1;
		simples_ict[i][j] = get_angle_index(atdex,bodex,andex)+1;
	    }
	    if(simples[j].get_type()==TORS_TYPE) {
		atdex = nuclear_ict[i][simples[j].get_atom()]-1;
		bodex = nuclear_ict[i][simples[j].get_bond()]-1;
		andex = nuclear_ict[i][simples[j].get_angle()]-1;
		todex = nuclear_ict[i][simples[j].get_tors()]-1;
		simples_ict[i][j] = get_torsion_index(
		    atdex,bodex,andex,todex)+1;
		if( (ct.sym_ops[i][0] == 'S') || (ct.sym_ops[i][0] == 'I'))
		    simples_ict[i][j] *= -1;
	    }
	}
    }

    /* determine representations of the coordinate vectors */
    int new_index;
    for(coord=0;coord<num_nonzero;++coord) {
	for(i=0;i<num_simples;++i)
	    result[i] = 0.0;
	for(i=0;i<ct.num_irreps;++i) {
	    for(j=0;j<num_simples;++j) {
		new_index = abs(simples_ict[i][j])-1;
		result[new_index] =
		    deloc_define[coord][j]*((double) simples_ict[i][j])/
		    ((double) abs(simples_ict[i][j])); 
	    }
	    for(j=0;j<num_simples;++j) 
		reps[coord][i] += deloc_define[coord][j]*result[j];
	}
    }
    free(result);


    /*--------------------------
      reduce the representations
      ------------------------*/

    initial_irrep = 
	math_tools::rep_reduce(point_group, reps, num_nonzero);
    free_matrix(reps,num_nonzero);


    /*----------------------------
      perform symmetry projections
      --------------------------*/


    nspan = nproj = 0;
    for(coord=0;coord<num_nonzero;++coord) {
	for(i=0;i<ct.num_irreps;++i) {
	    
	    /* this counts everything */
	    if( initial_irrep[coord][i] > 1.0e-14 ) 
		++nspan; 

	    /* this counts only what we want to keep */
	    if( initial_irrep[coord][i] > IRREP_TOL )
		++nproj; 
	}
    }       

    fprintf(outfile,"\n  %d coordinates will be symmetry projected\n",nproj); 
    if(print_lvl>=RIDICULOUS_PRINT) {
	fprintf(outfile,"\n  Representation reduction coefficients:\n");
	print_mat(initial_irrep,num_nonzero,ct.num_irreps,outfile);
    }

    results = init_matrix(ct.num_irreps,num_simples);
    project = init_int_array(ct.num_irreps);
    proj_coords = init_matrix(nproj,num_simples);
    proj_irreps = init_int_array(nproj);
    
    /* project out irreps for each coordinate */
    p=p2=-1;
    for(coord=0;coord<num_nonzero;++coord) {
	nproj_coord=0;
	for(i=0;i<ct.num_irreps;++i) 
	    project[i] = 0;
	for(i=0;i<ct.num_irreps;++i)
	    if(initial_irrep[coord][i] > IRREP_TOL) {
		    project[i] = 1;
		    ++nproj_coord;
	    }

	if(nproj_coord) {
	    
	    /* determine result of each symmetry operation */
	    for(i=0;i<ct.num_irreps;++i)
		    for(j=0;j<num_simples;++j) {
			new_index = abs(simples_ict[i][j])-1;
			results[i][new_index] = 
			    ((double)simples_ict[i][j])
			    / ((double)abs(simples_ict[i][j])) 
			    * deloc_define[coord][j];
		    }    
			
	    /* project */
	    
	    proj_temp = math_tools::rep_project(
		point_group,num_simples,results,project);
	    for(i=0;i<nproj_coord;++i) {
		++p;
		for(j=0;j<num_simples;++j)
		    proj_coords[p][j] = proj_temp[i][j];
	    }
	    for(i=0;i<ct.num_irreps;++i)
		if(project[i]) {
		    proj_irreps[++p2]=i;
		}
	    free_matrix(proj_temp,nproj_coord);
	}	
    }
    free_int_matrix(simples_ict);   
    free_matrix(results,ct.num_irreps);
    free(project);
    free_matrix(initial_irrep,num_nonzero);
    free_matrix(deloc_define,num_nonzero);

    if(print_lvl>=RIDICULOUS_PRINT) {
	p=0;
	fprintf(outfile,"\n  Projected coordinate irreps:\n");
	for(i=0;i<nproj;++i) {
	    if(p==8) {
		p=0;
		fprintf(outfile,"\n");
	    }
	    fprintf(outfile,"   %5s",ct.irrep_labels[proj_irreps[i]]);
	    ++p;
	}
	fprintf(outfile,"\n  Projected coordinates:\n");
	print_mat(proj_coords,nproj,num_simples,outfile);
    }
    

    /*--------------------------
      orthonormalize coordinates
      ------------------------*/
    
    n_indep = init_int_array(ct.num_irreps);
    ortho = (double***) malloc(ct.num_irreps*sizeof(double**));
    
    /* orthonormalize by irrep */
    for(i=0;i<ct.num_irreps;++i) {
	p=-1; ncoords_irrep = 0;
	
	for(j=0;j<nproj;++j) 
	    if(proj_irreps[j]==i) 
		++ncoords_irrep;
	
	if(ncoords_irrep) {
	    to_ortho = init_matrix(ncoords_irrep,num_simples);
	    for(j=0;j<nproj;++j)
		if(proj_irreps[j]==i) {
		    ++p;
		    for(k=0;k<num_simples;++k)
			to_ortho[p][k] = proj_coords[j][k];
		}
	 
	    ortho[i] = math_tools::orthogonalize(
		ncoords_irrep,num_simples,to_ortho,1,1.0e-2,&n_indep[i]);
	    free_matrix(to_ortho,ncoords_irrep);
	}

    }		
    free_matrix(proj_coords,nproj);
    free(proj_irreps);


    /*--------------------------------
      copy into global arrays,
      clean up, check number of coords
      ------------------------------*/

    num_coords = 0;
    for(i=0;i<ct.num_irreps;++i)
	num_coords += n_indep[i];

    deloc_define = init_matrix(num_coords,num_simples);
    deloc_irrep = init_int_array(num_coords);

    p=-1;
    for(i=0;i<ct.num_irreps;++i)
	for(j=0;j<n_indep[i];++j) {
	    deloc_define[++p] = ortho[i][j];
	    deloc_irrep[p] = i;
	}
    free(ortho); /* only the top level double*** is freed! */
    free(n_indep);
    
    fprintf(outfile,"\n  %d coordinates remaining after orthogonalization\n",
	   num_coords);  
    p=0;
    fprintf(outfile,"\n  Coordinate irreps:\n");
    for(i=0;i<num_coords;++i) {
	if(p==8) {
	    p=0;
	    fprintf(outfile,"\n");
	}
	fprintf(outfile,"   %5s",ct.irrep_labels[deloc_irrep[i]]);
	++p;
    }
    fprintf(outfile,"\n\n  Final delocalized coordinates:\n",num_coords);
    print_mat(deloc_define,num_coords,num_simples,outfile);

    return;
}

    
/*-------------------------------------------------------------------------*/
/*! \fn deloc::get_bond_index(int at, int bo);
  \brief Returns bond index. */
/*-------------------------------------------------------------------------*/
	   
int deloc::get_bond_index(int at, int bo) {

    int i, index;
    for(i=0;i<num_simples;++i) 
	if(simples[i].get_type() == BOND_TYPE) 
	    if( ((simples[i].get_atom() == at) && 
		 (simples[i].get_bond() == bo)) ||
		((simples[i].get_atom() == bo) && 
		 (simples[i].get_bond() == at)))
		index = i;
	
    
    return index;
}



/*-------------------------------------------------------------------------*/
/*! \fn deloc::get_angle_index(int at, int bo, int an);
  \brief Returns angle index. */
/*-------------------------------------------------------------------------*/
	   
int deloc::get_angle_index(int at, int bo, int an) {
    
    int i, index=0;
    for(i=0;i<num_simples;++i) 
	if(simples[i].get_type() == ANGLE_TYPE) 
	    if( ((simples[i].get_atom() == at) && 
		 (simples[i].get_bond() == bo) &&
		 (simples[i].get_angle() == an)) ||
		((simples[i].get_atom() == an) && 
		 (simples[i].get_bond() == bo) &&
		 (simples[i].get_angle() == at))) 
		index = i;
    
    return index;
}



/*-------------------------------------------------------------------------*/
/*! \fn deloc::get_torsion_index(int at, int bo, int an, int to);
  \brief Returns torsion index. */
/*-------------------------------------------------------------------------*/
	   
int deloc::get_torsion_index(int at, int bo, int an, int to) {
    
    int i, index;
    for(i=0;i<num_simples;++i) 
	if(simples[i].get_type() == TORS_TYPE) 
	    if( ((simples[i].get_atom() == at) && 
		 (simples[i].get_bond() == bo) &&
		 (simples[i].get_angle() == an) && 
		 (simples[i].get_tors() ==to)) ||
		((simples[i].get_atom() == to) && 
		 (simples[i].get_bond() == an) &&
		 (simples[i].get_angle() == bo) && 
		 (simples[i].get_tors() == at)) )
	    index = i;
    
    return index;
}

