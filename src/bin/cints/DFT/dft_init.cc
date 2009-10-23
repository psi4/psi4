/*! \file dft_init.cc
    \ingroup CINTS
    \author Shawn Brown
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"functional.h"
#include"functional_u.h"
#include"calc_den_fast.h"
#include"calc_den_u.h"
#include"calc_grad_fast.h"
#include"calc_close_basis.h"
#include"calc_close_basis_u.h"
#include"lebedev_init.h"
#include"grid_init.h"
#include"free_grid_structs.h"


namespace psi { namespace cints {
  void dft_init(void){
    
    int errcod;
    int i,j;
    int si,sj; 
    int si_n_ao;
    int sj_n_ao;
    int rpointstmp;
    int angpointstmp;
    int n_shells;
    int elem_in_func;
    char *exchstring;
    char *corrstring;
    char *gridstring;
    char *funcstring;
    double tmpa, tmpb;
    
    /* DFT input parsing */
    elem_in_func = 0;
    errcod = ip_count("FUNCTIONAL",&elem_in_func,0);
    if(elem_in_func == 0){
	errcod=ip_string("FUNCTIONAL",&funcstring,0);
	if(errcod == IPE_OK){
	    if(!strcmp(funcstring,"XALPHA")){
		UserOptions.hf_exch = 0.0;
		if(UserOptions.reftype == rhf){
		DFT_options.exchange_func = slater_ed;
		DFT_options.den_calc = calc_density_fast;
		DFT_options.correlation_func = no_funct;
		}
		else{
		    DFT_options.exchange_func = slater_ed;
		    DFT_options.den_calc = calc_density_u;
		    DFT_options.correlation_func = no_funct;
		}
	    }
	    else
		throw std::domain_error("Unrecognized fucntional specified with keyword FUNCTIONALu");
	}
	else
	    throw std::domain_error("Must define a functional with keyword FUNCTIONAL");
    }
    
    else if(elem_in_func == 2){
	errcod = ip_string("FUNCTIONAL",&exchstring,1,0);
	
        if(!strcmp(exchstring,"SLATER")||!strcmp(exchstring,"S")){
	    UserOptions.hf_exch = 0.0;
	    if(UserOptions.reftype == rhf){
		DFT_options.exchange_func = slater_ed;
		DFT_options.den_calc = calc_density_fast;
	    }
	    else{
	      DFT_options.exchange_func = slater_u_ed;     
	      DFT_options.den_calc = calc_density_u;
	    }
	}
	else if(!strcmp(exchstring,"B88")||!strcmp(exchstring,"B")){
	    UserOptions.hf_exch = 0.0;
	    if(UserOptions.reftype == rhf){
		DFT_options.exchange_func = Becke88_ed;
		DFT_options.den_calc = calc_grad_fast;
		exchstring[0] = '\0';
	    }
	    else{
		throw std::domain_error("UHF functional not implemented");
	    }
	}
	else if(!strcmp(exchstring,"NONE")){
	    UserOptions.hf_exch = 0.0;
	    if(UserOptions.reftype == rhf){
		DFT_options.exchange_func = no_funct;
		DFT_options.den_calc = calc_density_fast;
		exchstring[0] = '\0';
	    }
	    else{
		DFT_options.exchange_func = no_funct_u;
		DFT_options.den_calc = calc_density_u;
		exchstring[0] = '\0';
	    }
	}
	else
	    throw std::domain_error("Unrecognized or nonimplemented exchange functional specified");
	
	errcod = ip_string("FUNCTIONAL",&corrstring,1,1);
	
	if(!strcmp(corrstring,"NONE")){
	    UserOptions.hf_exch = 0.0;
	    if(UserOptions.reftype == rhf){
		DFT_options.correlation_func = no_funct;	
		corrstring[0] = '\0';
	    }
	    else{
		DFT_options.correlation_func = no_funct_u;		
		corrstring[0] = '\0';
	    }
	}
        else if(!strcmp(corrstring,"VWN5")){
	    if(UserOptions.reftype == rhf){
		DFT_options.correlation_func = VWN5_ed;
	    }
	    else{
		DFT_options.correlation_func = VWN5_u_ed;
	    }
	}           
	else if(!strcmp(corrstring,"VWN4")){
	    if(UserOptions.reftype == rhf){
		DFT_options.correlation_func = VWN4_ed;
	    }
	    else
		DFT_options.correlation_func = VWN4_u_ed;      
        }        
	else
	    throw std::domain_error("Unrecognized or nonimplemented correlation fuctional specified");

	i = strlen(exchstring) + strlen(corrstring);
	if (i == 0)
	  funcstring = strdup("NONE");
	else {
	  funcstring = (char *) malloc(sizeof(char)* (i + 1));
	  sprintf(funcstring,"%s%s",exchstring,corrstring);
	}
	free(exchstring); free(corrstring);
    }
    else
	throw std::domain_error("Something wrong in the specification of the FUNCTIONAL keyword");    
  }
  
  void cleanup_dft_options(DFT_options_t DFT_options){
    
    free(DFT_options.basis);
    free(DFT_options.Bragg);
    cleanup_grid_type(DFT_options.grid);
  }
}

}
