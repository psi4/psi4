/*! \file lebedev_init.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

namespace psi { namespace CINTS {

int generate_points(int type, int start, double a, double b, double v);

double *x,*y,*z,*weight;

leb_sphere_t lebedev_init(int degree){

    int i;
    int start=0;
    double a,b,c,v;
    leb_sphere_t leb_tmp;
    
    x = init_array(degree);
    y = init_array(degree);
    z = init_array(degree);
    weight = init_array(degree);
    
    switch (degree){
	
    case 6:
	
	v = 0.1666666666666667;
	start = generate_points(1,start,a,b,v);
	break;
    
    case 14:
	
	v = 0.06666666666666667;
	start = generate_points(1,start,a,b,v);
	v = 0.07500000000000000;
	start = generate_points(3,start,a,b,v);
	break;
	
     case 26:
	
	v = 0.04761904761904762;
	start = generate_points(1,start,a,b,v);
	v = 0.03809523809523810;
	start = generate_points(2,start,a,b,v);
	v = 0.03214285714285714;
	start = generate_points(3,start,a,b,v);
	break;
	
    case 38:
	v = 0.009523809523809524;
	start = generate_points(1,start,a,b,v);
	v =0.03214285714285714;
	start = generate_points(3,start,a,b,v);
	a = 0.4597008433809831;
	v = 0.02857142857142857;
	start = generate_points(5,start,a,b,v);
	break;

    case 50:
	
	v = 0.01269841269841270;
	start = generate_points(1,start,a,b,v);
	v = 0.02257495590828924;
	start = generate_points(2,start,a,b,v);
	v = 0.02109375000000000;
	start = generate_points(3,start,a,b,v);
	a = 0.3015113445777636;
	v = 0.02017333553791887;
	start = generate_points(4,start,a,b,v);
	break;
	
    case 74:
	
	v = 0.0005130671797338464;
	start = generate_points(1,start,a,b,v);
	v = 0.01660406956574204;
	start = generate_points(2,start,a,b,v);
	v = -0.02958603896103896;
	start = generate_points(3,start,a,b,v);
	a = 0.4803844614152614;
	v = 0.02657620708215946;
	start = generate_points(4,start,a,b,v);
	a = 0.3207726489807764;
	v = 0.01652217099371571;
	start = generate_points(5,start,a,b,v);
	break;
	
    case 86:
	
	v = 0.01154401154401154;
	start = generate_points(1,start,a,b,v);
	v = 0.01194390908585628;
	start = generate_points(3,start,a,b,v);
	a = 0.3696028464541502;
	v = 0.01111055571060340;
	start = generate_points(4,start,a,b,v);
	a = 0.6943540066026664;
	v = 0.01187650129453714;
	start = generate_points(4,start,a,b,v);
	a = 0.3742430390903412;
	v = 0.01181230374690448;
	start = generate_points(5,start,a,b,v);
	break;
	
    case 110:
	
	v = 0.003828270494937162;
	start = generate_points(1,start,a,b,v);
	v = 0.009793737512487512;
	start = generate_points(3,start,a,b,v);
	a = 0.1851156353447362;
	v = 0.008211737283191111;
	start = generate_points(4,start,a,b,v);
	a = 0.6904210483822922;
	v = 0.009942814891178103;
	start = generate_points(4,start,a,b,v);
	a = 0.3956894730559419;
	v =0.009595471336070963;
	start = generate_points(4,start,a,b,v);
	a = 0.4783690288121502;
	v = 0.009694996361663028;
	start = generate_points(5,start,a,b,v);
	break;
	
    case 146:
	
	v = 0.0005996313688621381;
	start = generate_points(1,start,a,b,v);
	v = 0.007372999718620756;
	start = generate_points(2,start,a,b,v);
	v = 0.007210515360144488;
	start = generate_points(3,start,a,b,v);
	a = 0.6764410400114264;
	v = 0.007116355493117555;
	start = generate_points(4,start,a,b,v);
	a = 0.4174961227965453;
	v = 0.006753829486314477;
	start = generate_points(4,start,a,b,v);
	a = 0.1574676672039082;
	v = 0.007574394159054034;
	start = generate_points(4,start,a,b,v);
	a = 0.1403553811713183;
	b = 0.4493328323269557;
	v = 0.006991087353303262;
	start = generate_points(6,start,a,b,v);
	break;
	
    case 170:
	
	v = 0.005544842902037365; 
	start = generate_points(1,start,a,b,v);
	v = 0.006071332770670752;
	start = generate_points(2,start,a,b,v);
	v = 0.006383674773515093;
	start = generate_points(3,start,a,b,v);
	a = 0.2551252621114134;
	v = 0.005183387587747790;
	start = generate_points(4,start,a,b,v);
	a = 0.6743601460362766;
	v = 0.006317929009813725;
	start = generate_points(4,start,a,b,v);
	a = 0.4318910696719410;
	v = 0.006201670006589077;
	start = generate_points(4,start,a,b,v);
	a = 0.2613931360335988;
	v = 0.005477143385137348;
	start = generate_points(5,start,a,b,v);
	a = 0.4990453161796037;
	b = 0.1446630744325115;
	v = 0.005968383987681156;
	start = generate_points(6,start,a,b,v);
	break;
	
    case 194:
	
	v = 0.001782340447244611;
	start = generate_points(1,start,a,b,v);
	v = 0.005716905949977102;
	start = generate_points(2,start,a,b,v);
	v = 0.005573383178848738;
	start = generate_points(3,start,a,b,v);
	a = 0.6712973442695226;
	v = 0.005608704082587997;
	start = generate_points(4,start,a,b,v);
	a = 0.2892465627575439;
	v = 0.005158237711805383;
	start = generate_points(4,start,a,b,v);
	a = 0.4446933178717437;
	v = 0.005518771467273614;
	start = generate_points(4,start,a,b,v);
	a = 0.1299335447650067;
	v = 0.004106777028169394;
	start = generate_points(4,start,a,b,v);
	a = 0.3457702197611283;
	v = 0.005051846064614808;
	start = generate_points(5,start,a,b,v);
	a = 0.1590417105383530;
	b = 0.8360360154824589;
	v = 0.005530248916233094;
	start = generate_points(6,start,a,b,v);
	break;
	 
	/* need to add 230, 266 */
	
    case 302:
	
	v = 0.8545911725128148E-3; 
   	start = generate_points(1,start,a,b,v); 
   	v = 0.3599119285025571E-2; 
   	start = generate_points(3,start,a,b,v); 
   	a=0.3515640345570105E+0; 
   	v=0.3449788424305883E-2; 
   	start = generate_points(4,start,a,b,v); 
   	a=0.6566329410219612E+0; 
   	v=0.3604822601419882E-2; 
   	start = generate_points(4,start,a,b,v); 
   	a=0.4729054132581005E+0; 
   	v=0.3576729661743367E-2; 
   	start = generate_points(4,start,a,b,v); 
   	a=0.9618308522614784E-1; 
   	v=0.2352101413689164E-2; 
   	start = generate_points(4,start,a,b,v); 
   	a=0.2219645236294178E+0; 
   	v=0.3108953122413675E-2; 
   	start = generate_points(4,start,a,b,v); 
   	a=0.7011766416089545E+0; 
   	v=0.3650045807677255E-2; 
   	start = generate_points(4,start,a,b,v); 
   	a=0.2644152887060663E+0; 
   	v=0.2982344963171804E-2; 
   	start = generate_points(5,start,a,b,v); 
   	a=0.5718955891878961E+0; 
   	v=0.3600820932216460E-2; 
   	start = generate_points(5,start,a,b,v); 
   	a=0.2510034751770465E+0; 
   	b=0.8000727494073952E+0; 
   	v=0.3571540554273387E-2; 
   	start = generate_points(6,start,a,b,v); 
   	a=0.1233548532583327E+0; 
   	b=0.4127724083168531E+0; 
   	v=0.3392312205006170E-2; 
   	start = generate_points(6,start,a,b,v); 
	break;
	
	/* need to add 350, 434, 590, 770, 974, 1202, 1454, 1730,
	   2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294,
	   5810*/ 
    default:
	throw std::domain_error("\nAngular grid unrecognized");

    }
    leb_tmp.n_ang_points = degree;
    leb_tmp.points = (leb_point_t *)
	malloc(degree*sizeof(leb_point_t));
    for(i=0;i<degree;i++){
	leb_tmp.points[i].p_cart.x = x[i];
	leb_tmp.points[i].p_cart.y = y[i];
	leb_tmp.points[i].p_cart.z = z[i];
	leb_tmp.points[i].ang_weight = weight[i];
    }
    free(x);
    free(y);
    free(z);
    free(weight);
    return leb_tmp;
}





int generate_points(int type, int start, double a, double b, double v){
    
    int end;
    double c;

    switch (type){
	
    case 1:
	a = 1.0;
	
	x[start] = a;
	y[start] = 0.0;
	z[start] = 0.0;
	weight[start] = v;
	
	x[start+1] = -a;
	y[start+1] = 0.0;
	z[start+1] = 0.0;
	weight[start+1] = v;

	x[start+2] = 0.0;
	y[start+2] = a;
	z[start+2] = 0.0;
	weight[start+2] = v;

	x[start+3] = 0.0;
	y[start+3] = -a;
	z[start+3] = 0.0;
	weight[start+3] = v;

	x[start+4] = 0.0;
	y[start+4] = 0.0;
	z[start+4] = a;
	weight[start+4] = v;
 
	x[start+5] = 0.0;
	y[start+5] = 0.0;
	z[start+5] = -a;
	weight[start+5] = v;
	end = start+6;
	break;
	
    case 2:
	a = sqrt(0.5);
	x[start] = 0.0;
	y[start] = a;
	z[start] = a;
	weight[start] = v;
	
	x[start+1] = 0.0;
	y[start+1] = -a;
	z[start+1] = a;
	weight[start+1] = v;

	x[start+2] = 0.0;
	y[start+2] = a;
	z[start+2] = -a;
	weight[start+2] = v;

	x[start+3] = 0.0;
	y[start+3] = -a;
	z[start+3] = -a;
	weight[start+3] = v;
	
	x[start+4] = a;
	y[start+4] = 0.0;
	z[start+4] = a;
	weight[start+4] = v;

	x[start+5] = a;
	y[start+5] = 0.0;
	z[start+5] = -a;
	weight[start+5] = v;
	
	x[start+6] = -a;
	y[start+6] = 0.0;
	z[start+6] = a;
	weight[start+6] = v;
	
	x[start+7] = -a;
	y[start+7] = 0.0;
	z[start+7] = -a;
	weight[start+7] = v;
	
	x[start+8] = a;
	y[start+8] = a;
	z[start+8] = 0.0;
	weight[start+8] = v;

	x[start+9] = -a;
	y[start+9] = a;
	z[start+9] = 0.0;
	weight[start+9] = v;
	
	x[start+10] = a;
	y[start+10] = -a;
	z[start+10] = 0.0;
	weight[start+10] = v;

	x[start+11] = -a;
	y[start+11] = -a;
	z[start+11] = 0.0;
	weight[start+11] = v;
	end = start+12;
	break;
	
    case 3:
	a = sqrt(1.0/3.0);
	x[start] = a;
	y[start] = a;
	z[start] = a;
	weight[start] = v;
	
	x[start+1] = -a;
	y[start+1] = a;
	z[start+1] = a;
	weight[start+1] = v;
	
	x[start+2] = a;
	y[start+2] = -a;
	z[start+2] = a;
	weight[start+2] = v;

	x[start+3] = a;
	y[start+3] = a;
	z[start+3] = -a;
	weight[start+3] = v;
	
	x[start+4] = -a;
	y[start+4] = -a;
	z[start+4] = a;
	weight[start+4] = v;

	x[start+5] = a;
	y[start+5] = -a;
	z[start+5] = -a;
	weight[start+5] = v;
	
	x[start+6] = -a;
	y[start+6] = a;
	z[start+6] = -a;
	weight[start+6] = v;
	
	x[start+7] = -a;
	y[start+7] = -a;
	z[start+7] = -a;
	weight[start+7] = v;
	end = start+8;
	break;
	
    case 4:
    /* In this case A is inputed */
	b = sqrt(1.0 - 2.0*a*a);
	x[start] = a;
	y[start] = a;
	z[start] = b;
	weight[start] = v;
	
	x[start+1] = -a;
	y[start+1] = a;
	z[start+1] = b;
	weight[start+1] = v;
	
	x[start+2] = a;
	y[start+2] = -a;
	z[start+2] = b;
	weight[start+2] = v;

	x[start+3] = a;
	y[start+3] = a;
	z[start+3] = -b;
	weight[start+3] = v;
	
	x[start+4] = -a;
	y[start+4] = -a;
	z[start+4] = b;
	weight[start+4] = v;

	x[start+5] = -a;
	y[start+5] = a;
	z[start+5] = -b;
	weight[start+5] = v;
	
	x[start+6] = a;
	y[start+6] = -a;
	z[start+6] = -b;
	weight[start+6] = v;
	
	x[start+7] = -a;
	y[start+7] = -a;
	z[start+7] = -b;
	weight[start+7] = v;
	
	x[start+8] = -a;
	y[start+8] = b;
	z[start+8] = a;
	weight[start+8] = v;

	x[start+9] = a;
	y[start+9] = -b;
	z[start+9] = a;
	weight[start+9] = v;
	
	x[start+10] = a;
	y[start+10] = b;
	z[start+10] = -a;
	weight[start+10] = v;

	x[start+11] = -a;
	y[start+11] = -b;
	z[start+11] = a;
	weight[start+11] = v;
	
	x[start+12] = -a;
	y[start+12] = b;
	z[start+12] = -a;
	weight[start+12] = v;
	
	x[start+13] = a;
	y[start+13] = -b;
	z[start+13] = -a;
	weight[start+13] = v;
	
	x[start+14] = -a;
	y[start+14] = -b;
	z[start+14] = -a;
	weight[start+14] = v;

	x[start+15] = a;
	y[start+15] = b;
	z[start+15] = a;
	weight[start+15] = v;
	
	x[start+16] = b;
	y[start+16] = a;
	z[start+16] = a;
	weight[start+16] = v;

	x[start+17] = -b;
	y[start+17] = a;
	z[start+17] = a;
	weight[start+17] = v;
	
	x[start+18] = b;
	y[start+18] = -a;
	z[start+18] = a;
	weight[start+18] = v;
	
	x[start+19] = b;
	y[start+19] = a;
	z[start+19] = -a;
	weight[start+19] = v;
	
	x[start+20] = -b;
	y[start+20] = -a;
	z[start+20] = a;
	weight[start+20] = v;

	x[start+21] = -b;
	y[start+21] = a;
	z[start+21] = -a;
	weight[start+21] = v;
	
	x[start+22] = b;
	y[start+22] = -a;
	z[start+22] = -a;
	weight[start+22] = v;

	x[start+23] = -b;
	y[start+23] = -a;
	z[start+23] = -a;
	weight[start+23] = v;
	end = start + 24;
	break;
	
    case 5:
  	/* A is inputed in this case as well*/ 
	b=sqrt(1-a*a);
	x[start] = a;
	y[start] = b;
	z[start] = 0.0;
	weight[start] = v;
	
	x[start+1] = -a;
	y[start+1] = b;
	z[start+1] = 0.0;
	weight[start+1] = v;
	
	x[start+2] = a;
	y[start+2] = -b;
	z[start+2] = 0.0;
	weight[start+2] = v;

	x[start+3] = -a;
	y[start+3] = -b;
	z[start+3] = 0.0;
	weight[start+3] = v;
	
	x[start+4] = b;
	y[start+4] = a;
	z[start+4] = 0.0;
	weight[start+4] = v;

	x[start+5] = -b;
	y[start+5] = a;
	z[start+5] = 0.0;
	weight[start+5] = v;
	
	x[start+6] = b;
	y[start+6] = -a;
	z[start+6] = 0.0;
	weight[start+6] = v;
	
	x[start+7] = -b;
	y[start+7] = -a;
	z[start+7] = 0.0;
	weight[start+7] = v;
	
	x[start+8] = a;
	y[start+8] = 0.0;
	z[start+8] = b;
	weight[start+8] = v;

	x[start+9] = -a;
	y[start+9] = 0.0;
	z[start+9] = b;
	weight[start+9] = v;
	
	x[start+10] = a;
	y[start+10] = 0.0;
	z[start+10] = -b;
	weight[start+10] = v;
	
	x[start+11] = -a;
	y[start+11] = 0.0;
	z[start+11] = -b;
	weight[start+11] = v;
	
	x[start+12] = b;
	y[start+12] = 0.0;
	z[start+12] = a;
	weight[start+12] = v;
       
	x[start+13] = -b;
	y[start+13] = 0.0;
	z[start+13] = a;
	weight[start+13] = v;
	
	x[start+14] = b;
	y[start+14] = 0.0;
	z[start+14] = -a;
	weight[start+14] = v;

	x[start+15] = -b;
	y[start+15] = 0.0;
	z[start+15] = -a;
	weight[start+15] = v;
	
	x[start+16] = 0.0;
	y[start+16] = a;
	z[start+16] = b;
	weight[start+16] = v;

	x[start+17] = 0.0;
	y[start+17] = -a;
	z[start+17] = b;
	weight[start+17] = v;
	
	x[start+18] = 0.0;
	y[start+18] = a;
	z[start+18] = -b;
	weight[start+18] = v;
	
	x[start+19] = 0.0;
	y[start+19] = -a;
	z[start+19] = -b;
	weight[start+19] = v;

	x[start+20] = 0.0;
	y[start+20] = b;
	z[start+20] = a;
	weight[start+20] = v;

	x[start+21] = 0.0;
	y[start+21] = -b;
	z[start+21] = a;
	weight[start+21] = v;

	x[start+22] = 0.0;
	y[start+22] = b;
	z[start+22] = -a;
	weight[start+22] = v;

	x[start+23] = 0.0;
	y[start+23] = -b;
	z[start+23] = -a;
	weight[start+23] = v;
	end = start + 24;
	break;
	
    case 6:
       /* both A and B are inputed in this case */ 
	c=sqrt(1.0 - a*a - b*b);
	x[start] = a;
	y[start] = b;
	z[start] = c;
	weight[start] = v;

	x[start+1] = -a;
	y[start+1] = b;
	z[start+1] = c;
	weight[start+1] = v;

	x[start+2] = a;
	y[start+2] = -b;
	z[start+2] = c;
	weight[start+2] = v;

	x[start+3] = a;
	y[start+3] = b;
	z[start+3] = -c;
	weight[start+3] = v;

	x[start+4] = -a;
	y[start+4] = -b;
	z[start+4] = c;
	weight[start+4] = v;

	x[start+5] = a;
	y[start+5] = -b;
	z[start+5] = -c;
	weight[start+5] = v;

	x[start+6] = -a;
	y[start+6] = b;
	z[start+6] = -c;
	weight[start+6] = v;

	x[start+7] = -a;
	y[start+7] = -b;
	z[start+7] = -c;
	weight[start+7] = v;

	x[start+8] = b;
	y[start+8] = a;
	z[start+8] = c;
	weight[start+8] = v;

	x[start+9] = -b;
	y[start+9] = a;
	z[start+9] = c;
	weight[start+9] = v;

	x[start+10] = b;
	y[start+10] = -a;
	z[start+10] = c;
	weight[start+10] = v;

	x[start+11] = b;
	y[start+11] = a;
	z[start+11] = -c;
	weight[start+11] = v;

	x[start+12] = -b;
	y[start+12] = -a;
	z[start+12] = c;
	weight[start+12] = v;

	x[start+13] = b;
	y[start+13] = -a;
	z[start+13] = -c;
	weight[start+13] = v;

	x[start+14] = -b;
	y[start+14] = a;
	z[start+14] = -c;
	weight[start+14] = v;

	x[start+15] = -b;
	y[start+15] = -a;
	z[start+15] = -c;
	weight[start+15] = v;

	x[start+16] = c;
	y[start+16] = a;
	z[start+16] = b;
	weight[start+16] = v;

	x[start+17] = -c;
	y[start+17] = a;
	z[start+17] = b;
	weight[start+17] = v;

	x[start+18] = c;
	y[start+18] = -a;
	z[start+18] = b;
	weight[start+18] = v;

	x[start+19] = c;
	y[start+19] = a;
	z[start+19] = -b;
	weight[start+19] = v;

	x[start+20] = -c;
	y[start+20] = -a;
	z[start+20] = b;
	weight[start+20] = v;

	x[start+21] = c;
	y[start+21] = -a;
	z[start+21] = -b;
	weight[start+21] = v;

	x[start+22] = -c;
	y[start+22] = a;
	z[start+22] = -b;
	weight[start+22] = v;

	x[start+23] = -c;
	y[start+23] = -a;
	z[start+23] = -b;
	weight[start+23] = v;

	x[start+24] = c;
	y[start+24] = b;
	z[start+24] = a;
	weight[start+24] = v;

	x[start+25] = -c;
	y[start+25] = b;
	z[start+25] = a;
	weight[start+25] = v;

	x[start+26] = c;
	y[start+26] = -b;
	z[start+26] = a;
	weight[start+26] = v;

	x[start+27] = c;
	y[start+27] = b;
	z[start+27] = -a;
	weight[start+27] = v;

	x[start+28] = -c;
	y[start+28] = -b;
	z[start+28] = a;
	weight[start+28] = v;

	x[start+29] = c;
	y[start+29] = -b;
	z[start+29] = -a;
	weight[start+29] = v;

	x[start+30] = -c;
	y[start+30] = b;
	z[start+30] = -a;
	weight[start+30] = v;

	x[start+31] = -c;
	y[start+31] = -b;
	z[start+31] = -a;
	weight[start+31] = v;

	x[start+32] = a;
	y[start+32] = c;
	z[start+32] = b;
	weight[start+32] = v;

	x[start+33] = -a;
	y[start+33] = c;
	z[start+33] = b;
	weight[start+33] = v;

	x[start+34] = a;
	y[start+34] = -c;
	z[start+34] = b;
	weight[start+34] = v;

	x[start+35] = a;
	y[start+35] = c;
	z[start+35] = -b;
	weight[start+35] = v;

	x[start+36] = -a;
	y[start+36] = -c;
	z[start+36] = b;
	weight[start+36] = v;

	x[start+37] = a;
	y[start+37] = -c;
	z[start+37] = -b;
	weight[start+37] = v;

	x[start+38] = -a;
	y[start+38] = c;
	z[start+38] = -b;
	weight[start+38] = v;

	x[start+39] = -a;
	y[start+39] = -c;
	z[start+39] = -b;
	weight[start+39] = v;

	x[start+40] = b;
	y[start+40] = c;
	z[start+40] = a;
	weight[start+40] = v;

	x[start+41] = -b;
	y[start+41] = c;
	z[start+41] = a;
	weight[start+41] = v;

	x[start+42] = b;
	y[start+42] = -c;
	z[start+42] = a;
	weight[start+42] = v;

	x[start+43] = b;
	y[start+43] = c;
	z[start+43] = -a;
	weight[start+43] = v;

	x[start+44] = -b;
	y[start+44] = -c;
	z[start+44] = a;
	weight[start+44] = v;

	x[start+45] = b;
	y[start+45] = -c;
	z[start+45] = -a;
	weight[start+45] = v;

	x[start+46] = -b;
	y[start+46] = c;
	z[start+46] = -a;
	weight[start+46] = v;
 
	x[start+47] = -b;
	y[start+47] = -c;
	z[start+47] = -a;
	weight[start+47] = v;
	end = start + 48;
	break;

    default:
	throw std::domain_error("\n Why in the hell did you specify that case, didn't you read the code dipshit?");
    }
    return end;
}
};}
