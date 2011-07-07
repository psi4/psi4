/// \file testdiff1Db.cc
/// \brief testing for diff() in 1D

#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/mraimpl.h>

using namespace madness;

typedef Vector<double,1> coordT;
typedef Function<double,1> functionT;
typedef FunctionFactory<double,1> factoryT;
typedef Tensor<double> tensorT;

static const double PI = 3.1415926535897932384;
static const int k = 9 ; // Wavelet order (usually precision + 2)
static const double thresh = 1.e-7 ; // Precision
static const int init_lev = 4; 
static const int test_axis = 0;

void compare(World& world, functionT test, functionT exact, const char *str)
{
   double error = (exact - test).norm2() ;

   if (world.rank() == 0) {
       std::cerr << "Error in " << str << ": " << error ;

       if (error < thresh) {
           std::cerr << " PASSED " << std::endl ;
       }
       else {
           std::cerr << " FAILED " << std::endl ;
       }
   }
   return ;
}

double sech(double val)
{
  double tmp = exp(val) ;
  return(2. / (tmp + 1. / tmp) ) ;

}


// Testing the derivatives of a function 1+x
// for a variety of boundary conditions
static double u_exact(const coordT &pt) {
  double x = pt[0] ;
  return (cos(x)*(1. + tanh(4*x)) ) ;
}

static double du_exact(const coordT &pt) {
  double x = pt[0] ;
  return (4.*cos(x)*sech(4.*x)*sech(4.*x) - sin(x)*(1.+tanh(4.*x))) ;
}

static double left_dirichlet(const coordT &pt) {
  return (1.0) ;
}
static double right_dirichlet(const coordT &pt) {
  return (cos(1.0)*(1.+tanh(4.0))) ;
}
static double left_neumann  (const coordT &pt) {
  return (4.0) ;
}
static double right_neumann  (const coordT &pt) {
  return (4.*cos(1.)*sech(4.)*sech(4.0) - sin(1.)*(1.+tanh(4.))) ;
}


int main(int argc, char** argv) {
    
	MPI::Init(argc, argv);
	World world(MPI::COMM_WORLD);
        startup(world,argc,argv);

        std::cout.precision(6);

       // Function defaults
        FunctionDefaults<1>::set_k(k);
        FunctionDefaults<1>::set_thresh(thresh);
        FunctionDefaults<1>::set_refine(true );
        FunctionDefaults<1>::set_autorefine(true );
        FunctionDefaults<1>::set_initial_level(init_lev);
        FunctionDefaults<1>::set_cubic_cell( 0. , 1.);

        BoundaryConditions<1> bc;

        functionT  u      = factoryT(world).f( u_exact );
        functionT due     = factoryT(world).f(du_exact );

        functionT  left_d = factoryT(world).f( left_dirichlet) ;
        functionT right_d = factoryT(world).f(right_dirichlet) ;
        functionT  left_n = factoryT(world).f( left_neumann  ) ;
        functionT right_n = factoryT(world).f(right_neumann  ) ;

        // Right B.C.: Dirichlet
        // Left  B.C.: Free
        bc(0,0) = BC_DIRICHLET ;
        bc(0,1) = BC_FREE ;

        Derivative<double,1> dx1(world, test_axis, bc, left_d, right_d, k) ;
        functionT du1 = dx1(u) ;
        compare(world, du1, due, "du1") ;

        // Right B.C.: Free
        // Left  B.C.: Dirichlet
        bc(0,0) = BC_FREE ;
        bc(0,1) = BC_DIRICHLET ;

        Derivative<double,1> dx2(world, test_axis, bc, left_d, right_d, k) ;
        functionT du2 = dx2(u) ;
        compare(world, du2, due, "du2") ;

        // Right B.C.: Neumann
        // Left  B.C.: Free
        bc(0,0) = BC_NEUMANN ;
        bc(0,1) = BC_FREE ;

        Derivative<double,1> dx3(world, test_axis, bc, left_n, right_n, k) ;
        functionT du3 = dx3(u) ;
        compare(world, du3, due, "du3") ;

        // Right B.C.: Free
        // Left  B.C.: Neumann
        bc(0,0) = BC_FREE ;
        bc(0,1) = BC_NEUMANN ;

        Derivative<double,1> dx4(world, test_axis, bc, left_n, right_n, k) ;
        functionT du4 = dx4(u) ;
        compare(world, du4, due, "du4") ;

         world.gop.fence();

    finalize();
    
    return 0;
}

        
        
