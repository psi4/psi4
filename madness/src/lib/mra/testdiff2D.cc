/// \file testdiff2D.cc
/// \brief testing for diff() in 2D

#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/mraimpl.h>
#include <string>

using namespace madness;

typedef Vector<double,2> coordT;
typedef Function<double,2> functionT;
typedef FunctionFactory<double,2> factoryT;
typedef Tensor<double> tensorT;

static const double PI = 3.1415926535897932384;
static const int k = 9 ; // Wavelet order (usually precision + 2)
static const double thresh = 1.e-7 ; // Precision
static const int init_lev = 2; 
static int test_axis = 0; 
static const double Length = 2.;

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


// Testing the derivatives of a function 1+x
// for a variety of boundary conditions
static double u_exact(const coordT &pt) {
  return (1.0+pt[0]*pt[1]) ;
}

static double dudx_exact(const coordT &pt) {
  return (pt[1]) ;
}

static double dudy_exact(const coordT &pt) {
  return (pt[0]) ;
}

static double xleft_dirichlet(const coordT &pt) {
  return (1.) ;
}
static double xright_dirichlet(const coordT &pt) {
  double x = Length, y=pt[1];
  return (1.0+x*y) ;
}
static double xleft_neumann  (const coordT &pt) {
  return (pt[1]) ;
}
static double xright_neumann  (const coordT &pt) {
  return (pt[1]) ;
}

static double yleft_dirichlet(const coordT &pt) {
  return (1.) ;
}
static double yright_dirichlet(const coordT &pt) {
  double x=pt[0], y = Length;
  return (1.0+x*y) ;
}
static double yleft_neumann  (const coordT &pt) {
  return (pt[0]) ;
}
static double yright_neumann  (const coordT &pt) {
  return (pt[0]) ;
}


int main(int argc, char** argv) {
    initialize(argc, argv);
	World world(MPI::COMM_WORLD);
        startup(world,argc,argv);

        std::cout.precision(6);

       // Function defaults
        FunctionDefaults<2>::set_k(k);
        FunctionDefaults<2>::set_thresh(thresh);
        FunctionDefaults<2>::set_refine(true );
        FunctionDefaults<2>::set_autorefine(true );
        FunctionDefaults<2>::set_initial_level(init_lev);
        FunctionDefaults<2>::set_cubic_cell( 0. , Length);

        BoundaryConditions<2> bc;

        functionT  u      = factoryT(world).f(   u_exact );
        functionT dudxe   = factoryT(world).f(dudx_exact );
        functionT dudye   = factoryT(world).f(dudy_exact );

        functionT  xleft_d = factoryT(world).f( xleft_dirichlet) ;
        functionT xright_d = factoryT(world).f(xright_dirichlet) ;
        functionT  xleft_n = factoryT(world).f( xleft_neumann  ) ;
        functionT xright_n = factoryT(world).f(xright_neumann  ) ;

        functionT  yleft_d = factoryT(world).f( yleft_dirichlet) ;
        functionT yright_d = factoryT(world).f(yright_dirichlet) ;
        functionT  yleft_n = factoryT(world).f( yleft_neumann  ) ;
        functionT yright_n = factoryT(world).f(yright_neumann  ) ;

        // Derivative in the x-direction
        test_axis = 0;

        // X Right B.C.: Dirichlet
        // X Left  B.C.: Free
        bc(0,0) = BC_DIRICHLET ;
        bc(0,1) = BC_FREE ;
        bc(1,0) = BC_FREE ;
        bc(1,1) = BC_FREE ;
        bc = bc ;

/*
//Create a derived class FreeSpaceDerivative, PeriodicDerivative
        FreeSpaceDerivative<double,2> dx1(world, k, test_axis) ; // Free Space
        FreeSpaceDerivative<double,2> dx1(world, k, test_axis, 2) ; // Free Space
        PeriodicDerivative<double,2> dx1(world, k, test_axis, 1) ; // Periodic
*/
        Derivative<double,2> dx1(world, test_axis, bc, xleft_d, xright_d, k) ;
        functionT dudx1 = dx1(u) ;
        compare(world, dudx1, dudxe, "dudx1") ;

        // X Right B.C.: Free
        // X Left  B.C.: Dirichlet
        bc(0,0) = BC_FREE ;
        bc(0,1) = BC_DIRICHLET ;
        bc = bc ;

        Derivative<double,2> dx2(world, test_axis, bc, xleft_d, xright_d, k) ;
        functionT dudx2 = dx2(u) ;
        compare(world, dudx2, dudxe, "dudx2") ;

        // X Right B.C.: Neumann
        // X Left  B.C.: Free
        bc(0,0) = BC_NEUMANN ;
        bc(0,1) = BC_FREE ;
        bc = bc ;

        Derivative<double,2> dx3(world, test_axis, bc, xleft_n, xright_n, k) ;
        functionT dudx3 = dx3(u) ;
        compare(world, dudx3, dudxe, "dudx3") ;

        // X Right B.C.: Free
        // X Left  B.C.: Neumann
        bc(0,0) = BC_FREE ;
        bc(0,1) = BC_NEUMANN ;
        bc = bc ;

        Derivative<double,2> dx4(world, test_axis, bc, xleft_n, xright_n, k) ;
        functionT dudx4 = dx4(u) ;
        compare(world, dudx4, dudxe, "dudx4") ;


        // Derivative in the y-direction
        test_axis = 1;

        // Y Right B.C.: Dirichlet
        // Y Left  B.C.: Free
        bc(0,0) = BC_FREE ;
        bc(0,1) = BC_FREE ;
        bc(1,0) = BC_DIRICHLET ;
        bc(1,1) = BC_FREE ;
        bc = bc ;

        Derivative<double,2> dy1(world, test_axis, bc, yleft_d, yright_d, k) ;
        functionT dudy1 = dy1(u) ;
        compare(world, dudy1, dudye, "dudy1") ;

        // Y Right B.C.: Free
        // Y Left  B.C.: Dirichlet
        bc(1,0) = BC_FREE ;
        bc(1,1) = BC_DIRICHLET ;
        bc = bc ;

        Derivative<double,2> dy2(world, test_axis, bc, yleft_d, yright_d, k) ;
        functionT dudy2 = dy2(u) ;
        compare(world, dudy2, dudye, "dudy2") ;

        // Y Right B.C.: Neumann
        // Y Left  B.C.: Free
        bc(1,0) = BC_NEUMANN ;
        bc(1,1) = BC_FREE ;
        bc = bc ;

        Derivative<double,2> dy3(world, test_axis, bc, yleft_n, yright_n, k) ;
        functionT dudy3 = dy3(u) ;
        compare(world, dudy3, dudye, "dudy3") ;

        // Y Right B.C.: Free
        // Y Left  B.C.: Neumann
        bc(1,0) = BC_FREE ;
        bc(1,1) = BC_NEUMANN ;
        bc = bc ;

        Derivative<double,2> dy4(world, test_axis, bc, yleft_n, yright_n, k) ;
        functionT dudy4 = dy4(u) ;
        compare(world, dudy4, dudye, "dudy4") ;

         world.gop.fence();

    finalize();
    
    return 0;
}

        
        
