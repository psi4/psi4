#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/sdf_shape_3D.h>
#include <mra/funcplot.h>
#include <examples/molecularmask.h>
#include <constants.h>
#include <vector>

using namespace madness;
using namespace std;

// Crude macro for timing
double XXstart;
#define TIME(MSG,X) XXstart=wall_time();         \
                    X; \
                    if (world.rank() == 0) print("timer:",MSG,"used",wall_time()-XXstart) \

// Area of two intersecting spheres separated by d
double area_two_spheres(double r1, double r2, double d) {
    if (d > r1+r2) d = r1 + r2;
    return constants::pi*(2*d*r1*r1+2*d*r2*r2+r1*r1*r1+r1*d*d-r1*r2*r2+r2*d*d-r2*r1*r1+r2*r2*r2)/d;
}

// Volume of two intersecting spheres separated by d
double volume_two_spheres(double r1, double r2, double d) {
    if (d > r1+r2) return 0.0;
    double overlap = constants::pi*(r1+r2-d)*(r1+r2-d)*(d*d + 2*d*r1 + 2*d*r2 - 3*r1*r1 - 3*r2*r2 + 6*r1*r2)/(12*d);
    return 4.0*constants::pi*(r1*r1*r1 + r2*r2*r2)/3.0 - overlap;
}


int main(int argc, char **argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);

    const int k = 6; // wavelet order
    const double thresh = 1e-4; // truncation threshold
    const double L = 5; // box is [-L,L]
    const int natom = 2; // number of atoms
    const double sigma = 0.1; // Surface width

    // Function defaults
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_cubic_cell(-L, L);
    FunctionDefaults<3>::set_initial_level(2);

    // Set up atomic coordinates and radii
    vector<double> atomic_radii(natom);
    vector<coord_3d> atomic_coords(natom);
    for (int i=0; i<natom; i++) {
        atomic_radii[i] = 1.0 + i*0.5;
        atomic_coords[i][0] = 0.0;
        atomic_coords[i][1] = 0.0;
        atomic_coords[i][2] = i*1.5;
        print("atom",i,atomic_coords[i],atomic_radii[i]);
    }

    real_functor_3d volume_functor(new MolecularVolumeMask(sigma, atomic_radii, atomic_coords));
    real_functor_3d surface_functor(new MolecularSurface(sigma, atomic_radii, atomic_coords));

    TIME("make volume ", real_function_3d volume = real_factory_3d(world).functor(volume_functor));
    TIME("make surface", real_function_3d surface = real_factory_3d(world).functor(surface_functor).truncate_on_project());

    print("the volume is", volume.trace());
    print("the area   is", surface.trace());

    if (natom == 2) {
        double d = distance(atomic_coords[0],atomic_coords[1]);
        double r1 = atomic_radii[0];
        double r2 = atomic_radii[1];
        print("d",d,"r1",r1,"r2",r2,"r1+r2-d",r1+r2-d);
        print("analytic volume intersecting spheres", volume_two_spheres(r1, r2, d));
        print("analytic area   intersecting spheres", area_two_spheres(r1, r2, d));
    }

    std::vector<long> npt(3,401);
    Tensor<double> cell(3,2);
    cell(_,0) = -5;
    cell(_,1) =  5;

    TIME("plot surface",plotdx(surface, "surface.dx"));
    TIME("plot volume ",plotdx(volume, "volume.dx", cell, npt));
    TIME("plot povray ",plotpovray(volume, "volume.df3", cell, npt));

    finalize();
    return 0;
}
