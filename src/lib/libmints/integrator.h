#ifndef libmints_solver_integrator_H
#define libmints_solver_integrator_H
/*
 *  integrator.h
 *  Declaration of class Integrator for use in KS-DFT
 *  and pseudospectral methods
 *
 *  Created by Robert Parrish on 02/15/10.
 *
 */

#include "integrator_defines.h"
#include <stdlib.h>
#include <vector>
#include <string>

namespace psi { 

class Molecule;
class Vector3;
class GridBlock;
class Options;
class PSIO;

/*! \ingroup MINTS */
//! Integration Points/Weights container class
/*! Generates integration grids and computes points and weights
    for common molecular integration schemes
    References: Becke JCP 1988, Treutler+Alrichs JCP 1995
*/

class Integrator {
public:
    Integrator(shared_ptr<Molecule> mol, shared_ptr<PSIO> psio);
    ~Integrator();
    static shared_ptr<Integrator>  build_integrator(shared_ptr<Molecule> molecule, \
        shared_ptr<PSIO> psio, Options& options);

    void setRadialQuadrature(RadialScheme type, int npoints);
    void setSphericalQuadrature(SphericalScheme type, int npoints);
    void setNuclearWeighting(NuclearScheme type);
    void setVoronoiWeighting(VoronoiScheme type, double alpha = 0.64);
    void setNamedIntegrator(SpecialScheme name);
    void setBoxParameters(double delta, double overage);   
    void setVoronoiParameters(double a1 , double a2);   
    void setPointGrouping(GroupingScheme s);
    void setVoronoiCoordinates(CoordinateScheme s);
    void setPruningMode(PruningScheme s);
 
    int buildGrid(int max_block_size);
     
    unsigned long int getNPoints() const { return npoints_; }
    unsigned long int getNBlocks() const { return nblocks_; }
    int getBlockSize() const { return block_size_; }
    
    shared_ptr<GridBlock> getBlock(int N);
    shared_ptr<GridBlock> getAddress(unsigned long int n);
 
    std::string getString();

protected:
    void partitionGrid();
    void free_grid();
    void groupPointsVoronoi();
    void groupPointsBoxes();
    void applyNuclearWeights();
 
    RadialQuadrature getBeckeRadial(int n, double xi);
    RadialQuadrature getMuraRadial(int n, double xi);
    RadialQuadrature getMultiExpRadial(int n, double xi);
    RadialQuadrature getEMRadial(int n, double xi);
    RadialQuadrature getTreutlerRadial(int n, double xi, double alpha = 0.6);

    SphericalQuadrature getEMSpherical(int n);
    SphericalQuadrature getLebedevSpherical(int n);
    int getLebedevReccurencePoints(int type, int start, double a, double b, double v, const SphericalQuadrature & leb);
    
    int getSphericalOrderConstant(double, int, int, int);
    int getSphericalOrderSG1(double, int, int, int);
    int getSphericalOrderSG0(double, int, int, int);
    int getSphericalOrderAutomatic(double, int, int, int);

protected:

    shared_ptr<Molecule> mol_;

    unsigned long int npoints_;
    unsigned long int nblocks_;
    int block_size_;

    std::vector<unsigned long int> block_starts_;
    std::vector<unsigned long int> block_sizes_;

    int nspherical_;
    int nradial_;

    RadialScheme radial_type_; 
    SphericalScheme spherical_type_; 
    SpecialScheme special_type_; 
    NuclearScheme nuclear_type_; 
    VoronoiScheme voronoi_type_; 
    GroupingScheme group_type_; 
    CoordinateScheme coord_type_;
    PruningScheme prune_type_;

    double stratmann_a_;
    double box_overage_;
    double box_delta_;
    double vor_a1_;
    double vor_a2_;

    int* Z_;
    double* x_;
    double* y_;
    double* z_;
    double* w_;

    shared_ptr<PSIO> psio_;

    bool built_; 
};

typedef shared_ptr<Integrator> SharedIntegrator;

}
#endif
