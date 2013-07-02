#include "extra.h"
#include <world/parar.h>
//using namespace madness;
std::ostream& operator<<(std::ostream& s, const InputParameters& p) {
    s << p.L<< " " << p.Lsmall<< " " << p.Llarge<< " " << p.F << " " << p.omega <<
        " " << p.ncycle << " " << p.Z << " " << p.R[0]<< " " << p.k<< " " <<
        p.thresh<< " " << p.cut<< " " << p.iState << " " << p.prefix<< " " << p.ndump<< " " <<
        p.nplot << " " << p.nprint << " "  << p.nloadbal << " " << p.nio << p.tScale << std::endl;
return s;
}

InputParameters param;


const char* wave_function_filename(int step) {
    static char fname[1024];
    sprintf(fname, "%s-%5.5d", prefix.c_str(), step);
    return fname;
}
bool wave_function_exists(World& world, int step) {
    return archive::ParallelInputArchive::exists(world, wave_function_filename(step));
}
void wave_function_store(World& world, int step, const complex_functionT& psi) {
    archive::ParallelOutputArchive ar(world, wave_function_filename(step), nIOProcessors);
    ar & psi;
}
complex_functionT wave_function_load(World& world, int step) {
    complex_functionT psi;
    archive::ParallelInputArchive ar(world, wave_function_filename(step));
    ar & psi;
    return psi;
}
template<class T>
std::string toString( const T& a ) {
    std::ostringstream o;
    o << a[0] << ", " << a[1] << ", " << a[2];
    return o.str();
}
void loadDefaultBasis(World& world, std::vector<WF>& boundList, double Z) {
    PRINT("Loading the default basis");
    const int NBSt = 3;    //Number of Bound States
    const int bSt[][3] = { {1,0,0},
                           {2,0,0},
                           {2,1,0} };
    for( int i=0; i<NBSt; i++ ) {
       boundList.push_back( WF(toString(bSt[i]), 
                 FunctionFactory<complexd,NDIM>(world).functor(functorT(
                 new BoundWF(Z , bSt[i][0], bSt[i][1], bSt[i][2]) ))));
    }
    PRINTLINE("Done loading the standard basis");
}

void loadList(World& world, std::vector<std::string>& boundList, std::vector<std::string>& unboundList) {
    std::ifstream bound("bound.num");
    std::ifstream unbound("unbound.num");
    if( ! bound.is_open() && ! unbound.is_open() ) {
        PRINTLINE("bound.num and unbound.num not found");
        boundList.push_back("1 0 0");
        boundList.push_back("2 1 0");
    } else {
        if(bound.is_open()) {
            int n,l,m;
            std::ostringstream nlm;
            while(bound >> n) {
                bound >> l;
                bound >> m;
                nlm << n << " " << l << " " << m;
                boundList.push_back(nlm.str());
                nlm.str("");
            }
            bound.close();
        } else PRINTLINE("bound.num not found");
        if(unbound.is_open()) {
            double kx, ky, kz;
            std::ostringstream kxyz;
            while(unbound >> kx) {
                unbound >> ky;
                unbound >> kz;
                kxyz.precision( 8 );
                kxyz << std::fixed;
                kxyz << kx << " " << ky << " " << kz;
                unboundList.push_back(kxyz.str());
                kxyz.str("");
            }
            unbound.close();
        } else PRINTLINE("unbound.num not found");
    }
}

// void loadBasis(World& world, std::vector<WF>& boundList, std::vector<WF>& unboundList, const double Z, const double cutoff) {
//     std::ifstream bound("bound.num");
//     std::ifstream unbound("unbound.num");
//     if( ! bound.is_open() && ! unbound.is_open() ) {
//         PRINTLINE("bound.num and unbound.num not found");
//         loadDefaultBasis(world,boundList,Z);
//     } else {
//         if(bound.is_open()) {
//             PRINTLINE("Calculating bound quantum states");
//             int n,l,m;
//             std::ostringstream nlm;
//             while(bound >> n) {
//                 bound >> l;
//                 bound >> m;
//                 nlm << n << l << m;
//                 double start = wall_time();
//                 PRINT(nlm.str());                
//                 boundList.push_back(WF(nlm.str()+"           ",
//                                        FunctionFactory<complexd,NDIM>(world).
//                                        functor(functorT(new BoundWF(Z,n,l,m)))));
//                 double used = wall_time() - start;
//                 PRINTLINE("\t" << used << " sec" );
//                 nlm.str("");
//             }
//             bound.close();
//         } else PRINTLINE("bound.num not found");
//         if(unbound.is_open()) {
//             PRINTLINE("Calculating unbound quantum states");
//             double kx, ky, kz;
//             std::ostringstream kxyz;
//             while(unbound >> kx) {
//                 unbound >> ky;
//                 unbound >> kz;
//                 kxyz.precision( 2 );
//                 kxyz << std::fixed;
//                 kxyz << kx << " " << ky << " " << kz;
//                 PRINT(kxyz.str());
//                 double start = wall_time();
//                 const double kvec[] = {kx, ky, kz};
//                 unboundList.push_back(WF(kxyz.str(),
//                                        FunctionFactory<complexd,NDIM>(world).
//                                        functor(functorT(new PhiK(Z,kvec,cutoff)))));
//                 double used = wall_time() - start;
//                 PRINTLINE("\t" << used << " sec");
//                 kxyz.str("");
//             }
//             unbound.close();
//         } else PRINTLINE("unbound.num not found");
//     }
// }
complexd zdipole( const vector3D& r) {
    return complexd(r[2],0.0);
}

/*****************************************************************
 * Dipole matrix elements available for perturbation calculations
 * |<phi_A|z|phi_B>|^2
 *****************************************************************/
void projectZdip(World& world, std::vector<WF> stateList) {
    std::cout.precision(8);
    complex_functionT z = complex_factoryT(world).f(zdipole);
    complexd output;
    std::vector<WF>::iterator basisI;
    std::vector<WF>::iterator basisII;
    PRINT(std::endl << "\t\t|<basis_m|z|basis_n>|^2 " << std::endl << "\t\t");
    for(basisII=stateList.begin(); basisII != stateList.end(); basisII++) {
        PRINT("|" << basisII->str << ">");
    }
    PRINT("\n");
    for(basisI = stateList.begin(); basisI !=stateList.end(); basisI++ ) {
        PRINT("<" << basisI->str << "|" );
        for(basisII = stateList.begin(); basisII != stateList.end(); basisII++) {
            output = inner(basisII->func,z*basisI->func); 
            PRINT(" " << real(output*conj(output)) <<"\t");
        }
        PRINT("\n");
    }
    PRINT("\n");
}


void compare1F1(World& world, double cutoff) {
    //load param
    std::string tag;
    double rMIN = 0.0;
    double rMAX = 10.0;
    double dr   = 1.0;
    double k    = 1.0; //NOT USED
    double Z    = 1.0;
    /***************************************
     *Load graphing parameters from the file: param
     * rMIN 0.0
     * rMAX 10.0
     * dr   1.0
     * TH   0.0
     ****************************************/
    std::ifstream f("param");
    if( f.is_open() ) {
        while(f >> tag) {
            if (tag[0] == '#') {
                char ch;
                PRINTLINE("    comment  " << tag.c_str());
                while (f.get(ch)) {
                    PRINTLINE(ch);
                    if (ch == '\n') break;
                }
            }
            else if (tag == "rMIN") {
                f >> rMIN;
                PRINTLINE("rMIN = " << rMIN);
            }
            else if (tag == "rMAX") {
                f >> rMAX;
                PRINTLINE("rMAX = " << rMAX);
            }
            else if (tag == "dr") {
                f >> dr;
                PRINTLINE("dr = " << dr);
            }
            else if (tag == "k") {
                f >> k;
                PRINTLINE("k = " << k);
            }
        }
    }
    //make functor
    const madness::Vector<double,3> kvec = vec(0.0, 0.0, 1.0);
    PhiK phi_k =  PhiK(Z, kvec, cutoff);
    complexd ONE(1.0,0.0);
    complexd I(0.0,1.0);
    std::cout << std::fixed;
    for(double r=rMIN; r<rMAX; r+=dr) {
        complexd ZZ(0.0,-r);
        std::cout.precision(2);
        PRINT(r                         << "\t");
        std::cout.precision(8);
        PRINT(real(conhyp(-I/k,ONE,ZZ)) << "\t");
        PRINT(imag(conhyp(-I/k,ONE,ZZ)) << "\t");
        PRINT(real(phi_k.aForm(ZZ))     << "\t");
        PRINT(imag(phi_k.aForm(ZZ))     << "\n");
    }
}

double myreal(double t) {return t;}
double myreal(const double_complex& t) {return real(t);}

// Invoke as \c u(r/c)/c where \c c is the radius of the smoothed volume.  
static double smoothed_potential(double r) {
    double r2 = r*r;
    double pot;
    if (r > 6.5){
        pot = 1.0/r;
    } else if (r > 1e-2) {
        pot = erf(r)/r + exp(-r2)*0.56418958354775630;
    } else{
        pot = 1.6925687506432689-r2*(0.94031597257959381-r2*(0.39493270848342941-0.12089776790309064*r2));
    }
    return pot;
}


// Nuclear attraction potential
complexd V(const vector3D& r) {
    const double x=r[0], y=r[1], z=r[2];
    double sum = 0.0;
    for (int i=0; i<param.natom; i++) {
      double xx = x-param.R[i][0];
      double yy = y-param.R[i][1];
      double zz = z-param.R[i][2];
      double rr = sqrt(xx*xx+yy*yy+zz*zz);
      sum +=  -param.Z[i]*smoothed_potential(rr/param.cut)/param.cut;
    }
    return complexd(sum,0.0);
}

// Given psi and V evaluate the energy ... leaves psi compressed, potn reconstructed
template <typename T>
double energy(World& world, const Function<T,3>& psi, const Function<T,3>& potn) {
    // First do all work in the scaling function basis
    psi.reconstruct();
    bool DOFENCE = false;
    Derivative<T,3> Dx(world,0), Dy(world,1), Dz(world,2);
    Function<T,3> dx = Dx(psi,DOFENCE);
    Function<T,3> dy = Dy(psi,DOFENCE);
    Function<T,3> dz = Dz(psi,DOFENCE);
    Function<T,3> Vpsi = psi*potn;
    // Now do all work in the wavelet basis
    psi.compress(DOFENCE); Vpsi.compress(DOFENCE); dx.compress(DOFENCE); dy.compress(DOFENCE); dz.compress(true);
    T S = psi.inner(psi);
    T PE = psi.inner(Vpsi);
    T KE = 0.5*(inner(dx,dx) + inner(dy,dy) + inner(dz,dz));
    T E = (KE+PE)/S;
    dx.clear(); dy.clear(); dz.clear(); Vpsi.clear(); // To free memory on return
    world.gop.fence();
    return myreal(E);
}
void converge(World& world, complex_functionT& potn, complex_functionT& psi, double& eps) {
     for (int iter=0; iter<30; iter++) {
         SeparatedConvolution<double,NDIM> op =
             BSHOperator3D(world, sqrt(-2*eps), param.cut, param.thresh);
         complex_functionT Vpsi = (potn*psi);
         Vpsi.scale(-2.0).truncate();
         complex_functionT tmp = apply(op,Vpsi).truncate(param.thresh);
         double norm = tmp.norm2();
         complex_functionT r = tmp-psi;
         double rnorm = r.norm2();
         double eps_new = eps - 0.5*real(inner(Vpsi,r))/(norm*norm);
         // if (world.rank() == 0) {
         //     print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
         // }
         psi = tmp.scale(1.0/norm);
         eps = eps_new;
         if (rnorm < std::max(1e-5,param.thresh)) break;
     }
     psi.truncate(param.thresh);     
 }


void compareGroundState(World& world, double Z) {
    //Create the softened ground state for a list of cuts
    const int nCUT = 5;
    double cut[nCUT] = { 0.3, 0.2, 0.1, 0.05, 0.03 };
    complexd output;
    if (world.rank() == 0) param.read("input");
    PRINTLINE("            |<1s|psi0>|^2 \t |<2s|psi0>|^2 \t |<3s|psi0>|^2 \t |<4s|psi0>|^2");
    //make |ns> states
    complex_functionT oneS = FunctionFactory<complexd,NDIM>(world).
        functor(functorT(new BoundWF(Z, 1, 0, 0)));
    complex_functionT twoS = FunctionFactory<complexd,NDIM>(world).
        functor(functorT(new BoundWF(Z, 2, 0, 0)));
    complex_functionT threeS = FunctionFactory<complexd,NDIM>(world).
        functor(functorT(new BoundWF(Z, 3, 0, 0)));
    complex_functionT fourS = FunctionFactory<complexd,NDIM>(world).
        functor(functorT(new BoundWF(Z, 4, 0, 0)));
    for(int i=0; i<nCUT; i++) {
        //Generate psi0 cut
        complex_functionT psi0 = FunctionFactory<complexd,NDIM>(world).
            functor(functorT(new BoundWF(Z, 1, 0, 0)));
        psi0.scale(1/psi0.norm2());
        psi0.truncate();//Does this throw away small coefficents?
        psi0.scale(1/psi0.norm2());
        //Make potn
        param.cut = cut[i]; //The function V calls param.cut
        complex_functionT potn = complex_factoryT(world).f(V);
        potn.truncate(param.thresh);
        double eps = energy(world, psi0, potn);
        converge(world, potn, psi0, eps);
        output = inner(oneS,psi0);
        PRINT(std::setprecision(2) << "cut = " << cut[i] << "  " <<std::setprecision(12) 
              << real(conj(output)*output) << "\t");
        output = inner(twoS,psi0);
        PRINT(real(conj(output)*output) << "\t");
        output = inner(threeS,psi0);
        PRINT(real(conj(output)*output) << "\t");
        output = inner(fourS,psi0);
        PRINTLINE(real(conj(output)*output) << "\t");
    }
}


/*************************************************
 * If you're curious about a wave function's value
 *************************************************/
void printBasis(World& world, double Z, double cutoff) {
    complexd output, output2;
    double sinTH, cosTH, sinPHI, cosPHI;
    std::string tag;
    double rMIN = 0.0;
    double rMAX = 10.0;
    double dr = 1.0;
    double TH = 0.0;
    double PHI = 0.0;
    double k = 1.0;
    /***************************************
     *Load graphing parameters from the file: param
     * rMIN 0.0
     * rMAX 10.0
     * dr   1.0
     * TH   0.0
     ****************************************/
    std::ifstream f("param");
    std::cout << std::fixed;
    std::cout.precision(1);
    if( f.is_open() ) {
        while(f >> tag) {
            if (tag[0] == '#') {
                char ch;
                PRINTLINE("    comment  " << tag.c_str());
                while (f.get(ch)) {
                    PRINTLINE(ch);
                    if (ch == '\n') break;
                }
            }
            else if (tag == "rMIN") {
                f >> rMIN;
                PRINTLINE("rMIN   = " << rMIN);
            }
            else if (tag == "rMAX") {
                f >> rMAX;
                PRINTLINE("rMAX   = " << rMAX);
            }
            else if (tag == "dr") {
                f >> dr;
                PRINTLINE("dr     = " << dr);
            }
            else if (tag == "TH") {
                f >> TH;
                PRINTLINE("TH     = " << TH);
            }
            else if (tag == "k") {
                f >> k;
                PRINTLINE("k      = " << k);
            }
            else if (tag == "cutoff") {
                f >> cutoff;
                PRINTLINE("cutoff = " << cutoff);
            }
        }
    }
    //make functions
    Vector<double,3> kvec = vec(0.0, 0.0, k);
    PhiK phi_k(Z, kvec, cutoff);
    //for(double TH=0; TH<3.14; TH+=0.3 ) {
    //    for(double r=0; r<sqrt(3)*phi_k.domain*phi_k.k; r+=1.0 ) {
    PRINT("r \tRe:phi(rVec) \t Im:phi(rVec) \t");
    PRINTLINE("Re:diff(rVec) \t Im:diff(rVec) \t");
    for(double r=rMIN; r<rMAX; r+=dr ) {
        std::cout.precision(3);
        std::cout << std::fixed;
        cosTH =  std::cos(TH);
        sinTH =  std::sin(TH);
        cosPHI = std::cos(PHI);
        sinPHI = std::sin(PHI);
        Vector<double,3> rvec = vec(r*sinTH*cosPHI, r*sinTH*sinPHI, r*cosTH);
        output = phi_k(rvec);
        PRINT(r);
        std::cout.precision(7);
        PRINT("\t" << real(output) << "\t" << imag(output));
//         std::cout << std::scientific;
//         PRINT("\t" << phi_k.diffR(2*k*r) << "\t" << phi_k.diffI(2*k*r));
//      PRINT("\t" << real(phi_k.f11(2*k*r)) << "\t" << imag(phi_k.f11(2*k*r)));
//      PRINT("\t" << real(conhyp(-I/k,ONE,-I*k*r)) << "\t" << imag(conhyp(-I/k,ONE,-I*k*r)));
        PRINTLINE(" ");
    }
    //    use sed to make the complexd output standard
    //    system("sed -i '' -e's/\\+/, /' -e's/j//' f11.out");
}

/****************************************************************************
 * Reproduces atomic form integrals given by Dz Belkic's analytic expression
 ****************************************************************************/
void belkic(World& world, double cutoff) {
    /************************************************************************
     * qVec is the momentum transfered from the laser field
     * kVec is the momentum of the ejected electron
     * I will take these to be colinear
     ************************************************************************/
    PRINTLINE("1 0 0");
    double Z = 1.0;
    complex_functionT b1s = complex_factoryT(world).functor(functorT( 
                                                      new BoundWF(Z, 1, 0, 0) ));
    double dARR[3] = {0, 0, 0.5};
    const vector3D kVec(dARR);
    PRINTLINE("|" << kVec << ">");
    PhiK phik = PhiK(Z, kVec, cutoff);
    phik.Init(world);
    complex_functionT phi_k = complex_factoryT(world).functor(functorT( &phik ));
    dARR[2] =  1.5; // {0, 0, 1.5}
    const vector3D qVec(dARR);
    PRINTLINE("Exp[I" << qVec << ".r>");
    complex_functionT expikDOTr = complex_factoryT(world).functor(functorT(
                                                           new Expikr(qVec) ));
    PRINTLINE("<k=0.5| Exp[iqVec.r] |100>");
    complexd output = inner(phi_k, expikDOTr*b1s);
    PRINTLINE(output);
}
