#ifndef EXTRA_H
#define EXTRA_H
#include "wavef.h"

struct InputParameters {
  static const int MAXNATOM=99;

    // IF YOU ADD A NEW PARAMETER DON'T FORGET TO INCLUDE IT IN
    // a) read()
    // b) serialize()
    // c) operator<<()
  
  double L;           // Box size for the simulation
  double Lsmall;      // Box size for small (near nucleus) plots
  double Llarge;      // Box size for large (far from nucleus) plots
  double F;           // Laser field strength
  double omega;       // Laser frequency
  double ncycle;      // Number of laser cycles in envelope
  int natom;          // Number of atoms
  double Z[MAXNATOM]; // Nuclear charge of atoms
  double R[MAXNATOM][3]; // Coordinates of atoms
  int k;              // wavelet order
  double thresh;      // precision for truncating wave function
  double safety;      // additional precision (thresh*safety) for operators and potential
  double cut;         // smoothing parameter for 1/r (same for all atoms for now)
  std::string iState ; // initial state = "1s" or "2s"
  std::string prefix; // Prefix for filenames
  int ndump;          // dump wave function to disk every ndump steps
  int nplot;          // dump opendx plot to disk every nplot steps
  int nprint;         // print stats every nprint steps
  int nloadbal;       // load balance every nloadbal steps
  int nio;            // Number of IO nodes 
  double tScale;      // Scaling parameter for optimization
  double target_time; // Target end-time for the simulation
  
  void read(const char* filename) {
    std::ifstream f(filename);
    std::string tag;
    iState = "1s";
    printf("\n");
    printf("       Simulation parameters\n");
    printf("       ---------------------\n");
    while(f >> tag) {
        if (tag[0] == '#') {
            char ch;
            printf("    comment  %s ",tag.c_str());
            while (f.get(ch)) {
                printf("%c",ch);
                if (ch == '\n') break;
            }
        }
        else if (tag == "L") {
            f >> L;
            printf("             L = %.1f\n", L);
        }
        else if (tag == "Lsmall") {
            f >> Lsmall;
            printf("        Lsmall = %.1f\n", Lsmall);
        }
        else if (tag == "Llarge") {
            f >> Llarge;
            printf("        Llarge = %.1f\n", Llarge);
        }
        else if (tag == "F") {
            f >> F;
            printf("             F = %.6f\n", F);
        }
        else if (tag == "omega") {
            f >> omega;
            printf("         omega = %.6f\n", omega);
        }
        else if (tag == "ncycle") {
            f >> ncycle;
            printf("         ncycle = %.6f\n", ncycle);
        }
        else if (tag == "natom") {
            f >> natom;
            printf("         natom = %d\n", natom);
            for (int i=0; i<natom; i++) {
                f >> Z[i] >> R[i][0] >> R[i][1] >> R[i][2];
                printf("           atom %2d   %.1f  %10.6f  %10.6f  %10.6f\n", i, Z[i], R[i][0], R[i][1], R[i][2]);
            }
        }
        else if (tag == "k") {
            f >> k;
            printf("             k = %d\n", k);
        }
        else if (tag == "thresh") {
            f >> thresh;
            printf("        thresh = %.1e\n", thresh);
        }
        else if (tag == "safety") {
            f >> safety;
            printf("        safety = %.1e\n", safety);
        }
        else if (tag == "cut") {
            f >> cut;
            printf("           cut = %.2f\n", cut);
        }
        else if (tag == "iState") {
            f >> iState;
            printf("        iState = %s\n", iState.c_str());
        }
        else if (tag == "prefix") {
            f >> prefix;
            printf("        prefix = %s\n", prefix.c_str());
        }
        else if (tag == "ndump") {
            f >> ndump;
            printf("         ndump = %d\n", ndump);
        }
        else if (tag == "nplot") {
            f >> nplot;
            printf("         nplot = %d\n", nplot);
        }
        else if (tag == "nprint") {
            f >> nprint;
            printf("         nprint = %d\n", nprint);
        }
        else if (tag == "nloadbal") {
            f >> nloadbal;
            printf("       nloadbal = %d\n", nloadbal);
        }
        else if (tag == "nio") {
            f >> nio;
            printf("            nio = %d\n", nio);
        }
        else if (tag == "target_time") {
            f >> target_time;
            printf("    target_time = %.3f\n", target_time);
        }
        else if (tag == "tScale") {
            f >> tScale;
            printf("         tScale = %.5f\n", tScale);
        }
        else {
            MADNESS_EXCEPTION("unknown input option", 0);
        }
    }
  }
    
  template <typename Archive>
  void serialize(Archive & ar) {
    ar & L & Lsmall & Llarge & F & omega & ncycle & natom & Z;
    ar & archive::wrap(&(R[0][0]), 3*MAXNATOM);
    ar & k & thresh & safety & cut & iState & prefix & ndump & nplot & nprint & nloadbal & nio;
    ar & target_time & tScale;
  }
};

std::ostream& operator<<(std::ostream& s, const InputParameters& p);

struct WF {
    std::string str;
    complex_functionT func;
    WF(const std::string& STR, const complex_functionT& FUNC) 
        : str(STR)
        , func(FUNC) 
    {
        func.truncate();
    }
};


const char* wave_function_filename(int step);
bool wave_function_exists(World& world, int step);
void wave_function_store(World& world, int step, const complex_functionT& psi);
complex_functionT wave_function_load(World& world, int step);

// This controls the distribution of data across the machine
class LevelPmap : public WorldDCPmapInterface< Key<3> > {
private:
    const int nproc;
public:
    LevelPmap() : nproc(0) {};
    
    LevelPmap(World& world) : nproc(world.nproc()) {}
    
    // Find the owner of a given key
    ProcessID owner(const Key<3>& key) const {
        Level n = key.level();
        if (n == 0) return 0;
        hashT hash;

        // This randomly hashes levels 0-2 and then
        // hashes nodes by their grand-parent key so as
        // to increase locality separately on each level.
        //if (n <= 2) hash = key.hash();
        //else hash = key.parent(2).hash();

        // This randomly hashes levels 0-3 and then 
        // maps nodes on even levels to the same
        // random node as their parent.
        // if (n <= 3 || (n&0x1)) hash = key.hash();
        // else hash = key.parent().hash();

        // This randomly hashes each key
        hash = key.hash();

        return hash%nproc;
    }
};


template<class T>
std::string toString( const T& a );
void loadDefaultBasis(World& world, std::vector<WF>& boundList, double Z);
void loadList(World& world, std::vector<std::string>& boundList, std::vector<std::string>& unboundList);
complexd zdipole( const vector3D& r);
void projectZdip(World& world, std::vector<WF> stateList);
void compare1F1(World& world, double cutoff);


complexd V(const vector3D& r);
double myreal(double t);
double myreal(const double_complex& t);
// Given psi and V evaluate the energy ... leaves psi compressed, potn reconstructed
template <typename T>
double energy(World& world, const Function<T,3>& psi, const Function<T,3>& potn);
void converge(World& world, complex_functionT& potn, complex_functionT& psi, double& eps);
void compareGroundState(World& world, double Z);
void printBasis(World& world, double Z, double cutoff = 10.0);
void belkic(World& world, double cutoff);
#endif
