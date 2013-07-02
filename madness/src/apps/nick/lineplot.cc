#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <constants.h>
#include <tensor/vmath.h>
#include <mra/lbdeux.h>
#include <complex>
#define PRINT(str) if(world.rank()==0) std::cout << str
#define PRINTLINE(str) if(world.rank()==0) std::cout << str << std::endl

using std::cout;
using std::cin;
using std::endl;

using namespace madness;
typedef std::complex<double> double_complex;
typedef Function<double_complex,3> complex_functionT;
typedef Vector<double,3> coordT;

const char* wave_function_filename(int step) {
    static char fname[1024];
    sprintf(fname, "%s-%5.5d", "data", step);
    return fname;
}

bool wave_function_exists(World& world, int step) {
    return ParallelInputArchive::exists(world, wave_function_filename(step));
}

complex_functionT wave_function_load(World& world, int step) {
    complex_functionT psi;
    ParallelInputArchive ar(world, wave_function_filename(step));
    ar & psi;
    return psi;
}

void line_plot(World& world, int step, complex_functionT& psi) {
    static const int npt = 10001;
    const double L = 300.0;
    double_complex v[10001];
    psi.reconstruct();
    for (int i=0; i<npt; i++)
        v[i] = 0.0;
    for (int i=world.rank(); i<npt; i+=world.size()) {
        double z = -L + 2.0*i*L/(npt-1);
        coordT r(0.0);
        r[2] = z;
        v[i] = psi.eval(r);
    }
    world.gop.fence();
    world.gop.sum(v, npt);
    if (world.rank() == 0) {
        char buf[256];
        sprintf(buf, "%s.lineplot", wave_function_filename(step));
        std::ofstream f(buf);
        f.precision(10);
        for (int i=0; i<npt; i++) {
            double z = -L + 2.0*i*L/(npt-1);
            f << z << " " << v[i] << "\n";
        }
    }
    world.gop.fence();
}

int main(int argc, char** argv) {
    initialize(argc,argv);
    World world(MPI::COMM_WORLD);
    startup(world,argc,argv);
    //Extract number from the first argument list
    try {
        int tMAX, tMIN, dt;
        PRINTLINE("Inside the try loop");
        switch (argc) {
        case 1:
            PRINTLINE("one");
            tMAX = atoi(argv[0]);
            tMIN = 0;
            dt = 10;
            PRINTLINE("tMAX = " << tMAX << "\t tMIN = " << tMIN << "\t dt = " << dt);
            break;
        case 2:
            PRINTLINE("two");
            tMAX = atoi(argv[1]);
            tMIN = atoi(argv[0]);
            dt = 10;
            PRINTLINE("tMAX = " << tMAX << "\t tMIN = " << tMIN << "\t dt = " << dt);
            break;
        case 3:
            PRINTLINE("three");
            tMAX = atoi(argv[2]);
            tMIN = atoi(argv[0]);
            dt   = atoi(argv[1]);
            PRINTLINE("tMAX = " << tMAX << "\t tMIN = " << tMIN << "\t dt = " << dt);
            break;
        default:
            PRINTLINE("   One argument: tMAX       DEFAULT tMIN=0 DEFAULT dt=10");
            PRINTLINE("  Two arguments: tMAX tMIN                 DEFAULT dt=10");
            PRINTLINE("Three arguments: tMAX tMIN dt");
            exit(1);
        }
        // Add MPI functionality ... maybe
        //        for (int i=world.rank(); i<npt; i+=world.size()) {
        PRINTLINE("Before for loop");
        for(int i=tMIN + dt*world.rank(); i<=tMAX; i+=world.size()) {
            if(wave_function_exists(world,i)) {
                complex_functionT psi = wave_function_load(world,i);
                line_plot(world, i, psi);
            }
        }
    }
    catch (const MPI::Exception& e) {
        error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
        print(e); std::cout.flush();
        error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
        print(e); std::cout.flush();
        error("caught a Tensor exception");
    } catch (const char* s) {
        print(s); std::cout.flush();
        error("caught a c-string exception");
    } catch (char* s) {
        print(s); std::cout.flush();
        error("caught a c-string exception");
    } catch (const std::string& s) {
        print(s); std::cout.flush();
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what()); std::cout.flush();
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }


    world.gop.fence();

    ThreadPool::end();
    print_stats(world);
    finalize();
    return 0;
}
