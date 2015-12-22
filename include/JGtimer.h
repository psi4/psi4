#include<sys/time.h>
#include <psi4-dec.h>
#include <string>

using namespace boost;
using namespace std;
using namespace psi;

namespace psi {

class timJG
{
  private:
    struct timeval tv1;
    struct timeval tv2;
//    string message;
    double timespent = 0.0;
    char* message;
  public:
//    timJG(const string& inmsg) : message(inmsg) {}
    timJG() {}
    void start() {
      gettimeofday(&tv1, NULL);
    }
//    void stop() {
//      gettimeofday(&tv2, NULL);
//      timespent = (tv2.tv_sec - tv1.tv_sec) + 
//      (tv2.tv_usec - tv1.tv_usec) / 1000.0 / 1000.0;
//      outfile-><< message;
//      outfile->Printf(" took %f seconds.\n", timespent);
//    }
    void stop(const char* msg) {
      gettimeofday(&tv2, NULL);
      timespent = (tv2.tv_sec - tv1.tv_sec) + 
      (tv2.tv_usec - tv1.tv_usec) / 1000.0 / 1000.0;
      outfile->Printf(msg);
      outfile->Printf(" took %f seconds.\n", timespent);
      timespent = 0.0;
    }
    void cumulate() {
      gettimeofday(&tv2, NULL);
      timespent += (tv2.tv_sec - tv1.tv_sec) + 
      (tv2.tv_usec - tv1.tv_usec) / 1000.0 / 1000.0;
    }
    void print(const char* msg) {
      outfile->Printf(msg);
      outfile->Printf(" took %f seconds.\n", timespent);
      timespent = 0.0;
    }
};

}
