#include <vector>
#include "2dym.h"
#include "2dym-analytic.h"
#include <map>

using namespace std;

int main(int argc, char *argv[]) {
    if(argc == 4){
        // Define parameters
        double A = atof(argv[1]);
        double chi = 2;
        int nF_min = atof(argv[2]);
        int nF_max = atof(argv[3]);

        // Start loop over nF
        for(int nF=nF_min; nF <= nF_max; nF++){
            int N = 2*nF + 1;

            //Fetch saddlepoint diagrams. First one is basically just the first approximation, second is better
            diagram d = find_saddlepoint(A, nF, chi, 1);
            diagram d2 = find_saddlepoint(A, nF, chi, 1000000);

            //Get energies of the diagrams, and also the analytical energy prediction E0
            double E0 = chi/4.0 * log(2 * A / chi) - A / 24.0 - 3 * chi / 8.0;
            double E = (A * d.C2 / double(2*N) - chi*d.logdim) / double(N * N);
            double E2 = (A * d2.C2 / double(2*N) - chi*d2.logdim) / double(N * N);
            // if(d.h == d2.h){
            //     cerr << "At N = " + to_string(N) << " you picked the same configs!" << endl;
            // }

            cerr << "At N = " + to_string(N) << " the difference is " << E0 - E2 << endl;
            cerr << "the analytic energy is " << E0 << endl;
            cerr << "the numerical energy is " << E2 << endl;
            cerr << endl;
        }
    }

    else{
        return 0;
    }
}