#include "2dym.h"
#include "2dym-analytic.h"

int main(int argc,char *argv[] ) {

  if(argc == 4) {
    double A = atof(argv[1]);
    double chi = atof(argv[2]);
    int nF = atof(argv[3]);

    stats s = run_simulation(A, nF, chi);
    plot_data_exp(s);

    params p = find_params(A/chi);
    plot_model(p, nF);
  } else {
    cerr << "fermion-config requires command line arguments!" << endl;
  }
  
}
