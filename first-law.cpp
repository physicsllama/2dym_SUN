#include "2dym.h"
#include "2dym-analytic.h"

using namespace std;

int main(int argc, char *argv[]) {

  double A_max = 20.0;
  double A_min = 4.0;
  double A_step = 1.0;
  double chi = 2.0;
  int nF_max = 10;
  int nF_min = 2;
  int nF_step = 1;
  double initial = 0.01;

  for(double A = A_max; A > A_min; A -= A_step) {
    stringstream strm;
    strm << "./exp/S_first/Sfirst_A_" << fixed << setprecision(2) << A << ".dat";
    ofstream output(strm.str());
    for(int n = nF_min; n <= nF_max; n += nF_step) {
      output << 2*n + 1 << " " << S_shannon_E(A, chi, n, initial) - S_shannon_E(A_max, chi, n ,initial) << endl;
    }
    output.close();
  }
}

double S_shannon_E(double A0, double chi0, int nF, double initial) {
  double dTau = 0.005;
  int num_its = 500000;
  double total = 0;
  int N = 2*nF + 1;
  // ofstream output("./exp/S_shan_anneal/Sshan_E_N_" + to_string(2*nF+1) + "_A_" + to_string(A0) + ".dat");
  // ofstream output2("./exp/E_tau_anneal/E_tau_N_" + to_string(2*nF+1) + "_A_" + to_string(A0) + ".dat");

  diagram d = initial_config(A0 / chi0, nF);
  // Holds the energy of the fermions calculated at the current value of tau + dTau
  double E_tauP = perform_iterations(&d, A0, nF, chi0, 1, num_its, metropolis).expE;
  // " " calculated at tau
  double E_tau = perform_iterations(&d, A0/double(1 - dTau), nF, chi0/double(1 - dTau), 1.0 - dTau, num_its, metropolis).expE;
  // " " calculated at tau - dTau
  double E_tauM;

  for(double tau = 1-dTau; tau >= initial; tau-=dTau){
    //output result of integration
    // output << tau << " " << total << endl;
    // output2 << tau << " " << E_tau << endl;

    // Performs the integration step
    E_tauM = perform_iterations(&d, A0/double(tau - dTau), nF, chi0/double(tau - dTau), tau - dTau, num_its, metropolis).expE;
    total += (E_tauP - E_tauM)/(2.0*tau);

    //update E's for next step of loop
    E_tauP = E_tau;
    E_tau = E_tauM;
  }
  cerr << "S_shannon for nF = " << nF << " and A = " << A0 << " is " << total << "!" << endl;
  // output.close();

  return total;
}