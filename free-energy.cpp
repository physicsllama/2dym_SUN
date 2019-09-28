#include "2dym-analytic.h"
#include "2dym.h"

int main(int argc, char *argv[]) {

  double chi = 2.0;
  int N = 81;

  //Range of A to plot the free energy over
  double A_max = pi*pi;
  double A_min = 2.0;
  double Astep = 0.01;

  //option to plot the free energy calculated from the saddle point
  //and the free energy from the saddle point + the determinant correction
  bool plot_saddle = false;

  //option to plot the experimental free energy
  bool plot_exp = true;

  //output files
  ofstream output("./exp/F/F_DK.dat");
  ofstream output2("./exp/F/F_DTV.dat");
  ofstream output3("./exp/F/F_GS.dat");
  ofstream output4("./exp/F/F_GM.dat");

  //perform the integration and plotting
  double totalF = 0;
  double dF;
  for(double A = A_max; A>=A_min; A-= Astep) {
    output << A << " " << totalF << endl;
    // Derivative of F w.r.t. A
    dF = Fprime_DK(A);
    //Numerical integration of F_DK
    // output << A << " " << totalF << endl;
    totalF -= dF * Astep;
    //Output of F_DTV
    output2 << A << " " << F_DTV(A, chi) - F_DTV(A_max, chi) << endl;
    //Output of F_GS
    output3 << A << " " << F_GS(A, chi) - F_GS(A_max, chi) << endl;
    //Output of F_GM
    output4 << A << " " << F_GM(A) - F_GM(A_max) << endl;
    
  }

  //plot the free energy calculated from the saddle point (if that option is on)
  if(plot_saddle) {
    //output data files
    ofstream output5("../exp/F/F_saddle.dat");
    ofstream output6("./exp/F/F_det.dat");

    //collect reference points to subtract off
    pair<double,double> F = F_corrected(A_max, chi, (N-1)/2);
    double F_saddle_ref = F.first;
    double F_det_ref = F.second;
    
    //fill data files
    for(double A = A_max; A>=A_min; A-= Astep) {
      //output of F_corrected
      F = F_corrected(A, chi, (N-1)/2);
      output5 << A << " " << F.first - F_saddle_ref << endl;
      output6 << A << " " << F.second - F_det_ref << endl;
    }

    output5.close();
    output6.close();
  }

  //plots the experimental free energy calculated from the saddle point (if that option is on)
  if(plot_exp) {
    vector<stats> s = run_range_of_simulations_exp(A_max, A_min, Astep, (N-1)/2);
    F_total_exp(s, Astep);
  }

  output.close();
  output2.close();
  output3.close();
  output4.close();
}

pair<double,double> F_determinant(double A, double chi, int nF){
  int N = 2*nF + 1;
  //get diagram in correct configuration
  diagram d = find_saddlepoint(A, nF, chi, 100000);
  perform_iterations(&d, A, nF, chi, 1.0, 10000, gradientdescent);

  double E = A * d.C2 / double(2*N) - chi*d.logdim;
  double det = find_determinant(A, chi, nF, d.h);
  return pair<double,double>((- E)/double(N*N), (- E - 0.5*log(det))/double(N*N));
  //- (1.0/24.0)*A + 0.4
}

double F_DTV(double A, double chi) {
  return -chi/4.0 *log(2*A/chi) + (A/24.0) + 3.0*chi/8.0;
}

double F_GS(double A, double chi) {
  return +chi/2.0 * log(chi/(2.0*A)) + chi + A/24.0;
}

double F_GM(double A) {
  return -A/24.0 - 0.5*log(A);
}

void F_total_exp(vector<stats> data, double incr) {
  double totalF = 0;
  double dF1;
  double dF2;
  int N = data[0].N;
  string filename = "./exp/F/F_exp.dat";
  ofstream output(filename, ofstream::out);
  for(int i=0; i < data.size() -1; i++) {
    //printing the value of the integral
    output << data[i].A << " " << totalF << endl;
    // Derivative of F w.r.t. A
    dF1 = Fprime_exp(data[i]);
    dF2 = Fprime_exp(data[i+1]);
    //Numerical integration of F
    totalF -= ((dF1 + dF2)/2.0)*incr;

  }
  output.close();
}