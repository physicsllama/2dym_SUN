#include "2dym.h"
#include "2dym-analytic.h"

using namespace std;

int main(int argc, char *argv[]) {

  //choose whether to plot code as a function of N or A
  bool plotAsFunctionOfN = true;

  //set parameters
  double A_max = 20.0;
  double A_min = 4.0;
  double Astep = 0.01;
  int nF_max = 10;
  int nF_min = 2;
  int nFstep = 1;
  double chi = 2.0;

  
  //containers to store the entropies over the range of Ns AND the range of As
  vector<vector<double>> S_shan;
  vector<vector<double>> S_tot;
  vector<vector<double>> S_boltz;
  vector<vector<double>> S_sub;

  //RUNS SIMULATIONS AND COMPUTES ENTROPIES
  vector<double> Ns;
  for(int n = nF_min; n <= nF_max; n += nFstep) {
    vector<stats> data = run_range_of_simulations_exp(A_max, A_min, Astep, n, chi);

    //compute shannon entropy via the variance method and subtraction method 
    //and adds the data from this value of nF to the data collected at all other values of nF
    S_shan.push_back(S_shannon_exp(data, Astep));
    S_tot.push_back(S_total_exp(data, Astep));
    S_boltz.push_back(S_boltz_exp(data));
    S_sub.push_back(S_subtract(S_total_exp(data, Astep), S_boltz_exp(data)));

    //keep track of the Ns
    Ns.push_back(2*n +1);
  }

  //collect range of As in a vector
  vector<double> As;
  for(double A = A_max; A > A_min; A -= Astep){
    As.push_back(A);
  }
  
  //PLOTS ENTROPIES
  if(plotAsFunctionOfN) {
    print_functionN("./exp/S_shan/Sshan_A_", S_shan, As, Ns);
    print_functionN("./exp/S_tot/Stot_A_", S_tot, As, Ns);
    print_functionN("./exp/S_boltz/Sboltz_A_", S_boltz, As, Ns);
    print_functionN("./exp/S_sub/Ssub_A_", S_sub, As, Ns);
  } else {
    for(int i = 0; i < S_shan.size(); i++) {
      stringstream strm;
      strm << fixed << setprecision(0) << Ns[i] << ".dat";
      print("./exp/S_shan/Sshan_N_" + strm.str(), As, S_shan[i]);
      print("./exp/S_tot/Stot_N_" + strm.str(), As, S_tot[i]);
      print("./exp/S_boltz/Sboltz_N_" + strm.str(), As, S_boltz[i]);
      print("./exp/S_sub/Ssub_N_" + strm.str(), As, S_sub[i]);
    }
  }

}

double dA_S_exp(stats s) {
  return -s.A*((1/double(4*s.N*s.N))*s.varC2) - s.chi*((-1/double(2*s.N))*s.covar);
}

vector<double> S_shannon_exp(vector<stats> data, double incr) {
  double total = 0;
  vector<double> S_shan;
  double dS;
  for(stats s : data) {
    S_shan.push_back(total);
    // Derivative of S_shannon w.r.t. A
    dS = dA_S_exp(s);
    //Numerical integration
    total -= dS * incr;
  }
  return S_shan;
}

vector<double> S_total_exp(vector<stats> data, double incr) {
  double totalF = 0;
  double totalS = 0;
  double dF1;
  double dF2;
  int N = data[0].N;
  double F3prime_max;
  double range;

  vector<double> S_tot;
  vector<double> error;

  for(int i=0; i < data.size() -1; i++){
    dF1 = Fprime_exp(data[i]);
    dF2 = Fprime_exp(data[i+1]);
    totalF -= (dF1 + dF2) * incr/2.0;
    totalS = (N*N)*(totalF - data[i].A*dF1 + (data.front().A)*Fprime_exp(data.front()));
    S_tot.push_back(totalS);

    //find error using analytic expression for F(A)
    F3prime_max = 1.0/(2.0 * data[i+1].A);
    range = data[0].A - data[i+1].A;
    error.push_back((F3prime_max*range*range*range*N*N)/double(12*(i+1)*(i+1)));
  }

  return S_tot;
}

vector<double> S_boltz_exp(vector<stats> data) {
  vector<double> S_boltz;
  double initial = data.front().chi * data.front().explogdim;
  for(stats s : data) {
    S_boltz.push_back(s.chi*s.explogdim - initial);
  }
  return S_boltz;
}

vector<double> S_subtract(vector<double> S_tot, vector<double> S_boltz) {
  vector<double> S_sub;
  for(int i = 0; i < S_tot.size(); i++) {
    S_sub.push_back(S_tot[i] - S_boltz[i]);
  }
  return S_sub;
}

void print(string filename, vector<double> x, vector<double> y) {
  ofstream output(filename);
  for(int i = 0; i < x.size(); i++) {
    output << x[i] << " " << y[i] << endl;
  }
  output.close();
}

void print_functionN(string partial_filename, vector<vector<double>> data, vector<double> As, vector<double> Ns) {

  double A;
  for(int i = 0; i < As.size(); i++) {
    
    A = As[i];
    if(i % 100 == 0) {
      vector<double> data_oneA;

      //collects all data that have this values of A (over the range of N)
      for(int j = 0; j < Ns.size(); j++) {
        data_oneA.push_back(data[j][i]);
      }

      //prints the data
      stringstream strm;
      strm << partial_filename << fixed << setprecision(2) << A << ".dat";  
      string filename = strm.str();
      print(filename, Ns, data_oneA);
    }
    
  }

}