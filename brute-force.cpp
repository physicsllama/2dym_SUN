#include <vector>
#include <map>
#include "2dym.h"
#include "2dym-analytic.h"

using namespace std;

int main(int argc, char *argv[]) {

  //set parameters
  double A_max = 20.0;
  double A_min = 9.0;  //don't make lower than pi^2
  double A_step = 1.0;
  double A_ref = 20.0;
  double chi = 2.0;
  int nF_max = 20;
  int nF_min = 2;
  int nF_step = 1;
  int max_terms = 8000000;

  //calculate S_shannon using the brute force method
  bool plotSshannon = false;

  //calculate S_boltzmann using the brute force method and subtract off its value
  //calculated from 
  bool plotSboltz = true;

  if(plotSshannon) {
    for(double A = A_max; A > A_min; A -= A_step) {
      stringstream strm;
      strm << "./exp/S_brute/Sbrute_A_" << fixed << setprecision(2) << A << ".dat";
      ofstream output(strm.str());
      for(int n = nF_min; n <= nF_max; n += nF_step) {
        output << 2*n + 1 << " " << S_shan_an(A, chi, n, 0.001, brute_logZ) - S_shan_an(A_ref, chi, n, 0.001, brute_logZ) << endl;
      }
      output.close();
    }
  }
  
  if(plotSboltz) {
    for(double A = A_max; A > A_min; A -= A_step) {
      stringstream strm, strm2;
      strm << "./exp/S_boltz_lin/Sboltz_A_" << fixed << setprecision(2) << A << ".dat";
      strm2 << "./exp/S_boltz_quad/Sboltz_A_" << fixed << setprecision(2) << A << ".dat";
      ofstream output(strm.str());
      // ofstream output2(strm2.str());
      for(int n = nF_min; n <= nF_max; n += nF_step) {

        //calculate S_boltzmann with the brute force partition function /except/ for the saddle point term
        output << 2*n + 1 << " " << (S_boltz_5(A, chi, n, 0.001, justbrute_logZ) - S_boltz_5(A_ref, chi, n, 0.001, justbrute_logZ)) << endl;

        //calulate S_boltzmann with entire partition function
        // output2 << 2*n + 1 << " " << (S_boltz(A, chi, n, 0.001, brute_logZ) - S_boltz(A_ref, chi, n, 0.001, brute_logZ)) << endl;
      }
      output.close();
      // output2.close();
    }
  }
  
}

double brute_force_Z(double A, int nF, double chi, int max_terms) {
  //get initial config
  int N = 2*nF+1;
  diagram d = find_saddlepoint(A, nF, chi, 100000);
  perform_iterations(&d, A, nF, chi, 1.0, 10000, gradientdescent);

  //create priority queue
  vector<pair<diagram, double>> queue;
  queue.push_back(pair<diagram, double>(d,d.C2*A/double(2*N) - chi*d.logdim));
  make_heap(queue.begin(), queue.end(), greater1());

  //create map to record configurations previously seen
  map<int, vector<vector<int>>> prev;
  vector<vector<int>> vec;
  vec.push_back(d.h);
  prev[hash_func(d.h)] = vec;

  double Z = exp(-d.C2*A/double(2*N) + chi*d.logdim);
  double Zprev = Z;
  cerr << "A = " << A << ", nF = " << nF << ", chi = " << chi << " done!" << endl;
  cerr << "    Initial partition function : " << Z << endl;

  //every time the number of configurations added to the partition
  //function doubles, the while loop checks how much the partition function
  //has changed; if less than "delta", the while loop breaks
  long counter = 0;
  long base = 2;
  int f = 13;
  double delta = 0.01;
  pair<diagram,double> next;
  diagram temp;
  double newE;
  int hash_val;

  //file to track building up of partition function
  ofstream output("./E_brute.dat");
  ofstream output2("./Z_brute.dat");

  //storage for the Boltzmann factors
  vector<double> factors;

  int collisions;
  while((counter != pow(base, f) || (Z - Zprev)/Zprev > delta) && !queue.empty() && counter < max_terms) {

    //update the exponent in 2^f so that we run for double the number of configs
    //we have so far if necessary
    if( counter == pow(base,f)) {
      Zprev = Z;
      f++;
    }

    //get the lowest energy config
    next = queue.front();
    d = next.first;
    pop_heap(queue.begin(), queue.end(), greater1());
    queue.pop_back();

    //generate neighbors of lowest energy config by iterating over edges
    int dh = 0;
    int j = 0;
    for(int e0 = 0; e0 <= d.edges -1; e0++) {

      //find e0
      while(j < N) {
        // Check if this is a "rising edge"
        if(e0 % 2 == 0 && ((e0 == 0 && j == 0) || d.h[j] > d.h[j-1]+1)) {
          dh = -1;
          break;
        }
        // Check if this is a "falling edge"
        if(e0 % 2 != 0 && ((e0 == d.edges - 1 && j == N-1) || d.h[j+1] > d.h[j] + 1)) {
          dh = 1;
          break;
        }
        j++;
      }

      //create new config
      temp = d;
      temp.h[j] += dh;
      hash_val = hash_func(temp.h);

      //check if this config has been seen before
      //and insert if not
      if(prev.count(hash_val) != 0) {
        //collect the vector of configs stored at that hash
        vec = prev.at(hash_val);

        //see if any of those configs are the one we just generated
        bool match;
        for(vector<int> h : vec) {
          match = configs_match(h,temp.h);
          if(match) {
            //stop checking if we find a match
            break;
          } else {
            collisions++;
          }
        }

        //if we found a match, move on to next neighbor, i.e. next iteration of for loop
        if(match) {
          continue;
        }
      }

      //add the new configuration into the map
      vec.push_back(temp.h);
      prev[hash_val] = vec;

      //calculate dC2 and dlogdim
      double dC2 = 2*d.h[j]*dh + 1;
      double dlogdim = 0;
      for(int k = 0; k < N; k++) {
        if(j == k) continue;
        dlogdim += log(1 + dh / double(d.h[j] - d.h[k]));
      }

      // Find the number of edges of the new diagram
      int new_edges = d.edges;
      if(d.h[j+dh] == d.h[j]+2*dh) new_edges -= 2;
      if(d.h[j-dh] == d.h[j]-dh) new_edges += 2;

      //finalize the calculation of temp and its energy
      temp.edges = new_edges;
      temp.C2 = d.C2 + dC2;
      temp.logdim = d.logdim + dlogdim;
      newE = A*temp.C2/double(2*N) - chi*temp.logdim;

      //add to the priority queue
      queue.push_back(pair<diagram,double>(temp, newE));
      push_heap(queue.begin(), queue.end(), greater1());

      //move on to the next fermion if necessary
      if(e0 % 2 == 1) {
        j++;
      }
    }
    // add next to partition function
    if(counter < 100) {
      output << counter << " " << next.second << endl;
    }
    
    if(counter % 10 == 0 && counter < 2000) {
      output2 << counter << " " << Z << endl;
    }

    Z += exp(-next.second);
    factors.push_back(exp(-next.second));
    counter++;
  }
  output.close();
  output2.close();

  //calculate Shannon entropy
  double S_shan = 0;
  double p;
  for(double d : factors) {
    p = d/Z;
    S_shan -= (p)*log(p);
  }
  cerr << "    S_shan = " << S_shan << endl;

  //output statistics
  if(max_terms != 0 && counter == max_terms) {
    cerr << "    MAY NOT HAVE CONVERGED!" << endl;
  }
  cerr << "    Final partition function : " << Z << endl;
  cerr << "    Counter is " << counter << ", 2^f is " << pow(base, f) << endl;
  cerr << "    Change : " << (Z - Zprev)/Zprev << endl;
  cerr << "    Num collisions: " << collisions << endl;
  cerr << "    Final queue size : " << queue.size() << endl << endl;

  return Z;
}

bool configs_match(vector<int> h1, vector<int> h2) {
  for(int i = 0; i < h1.size(); i++) {
    if(h1[i] != h2[i]){
      return false;
    }
  }
  return true;
}

double S_shan_an(double A, double chi, int nF, double incr, double (*logZ)(double A0, double chi0, int nF0)) {
  double dlogZ_A = (logZ(A + incr, chi, nF) - logZ(A -incr, chi, nF))/double(2*incr);
  double dlogZ_chi = (logZ(A, chi + incr, nF) - logZ(A, chi - incr, nF))/double(2*incr);
  return logZ(A, chi, nF) - A*dlogZ_A - chi*dlogZ_chi;
}

double S_shan_an_5(double A, double chi, int nF, double incr, double (*logZ)(double A0, double chi0, int nF0)) {
  double dlogZ_A = (-logZ(A + 2*incr, chi, nF) + 8.0*logZ(A + incr, chi, nF) - 8.0*logZ(A - incr, chi, nF) + logZ(A - 2*incr, chi, nF))/double(12*incr);
  double dlogZ_chi = (-logZ(A, chi + 2*incr, nF) + 8.0*logZ(A, chi + incr, nF) - 8.0*logZ(A, chi - incr, nF) + logZ(A, chi - 2*incr, nF))/double(12*incr);
  return logZ(A, chi, nF) - A*dlogZ_A - chi*dlogZ_chi;
}

double S_boltz(double A, double chi, int nF, double incr, double (*logZ)(double A0, double chi0, int nF0)){
  double dlogZ_chi = (logZ(A, chi + incr, nF) - logZ(A, chi - incr, nF))/double(2*incr);
  return 2*dlogZ_chi;
}

double S_boltz_5(double A, double chi, int nF, double incr, double (*logZ)(double A0, double chi0, int nF0)) {
  double dlogZ_chi = (-logZ(A, chi + 2*incr, nF) + 8.0*logZ(A, chi + incr, nF) - 8.0*logZ(A, chi - incr, nF) + logZ(A, chi - 2*incr, nF))/double(12*incr);
  return 2*dlogZ_chi;
}

double justbrute_logZ(double A, double chi, int nF) {
  return log(brute_force_Z(A, nF, chi, 1100000)) - log(brute_force_Z(A, nF, chi, 0));
}

double brute_logZ(double A, double chi, int nF) {
  return log(brute_force_Z(A, nF, chi, 1100000));
}

double saddle_point_logZ(double A, double chi, int nF) {
  return log(brute_force_Z(A,nF,chi,1));
}
/*
double perturbed_logZ(double A, double chi, int nF) {
  //get diagram in correct configuration
  diagram d = initial_config(A / chi, nF);
  perform_iterations(&d, A, nF, chi, 1.0, 10000, gradientdescent);

  return saddle_point_logZ(A, chi, nF) - 0.5*log(find_determinant(A, chi, nF, d.h));
}

double justdet_logZ(double A, double chi, int nF) {
  //get diagram in correct configuration
  diagram d = initial_config(A / chi, nF);
  perform_iterations(&d, A, nF, chi, 1.0, 10000, gradientdescent);

  return - 0.5*log(find_determinant(A, chi, nF, d.h));
}

double metro_logZ(double A, double chi, int nF) {
  diagram d = initial_config(A / chi, nF);
  return log(S_shan_metro(&d, A, nF, chi));
}
*/