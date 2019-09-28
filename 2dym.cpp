#include <vector>
#include <random>
#include <functional>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <assert.h>
#include <algorithm>
#include <set>
#include <map>
#include "2dym.h"
#include "2dym-analytic.h"

using namespace std;

/*
  Initialize random number generators.
  Later on we can call uniform_dist() and it will return a uniform random number between 0 and 1.
*/
  default_random_engine generator;
  uniform_real_distribution<double> uniform_dist(0.0,1.0);
  auto random_real = bind(uniform_dist, generator);

/*
  To mutate we generate a random number r between 1 and the number of edges.
  We count from the beginning of the diagram to the rth edge - depending on
  if it's a rising or falling edge, we either move a fermion down one level or
  up one level, h'_i = h_i + dh, where dh = +- 1.
  Update C2 (which just means adding h'_i^2 - h_i^2)
  Update log dim R (which involves looping over the other fermions)
  Update the number of edges:
    It goes up by two if h_i-dh is filled.
    It goes down by two if h'_i+dh is filled.
  Accept or reject the new configuration based on exp(-C_2 + 2 log dim R)
  but also the ratio of number of edges in the old configuration versus the
  new configuration.

  Now we can collect statistics.
  - The first thing we can do is collect a count of
  how many times each site was occupied - this should give us a semicircle.
  - We can determine the expectation value of C2 and log dim R.
  - Once we have <C_2> = E(T) as a function of T, we can determine S(T)
  by integrating the first law, dS = dE/T. Given a set of samples E_i = E(T_i),
    S(T_i) = sum_{0 < j < i} (E(T_j) - E(T_{j-1})) / ((T_j+T{j-1})/2)
  where T_0 = 0 and S(T_0) = 0.

  What else can we determine?
  - Density of states as a function of energy?
*****/
#if 0
int main(int argc,char *argv[] ) {

  /*
    Introduce the constant a, which is given by a = g^2 A N / 2 = lambda A/2. Command line arguments
    can be used to set the value of a and chi, and compute Sshannon for these values (varying them together)
    over a range of values of N.
  */
  if( argc > 2 ) {
    //first command line argument
    double A = atof(argv[1]);
    //second
    double chi = atof(argv[2]);
    // third
    int nF = atof(argv[3]);
    // //fourth
    int num_its = atof(argv[4]);

    // ofstream output("./F_brute.dat");
    // double initial = log(brute_force_Z(10.0, nF, chi, num_its));
    // int N = 2*nF +1;
    // for(double a = 10.0; a >= 6.6; a -= 0.2) {
    //   output << a << " " << (1.0/double(N*N))*(log(brute_force_Z(a, nF, chi, num_its)) - initial)  << endl;
    // }
    // output.close();

    // for(double A = 19.0; A >= 7.0; A -= 1.0) {
    //   ofstream output1("./exp/S_shan_metro/Sshan_A_" + to_string(A) + ".dat");
    //   double incr = 0.01;
    //   for(int n = 2; n <= 10; n++) {
    //     cerr << "A: " << A << ", nf: " << n << endl;
    //     output1 << 2*n +1 << " " << S_shan_an(A, chi, n, incr,metro_logZ) - S_shan_an(20.0, chi, n, incr,metro_logZ) << endl;
    //   }
    //   output1.close();
    // }

    cerr << "Partition function: " << Zf << ", saddle : " << Z << endl;
    cerr << log(Zf) - log(Z) << endl;

    // //calls the main function from the analytic code with a acting as a 'command line argument'.  This
    // //finds the parameters and plots the model for this value of a.
    // main_analytic(2, A, chi);
  } else {

    double A_start = 10.0;
    double A_end = 3.9;
    double incr = 0.01;
    double chi = 2.0;

    //Computes various entropies and F'(A) for a range of ~coupling constant a over a range of N
    for(int nF = 2; nF <= 14; nF++){

      // runs simulations
      vector<stats> data = run_range_of_simulations_exp(A_start, A_end, incr, nF);

      //computes entropies
      S_total_exp_rec(data, incr);
      S_total_exp(data, incr);
      S_total_exp_under(data, incr);
      S_boltz_exp(data);
      S_total(A_start, A_end, incr, 2*nF +1);
      plot_F(A_start, A_end, incr, 2*nF +1);
      F_total_exp(data, incr);

    }
  }
}
#endif

stats run_simulation(double A, int nF, double chi, double tau, int num_its) {
  /* Initialize the diagram to the trivial irrep
     The trivial representation has h_1 = -n_F, ..., h_N = n_F, where N = 2 n_F + 1.
     C_2 = 0, log dim R = 0, and two edges.
  */
  diagram d = find_saddlepoint(A, nF, chi, num_its);

  return perform_iterations(&d, A, nF, chi, tau, num_its, metropolis);
}

//needs update to correct use of hash function
stats perform_iterations(diagram* d, double A, int nF, double chi, double tau, int iterations, bool (*compare)(double rand, double A, double chi, int N, double dC2, double dlogdim, int edges, int new_edges)) {

  /* Initially take N constant, but we may want to do some extrapolation later */
  int N = 2*nF + 1;

  int burnin = 100000;

  /*Variables we will use to collect statistics.
    Extend by 3*nF in either direction - this is somewhat arbitrary. */
  vector<double> counts(6*nF+1, 0);
  double totalC2 = d->C2;
  double totallogdim = d->logdim;
  double totalC2_sq = d->C2*d->C2;
  double totalC2_cub = d->C2*d->C2*d->C2;
  double totalC2logdim = d->C2*d->logdim;
  double totalE = -chi*tau*d->logdim + A*tau*d->C2/double(2*N);
  double totalE_sq = totalE*totalE;
  int newh;
  int accepts = 0;
  // ofstream output("./exp/E/E_A_" + to_string(A) + ".dat");

  //create a set to track unique configs visited
  set<double> prev;
  int hash_val;
  prev.insert(hash_func(d->h));

  //create a heap to track min energy configs visited
  // vector<double> queue;
  // queue.push_back(d->C2*A/double(2*N) - chi*d->logdim);
  // make_heap(queue.begin(), queue.end(), greater2());

  // Number of configurations we want to sample 
  int real_iterations = 0;

  while(real_iterations < 750000) {

    // We want to finish the burnin period at some point
    if(burnin>=0){
      burnin -= 1;
    }
    

    // output << i << " " << -chi*tau*d->logdim + A*tau*d->C2/double(2*N) << endl;

    /* Choose a random edge of the distribution to perturb, e_0. */
    int e0 = uniform_int_distribution<int>(0,d->edges-1)(generator);

    /* Iterate over the diagram to find the e_0'th edge.
       Perturbing this edge will correspond to changing h[j] -> h[j] + dh. */
    int e = 0;
    int dh = 0;
    int j = 0;
    for(j = 0; j < N; j++) {
      // Check if this is a "rising edge"
      if(j == 0 || d->h[j] > d->h[j-1]+1) {
        if(e == e0) {
          dh = -1;
          break;
        }
        e++;
      }
      // Check if this is a "falling edge"
      if(j == N-1 || d->h[j+1] > d->h[j] + 1) {
        if(e == e0) {
          dh = 1;
          break;
        }
        e++;
      }
    }

    // We don't want to move the rightmost fermion
    if(j==N-1) continue;

    assert(d->h[N-1] > d->h[N-2]);


    /* Now we want to consider mutating d->h[j] to d->h[j]+dh.
       To decide whether to keep the new configuration we have to find the difference in effective action.
       This means finding how much C_2(R) and log dim(R) change */

    // First we will define a new vector which shifts the fermions to have the rightmost one as reference.

    // vector<int> H;

    // for(int k=0; k<N; k++){
    //   int new_h = d->h[k] - d->h[N-1] + (N-1) / 2;
    //   H.push_back(new_h);
    // }

    // Make sure rightmost fermion is where it should be.
    assert(d->h[N-1] == (N-1)/2);
    

    // Change in number of boxes dbox always has a 1.
    double dbox = 1;
    double dlogdim = 0;
    for(int k = 0; k < N; k++) {
      dbox += 2*d->h[k]*dh;
      if(j == k) continue;
      dlogdim += log(1 + dh / double(d->h[j] - d->h[k]));
    }

    //Include change in number of boxes in addition to typical U(N) case.
    double dC2 = 2*d->h[j]*dh + 1 - dbox / double(N);

    // Find the number of edges of the new diagram
    int new_edges = d->edges;
    if(d->h[j+dh] == d->h[j]+2*dh) new_edges -= 2;
    if(d->h[j-dh] == d->h[j]-dh) new_edges += 2;

    /* Find the probability p to jump to the other configuration (p > 1 means we definitely do it)
       We have to correct for the change in the number of edges - configurations
       with fewer edges are more likely to stay put.
    */
    // double p = exp(-A * dC2/double(2 * N) + chi*dlogdim) * (d->edges / double(new_edges));
    // double p = exp(-A * dC2/double(2 * N) + chi*dlogdim);

    if(compare(random_real(), A, chi, N, dC2, dlogdim, d->edges, new_edges)) {
      //Switch to new state
      newh = d->h[j] + dh;
      // for(int i : d->h) {
      //   assert(i != newh);
      // }
      d->h[j] += dh;
      d->C2 += dC2;
      d->logdim += dlogdim;
      d->edges = new_edges;
      accepts++;

      // hash_val = hash_func(d->h);
      // if(prev.count(hash_val) == 0) {
        //prev.insert(hash_val);

        // queue.push_back(d->C2*A/double(2*N) - chi*d->logdim);
        // push_heap(queue.begin(), queue.end(), greater2());
    }

    // Ensure rightmost fermion still hasn't moved.
    assert(d->h[N-1] == (N-1)/2);

    /* Once we are through the initial burn-in period, collect statistics.*/
    if(burnin <= 0) {
      for(int j = 0; j < N; j++) {
        if(d->h[j] >= -3*nF && d->h[j] <= 3*nF) {
          counts[d->h[j]+3*nF]++;
        }
      }
      totalC2 += d->C2;
      totallogdim += d->logdim;
      totalC2_sq += d->C2*d->C2;
      totalC2_cub += d->C2*d->C2*d->C2;
      totalC2logdim += d->logdim*d->C2;
      totalE += -chi*tau*d->logdim + A*tau*d->C2/double(2*N);
      totalE_sq += (-chi*tau*d->logdim + A*tau*d->C2/double(2*N)) * (-chi*tau*d->logdim + A*tau*d->C2/double(2*N));
      real_iterations += 1;
    }
  }

  // Normalize
  for(int j = 0; j < 6*nF+1; j++) {
    counts[j] = counts[j]/(double(real_iterations));
    // output2 << double(j - 3*nF)/double(N) << " " << counts_final[j] << endl;
  }
  // output2.close();

  //Calculate expectation values and variances
  double len = (double(real_iterations));
  double expC2 = totalC2/len;
  double expC2_sq = totalC2_sq/len;
  double expC2_cub = totalC2_cub/len;
  double explogdim = totallogdim/len;
  double expE = totalE/len;
  double varC2 = totalC2_sq/len - expC2*expC2;
  double covar = totalC2logdim/len - expC2*explogdim;
  double varE = totalE_sq/len - expE*expE;

  stats s = {counts, A, chi, N, expC2, expC2_sq, expC2_cub, explogdim, expE, varC2, covar, varE};

  // cerr << "unique configs: "<< prev.size() << endl;

  int counter = 0;
  double Z = 0;
  double E = 0;
  ofstream output("E_metro.dat");
  ofstream output2("Z_metro.dat");
  // while(!queue.empty()) {
  //   E = queue.front();
  //   output << counter << " " << E << endl;
  //   Z += exp(-E);
  //   output2 << counter << " " << Z << endl;

  //   pop_heap(queue.begin(), queue.end(), greater2());
  //   queue.pop_back();

  //   counter++;
  // }
  cerr << "VarC2 = " << varC2 << endl;
  cerr << "covar = " << covar << endl;
  // cerr << "Covar = " << covar << endl;
  // cerr << "len: " << real_iterations << endl;
  // cerr << "burnin: " << burnin << endl;
  output.close();
  output2.close();

  return s;

}

int hash_func(vector<int> h) {
  unsigned int b = 378551;
  unsigned int a = 63689;
  int tot = 0;
  unsigned int k;
  for(k = 0; k < h.size(); k++) {
    tot = (tot * a) + h[k];
    a = a * b;
  }
  return tot;
}

double Fprime_exp(stats s) {
  double Fprime = (-1.0/double(2.0*s.N*s.N*s.N))*s.expC2;
  // cerr << "Fprime = " << Fprime << endl;
  // cerr << "C2 = " << s.expC2 << endl;
  return Fprime;
}

void plot_data_exp(stats s) {
  stringstream strm;
  strm << "./exp/fermions/A" << fixed << setprecision(2) << s.A << "_" << "N" << s.N << ".dat";  string filename = strm.str();
  ofstream output(filename, ofstream::out);

  /* Dump average counts to console */
  for(int j = 0; j < 3*s.N -2; j++) {
    output << (j - 1.5*(s.N-1))/s.N << ' ' << s.counts[j] << endl;
  }

  output.close();
}

vector<stats> run_range_of_simulations_exp(double A_start, double A_end, double A_incr, int nF, double chi, int num_its) {
  vector<stats> s;
  int N = 2*nF +1;
  int counter = 0;
  for(double A = A_start; A > A_end; A -= A_incr) {
    s.push_back(run_simulation(A, nF, chi, 1, num_its));
    cerr << "Simulation with A = " << A << ", N = " << N << " done!" << endl;
    // if(counter>=1){
    //   double numb = (s[counter].expC2 - s[counter-1].expC2)/s[counter-1].expC2;
    //   if(numb > 100){
    //     assert(false);
    //   }
    // }
    }
  return s;
}

vector<stats> run_Nrange_of_simulations_exp(double nF_start, double nF_end, double nF_incr, double A, double chi) {
  vector<stats> s;
  for(double nF = nF_start; nF < nF_end; nF += nF_incr) {
    s.push_back(run_simulation(A, nF, chi));
    plot_data_exp(s.back());
    cerr << "Simulation with A = " << A << ", N = " << 2*nF +1 << " done!" << endl;
  }
  return s;
}

diagram initial_config(double ratio, int nF) {
  //Find parameters for analytic model rho
  params p = find_params(ratio);
  diagram d;

  // Initialize stuff
  double I = 0;
  int h = -3*nF;
  int N = 2*nF + 1;

  //Numerically integrate rho
  double total = 0;
  for(int i = -3*nF; i <= 3*nF; i++) {
    total += rho(p.a, p.b, double(i) / double(N)) / double(N);
  }

  //We integrate the model until the value of the integral
  //has reached that i/N
  I = -0.5*(total/double(N));
  int n = 0; //number of fermions
  for(int i = -3*nF; i <= 3*nF; i++) {
    I += rho(p.a, p.b, double(i) / double(N))/double(N);
    if(I > (total*n)/double(N)) {
      d.h.push_back(i);
      n++;
    }
  }

  // Check that we have placed the correct number of fermions
  assert(d.h.size() == N);
  // And satisfy Pauli exclusion principle
  for(int i = 0; i < d.h.size(); i++) {
    for(int j = 0; j < i; j++) {
      assert(d.h[i] != d.h[j]);
    }
  }

  // Define vector that shifts fermions so that rightmost one is in correct spot to ensure n_N=0.
  vector<int> H;
  for(int j=0; j<N; j++){
    int new_h = d.h[j] - d.h[N-1] + (N-1) / 2;
    H.push_back(new_h);
  }

  d.h = H;

  //Calls diagram_props to calculate the initial number of edges, and values of C2 and logdim
  return diagram_props(d, nF);
}

diagram diagram_props(diagram d, int nF) {

  int N = 2*nF +1;

  //initialize C2 and logdim with their respective constants
  d.C2 = -(1.0/12.0)*N*(N*N - 1);
  d.logdim = 0;
  for(int j = 1; j < N; j++){
    for(int i = 1; i < j; i ++){
      d.logdim -= log(j - i);
    }
  }

  // Define fermions wrt rightmost fermion.
  // vector<int> H;

  // for(int j=0; j<N; j++){
  //   int new_h = d.h[j] - d.h[N-1] + (N-1) / 2;
  //   H.push_back(new_h);
  // }

  //Compute C2 and logdim from the heights h
  int h;
  int box = 0;
  for(int h : d.h) {
    box -= h;
    d.C2 += h*h;
    for(int h2: d.h){
      if(h2 < h){
        d.logdim += log(h - h2);
      } else{
        break;
      }
    }
  }
  d.C2 -= double(box * box) / double(N);

  d.edges = count_edges(d);

  return d;
}

int count_edges(diagram d) {
  //Compute the number of edges
  int h;
  int edges = 0;
  for(int i= 0; i < d.h.size(); i++) {
    h = d.h[i];
    //We get one edge each for the first and last fermion being able to move "outwards"
    if(i==0 || i==d.h.size()-1){
      edges++;
    }
    //We get an edge for each fermion able to increase its height
    if(i != d.h.size()-1 && d.h[i+1] != h+1) {
      edges++;
    }
    //And for each fermion able to decrease its height
    if(i != 0 && d.h[i-1] != h-1) {
      edges++;
    }
  }
  return edges;
}

diagram find_saddlepoint(double A, int nF, double chi, int iterations) {

  /* Initially take N constant, but we may want to do some extrapolation later */
  int N = 2*nF + 1;

  diagram d = initial_config(A / chi, nF);
  diagram minE_config = d;
  double minE = A * d.C2 / double(2*N) - chi*d.logdim;
  double newE;

  for(int i = 0; i < iterations; i++) {

    /* Choose a random edge of the distribution to perturb, e_0. */
    int e0 = uniform_int_distribution<int>(0,d.edges-1)(generator);

    /* Iterate over the diagram to find the e_0'th edge.
       Perturbing this edge will correspond to changing h[j] -> h[j] + dh. */
    int e = 0;
    int dh = 0;
    int j = 0;
    for(j = 0; j < N; j++) {
      // Check if this is a "rising edge"
      if(j == 0 || d.h[j] > d.h[j-1]+1) {
        if(e == e0) {
          dh = -1;
          break;
        }
        e++;
      }
      // Check if this is a "falling edge"
      if(j == N-1 || d.h[j+1] > d.h[j] + 1) {
        if(e == e0) {
          dh = 1;
          break;
        }
        e++;
      }
    }

    if(j==N-1) continue;

    /* Now we want to consider mutating d->h[j] to d->h[j]+dh.
       To decide whether to keep the new configuration we have to find the difference in effective action.
       This means finding how much C_2(R) and log dim(R) change */

    // The energies are computed with H, which shifts fermions to ensure n_N = 0.
    // vector<int> H;
    

    // for(int j=0; j<N; j++){
    //   int new_h = d.h[j] - d.h[N-1] + (N-1) / 2;
    //   H.push_back(new_h);
    // }

    double dbox = 1;
    double dlogdim = 0;
    for(int k = 0; k < N; k++) {
      dbox += 2*d.h[k]*dh;
      if(j == k) continue;
      dlogdim += log(1 + dh / double(d.h[j] - d.h[k]));
    }
    double dC2 = 2*d.h[j]*dh + 1 - dbox / double(N);

    assert(d.h[N-1] == (N-1)/2);

    // Find the number of edges of the new diagram
    int new_edges = d.edges;
    if(d.h[j+dh] == d.h[j]+2*dh) new_edges -= 2;
    if(d.h[j-dh] == d.h[j]-dh) new_edges += 2;

    /* Find the probability p to jump to the other configuration (p > 1 means we definitely do it)
       We have to correct for the change in the number of edges - configurations
       with fewer edges are more likely to stay put.
    */
    double p = exp(-A * dC2/double(2 * N) + chi*dlogdim) * (d.edges / double(new_edges));

    if(p > random_real()) {
      //Switch to new state
      d.h[j] += dh;
      d.C2 += dC2;
      d.logdim += dlogdim;
      d.edges = new_edges;

      //check if we have a new minimum energy config
      newE = A * d.C2 / double(2*N) - chi*d.logdim;
      if(newE < minE) {
        minE = newE;
        minE_config = d;
      }
      assert(minE_config.h[N-1]== (N-1)/2);
    }

  }

  return minE_config;
}

bool metropolis(double rand, double A, double chi, int N, double dC2, double dlogdim, int edges, int new_edges) {
  double p = exp(-A * dC2/double(2 * N) + chi*dlogdim) * (edges / double(new_edges));
  return p > rand;
}

bool gradientdescent(double rand, double A, double chi, int N, double dC2, double dlogdim, int edges, int new_edges) {
  double p = exp(-A * dC2/double(2 * N) + chi*dlogdim);
  return p > 1;
}
