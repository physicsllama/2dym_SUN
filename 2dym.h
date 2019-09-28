#ifndef EXPERIMENTAL_H
#define EXPERIMENTAL_H

#include <vector>
#include <random>
#include <functional>
#include <limits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

/*
  Calculate observables in 2D Yang-Mills theory using the Metropolis algorithm

  In the free fermion description, a configuration is a sequence of integers:
    h_1 < h_2 < .... < h_N
    C_2(R) = sum_i h_i^2
    log dim(R) = sum_{i < j} log(h_j-h_i)
  Where both expressions are up to a constant we subtract to make the trivial
  representation have C_2(trivial) = log dim(trivial) = 0.
  We will also keep track of the number of "edges", defined as:
    # edges = 2 + 2 * #{i : h_{i+1} != h_i + 1}
  The number of edges ranges between 2 and N. 
*/
struct diagram {
  vector<int> h;
  double C2, logdim;
  int edges;
};

/*
  Structure to hold all of the statistics we are collecting from the simulation
  for a run at given values of the parameters a and nF.  Includes the expectation
  values of C2 and logdim, as well as their variances.
*/
struct stats {
  vector<double> counts;
  double A, chi;
  int N;
  // Expectation value of C2(R)
  double expC2;
  double expC2_sq;
  double expC2_cub;
  // Expectation value of logdim(R)
  double explogdim;
  // Expectation value of the fermion energy
  double expE;
  // The variance of C2 = <C2(R)^2> - <C2(R)>^2
  double varC2;
  // The covariance of C2 and logdim = <C2(R) logdim(R)> - <C2(R)><logdim(R)>
  double covar;
  // The variance of chi*logdim(R) - A*(lambda C2(R) / 2 N), i.e. the fermion energy
  double varE;
};

struct greater1{
  bool operator()(const pair<diagram,double>& a,const pair<diagram,double>& b) const{
    return a.second > b.second;
  }
};

struct greater2{
  bool operator()(const double& a,const double& b) const{
    return a > b;
  }
};

/*
  Runs the Metropolis Monte Carlo simulation for a given value of the parameters a = g^2 A N / 2 and

  the number of fermions nF.  Returns the statistics collected from the simulation.
*/
stats run_simulation(double A, int nF, double chi = 2, double tau = 1.0, int num_its = 500000);

/* Performs the iteration steps for the simulation.  Called by run_simuation and S_shannon_E. */
stats perform_iterations(diagram* d, double A, int nF, double chi, double tau, int iterations, bool (*compare)(double rand, double A, double chi, int N, double dC2, double dlogdim, int edges, int new_edges));

/*Performs the metropolis algorithim, starting from initial_config, storing and returning the
lowest energy config touched.  Running this, followed by perform_iterations set to gradientdescent,
is the current best method for finding the saddle point config.*/
diagram find_saddlepoint(double A, int nF, double chi, int iterations);

/* Finds the partition function by adding configurations on in increasing energy order until 
adding twice as many terms doesn't increase the partition function by a certain fraction (or a 
maximum number of terms is reached.*/
double brute_force_Z(double A, int nF, double chi, int max_terms);

/*Brute force computation of Z that doesn't catch hash collisions.*/
double naive_hash_Z(double A, int nF, double chi, int max_terms);

double S_shan_metro(diagram* d, double A, int nF, double chi, double tau = 1.0, int iterations = 500000);

/*Generates and integer hash of a fermion height vector. */
int hash_func(vector<int> h);

/*Checks if two configurations of fermions are the same. */
bool configs_match(vector<int> h1, vector<int> h2);

/*
  Function to plot normalized distribution of counts.
*/
void plot_data_exp(stats s);

/*
  Runs simulation for every value of alpha between a_start and a_end, in increments of a_incr.
  Note: a_start must be LARGER than a_end.
*/
vector<stats> run_range_of_simulations_exp(double A_start, double A_end, double A_incr, int nF=100, double chi = 2, int num_its = 500000);

/*
  Runs simulation for every value of nF between nF_start and nF_end, in increments of nF_incr.
  Note: nF_start must be smaller than nF_end.
*/
vector<stats> run_Nrange_of_simulations_exp(double nF_start, double nF_end, double nF_incr, double A, double chi = 2);

/*
  Computes the derivative of S_shannon with respect to the parameter alpha at a given value of
  alpha and N (which are wrapped up within the stats struct).
*/
double dA_S_exp(stats s);

/*
  Generates a plot of S_shannon(a) - S_shannon(a_start) for alpha in the range a_end to a_start.
*/
vector<double> S_shannon_exp(vector<stats> data, double incr);

/* Calulated the derivative of the "free energy" with respect to the parameter a using the 
expectation value of C2.
*/
double Fprime_exp(stats s);


void F_total_exp(vector<stats> data, double incr);

/* Calculates the total entropy by integrating the expression for F'(A) based on the expectation
value of C2 for a range of values of the parameter alpha, then outputs as a plot.
*/
vector<double> S_total_exp(vector<stats> data, double incr);

/* Calculates the total entropy, using the rectangle method of integration. */
double S_total_exp_rec(vector<stats> data, double incr);

void S_total_exp_under(vector<stats> data, double incr);

/* Calculates the Boltzmann entropy using the expectation value of logdim for data taken over
a range of values of the parameter a and creates a plot.
*/
vector<double> S_boltz_exp(vector<stats> data);

/*Calculates the Shannon entropy by subtracting the total and the Boltzmann entropies.*/
vector<double> S_subtract(vector<double> S_tot, vector<double> S_boltz);

/*Prints vectors x and y in two columns (for plotting with gnuplot) in a file named by
filename.*/
void print(string filename, vector<double> x, vector<double> y);

/* Prints (entropy) data as a function of N for every /integer/ value of A.  Here, "data"
is a vector of (vectors of entropy data taken over a range of A) for several different values
N.
*/
void print_functionN(string partial_filename, vector<vector<double>> data, vector<double> As, vector<double> Ns);

/* Calculates the inital configuration of fermions for a given value of the coupling constant
and number of fermions.
*/
diagram initial_config(double ratio, int nF);

/* Computes the number of edges, as well as the values of C2 and logdim for the configuration
of fermions contained in the height vector of a diagram d.
*/
diagram diagram_props(diagram d, int nF);

/* Counts edges of a diagram. */
int count_edges(diagram d);

/* Computes the derivative of G(tau) := Sshannon(A/tau, chi/tau) with respect to tau.
*/
double Gprime(double tau, double A, double chi, int nF);

/* Computes S_shannon by integrating Gprime from 0 to 1.
*/
void S_shannon_tau(double A, double chi, int nF);

/* Computes S_shannon by integrating the first law.
*/
double S_shannon_E(double A0, double chi0, int nF, double initial = 0.01);

/* Finds and prints the height vector for all saddle points of 3 fermions.
*/
void find_saddle_points(double A, double chi);

/* Computes fermion energy given a height vector.
*/
double compute_energy(vector<int> h, double A, double chi);

/* Comparison function that can be fed to perform_iterations to run Metropolis Monte Carlo. */
bool metropolis(double rand, double A, double chi, int N, double dC2, double dlogdim, int edges, int new_edges);

/* Comparison function that can be fed to perform_iterations to run gradient descent 
Monte Carlo, i.e. accepting only lower energy configs.*/
bool gradientdescent(double rand, double A, double chi, int N, double dC2, double dlogdim, int edges, int new_edges);

#endif