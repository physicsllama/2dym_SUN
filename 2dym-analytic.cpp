#include <iostream>
#include <fstream>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/tools/roots.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include "2dym-analytic.h"
#include "2dym.h"

using namespace std;
using namespace boost::math;
using namespace boost::math::double_constants;
using namespace boost::math::tools;
using namespace Eigen;

/* Define a functor that will allow us to solve for k = b/a given alpha.
*/
struct k_functor {
  k_functor(double const& ap) : alpha(ap) {}

  double operator()(double const& k) {
    if(k >= 1.0) return numeric_limits<double>::max();
    double K = ellint_1(k);
    double E = ellint_2(k);
    return 2*K*(2*E + (k*k-1)*K) - alpha;
  }
private:
  double alpha;
};


/* Evaluate the strong coupling density rho(h) given a and b */
double rho(double a, double b, double h) {
  if(abs(h) >= a) return 0.0;
  if(b == 0.0) {
    return 2.0*sqrt(a*a - h*h)/(pi*a*a);
  }
  if(abs(h) <= b) return 1.0;
  return (2*sqrt(a*a -h*h)*sqrt(h*h-b*b))/(pi*a * abs(h)) * ellint_3(b/a,b*b/(h*h));
}


void main_analytic(int argc, double A, double chi) {
  //If there is a command line argument alpha, finds a,b for that alpha and plots the model.
  if( argc > 1 ) {
    params p = find_params(A / chi);
    cerr << "b = " << p.b << ", a = " <<  p.a << ", A = " << A << endl;
    plot_model(p, 100);

  // cerr << Fprime(A / chi) << endl;

  //Otherwise finds and plots for every alpha from 3.0 to 7.0 in 0.1 intervals
  } else {
    S_total((pi*pi), 0.1, 0.1, 400);
  }
}

params find_params(double ratio) {
  double a,b;
  if(ratio > pi*pi / 2.0) {

    /* Solve for k. */
    long unsigned int its = 20;
    pair<double,double> r = bracket_and_solve_root(k_functor(ratio), 0.5, 2.0, true, eps_tolerance<double>(31), its);
    //cerr << "Converged in " << its << " iterations" << endl;
    double k = 0.5*(r.first + r.second);

    // Plug back in for a and b.
    a = (2.0/ratio) * ellint_1(k);
    b = a * k;
  } else {
    b = 0.0;
    a = sqrt(2.0/ratio);
  }

  params p = {ratio, a, b};
  return p;
}

void plot_model(params p, int nF) {

  stringstream strm;
  strm << "./exp/models/ratio_" << fixed << setprecision(2) << p.ratio << ".dat";
  ofstream output(strm.str());

  int N = 2*nF +1;

  // Plot the density and numerically integrate along the way.
  double total = 0.0;
  for(int j = 0; j < 6*nF+1; j++) {
    double h = (j -3*nF)/double(N);
    double rhoh = rho(p.a,p.b,h);
    output << h << ' ' << rhoh << endl;
    total += rhoh / double(N);
  }
  // output.close();

  // Should be close to 1
  cerr << "Numerical total: " << total << endl;
}

vector<params> find_range_of_params(double start, double finish, double increment, bool plot) {
  ofstream output;
  output.open("params.dat");
  vector<params> vec;

  for(double i = start; i >= finish; i -= increment) {
    params p = find_params(i);
    output << p.a << " " << p.b << endl;

    vec.push_back(p);
  }

  output.close();
  return vec;
}

double Fprime_DK(double A) {
  params p = find_params(A/2.0);
  if(A/2.0 > pi*pi/2.0) {
    double k = p.b/p.a;
    return (-1.0/6.0)*p.a*p.a + (1.0/12.0)*p.a*p.a*(1-k*k) + (1.0/24.0) - (1.0/96.0)*p.a*p.a*p.a*p.a*(1-k*k)*(1-k*k)*A;
  } else {
    return -1.0/(2.0 * A) + (1.0/24.0);
  }
}

void S_total(double A_start, double A_end, double incr, int N) {
  double totalF = 0;
  double Stot = 0;
  double dF;
  ofstream output("./an/S_tot/Stot_an_N_" + to_string(N) + ".dat", ofstream::out);

  double chi = 2.0;
  for(double A = A_start; A>=A_end; A-= incr) {
    //prints the integral
    output << A << " " << Stot << endl;
    // Derivative of F w.r.t. A
    dF = Fprime_DK(A);
    totalF -= dF * incr;
    //compute Stotal at this value of alpha
    Stot = (N*N)*(totalF - A*Fprime_DK(A) + A_start*Fprime_DK(A_start));
  }
  output.close();
}

double find_determinant(double A, double chi, int nF, vector<int> h) {
  int N = 2*nF+1;
  MatrixXd m(N,N);
  
  //create matrix
  for(int i=0; i<N; i++){
    double sum = 0;
    for(int j=0; j<N; j++){
      if(i!=j){
        sum += 1.0 / double((h[j] - h[i])*(h[j] - h[i]));
      }
    }
    for(int j=0; j<N; j++){
      if(i!=j){
        m(i,j) = -chi*1.0/double((h[i]-h[j]) * (h[i] - h[j]));
      }
      if(i==j){
        m(i,j) = A / double(N) + chi*sum;
      }
    }
  }

  SelfAdjointEigenSolver<MatrixXd> eigensolver(N);
  eigensolver.compute(m, EigenvaluesOnly);
  if(eigensolver.info() != Success) abort();
  cerr << eigensolver.eigenvalues().transpose() << endl;
  double det = 1;
  for(int i=0; i<N; i++){
    det *= eigensolver.eigenvalues()[i];
  }
  return det;
}