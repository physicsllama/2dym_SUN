#ifndef ANALYTIC_H
#define ANALYTIC_H

#include <iostream>
#include <fstream>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/tools/roots.hpp>

using namespace std;
using namespace boost::math;
using namespace boost::math::double_constants;
using namespace boost::math::tools;

struct params {
  double ratio, a, b;
};

/*
*/
void main_analytic(int argc, double A, double chi);

/* Function to analytically solve for parameters a,b*/
params find_params(double ratio);

/* Plots rho(h) and numerically integrates it.  Value of the numerical integrate is
output to terminal. */
void plot_model(params p, int nF);

/* Calls find_params function to find a,b for every ratio at each increment from start to finish.
These parameters are output to "params.dat".  If the plot variable is toggled on (true by default),
this will also generate the plot of rho for each value of ratio and store in the 'models' directory. 
Note: start must be greater than finish.*/
vector<params> find_range_of_params(double start, double finish, double increment, bool plot = true);

/* Computes the first derivative of the free energy using Douglas and Kazakhov's formula.
*/
double Fprime_DK(double A);

/* Computes the total entropy relative to the reference point a_start by integrating the analytical
expression for F'(A).
*/
void S_total(double A_start, double A_end, double incr, int N);


/* Evaluate the strong coupling density rho(h) given a and b */
double rho(double a, double b, double h);

/* Computes the first derivative of the free energy using our formula. */
double F_DTV(double A, double chi);

/* Computes the first derivative of the free energy using Gromov and Santos' formula. */
double F_GS(double A, double chi);

/* Computes the first derivative of the free energy using Gross and Matytsin's formula. */
double F_GM(double A);

/* Integrated Fprime, Fanalytic, F_GS, and F_GM, and plots the result along with F_corrected.*/
void plot_F(double A_start, double A_end, double incr, int N);

/* Finds the determinant of S'' */
double find_determinant(double A, double chi, int nF, vector<int> h);

void trace_log_method(double A, double chi, int nF);

/* Computes the free energy by approximating log(Z) as the energy of the saddle point config
(first member of the returned pair) and as the energy of the saddle point config perturbed
by the determinant (second member). */
pair<double,double> F_corrected(double A, double chi, int nF);

/* Computes the first derivative of the free energy found in F_corrected. */
pair<double,double> F_prime_corrected(double A, double chi, int nF, double incr);

/* Takes derivatives of log(Z), supplied as a parameter, to compute S_shannon. */
double S_shan_an(double A, double chi, int nF, double incr, double (*logZ)(double A0, double chi0, int nF0));

/* Takes derivatives of log(Z), supplied as a parameter, using a more sophisticated
numerical derivative formula that requires 5 evaluations of log(Z) for every derivative
computed, to compute S_shannon. */
double S_shan_an_5(double A, double chi, int nF, double incr, double (*logZ)(double A0, double chi0, int nF0));

/* Calculates S_boltz by differentiating the log(Z) with respect to chi.*/
double S_boltz(double A, double chi, int nF, double incr, double (*logZ)(double A0, double chi0, int nF0));

double S_boltz_5(double A, double chi, int nF, double incr, double (*logZ)(double A0, double chi0, int nF0));

/* log(Z) from the saddle point approximation (parameter for S_shan_an). */
double saddle_point_logZ(double A, double chi, int nF);

/* log(Z) from the saddle point approximation perturbed by the determinant (parameter for S_shan_an). */
double perturbed_logZ(double A, double chi, int nF);

/* log(Z) with just the determinant term (parameter for S_shan_an). */
double justdet_logZ(double A, double chi, int nF);

/* */
double justbrute_logZ(double A, double chi, int nF);

/* */
double brute_logZ(double A, double chi, int nF);

double metro_logZ(double A, double chi, int nF);

#endif
