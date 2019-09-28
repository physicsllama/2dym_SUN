
#include "2dym-analytic.h"

/* Estimate the entropy by integrating the binary entropy of the distribution */
double binary_entropy(double p) {
  if(p <= 0.0 || p >= 1.0) { return 0.0; }
  return -p*log(p) - (1.0-p)*log(1.0-p);
}

int main(int argc,char *argv[] ) {

  // Range of A to integrate
  double Amin = 4.0;
  double Amax = 20.0;
  double Astep = 0.01;
  double chi = 2.0;

  // How many steps to do the h integral
  int hsteps = 1000;

  for(double A = Amin; A <= Amax; A += Astep) {
    params p = find_params(A / chi);

    // Integrate the binary entropy
    double S = 0;
    for(int i = 0; i < hsteps; i++) {
      double dh = (p.a-p.b)/hsteps;
      double h = p.b + (i+0.5)*dh;
      double r = rho(p.a, p.b, h);
      S += 2*binary_entropy(r)*dh;
//      cout << h << ' ' << r << ' ' << S << endl;
    }

    cout << A << ' ' << S << endl;
  }


}
