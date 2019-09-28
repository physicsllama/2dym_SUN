/* A storage place for old functions that have fallen out of favor
and likely won't be included in the final version of the code, but 
should be kept around for convinience until a final judgment day comes.*/

double S_total_exp_rec(vector<stats> data, double incr) {
  double totalF = 0;
  double Stot = 0;
  double initial = data.front().chi * data.front().explogdim;
  double dF;
  stats s;
  stats prev;
  string filename = "./exp/S_tot_rec/Stot_exp_N_" + to_string(data.front().N) + ".dat";
  ofstream output(filename, ofstream::out);
  ofstream output2("./exp/S_tot_rec/Ssub_exp_N_" + to_string(data.front().N) + ".dat");
  for(int i = 1; i < data.size(); i++) {
    s=data[i];
    prev=data[i-1];
    //printing the value of the integral
    output << prev.A << " " << Stot << endl;
    output2 << prev.A << " " << Stot - (prev.chi*prev.explogdim - initial) << endl;
    // Derivative of F w.r.t. A
    dF = Fprime_exp(data[i]);
    //Numerical integration of F
    totalF -= dF * incr;
    //compute Stotal at this value of alpha
    Stot = (s.N*s.N)*(totalF - s.A*Fprime_exp(s) + (data.front().A)*Fprime_exp(data.front()));
  }
  output.close();

  return Stot;
}

void S_total_exp_under(vector<stats> data, double incr) {
  double totalF = 0;
  double Stot = 0;
  double dF;
  string filename = "./exp/S_tot_under/Stot_exp_N_" + to_string(data.front().N);
  ofstream output(filename, ofstream::out);
  for(stats s : data) {
    //printing the value of the integral
    output << s.A << " " << Stot << endl;
    // Derivative of F w.r.t. A
    dF = Fprime_exp(s);
    //Numerical integration of F
    totalF -= dF * incr;
    //compute Stotal at this value of alpha
    Stot = (s.N*s.N)*(totalF - s.A*Fprime_exp(s) + (data.front().A)*Fprime_exp(data.front()));
  }
  output.close();
}

double Gprime(double tau, double A, double chi, int nF) {
  int num_its = 2000;
  stats s = run_simulation(A/tau, nF, chi/tau, num_its);

  //returns the unbiased variance of the fermion energy
  return (1.0/(tau))*(double(num_its)/double(num_its - 1))*s.varE;
}

void S_shannon_tau(double A, double chi, int nF) {
  double incr = 1.0/10000.0;
  double total = 0;
  double g;

  ofstream S_out("./tau/S_shan/S_shan_N_" + to_string(2*nF +1) + ".dat");
  ofstream G_out("./tau/Gprime/Gprime_N_" + to_string(2*nF +1) + ".dat");

  //Integrates Gprime over a range of values of tau from 0 and 1
  // Generates a plot of Gprime over these values, and its integral, Sshannon
  for(double i = incr; i < 1; i+= incr) {
    g = Gprime(i, A, chi, nF);
    S_out << i << " " << total << endl;
    G_out << i << " " << g << endl;
    total += g * incr;
  }
  cerr << "Calculation of S_shan with N = " << 2*nF +1 << " done!" << endl;
}

void find_saddle_points(double A, double chi) {
  vector<int> X;
  int lower = 10;
  for(int i = -lower; i <= lower; i++) {
    X.push_back(i);
  }

  int range = X.size();
  double E;

  vector<vector<int>> minima;

  for(int i = 0; i < range; i++) {
    for(int j = 0; j < i; j++) {
      for(int k = 0; k < j; k++) {
        vector<int> h = {X[i],X[j],X[k]};
        E = compute_energy(h, A, chi);

        vector<int> hp = {X[i]+1, X[j], X[k]};
        vector<int> hm = {X[i]-1, X[j], X[k]};
        if(E < compute_energy(hp,A,chi) && E < compute_energy(hm,A,chi)) {
          hp = {X[i], X[j]+1, X[k]};
          hm = {X[i], X[j]-1, X[k]};
          if(E < compute_energy(hp,A,chi) && E < compute_energy(hm,A,chi)) {
            hp = {X[i], X[j], X[k]+1};
            hm = {X[i], X[j], X[k]-1};
            if(E < compute_energy(hp,A,chi) && E < compute_energy(hm,A,chi)) {
              minima.push_back(h);
              cerr << h[0] << " " << h[1] << " " << h[2] << " " << E << endl;
              cerr << "    " << E << " " << compute_energy(hp,A,chi) << " " << compute_energy(hm,A,chi) << endl;
            }
          }
        }
      }
    }
  }
}

double compute_energy(vector<int> h, double A, double chi) {
  int N = h.size();
  double C2 = 0;
  double logdim = 0;

  sort(h.begin(), h.end());

  for(int i = 0; i < N; i++) {
    C2 += (A/double(2*N))*(h[i]^2);
    for(int j =0; j < i; j++) {
      if(h[i] == h[j]) {
        return 10000;
      }
      logdim -= chi*log(h[i] - h[j]);
    }
  }

  return C2 + logdim;
}

//NEEDS TO BE UPDATED TO USE THE HASH FUNCTION CORRECTLY
double S_shan_metro(diagram* d, double A, int nF, double chi, double tau, int iterations) {

  /* Initially take N constant, but we may want to do some extrapolation later */
  int N = 2*nF + 1;
  int burnin = iterations/5;
  int accepts = 0;
  // ofstream output("./exp/E/E_A_" + to_string(A) + ".dat");

  //create a set to track unique configs visited
  set<double> prev;
  int hash_val;
  prev.insert(hash_func(d->h));

  //create a heap to track min energy configs visited
  vector<double> queue;
  queue.push_back(d->C2*A/double(2*N) - chi*d->logdim);
  make_heap(queue.begin(), queue.end(), greater2());

  for(int i = 0; i < iterations; i++) {

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

    /* Now we want to consider mutating d->h[j] to d->h[j]+dh.
       To decide whether to keep the new configuration we have to find the difference in effective action.
       This means finding how much C_2(R) and log dim(R) change */
    double dC2 = 2*d->h[j]*dh + 1;
    double dlogdim = 0;
    for(int k = 0; k < N; k++) {
      if(j == k) continue;
      dlogdim += log(1 + dh / double(d->h[j] - d->h[k]));
    }

    // Find the number of edges of the new diagram
    int new_edges = d->edges;
    if(d->h[j+dh] == d->h[j]+2*dh) new_edges -= 2;
    if(d->h[j-dh] == d->h[j]-dh) new_edges += 2;

    /* Find the probability p to jump to the other configuration (p > 1 means we definitely do it)
       We have to correct for the change in the number of edges - configurations
       with fewer edges are more likely to stay put.
    */
    double p = exp(-A * dC2/double(2 * N) + chi*dlogdim) * (d->edges / double(new_edges));
    // double p = exp(-A * dC2/double(2 * N) + chi*dlogdim);

    if(random_real() < p) {
      //Switch to new state
      d->h[j] += dh;
      d->C2 += dC2;
      d->logdim += dlogdim;
      d->edges = new_edges;
      accepts++;

      hash_val = hash_func(d->h);
      if(prev.count(hash_val) == 0) {
        prev.insert(hash_val);

        queue.push_back(d->C2*A/double(2*N) - chi*d->logdim);
        push_heap(queue.begin(), queue.end(), greater2());
      }
    }
  }
    
  // cerr << "unique configs: "<< prev.size() << endl;

  int counter = 0;
  double Z = 0;
  double E = 0;
  vector<double> factors;
  ofstream output("E_metro.dat");
  ofstream output2("Z_metro.dat");
  while(!queue.empty()) {
    E = queue.front();
    Z += exp(-E);
    factors.push_back(exp(-E));

    output << counter << " " << E << endl;
    output2 << counter << " " << Z << endl;

    pop_heap(queue.begin(), queue.end(), greater2());
    queue.pop_back();

    counter++;
  }
  output.close();
  output2.close();

  //compute S_shan
  // double prob;
  // double S_shan = 0;
  // for(double d : factors) {
  //   prob = d/Z;
  //   S_shan -= prob*log(prob);
  // }

  return Z;

}

//NUMERICAL ERROR FROM S_TOT
/*
double F3prime_max = abs((-1.0/double(8*N*N*N*N*N))*(data[0].expC2_cub - 3*data[0].expC2*data[0].expC2_sq + 2*data[0].expC2*data[0].expC2*data[0].expC2));
double F3prime_max_num = abs((Fprime_exp(data[2]) - 2*Fprime_exp(data[1]) + Fprime_exp(data[0]))/(incr*incr));

//find maximum absolute value of second derivative
    temp = abs((-1.0/double(8*N*N*N*N*N))*(data[i].expC2_cub - 3*data[i].expC2*data[i].expC2_sq + 2*data[i].expC2*data[i].expC2*data[i].expC2));
    cerr << "temp: " << temp << ", F3prime_max: " << F3prime_max << endl;
    if(temp > F3prime_max) {
      F3prime_max = temp;
      cerr << "Change!" << endl;
    }

    //compute error
    double error = (F3prime_max*range*range*range*N*N)/double(12*(i+1)*(i+1));

    cerr << "At A = " << data[i].A << ", Stot = " << Stot << ", Sboltz = " << Sb << ", Sshan = " << Stot - Sb << ", error = " << error << endl;

    // numerically find error
    if(i > 0) {
      temp_num = abs((Fprime_exp(data[i+1]) - 2*Fprime_exp(data[i]) + Fprime_exp(data[i-1]))/(incr*incr));
      if(temp_num > F3prime_max_num) {
        F3prime_max_num = temp_num;
      }
    }
    double error_num = (F3prime_max_num*range*range*range*N*N)/double(12*(i+1)*(i+1));
    cerr << "    numerical error: " << error_num << ", F3prime_max_num: " << F3prime_max_num << ", i : " << i << ", range: " << range << endl;

*/
