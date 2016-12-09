
// recover mu_e from (nb, Ye), assuming degenerate electron gas
inline double electron_potential(double nb, double Ye)
{
    return M_ELECTRON*sqrt(1+pow(nb*Ye*1.705199692e9, 2./3.));
}

// fast electron capture simulation, using the parametrization by Langake
double electron_capture_fit(int A, int Z, double T, double mu_e = M_ELECTRON, double Q = 1e10)
{
    const double K = 6146, beta = 4.6, dE = 2.5;
    if(Q>1e6) Q = beta_decay_Q(A, Z);
    if(Q<-1e6) return 0;
    const double chi = (Q-dE)/T;
    const double eta = chi+mu_e/T;

    double rate = 0.693 * beta/K * pow(T/M_ELECTRON, 5.);
    
    // gsl fermi integrals include a normalization factor 1/Gamma(j) hence the multiplication by (j-1)!
    rate *= 24*gsl_sf_fermi_dirac_int (4,eta) - 2 * chi * 6*gsl_sf_fermi_dirac_int (3,eta) + chi*chi * 2*gsl_sf_fermi_dirac_int (2, eta);

    return rate;
}
