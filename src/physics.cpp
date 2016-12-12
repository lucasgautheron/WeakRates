#include "common.h"

// Fast electron capture simulation, using the parametrization by Langake, Martinez-Pinedo et al.
double electron_capture_fit(int A, int Z, double T, double mu_e, double Q)
{
    // K [s]
    // Characteristic time associated with the capture rate.
    const double K = 6146;
    
    // Q [1]
    // Typical GT+F matrix element
    const double beta = 4.6;

    // dE [MeV]
    // Adjusted param accounting for the difference between the energy levels of the daughter & parent nuclei excited states.
    const double dE = 2.5;

    if(Q>1e6) Q = beta_decay_Q(A, Z);
    if(Q<-1e6) return 0;
    const double chi = (Q-dE)/T;
    const double eta = chi+mu_e/T;

    double rate = 0.693 * beta/K * pow(T/M_ELECTRON, 5.);
    
    // gsl fermi integrals include a normalization factor 1/Gamma(j) hence the multiplication by (j-1)!
    rate *= 24*gsl_sf_fermi_dirac_int (4,eta) - 2 * chi * 6*gsl_sf_fermi_dirac_int (3,eta) + chi*chi * 2*gsl_sf_fermi_dirac_int (2, eta);

    return rate;
}

// Bruenn 1985 + Horowitz 1997
double nucleus_scattering_cross_section(int A, int Z, double eps_neutrino, double density)
{
    const double c_a = 0.5, c_v = 0.5 + 2.0 * 0.23119;

    double y_bruenn = 2.0 / 5.0 * pow(1.07 * A, 2./3.) * pow(eps_neutrino, 2.) / (HBARC_CM * 1.0e13);

    double sigma_scattering_nucleus = pow(eps_neutrino, 2) * (SIGMA_ZERO / 16.0) 
      / pow(M_ELECTRON, 2) * pow(A, 2) * pow(((c_a - c_v) + (2 - c_a - c_v) 
      * (2 * Z - A) / A), 2) * (y_bruenn - 1 + 
      (1 + y_bruenn) * exp(-2*y_bruenn)) / pow(y_bruenn, 3);

    // ion-ion correlation
    double a_ion = pow((4.0 * M_PI / 3.0) * density, -1. / 3.);

    double gamma_ion = pow(Z, 2) * 7.2973525376e-3 * HBARC_CM 
      / (a_ion * eps_neutrino);

    double s_ion = 1.0 / (1.0 + exp((-1.0 * a_ion * eps_neutrino 
      / (HBARC_CM)) * (log(0.3/(0.3 + 3 * gamma_ion)) + 20.0
      / 3.0  - 4.5948111 + 0.23463194 * pow(gamma_ion, 0.5) - 
      0.040086846 * gamma_ion + 0.0015845167
      * pow(gamma_ion, 1.5))));

    sigma_scattering_nucleus *= s_ion;
    return sigma_scattering_nucleus;
}
