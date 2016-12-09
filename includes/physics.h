
// recover mu_e from (nb, Ye), assuming degenerate electron gas
inline double electron_potential(double nb, double Ye)
{
    return M_ELECTRON*sqrt(1+pow(nb*Ye*1.705199692e9, 2./3.));
}

// fast electron capture simulation, using the parametrization by Langake
double electron_capture_fit(int A, int Z, double T, double mu_e = M_ELECTRON, double Q = 1e10);

// Bruenn 1985 + Horowitz 1997
double nucleus_scattering_cross_section(int A, int Z, double eps_neutrino, double density);
