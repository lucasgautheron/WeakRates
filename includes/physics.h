
inline double gas_density(double T, double mu, double mass, const int g = 2)
{
    const double prefactor = g * pow(mass * CELERITY_FM / HBAR,3.) / (2*M_PI*M_PI);
    
    double integral = 0;
    const int N = 1000;
    double umin = 1, umax = 30*T/mass;
    const double du = (umax-umin)/double(N);
    for(int i = 0; i < N; ++i)
    {
        const double u = umin + (umax-umin)*(double(i)+0.5)/double(N);
        integral += u*sqrt(u*u-1) / (exp((u-mu/mass)*mass/T)-1) * du;
    }
    return prefactor * integral;
}

// recover mu_e from (nb, Ye), assuming degenerate electron gas
inline double degenerate_potential(double mass, double density)
{
    return mass*sqrt(1+pow(density*1.705199692e9, 2./3.));
}

inline double gas_potential(double T, double density, double mass, const int g = 2)
{
    const double max_error = 1e-3;
    const int max_iterations = 1000;
    const double x0 = degenerate_potential(mass, density); // first guess
    
    double x = x0, x1 = x0, error = max_error+1;
    
    int i = 0;
    while (error > max_error && i < max_iterations)
    {
        double cdensity = gas_density(T, x, mass, g);
        double dd = (cdensity-gas_density(T, x-max_error*x/10, mass, g))/(x*max_error/10);
        x1 = x - (cdensity-density)/(dd);
        error = fabs(x1-x);
        x = x1;
        ++i;
    }
    return x;
}

inline double average_neutrino_energy(double T, double mu_nu)
{
    return T*(6*gsl_sf_fermi_dirac_int (3,mu_nu/T)) / (2*gsl_sf_fermi_dirac_int (2,mu_nu/T));
}

inline double fermi_dirac(double E, double mu, double T)
{
    return 1./(exp((E-mu)/T)+1);
}

// fast electron capture simulation, using the parametrization by Langake
double electron_capture_fit(int A, int Z, double T, double mu_e = M_ELECTRON, double Q = 1e10);

// Bruenn 1985
double electron_capture_proton(double T, double nb, double mu_e, double mu_nu, double n_p, double n_n);

// Bruenn 1985 + Horowitz 1997
double nucleus_scattering_cross_section(int A, int Z, double eps_neutrino, double density);
