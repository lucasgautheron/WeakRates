
inline double gas_density(double T, double x, double mass, const int g = 2) // x = mu/m
{
    const double prefactor = g * pow(mass / HBAR / CELERITY_FM, 3.) / (2*M_PI*M_PI);
    
    double integral = 0;
    const int N = 512;
    double umin = 1, umax = 20*(T/mass+x);
    const double du = (umax-umin)/double(N);
    for(int i = 0; i < N; ++i)
    {
        const double u = umin + (umax-umin)*(double(i)+0.5)/double(N);
        integral += u*sqrt(u*u-1) / (exp((u-x)*mass/T)+1) * du;
    }
    return prefactor * integral;
}

// recover mu from density, assuming degenerate gas (T << Tf)
inline double degenerate_potential(double mass, double density)
{
    return mass*sqrt(1+pow(density*1.705199692e9, 2./3.));
}

inline double gas_potential(double T, double density, double mass, const int g = 2)
{
    const double max_error = 1e-5;
    const int max_iterations = 256;
    const double x0 = degenerate_potential(mass, density)/mass; // first guess

    // if cold enough, use the degenerate expression. 
    // (avoids infinities)
    const double fermi_temperature = 0.5 * pow(HBARC_FM/mass, 2.) * pow(3*M_PI*M_PI*density, 2./3.);
    if(T/fermi_temperature < 0.001) return x0*mass;
    
    double x = 1.01, x1 = 1.01, error = max_error+1;
    
    int i = 0;
    while (error > max_error && i < max_iterations)
    {
        double cdensity = gas_density(T, x, mass, g);
        double dd = (gas_density(T, x+max_error/100, mass, g)-cdensity)/(max_error/100);
        x1 = x - (cdensity-density)/(dd);
        error = fabs(x1-x);
        x = x1;
        ++i;
    }
    return x*mass;
}

inline double average_neutrino_energy(double T, double mu_nu)
{
    return mu_nu/T > -40 ? T*(6*gsl_sf_fermi_dirac_int (3,mu_nu/T)) / (2*gsl_sf_fermi_dirac_int (2,mu_nu/T)) : 3*T;
}

inline double fermi_dirac(double E, double mu, double T)
{
    return 1./(exp((E-mu)/T)+1);
}

// fast electron capture simulation, using the parametrization by Langake
double electron_capture_fit(int A, int Z, double T, double mu_e = M_ELECTRON, double Q = 1e10);

// Bruenn 1985
double electron_capture_proton(double T, double nb, double mu_e, double mu_nu, double eta_pn);

// Bruenn 1985 + Horowitz 1997
double nucleus_scattering_cross_section(int A, int Z, double eta, double eps_neutrino, double density);

// Bruno Peres
double  elec_capt_proton_effective(double mu_e, double mu_nu, double t, double muneut, double mup, double y_p, double y_n, double eta_pn);
double  elec_capt_heavy_nuclei_effective(double mue, double mu_nu, double na, double t, double muneut, double mup, double zheavy, double aheavy);
double eta_np_v3(double mun, double mup, double t);
double eta_pn_v3(double mun, double mup, double t);

inline void gausslegendre(int n,double x1,double x2,double* x,double* w) 
{  
  int j,m;
  double p1, p2, p3, pp, xl, xm, z, z1;
  //const double eps = 1e-15;
  const double eps = 1e-6;
  m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  // debut boucle 
  for (int i=1;i<=m;i++)
  {
    z=cos(M_PI*(double(i) - 0.25)/(double(n) + 0.5));
    int count = 0;
    do {
    p1 = 1;
    p2 = 0;
    for (j=1;j<=n;j++)
    {
      p3 = p2;
      p2 = p1;
      p1 = ((2*j - 1)* z * p2 - (j - 1)*p3)/j;
    }
    pp = n * (z * p1 - p2)/(z*z - 1);
    z1 = z;
    z = z1 - p1/pp;}
    while(fabs(z-z1) > eps && ++count < 64);
    x[i] = xm - xl * z;
    x[n+1-i] = xm + xl * z;
    w[i] = 2 * xl / ((1 - z*z) * pp * pp);
    w[n+1-i] = w[i];
  }
}
