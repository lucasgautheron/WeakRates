#define SQUARE(x) ((x)*(x))
#define CUBIC(x) ((x)*(x)*(x))

double electron_capture_ps_int(double x, void *params);

extern gsl_integration_workspace *gsl_workspace;
extern gsl_function electron_capture_ps_func,
                    electron_capture_proton_func;

int gsl_init();

inline double average_neutrino_energy(double T, double mu_nu)
{
    return mu_nu > -10*T ? T*(6*gsl_sf_fermi_dirac_int (3,mu_nu/T)) / (2*gsl_sf_fermi_dirac_int (2,mu_nu/T)) : 3*T;
}

inline double fermi_dirac_dimless(double E, double mu)
{
    return 1./(exp(E-mu)+1);
}

inline double fermi_dirac(double E, double mu, double T)
{
    return 1./(exp((E-mu)/T)+1);
}

// fast electron capture simulation, using the parametrization by Langake
double electron_capture_fit(int A, int Z, double T, double mu_e = M_ELECTRON, double Q = 1e10);

// phase space factor for electron capture
double electron_capture_ps(double T, double mu_e, double mu_nu, double Q);

// Bruenn 1985
inline double eta_nucl(double T, double n1, double n2, double mu1, double mu2)
{
    return (n1-n2)/(exp((mu1-mu2)/T)-1);
}
double electron_capture_proton(double T, double nb, double mu_e, double mu_nu, double eta_pn);

// Bruenn 1985 + Horowitz 1997
#ifdef DEBUG
double nucleus_scattering_cross_section(double A, double Z, double eta, double eps_neutrino, double density);
#else
double nucleus_scattering_cross_section(int A, int Z, double eta, double eps_neutrino, double density);
#endif

// Bruno Peres
inline double shell_capt_factor(double aheavy, double zheavy)
{
  double nheavy = aheavy-zheavy;
  double nh, np;
  if (zheavy < 20)
    {
      np = 0;
    }
  if (zheavy >= 20 && zheavy < 28)
    {
      np = zheavy - 20;
    }
  if (zheavy >= 28)
    {
      np = 8;
    }
  if (nheavy < 34)
    {
      nh = 6;
    }
  if (nheavy >= 34 && nheavy < 40)
    {
      nh = 40 - nheavy;
    }
  if (nheavy >= 40)
    {
      nh = 0;
    }
    return nh*np;
}
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
