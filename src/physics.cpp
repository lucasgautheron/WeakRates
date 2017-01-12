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

    if(eta < -40) return 0;

    double rate = 0.693 * beta/K * pow(T/M_ELECTRON, 5.);
    
    // gsl fermi integrals include a normalization factor 1/Gamma(j) hence the multiplication by (j-1)!
    rate *= 24*gsl_sf_fermi_dirac_int (4,eta) - 2 * chi * 6*gsl_sf_fermi_dirac_int (3,eta) + chi*chi * 2*gsl_sf_fermi_dirac_int (2, eta);

    return rate;
}

double electron_capture_ps(double T, double mu_e, double mu_nu, double Q)
{
    double integral = 0;
    const int N = 256;
    double Emin = max(0, M_ELECTRON-Q), Emax = 30*T;
    const double dE = (Emax-Emin)/double(N);
    for(int i = 0; i < N; ++i)
    {
        const double E = Emin + (Emax-Emin)*(double(i)+0.5)/double(N);
        integral += E * E * (1-fermi_dirac(E, mu_nu, T)) * fermi_dirac(E+Q, mu_e, T) * (E+Q)*(E+Q) * dE;
    }
    return integral;

}

// Bruenn 1985
double electron_capture_proton(double T, double nb, double mu_e, double mu_nu, double eta_pn)
{
    const double Vud = 0.97427;
    const double gA = 1.24, gV = 1.;
    const double Q = M_NEUTRON-M_PROTON;

    double rate = 4*CELERITY_FM * (gV*gV+3*gA*gA) * Vud*Vud * pow(2*M_PI, -3.) * FERMI_COUPLING*FERMI_COUPLING;
    rate *= eta_pn;
    
    double integral = 0;
    const int N = 256;
    double Emin = max(0, M_ELECTRON-Q), Emax = 30*T;
    const double dE = (Emax-Emin)/double(N);
    for(int i = 0; i < N; ++i)
    {
        const double E = Emin + (Emax-Emin)*(double(i)+0.5)/double(N);
        integral += E * E * (1-fermi_dirac(E, mu_nu, T)) * fermi_dirac(E+Q, mu_e, T) * (E+Q)*(E+Q) * sqrt(1-(M_ELECTRON/(E+Q))*(M_ELECTRON/(E+Q))) * dE;
    }
    rate *= integral;

    return rate;
}

// Bruenn 1985 + Horowitz 1997
double nucleus_scattering_cross_section(int A, int Z, double eta, double eps_neutrino, double density)
{
    const double c_a = 0.5, c_v = 0.5 + 2.0 * 0.23119;
    density *= 1e39; // fm^{-3} -> cm^{-3}

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

    const double f5_over_f3 = eta > -40 ? 120 * gsl_sf_fermi_dirac_int (5,eta) / (6 * gsl_sf_fermi_dirac_int (3,eta)) : 20;
    sigma_scattering_nucleus *= f5_over_f3; 

    return sigma_scattering_nucleus;
}

// Bruno Peres

double  elec_capt_proton_effective(double mu_e, double mu_nu, double t, double muneut, double mup, double y_p, double y_n, double eta_pn) 
{  
  double g_squared, ga, q;
  q = 1.2935; // MeV
  ga = 1.2695;
  g_squared = 5.2971899e-5; // (MeV-2 cm2) * 1e39  
  double me_c2 = M_ELECTRON; // MeV
  double h_bar_c = HBARC_CM; // MeV cm
  double c = 29979245800; //cm s-1
  /* integration */
  double rproton = 0;
  double rprotonold;
  double cste;
  cste = g_squared*(3*ga*ga+1)/(2*pow(M_PI*h_bar_c,3));
  double first_zero, second_zero, borne_inf, resultat;
  first_zero = min(mu_nu,mu_e-q);
  second_zero = max(mu_nu,mu_e-q);
  borne_inf = 0;
  const int n = 64;
  double x1,x2;
  double* x = new double[n];
  double* w = new double[n];
  if (first_zero > borne_inf)
    {
      // integration de borne_inf a first_zero
      x1 = borne_inf;
      x2 = first_zero;
      gausslegendre(n,x1,x2,x,w);
      resultat=0;
      for (int l=1;l<=n;l++)
	{
	  resultat = resultat + w[l] * cste*x[l]*x[l]*(x[l]+q)*(x[l]+q)*sqrt(1-(me_c2*me_c2/((x[l]+q)*(x[l]+q))))*((eta_pn/(1+exp((x[l]+q-mu_e)/t)))*(1.0-1.0/(1.0+exp((x[l]-mu_nu)/t))) );//-(eta_np*(1.0-1.0/(1+exp((x[l]+q-mu_e)/t)))*(1.0/(1+exp((x[l]-mu_nu)/t)))));
	}
      rproton = resultat;
      borne_inf = first_zero;
    }
  if (second_zero > borne_inf)
    {
      // integration de borne_inf a second_zero
      x1 = borne_inf;
      x2 = second_zero;
      gausslegendre(n,x1,x2,x,w);
      resultat=0;
      for (int o=1;o<=n;o++)
	{
	  resultat = resultat + w[o] * cste*x[o]*x[o]*(x[o]+q)*(x[o]+q)*sqrt(1-(me_c2*me_c2/((x[o]+q)*(x[o]+q))))*((eta_pn/(1.0+exp((x[o]+q-mu_e)/t)))*(1.0-1.0/(1.0+exp((x[o]-mu_nu)/t))) );//-(eta_np*(1.0-1.0/(1.0+exp((x[o]+q-mu_e)/t)))*(1.0/(1.0+exp((x[o]-mu_nu)/t)))));
	}
      rproton = rproton + resultat;
      borne_inf = second_zero; 
    }
  // integration de borne_inf a l'infini
  do
    {
      x1 = borne_inf;
      x2 = x1 + 4;
      gausslegendre(n,x1,x2,x,w);
      resultat=0;
      for (int p=1;p<=n;p++)
	{
	  resultat = resultat + w[p] * cste*x[p]*x[p]*(x[p]+q)*(x[p]+q)*sqrt(1.0-(me_c2*me_c2/((x[p]+q)*(x[p]+q))))*((eta_pn/(1.0+exp((x[p]+q-mu_e)/t)))*(1.0-1.0/(1.0+exp((x[p]-mu_nu)/t))) );//-(eta_np*(1.0-1.0/(1.0+exp((x[p]+q-mu_e)/t)))*(1.0/(1.0+exp((x[p]-mu_nu)/t)))));
	}
      rprotonold = rproton;
      rproton = rproton + resultat;
      borne_inf = x2;
    }
  while (fabs((rproton-rprotonold)/rproton) > 1e-10 && fabs(x2) < second_zero+100);
  //while (fabs((rproton-rprotonold)/rproton) > 1e-10);
  //while (++);
  delete [] x;
  delete [] w;
  rproton = rproton * c;
  return 1e-39 * rproton;
}


double  elec_capt_heavy_nuclei_effective(double mue, double mu_nu, double na, double t, double muneut, double mup, double zheavy, double aheavy) {  
  double g_squared, ga, np, nh, q, mun;  
  ga = 1.2695;
  g_squared = 5.2971899e-5; // (MeV-2 cm2) *1e39
  double c = 29979245800; //cm s-1
  double me_c2 = M_ELECTRON; // MeV
  double h_bar_c = HBARC_CM; // MeV cm 
  /* calcul de Np et Nh */
  double nheavy;
  nheavy = aheavy - zheavy;
  mun = mup * zheavy + muneut *(aheavy - zheavy);
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

  /* calcul de q */
  double delta = 3; // niveau fils excite en MeV

  // with Bruenn's formula
  q = max(muneut - mup + delta, 0.0);

  // integration 
  double eta_nu = mu_nu / t;
  double eta_e = mue / t;
  double rheavy = 0;
  double rheavyold;
  double cste;
  cste = 4./(pow(2.0*M_PI*h_bar_c,3))*g_squared*ga*ga*2./7.*na*np*nh;
  double first_eta, second_eta, borne_inf, resultat;
  first_eta = min(eta_nu,eta_e);
  second_eta = max(eta_nu,eta_e);
  borne_inf = me_c2 - q;
  const int n = 64;
  double x1,x2;
  double* x = new double[n];
  double* w = new double[n];
  if (first_eta > borne_inf)
    {
      // integration de borne_inf a first_eta
      x1 = borne_inf;
      x2 = first_eta;
      gausslegendre(n,x1,x2,x,w);
      resultat=0;
      for (int l=1;l<=n;l++)
	{
	  resultat = resultat + w[l] * cste*x[l]*x[l]*(x[l]+q)*(x[l]+q)*sqrt(1.0-(me_c2*me_c2/((x[l]+q)*(x[l]+q))))*((1.0/(1.0+exp((x[l]+q-mue)/t)))*(1.0-1.0/(1.0+exp((x[l]-mu_nu)/t))) );//-(exp((muneut-mup-q)/t)*(1.0-1.0/(1+exp((x[l]+q-mue)/t)))*1.0/(1.0+exp((x[l]-mu_nu)/t))));
	}
      rheavy = resultat;
      borne_inf = first_eta;
    }
  if (second_eta > borne_inf)
    {
      // integration de borne_inf a second_eta
      x1 = borne_inf;
      x2 = second_eta;
      gausslegendre(n,x1,x2,x,w);
      resultat=0;
      for (int o=1;o<=n;o++)
	{
	  resultat = resultat + w[o] * cste*x[o]*x[o]*(x[o]+q)*(x[o]+q)*sqrt(1.0-(me_c2*me_c2/((x[o]+q)*(x[o]+q))))*((1.0/(1.0+exp((x[o]+q-mue)/t)))*(1.0-1.0/(1.0+exp((x[o]-mu_nu)/t))) );//-(exp((muneut-mup-q)/t)*(1.0-1.0/(1+exp((x[o]+q-mue)/t)))*1.0/(1.0+exp((x[o]-mu_nu)/t))));
	}
      rheavy = rheavy + resultat;
      borne_inf = second_eta; 
    }
  // integration de borne_inf a l'infini
  do
    {
      x1 = borne_inf;
      x2 = x1 + 5; // a voir
      gausslegendre(n,x1,x2,x,w);
      resultat=0;
      for (int p=1;p<=n;p++)
	{
	  resultat = resultat + w[p] * cste*x[p]*x[p]*(x[p]+q)*(x[p]+q)*sqrt(1.0-(me_c2*me_c2/((x[p]+q)*(x[p]+q))))*((1.0/(1.0+exp((x[p]+q-mue)/t)))*(1.0-1.0/(1.0+exp((x[p]-mu_nu)/t))) );//-(exp((muneut-mup-q)/t)*(1.0-1.0/(1+exp((x[p]+q-mue)/t)))*1.0/(1.0+exp((x[p]-mu_nu)/t))));
	}
      rheavyold = rheavy;
      rheavy = rheavy + resultat;
      borne_inf = x2;
    }
  //while (fabs(x2) < second_eta+100);
  while ((fabs((rheavy-rheavyold)/rheavy) > 1e-10) && (x2 < 200));
  //while (fabs(rheavy-rheavyold)/rheavy > 1e-10);
  delete [] x;
  delete [] w;
  rheavy = rheavy * c / aheavy;
  return 1e-39 * rheavy;
}

double eta_np_v3(double mun, double mup, double t)
{
  double mp_c2 = M_PROTON; // MeV
  double mn_c2 = M_NEUTRON; // MeV (booklet p.126)
  double eta_n = mun / t;
  double eta_p = mup / t;
  double h_bar_c = HBARC_FM; // MeV.fm
  double borne_inf, resultat, eta_np_tot, eta_np_old;
  borne_inf = 0;
  int n = 64;
  int count = 0;
  double x1,x2;
  double* x = new double[n];
  double* w = new double[n];
  eta_np_tot = 0;
  do
    {
      count++;
      x1 = borne_inf;
      x2 = x1 + 0.1; // a voir
      gausslegendre(n,x1,x2,x,w);
      resultat=0;
      for (int p=1;p<=n;p++)
	{
	  resultat = resultat + w[p] * x[p]*x[p]/(M_PI*M_PI)/(1+exp(x[p]*x[p]/(2*mn_c2*t)-eta_n))*(1-1/(1+exp(x[p]*x[p]/(2*mp_c2*t)-eta_p)));
	}
      eta_np_old = eta_np_tot;
      eta_np_tot = eta_np_tot + resultat;
      borne_inf = x2;
    }
  while((fabs((eta_np_tot-eta_np_old)/eta_np_tot) > 1.0e-10) or (count < 100));
  delete [] x;
  delete [] w;
  eta_np_tot = eta_np_tot / pow(h_bar_c,3); // converts MeV^3 into fm^-3
  return eta_np_tot; // fm^-3
}

double eta_pn_v3(double mun, double mup, double t)
{
  double mp_c2 = 938.272013; // MeV
  double mn_c2 = 939.56536; // MeV (booklet p.126)
  double eta_n = mun / t;
  double eta_p = mup / t;
  double h_bar_c = 197.3269631; // MeV.fm
  double borne_inf, resultat, eta_pn_tot, eta_pn_old;
  borne_inf = 0;
  //const int n = 64;
  const int n = 32;
  int count = 0;
  double x1,x2;
  double* x = new double[n];
  double* w = new double[n];
  eta_pn_tot = 0;
  do
    {
      count++;
      x1 = borne_inf;
      x2 = x1 + 0.1; // a voir
      gausslegendre(n,x1,x2,x,w);
      resultat=0;
      for (int p=1;p<=n;p++)
	{
	  resultat = resultat + w[p] * x[p]*x[p]/(M_PI*M_PI)/(1+exp(x[p]*x[p]/(2*mp_c2*t)-eta_p))*(1-1/(1+exp(x[p]*x[p]/(2*mn_c2*t)-eta_n)));
	}
      eta_pn_old = eta_pn_tot;
      eta_pn_tot = eta_pn_tot + resultat;
      borne_inf = x2;
    }
  while((fabs((eta_pn_tot-eta_pn_old)/eta_pn_tot) > 1.0e-10) or (count < 100));
  delete [] x;
  delete [] w;
  eta_pn_tot = eta_pn_tot / pow(h_bar_c,3); // converts MeV^3 into fm^-3
  return eta_pn_tot; // fm^-3 
}
