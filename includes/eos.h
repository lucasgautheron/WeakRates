struct EOS_table
{
    int m_ln_rho, n_ln_t, o_y_e;
    double *ln_rho_eos,*ln_t_eos, *y_e_eos;
    virtual void dump() {}
    virtual const int size() { return 0; }
    virtual void allocate() {}

    virtual void read(const char *path, int *error) {}
    virtual void write(const char *path, int *error) {}
};


struct short_EOS_table : EOS_table
{
    int p_mu;
    double *mu_nu_eos, *elec_rate_tab_eos, *elec_rate_single_eos;
    double *scattering_xs_eos, *elec_rate_fast_eos;

    void read(const char *path, int *error);
    void write(const char *path, int *error);

    void dump()
    {
        printf("ln(rho)\n");
        for(int m = 200; m < m_ln_rho; ++m) printf("%e\n", ln_rho_eos[m]);
        printf("ln(T)\n");
        for(int n = 40; n < n_ln_t; ++n) printf("%e\n", ln_t_eos[n]);
        printf("Y_e\n");
        for(int o = 0; o < o_y_e; ++o) printf("%e\n", y_e_eos[o]);
        printf("mu_nu\n");
        for(int p = 0; p < p_mu; ++p) printf("%e\n", mu_nu_eos[p]);
    }

    const int size()
    {
        return m_ln_rho*n_ln_t*o_y_e*p_mu;
    }

    void allocate()
    {
        int s = size();
        ln_rho_eos = new double[m_ln_rho];
        ln_t_eos = new double[n_ln_t];
        y_e_eos = new double[o_y_e];
        mu_nu_eos = new double[p_mu];
        elec_rate_tab_eos = new double[s];
        elec_rate_fast_eos = new double[s];
        elec_rate_single_eos = new double[s];
        scattering_xs_eos = new double[s];
    }
};

struct full_EOS_table : EOS_table
{
    double *p_eos; double *e_eos;
    double *entropy_eos, *dp_drho_eos, *dp_deps_eos; 
    double *c_sound_squared_newtonian_eos, *mub_eos, *muq_eos; 
    double *mul_eos, *xp_eos, *xn_eos, *xa_eos, *xheavy_eos; 
    double *aheavy_eos, *zheavy_eos;
    int *eflg_eos;
    double *mu_nu_effective_frac, *thermo_eos, *yi_eos;

    void read(const char *path, int *error);

    const int size()
    {
        return m_ln_rho*n_ln_t*o_y_e;
    }

    void dump()
    {
        printf("ln(rho)\n");
        for(int m = 200; m < m_ln_rho; ++m) printf("%e\n", ln_rho_eos[m]);
        printf("ln(T)\n");
        for(int n = 40; n < n_ln_t; ++n) printf("%e\n", ln_t_eos[n]);
        printf("Y_e\n");
        for(int o = 0; o < o_y_e; ++o) printf("%e\n", y_e_eos[o]);
    }

    void allocate()
    {
        int s = size();

        ln_rho_eos = new double[m_ln_rho];
        ln_t_eos = new double[n_ln_t];
        y_e_eos = new double[o_y_e];
        p_eos = new double[s];
        e_eos = new double[s];
        entropy_eos = new double[s];
        dp_drho_eos = new double[s];
        dp_deps_eos = new double[s];
        c_sound_squared_newtonian_eos = new double[s];
        mub_eos = new double[s];
        muq_eos = new double[s];
        mul_eos = new double[s];
        xp_eos = new double[s];
        xn_eos = new double[s];
        xa_eos = new double[s];
        xheavy_eos = new double[s];
        aheavy_eos = new double[s];
        zheavy_eos = new double[s];
        eflg_eos = new int[s];

        thermo_eos = new double[s * 8];
        yi_eos = new double[s * 3];
    }
};
