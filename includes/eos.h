#define ALLOC(arr, type, size) arr = new type[size]; memset(arr, 0, size*sizeof(type));

enum { EOS_TYPE_LOW = 0, EOS_TYPE_COMPOSE };

struct EOS_table
{
    int m_ln_rho, n_ln_t, o_y_e;
    double *ln_rho_eos,*ln_t_eos, *y_e_eos;
    virtual void dump() {}
    virtual const int size() { return 0; }
    virtual void allocate() {}

    virtual int read(const char *path, int type) {}
    virtual int write(const char *path) {}
};


struct full_EOS_table;

struct short_EOS_table : EOS_table
{
    int p_mu;
    double *mu_nu_eos, *elec_rate_tab_eos, *elec_rate_single_eos, *elec_rate_fast_eos;
    double *scattering_xs_nu_eos, *scattering_xs_nu_bar_eos, *scattering_xs_nu_x_eos;
#ifdef DEBUG
    double *scattering_xs_nu_sna_eos;
#endif

    int read(const char *path, int type);
    int write(const char *path);

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

    void bind_full_table(full_EOS_table &t);

    const int size()
    {
        return m_ln_rho*n_ln_t*o_y_e*p_mu;
    }

    void allocate()
    {
        int s = size();
        ALLOC(ln_rho_eos, double, m_ln_rho);
        ALLOC(ln_t_eos, double, n_ln_t);
        ALLOC(y_e_eos, double, o_y_e);
        ALLOC(mu_nu_eos, double, p_mu);
        ALLOC(elec_rate_tab_eos, double, s);
        ALLOC(elec_rate_fast_eos, double, s);
        ALLOC(elec_rate_single_eos, double, s);
        ALLOC(scattering_xs_nu_eos, double, s);
        ALLOC(scattering_xs_nu_bar_eos, double, s);
        ALLOC(scattering_xs_nu_x_eos, double, s);
#ifdef DEBUG
        ALLOC(scattering_xs_nu_sna_eos, double, s);
#endif
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

    int read(const char *path, int type = EOS_TYPE_COMPOSE);

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

        ALLOC(ln_rho_eos, double, m_ln_rho);
        ALLOC(ln_t_eos, double, n_ln_t);
        ALLOC(y_e_eos, double, o_y_e);
        ALLOC(p_eos, double, s);
        ALLOC(e_eos, double, s);
        ALLOC(entropy_eos, double, s);
        ALLOC(dp_drho_eos, double, s);
        ALLOC(dp_deps_eos, double, s);
        ALLOC(c_sound_squared_newtonian_eos, double, s);
        ALLOC(mub_eos, double, s);
        ALLOC(muq_eos, double, s);
        ALLOC(mul_eos, double, s);
        ALLOC(xp_eos, double, s);
        ALLOC(xn_eos, double, s);
        ALLOC(xa_eos, double, s);
        ALLOC(xheavy_eos, double, s);
        ALLOC(aheavy_eos, double, s);
        ALLOC(zheavy_eos, double, s);
        ALLOC(eflg_eos, int, s);

        ALLOC(thermo_eos, double, s * 8);
        ALLOC(yi_eos, double, s * 3);
    }
};

void copyfile(const char *source, const char *destination);
void compile_compose_data(const char *path);
