struct EOS_table
{
    int m_ln_rho, n_ln_t, o_y_e, p_mu;
    double *ln_rho_eos,*ln_t_eos, *y_e_eos, *mu_nu_eos, *elec_rate_tab_eos;
    double *scattering_xs_eos, *elec_rate_eos;

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
};

void read_EOS_table(const char *path, EOS_table &table, int *error);
