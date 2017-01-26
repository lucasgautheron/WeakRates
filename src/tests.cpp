#include "common.h"

int main(int argc, const char *argv[])
{
    if(argc < 2)
    {
        std::cout << "Missing compose model! Can't run any further.\n";
        return 1;
    }
    
    gsl_init();

    // read EOS table
    // read abundances
    // browse EOS table and compute rates

    int entries = 0;

    // Read nuclear data (masses)
    entries = read_nuclear_data("data/mass.mas12");
    std::cout << "Read " << entries << " nuclear data entries\n";
    
    std::cout << "Computing Q values...\n";
    
    FILE *fp_nuclear = fopen("output/compare_qvalues.res", "w+");
    int prev_A = 1;
    for(auto it = nuclear_table.begin(); it != nuclear_table.end(); it++) {
        if(prev_A != it->second->A)
        {
            fprintf(fp_nuclear, "\n");
            prev_A = it->second->A;
        }
        std::array<int, 2> daughter = { it->second->A, it->second->Z+1 };
        fprintf(fp_nuclear, "%d %d %e %e %e %e %e %e %e %e %e %e\n", it->second->A, it->second->Z, it->second->Z * M_PROTON + (it->second->A-it->second->Z)*M_NEUTRON,
                it->second->m, mass_formula(it->second->A, it->second->Z), mass_dz10(it->second->A, it->second->Z), mass_SEMF(it->second->A, it->second->Z),
                it->second->beta_q, nuclear_table.count(daughter) ? it->second->m - nuclear_table[daughter]->m : 0,
                mass_formula(it->second->A, it->second->Z)-mass_formula(it->second->A, it->second->Z+1),
                mass_dz10(it->second->A, it->second->Z)-mass_dz10(it->second->A, it->second->Z+1),
                mass_SEMF(it->second->A, it->second->Z)-mass_SEMF(it->second->A, it->second->Z+1));
    }
    fclose(fp_nuclear);
    
    std::cout << "Computing test integrals...\n";
    
    FILE *fp_integrals =  fopen("output/compare_integrals.res", "w+");
    for(int i = 0; i < 100; ++i) for(int j = 0; j < 100; ++j) for(int k = 0; k < 50; ++k)
    {
        double mu_e = 10*double(i)/double(100), mu_nu = 100*double(j)/double(100), Q = -5+10*double(k)/50.;
        double rate = 24*gsl_sf_fermi_dirac_int (4,Q+mu_e) - 2 * Q * 6*gsl_sf_fermi_dirac_int (3,Q+mu_e) + Q*Q * 2*gsl_sf_fermi_dirac_int (2, Q+mu_e);
    
        fprintf(fp_integrals, "%.3f %.3f %.3f %e %e\n", mu_e, mu_nu, Q, electron_capture_ps(1., mu_e, mu_nu, Q), rate);
    }
    fclose(fp_integrals);

    // Read abundance data
    std::cout << "Reading compose data from " << argv[1] << "...\n";
    entries = read_abundance_data(argv[1]);
    std::cout << "Read " << entries << " nuclear abundance data entries\n";

    // Read full EOS data
    full_EOS_table full_table;
    full_table.read("data/eosls220_low.h5", EOS_TYPE_LOW);
    //full_table.read("data/eoscompose.h5", EOS_TYPE_COMPOSE);

    // Read rate EOS data
    short_EOS_table rates_table;
    rates_table.read("data/elec_capt_rate_ls220.h5", 0);

    // Read output rates
    short_EOS_table output_table;
    output_table.read("output/sigma_scattering_rate.h5", 0);
    

    entries = rates_table.size();
    std::cout << "Read " << entries << " EOS entries (" 
        << rates_table.m_ln_rho << "x" << rates_table.n_ln_t << "x" << rates_table.o_y_e << ")\n";

    rates_table.dump();

    entries = full_table.size();
    std::cout << "Read " << entries << " EOS entries\n";
    full_table.dump();
    

    FILE *fp_scattering = fopen("output/compare_neutrinos.res", "w+");
    FILE *fp_capture = fopen("output/compare_capture.res", "w+");
    FILE *fp_nuclei = fopen("output/compare_nuclei.res", "w+");
    for(int m = 0; m < rates_table.m_ln_rho; ++m)
    for(int n = 0; n < rates_table.n_ln_t; ++n)
    for(int o = 0; o < rates_table.o_y_e; ++o)
    for(int p = 0; p < rates_table.p_mu; ++p)
    {
        int i = o*(rates_table.n_ln_t*rates_table.m_ln_rho*rates_table.p_mu) + n*(rates_table.m_ln_rho*rates_table.p_mu) + m*rates_table.p_mu + p;
        int ii = m + n * rates_table.m_ln_rho + o * rates_table.m_ln_rho * rates_table.n_ln_t;
	
        const double mu_e = -full_table.muq_eos[ii] + full_table.mul_eos[ii], //- 0.51099891;
	     mu_nu = full_table.mul_eos[ii],
	     mu_neut = full_table.mub_eos[ii], //- 939.56536; // non relativistic mu
	     mu_p = full_table.mub_eos[ii] + full_table.muq_eos[ii];

        const double nb = exp(rates_table.ln_rho_eos[m]), T = exp(rates_table.ln_t_eos[n]), Y_e = rates_table.y_e_eos[o], ec_tab = rates_table.elec_rate_tab_eos[i];

        const double mu_nu_eff = mu_nu * rates_table.mu_nu_eos[p],
                     eta = -mu_nu_eff / T;

        //const double eta_pn = eta_pn_v3(mu_neut, mu_p, T);

        double conditions[3] = {T, nb, Y_e};
        
        double total_abundance = 0;
        std::vector<element> elements;
        get_abundances(conditions, elements);
        
        double aheavy = 0, zheavy = 0;
        for(int j = 0; j < elements.size(); ++j)
	{
            int A = elements[j].A, Z = elements[j].Z;
	    double abundance = elements[j].abundance;

            
            if(A < 4 || Z < 2 || abundance <= 0) continue;
            aheavy += A*abundance;
            zheavy += Z*abundance;

            total_abundance += abundance;
            
	}
	if(total_abundance > 1e-30)
        {
            aheavy /= total_abundance;
            zheavy /= total_abundance;
        }

        double /*mu_e = degenerate_potential(M_ELECTRON, nb*Y_e),*/ eps_mu = average_neutrino_energy(T, mu_nu_eff);

        if(T > 0.01 && T < 5 && nb > 1e-8 && nb < 1e-2 && Y_e > 0.1)
        {
            fprintf(fp_capture, "%.3f %e %.3f %e %e %e %e %.3f\n", T, nb, Y_e, mu_nu_eff, rates_table.elec_rate_tab_eos[i], output_table.elec_rate_fast_eos[i], output_table.elec_rate_tab_eos[i], (float)shell_capt_factor(aheavy, zheavy));
            fprintf(fp_scattering, "%.3f %e %.3f %e %e %e %e %e %e %e\n", T, nb, Y_e, mu_nu_eff, aheavy, zheavy, full_table.aheavy_eos[ii], full_table.zheavy_eos[ii], output_table.scattering_xs_nu_eos[i], output_table.scattering_xs_nu_sna_eos[i]);
            fprintf(fp_nuclei, "%.3f %.3f %.3f %.3f %e %e %d\n", aheavy, zheavy, full_table.aheavy_eos[ii], full_table.zheavy_eos[ii], total_abundance, full_table.xheavy_eos[ii], elements.size());
        }

        //printf("%e %e %e\n", T, mu_e, degenerate_potential(M_ELECTRON, nb*Y_e));
        //continue;

    }
    fclose(fp_scattering);
    fclose(fp_capture);
    fclose(fp_nuclei);

    return 0;
}
