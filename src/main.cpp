#include "common.h"

int main(int argc, const char *argv[])
{
    // read EOS table
    // read abundances
    // browse EOS table and compute rates

    int entries = 0;

    // Read nuclear data (masses)
    entries = read_nuclear_data("data/mass.mas12");
    std::cout << "Read " << entries << " nuclear data entries\n";

    // Read abundance data
    entries = read_abundance_data("data/abundances/1/eos.compo", "data/abundances/1/eos.t", "data/abundances/1/eos.nb", "data/abundances/1/eos.yq");
    std::cout << "Read " << entries << " nuclear abundance data entries\n";

    // Read rate EOS data
    short_EOS_table rates_table;
    int error;
    rates_table.read("data/elec_capt_rate_ls220.h5", &error);
    rates_table.write("output/sigma_scattering_rate.h5", &error);

    entries = rates_table.size();
    std::cout << "Read " << entries << " EOS entries\n";

    
    rates_table.dump();

    // Read full EOS data
    full_EOS_table full_table;
    full_table.read("data/eosls220_low.h5", &error);

    entries = full_table.size();
    std::cout << "Read " << entries << " EOS entries\n";
    full_table.dump();

    // Perform calculations

    int processed = 0;
    // poor man's multithreading
    #pragma omp parallel for
    for(int m = 0; m < rates_table.m_ln_rho; ++m)
    for(int n = 0; n < rates_table.n_ln_t; ++n)
    for(int o = 0; o < rates_table.o_y_e; ++o)
    for(int p = 0; p < rates_table.p_mu; ++p)
    {
        //int i = m*(table.n_ln_t*table.o_y_e*table.p_mu) + n*(table.o_y_e*table.p_mu) + o*table.p_mu + p;
        int i = o*(rates_table.n_ln_t*rates_table.m_ln_rho*rates_table.p_mu) + n*(rates_table.m_ln_rho*rates_table.p_mu) + m*rates_table.p_mu + p;
        int ii = m + n * rates_table.m_ln_rho + o * rates_table.m_ln_rho * rates_table.n_ln_t;
	
        const double mu_e = -full_table.muq_eos[ii] + full_table.mul_eos[ii], //- 0.51099891;
	     mu_nu = full_table.mul_eos[ii],
	     mu_neut = full_table.mub_eos[ii], //- 939.56536; // non relativistic mu
	     mu_p = full_table.mub_eos[ii] + full_table.muq_eos[ii];

        const double mu_nu_eff = mu_nu * rates_table.mu_nu_eos[p];

        const double nb = exp(rates_table.ln_rho_eos[m]), T = exp(rates_table.ln_t_eos[n]), Y_e = rates_table.y_e_eos[o], ec_tab = rates_table.elec_rate_tab_eos[i];

        const double eta_pn = eta_pn_v3(mu_neut, mu_p, T);

        double conditions[3] = {T, nb, Y_e};

        double /*mu_e = degenerate_potential(M_ELECTRON, nb*Y_e),*/ eps_mu = average_neutrino_energy(T, mu_nu_eff);

        //printf("%e %e %e\n", T, mu_e, degenerate_potential(M_ELECTRON, nb*Y_e));
        //continue;

        rates_table.elec_rate_single_eos[i] = 0;
        rates_table.elec_rate_fast_eos[i] = rates_table.elec_rate_tab_eos[i] = rates_table.scattering_xs_eos[i] = 0;

        double total_abundance = 0;
        std::vector<element> elements;
        get_abundances(conditions, elements);

        for(int j = 0; j < elements.size(); ++j)
	{
            int A = elements[j].A, Z = elements[j].Z;
	    double abundance = elements[j].abundance;
            
            if(A < 2 || Z < 2) continue;

            total_abundance += abundance;
            
            rates_table.elec_rate_fast_eos[i] += abundance * electron_capture_fit(A, Z, T, degenerate_potential(M_ELECTRON, nb*Y_e));
            rates_table.scattering_xs_eos[i] += abundance * nucleus_scattering_cross_section(A, Z, eps_mu, abundance*nb);
            rates_table.elec_rate_tab_eos[i] += elec_capt_heavy_nuclei_effective(mu_e, mu_nu_eff, abundance, T, mu_neut, mu_p, Z, A);
	}

        const double rb = elec_capt_proton_effective(mu_e, mu_nu_eff, T, mu_neut, mu_p, full_table.xp_eos[ii], full_table.xn_eos[ii], eta_pn);
        rates_table.elec_rate_tab_eos[i] += rb/nb;
        rates_table.elec_rate_single_eos[i] += rb/nb;
        double rheavy = 0; 
		 
	if ((full_table.xheavy_eos[ii] > 0.) && (full_table.aheavy_eos[ii] > 0.))
	{
	    // ott table
	    //xheavy_eos[index] = xheavy_eos[index] / aheavy_eos[index];
	    rates_table.elec_rate_single_eos[i] += elec_capt_heavy_nuclei_effective(mu_e, mu_nu_eff, full_table.xheavy_eos[ii], T, mu_neut, mu_p, full_table.zheavy_eos[ii], full_table.aheavy_eos[ii]);
	} 
	 /* eta_np en fm-3 compens\E9 par rho en fm-3
	  facteur 1e39 dans calc_rate compens\E9 ici
	  tout le reste en MeV, cm, s jusqu'ici */
	rates_table.elec_rate_tab_eos[i] *= 1e-39 * INV_COCO_TIME; 
        rates_table.elec_rate_single_eos[i] *= 1e-39 * INV_COCO_TIME;
        
        rates_table.elec_rate_fast_eos[i] += electron_capture_proton(T, nb, degenerate_potential(M_ELECTRON, nb*Y_e), mu_nu, eta_pn);
        rates_table.elec_rate_fast_eos[i] *= INV_COCO_TIME;

        ++processed;
        printf("%d / %d\n", processed, rates_table.size());
    }

    rates_table.write("output/sigma_scattering_rate.h5", &error);

    return 0;
}
