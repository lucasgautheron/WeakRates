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

    // Read EOS data
    EOS_table table;
    int error;
    read_EOS_table("data/elec_capt_rate_ls220.h5", table, &error);
    write_EOS_table("output/sigma_scattering_rate.h5", table, &error);

    entries = table.size();
    std::cout << "Read " << entries << " EOS entries\n";
    table.dump();

    // Perform calculations

    int processed = 0;
    // poor man's multithreading
    #pragma omp parallel for
    for(int m = 0; m < table.m_ln_rho; ++m)
    for(int n = 0; n < table.n_ln_t; ++n)
    for(int o = 0; o < table.o_y_e; ++o)
    for(int p = 0; p < table.p_mu; ++p)
    {
        //int i = m*(table.n_ln_t*table.o_y_e*table.p_mu) + n*(table.o_y_e*table.p_mu) + o*table.p_mu + p;
        int i = o*(table.n_ln_t*table.m_ln_rho*table.p_mu) + n*(table.m_ln_rho*table.p_mu) + m*table.p_mu + p;
        const double nb = exp(table.ln_rho_eos[m]), T = exp(table.ln_t_eos[n]), Y_e = table.y_e_eos[o], mu_nu = table.mu_nu_eos[p], ec_tab = table.elec_rate_tab_eos[i];

        double conditions[3] = {T, nb, Y_e};

        double mu_e = degenerate_potential(M_ELECTRON, nb*Y_e), eps_mu = average_neutrino_energy(T, mu_nu);

        table.elec_rate_tab_eos[i] = table.scattering_xs_eos[i] = 0;

        double total_abundance = 0;
        std::vector<element> elements;
        get_abundances(conditions, elements);

        double n_n = 0, n_p = 0;

        for(int j = 0; j < elements.size(); ++j)
	{
            int A = elements[j].A, Z = elements[j].Z;
	    double abundance = elements[j].abundance;
            if(A==1)
            {
                if(Z==0) n_n = abundance;
                else n_p = abundance;
            }
            
            if(A < 2 || Z < 2) continue;

            total_abundance += abundance;
            
            table.elec_rate_tab_eos[i] += abundance * electron_capture_fit(A, Z, T, mu_e);
            table.scattering_xs_eos[i] += abundance * nucleus_scattering_cross_section(A, Z, eps_mu, abundance*nb);
	}
        
        table.elec_rate_tab_eos[i] += electron_capture_proton(T, nb, mu_e, mu_nu, n_p, n_n);
    }

    write_EOS_table("output/sigma_scattering_rate.h5", table, &error);

    return 0;
}
