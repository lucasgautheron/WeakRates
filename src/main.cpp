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
    entries = read_abundance_data("data/EOS.compo", "data/EOS.t", "data/EOS.nb", "data/EOS.yq");
    std::cout << "Read " << entries << " nuclear abundance data entries\n";

    // Read EOS data
    EOS_table table;
    int error;
    read_EOS_table("data/elec_capt_rate_ls220.h5", table, &error);
    write_EOS_table("output/table.h5", table, &error);

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
        int i = m*(table.n_ln_t*table.o_y_e*table.p_mu) + n*(table.o_y_e*table.p_mu) + o*table.p_mu + p;
        const double nb = exp(table.ln_rho_eos[m]), T = exp(table.ln_t_eos[n]), Y_e = table.y_e_eos[o], mu_nu = table.mu_nu_eos[p], ec_tab = table.elec_rate_tab_eos[i];
        //if(nb >= 0.1) continue;
        double conditions[3] = {T, nb, Y_e};

        double mu_e = electron_potential(nb, Y_e), eps_mu = average_neutrino_energy(T, mu_nu);


        double rate = 0;
        double total_abundance = 0;
        table.scattering_xs_eos[i] = 0;
        std::vector<element> elements;
        get_abundances(conditions, elements);
        table.elec_rate_tab_eos[i] = table.scattering_xs_eos[i] = 0;

        for(int j = 0; j < elements.size(); ++j)
	{
            int A = elements[j].A, Z = elements[j].Z;
	    double abundance = elements[j].abundance;
            total_abundance += abundance;
            if(A < 2 || Z < 2) continue;
            
	    if(abundance > 1e-30)
            {
                //printf("%e (%e %e %e %e)\n", abundance, elements[j].v[0], elements[j].v[1], elements[j].v[2], elements[j].v[3]);
                table.elec_rate_tab_eos[i] += abundance * electron_capture_fit(A, Z, T, mu_e);
                table.scattering_xs_eos[i] += abundance * nucleus_scattering_cross_section(A, Z, eps_mu, abundance*nb);
            }
	}

        //if(total_abundance > 1e-15 && nb > 1e-4 && nb < 1e-2) printf("%d %d %e %e %e %e %e %e %e %e\n", i, elements.size(), T, nb, Y_e, mu_nu, ec_tab, table.elec_rate_tab_eos[i], table.scattering_xs_eos[i]*1e36, total_abundance);
        ++processed;
        if(processed % 10000 == 0) std::cout << "Processed " << processed << " states out of " << table.size() << std::endl;
        
    }

    write_EOS_table("output/table.h5", table, &error);

    return 0;
}
