#include "common.h"

struct nucleus_capture_data
{
    int A, Z;
    double Q;
    double T;
    double nb_Y; // nb*Y
    double mu_e;
    double lambda;
};

std::vector<nucleus_capture_data *> nucleus_capture;

int read_nucleus_capture(const char *path, std::vector<nucleus_capture_data *> &nucleus_capture)
{
    std::ifstream infile(path);
    double val;
    int count = 0;
    std::string line;
    std::istringstream iss;
    while (std::getline(infile, line)) {
        iss.clear();
        iss.str(line);
        nucleus_capture_data *nc = new nucleus_capture_data();
        iss >> nc->A;
        iss >> nc->Z;
        iss >> nc->Q;
        iss >> nc->T;
        nc->T /= 11.6594202899;
        iss >> nc->nb_Y;
        nc->nb_Y = pow(10., nc->nb_Y)/1.674e15;
        iss >> nc->mu_e;
        nc->mu_e += M_ELECTRON;
        for(int k = 0; k < 1; ++k) iss >> val; 
        iss >> nc->lambda;
        nc->lambda = pow(10., nc->lambda);
        ++count;
        nucleus_capture.push_back(nc);
    }
    return count;
}

struct thermo_state
{
    double T, nb, Y, mu_e, lambda;

};

std::vector<thermo_state *> thermo_states[2];

int read_thermo_states(const char *path, std::vector<thermo_state *> &thermo_states)
{
    std::ifstream infile(path);
    double val;
    int count = 0;
    std::string line;
    std::istringstream iss;
    while (std::getline(infile, line)) {
        iss.clear();
        iss.str(line);
        thermo_state *ts = new thermo_state();
        iss >> ts->nb;
        ts->nb /= (1.674e15);
        iss >> ts->Y;
        iss >> val;
        iss >> ts->T;
        iss >> ts->mu_e;
        for(int k = 0; k < 8; ++k) iss >> val; 
        iss >> ts->lambda;
        ++count;
        thermo_states.push_back(ts);
    }
    return count;
}

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
    entries = table.size();
    std::cout << "Read " << entries << " EOS entries\n";
    table.dump();

    // Perform tests

    // electron potential
    FILE *fp = fopen("potential.res", "w+");
    for(int i = 0; i < 1000; ++i)
    {
        double n = pow(10., -8+7.*double(i)/1000.);
        fprintf(fp, "%e %e\n", n*1.674e-24*1e39, electron_potential(n, 1.));
    }
    fclose(fp);

    /* single nucleus rates */
 
    // read single rates
    read_nucleus_capture("data/single_rates", nucleus_capture);

    fp = fopen("output/single_rates.res", "w+");
    for(int i = 0; i < nucleus_capture.size(); ++i)
    {
        nucleus_capture_data *nc = nucleus_capture[i];
        if(nc->A < 2) continue;
        if(nc->T < 0.1) continue;
        if(nc->nb_Y < 1e-7) continue;

        fprintf(fp, "%d %d %e %e %e %e %e\n", nc->A, nc->Z, nc->T, nc->nb_Y, nc->lambda,
        electron_capture_fit(nc->A, nc->Z, nc->T, nc->mu_e), electron_capture_fit(nc->A, nc->Z, nc->T, electron_potential(nc->nb_Y, 1.)));
    }
    fclose(fp);

    fp = fopen("output/single_rates_Q.res", "w+");

    const double Q[2] = { -18, 16 };
    for(int i = 0; i < 1000; ++i)
    {
        double _Q = Q[0]+(Q[1]-Q[0])*double(i)/1000.;
        double mu_e[2] = { electron_potential(1.32e-6, 0.447), electron_potential(1.12e-4, 0.361) };
        fprintf(fp, "%f %e %e\n", _Q, electron_capture_fit(39, 21, 0.68, mu_e[0], _Q), electron_capture_fit(39, 21, 1.3, mu_e[1], _Q));
    }
    fclose(fp);

    /* total rates check */

    // read CCSN trajectories
    read_thermo_states("data/trajectory_15", thermo_states[0]);
    read_thermo_states("data/trajectory_25", thermo_states[1]);

    for(int k = 0; k < 2; ++k) {
        fp = fopen(k == 0 ? "output/total_rates_15.res" : "output/total_rates_25.res", "w+");
        for(int i = 0; i < thermo_states[k].size(); ++i)
        {
            thermo_state *ts = thermo_states[k][i];
            double rate = 0, rate_corr = 0;
            double mu_e = electron_potential(ts->nb, ts->Y);
            double total_abundance = 0, abundance_error = 0;
            double conditions[3] = {ts->T, ts->nb, ts->Y};

            for(int A = 2; A < 250; ++A) for(int Z = 0; Z <= A; ++Z) if(A >= 20)
            {
                //if(abs(A-Z) > 50) continue;
                double vl = 0, vh = 0;
                double abundance = element_abundance(A, Z, conditions, &vl, &vh);
                total_abundance += abundance;
                abundance_error += (vh-vl)*(vh-vl);
                rate += abundance * electron_capture_fit(A, Z, ts->T, mu_e);
                rate_corr += abundance * electron_capture_fit(A, Z, ts->T, ts->mu_e);
            }
            fprintf(fp, "%e %e %e %e %e %e %e %e %e\n", ts->T, ts->nb, ts->Y, ts->lambda, rate, rate/total_abundance, rate * sqrt(abundance_error)/total_abundance / total_abundance);
        }
        fclose(fp);
    }


    // poor man's multithreading
    #pragma omp parallel for
    for(int m = 200; m < table.m_ln_rho; ++m)
    for(int n = 40; n < table.n_ln_t; ++n)
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

    for(int A = 2; A < 295; ++A) for(int Z = 0; Z <= A; ++Z)
    {
        double vl = 0, vh = 0;
        double abundance = element_abundance(A, Z, conditions);
            total_abundance += abundance;
        if(abundance > 1e-30)
            {
                rate += abundance * electron_capture_fit(A, Z, T, mu_e);
                table.scattering_xs_eos[i] += abundance * nucleus_scattering_cross_section(A, Z, eps_mu, abundance*nb);
            }
    }

        printf("%d %e %e %e %e %e %e %e %e\n", i, T, nb, Y_e, mu_nu, ec_tab, rate, table.scattering_xs_eos[i]*1e36, total_abundance);
    }

    return 0;
}
