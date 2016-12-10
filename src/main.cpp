#include "common.h"

struct EOS_table
{
    int m_ln_rho, n_ln_t, o_y_e, p_mu;
    double *ln_rho_eos,*ln_t_eos, *y_e_eos, *mu_nu_eos, *elec_rate_tab_eos;
    double *scattering_xs_eos, *elec_rate_eos;

    const int size()
    {
        return m_ln_rho*n_ln_t*o_y_e*p_mu;
    }
};

void read_EOS_table(const char *path, EOS_table &table, int *error);

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

    // Perform calculations

    printf("%e %e %e %e\n", table.ln_rho_eos[0], table.ln_rho_eos[1], table.ln_t_eos[0], table.ln_t_eos[1]);

    // poor man's multithreading
    #pragma omp parallel for
    for(int i = 0; i < table.size(); ++i)
    {
        const double nb = exp(table.ln_rho_eos[i]), T = exp(table.ln_t_eos[i]), Y_e = table.y_e_eos[i], mu_nu = table.mu_nu_eos[i], ec_tab = table.elec_rate_tab_eos[i];
        if(nb >= 0.1) continue;
        double conditions[3] = {T, nb, Y_e};
        double mu_e = electron_potential(nb, Y_e);
        double rate = 0;

	for(int A = 2; A < 295; ++A) for(int Z = 0; Z <= A; ++Z)
	{
	    double vl = 0, vh = 0;
	    double abundance = element_abundance(A, Z, conditions);
	    rate += abundance * electron_capture_fit(A, Z, T, mu_e);
	}

        printf("%d %e %e %e %e %e %e\n", i, T, nb, Y_e, mu_nu, ec_tab, rate);
    }

    nucleus_scattering_cross_section(1, 1, 10, 0);

    return 0;
}

void read_EOS_table(const char *path, EOS_table &table, int *error)
{
  hid_t file;
  herr_t status;
  hid_t dataset;

  // open file

  file = H5Fopen (path, H5F_ACC_RDONLY, H5P_DEFAULT);


  // read data sets
  dataset = H5Dopen (file, "pointsrho");
  status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &table.m_ln_rho);
  H5Dclose (dataset);  
  
  dataset = H5Dopen (file, "pointst");
  status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &table.n_ln_t);
  H5Dclose (dataset);  

  dataset = H5Dopen (file, "pointsye");
  status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &table.o_y_e);
  H5Dclose (dataset);  

  dataset = H5Dopen (file, "pointsmu");
  status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &table.p_mu);
  H5Dclose (dataset);

  int size = table.size();
  table.ln_rho_eos = new double[size+1];
  table.ln_t_eos = new double[size+1];
  table.y_e_eos = new double[size+1];
  table.mu_nu_eos = new double[size+1];
  table.elec_rate_tab_eos = new double[size+1];
  table.scattering_xs_eos = new double[size+1];

  dataset = H5Dopen (file, "logrho");
  status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, table.ln_rho_eos);
  H5Dclose (dataset);

  dataset = H5Dopen (file, "logt");
  status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, table.ln_t_eos);
  H5Dclose (dataset);

  dataset = H5Dopen (file, "ye");
  status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, table.y_e_eos);
  H5Dclose (dataset);

  dataset = H5Dopen (file, "mu_nu");
  status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, table.mu_nu_eos);
  H5Dclose (dataset);

  dataset = H5Dopen (file, "elec_capt_rate");
  status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, table.elec_rate_tab_eos);
  H5Dclose (dataset);

  // close file
  status = H5Fclose (file);
}
