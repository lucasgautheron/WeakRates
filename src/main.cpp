#include "common.h"

struct EOS_table
{
    int m_ln_rho, n_ln_t, o_y_e, p_mu;
    double *ln_rho_eos,*ln_t_eos, *y_e_eos, *mu_nu_eos, *elec_rate_tab_eos;
    double *scattering_xs_eos;
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

    printf("%e %e %e %e\n", table.ln_rho_eos[0], table.ln_rho_eos[1], table.ln_t_eos[0], table.ln_t_eos[1]);

    /*for(int i = 0; i < table.m_ln_rho; ++i)
    for(int j = 0; j < table.n_ln_t; ++j)
    for(int k = 0; k < table.o_y_e; ++k)
    for(int l = 0; l < table.p_mu; ++l)
    {
        //int n = i*(table.n_ln_t*table.o_y_e*table.p_mu) + j*(table.o_y_e*table.p_mu) + k*table.p_mu + l;
        int n = i + j*(table.m_ln_rho) + k*table.m_ln_rho*table.n_ln_t + l*table.m_ln_rho*table.n_ln_t*table.o_y_e;
        printf("%e\n", table.ln_rho_eos[n]);
        break;
    }*/

    nucleus_scattering_cross_section(1, 1, 10, 0);

    return 0;
}

void read_EOS_table(const char *path, EOS_table &table, int *error)
{
  hid_t file;
  herr_t status;
  hid_t dataset;

  int m_ln_rho_check, n_ln_t_check, o_y_e_check, p_mu_check;
  
  std::cout << std::endl << table.m_ln_rho << "  " << table.n_ln_t << "  " << table.o_y_e << "  " << table.p_mu << std::endl;

  // open file

  file = H5Fopen (path, H5F_ACC_RDONLY, H5P_DEFAULT);


  // read data sets
  dataset = H5Dopen (file, "pointsrho");
  status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_ln_rho_check);
  H5Dclose (dataset);  
  
  dataset = H5Dopen (file, "pointst");
  status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &n_ln_t_check);
  H5Dclose (dataset);  

  dataset = H5Dopen (file, "pointsye");
  status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &o_y_e_check);
  H5Dclose (dataset);  

  dataset = H5Dopen (file, "pointsmu");
  status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &p_mu_check);
  H5Dclose (dataset);

  std::cout << std::endl << m_ln_rho_check << "  " << n_ln_t_check << "  " << o_y_e_check << "  " << p_mu_check << std::endl;

  table.m_ln_rho = m_ln_rho_check;
  table.n_ln_t = n_ln_t_check;
  table.o_y_e = o_y_e_check;
  table.p_mu = p_mu_check;

  int size = m_ln_rho_check * n_ln_t_check * o_y_e_check * p_mu_check;
  table.ln_rho_eos = new double[size+1];
  table.ln_t_eos = new double[size+1];
  table.y_e_eos = new double[size+1];
  table.mu_nu_eos = new double[size+1];
  table.elec_rate_tab_eos = new double[size+1];
  table.scattering_xs_eos = new double[size+1];

    if (table.m_ln_rho == m_ln_rho_check && table.n_ln_t == n_ln_t_check && table.o_y_e == o_y_e_check && table.p_mu == p_mu_check)
    {
      *error = 0;
      
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

    }
  else
    {
      std::cout << "Error while reading EOS table: array size mismatch" << std::endl;
      *error = 1;
    }
  
  // close file
  status = H5Fclose (file);
}
