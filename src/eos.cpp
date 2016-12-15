#include "common.h"

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
  table.ln_rho_eos = new double[table.m_ln_rho];
  table.ln_t_eos = new double[table.n_ln_t];
  table.y_e_eos = new double[table.o_y_e];
  table.mu_nu_eos = new double[table.p_mu];
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
