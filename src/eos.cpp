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

void write_EOS_table(const char *path, EOS_table &table, int *error)
{
  hid_t file;
  herr_t status;
  hid_t dataset, dataspace;

  // open file

  file = H5Fcreate(path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


  // write data sets

  hsize_t one = 1;
  dataspace = H5Screate_simple(1, &one, NULL);
  dataset = H5Dcreate2(file, "pointsrho", H5T_NATIVE_INT, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &table.m_ln_rho);
  H5Dclose (dataset);
  H5Sclose(dataspace);

  dataspace = H5Screate_simple(1, &one, NULL);
  dataset = H5Dcreate2(file, "pointst", H5T_NATIVE_INT, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &table.n_ln_t);
  H5Dclose (dataset);
  H5Sclose(dataspace);


  dataspace = H5Screate_simple(1, &one, NULL);
  dataset = H5Dcreate2(file, "pointsye", H5T_NATIVE_INT, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &table.o_y_e);
  H5Dclose (dataset);
  H5Sclose(dataspace);

  dataspace = H5Screate_simple(1, &one, NULL);
  dataset = H5Dcreate2(file, "pointsmu", H5T_NATIVE_INT, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &table.p_mu);
  H5Dclose (dataset);
  H5Sclose(dataspace);


  hsize_t sz[1] = { table.m_ln_rho };
  dataspace = H5Screate_simple(1, sz, NULL);
  dataset = H5Dcreate2(file, "logrho", H5T_NATIVE_DOUBLE, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, table.ln_rho_eos);
  H5Dclose (dataset);
  H5Sclose(dataspace);


  sz[0] = table.n_ln_t;
  dataspace = H5Screate_simple(1, sz, NULL);
  dataset = H5Dcreate2(file, "logt", H5T_NATIVE_DOUBLE, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, table.ln_t_eos);
  H5Dclose (dataset);
  H5Sclose(dataspace);

  sz[0] = table.o_y_e;
  dataspace = H5Screate_simple(1, sz, NULL);
  dataset = H5Dcreate2(file, "ye", H5T_NATIVE_DOUBLE, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, table.y_e_eos);
  H5Dclose (dataset);
  H5Sclose(dataspace);

  sz[0] = table.p_mu;
  dataspace = H5Screate_simple(1, sz, NULL);
  dataset = H5Dcreate2(file, "mu_nu", H5T_NATIVE_DOUBLE, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, table.mu_nu_eos);
  H5Dclose (dataset);
  H5Sclose(dataspace);

  hsize_t dims[4] = { table.m_ln_rho, table.n_ln_t, table.o_y_e, table.p_mu };

  sz[0] = table.size();
  dataspace = H5Screate_simple(1, sz, NULL);
  dataset = H5Dcreate2(file, "elec_capt_rate", H5T_NATIVE_DOUBLE, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, table.elec_rate_tab_eos);
  H5Dclose (dataset);


  dataset = H5Dcreate2(file, "scattering_xs_eos", H5T_NATIVE_DOUBLE, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, table.scattering_xs_eos);
  H5Dclose (dataset);
  H5Sclose(dataspace);

  // close file
  status = H5Fclose (file);
}
