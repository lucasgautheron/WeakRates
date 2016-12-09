#include "common.h"

void read_ls_eos_table_hdf_cpp (const char *table, int *m_ln_rho, int *n_ln_t, int *o_y_e,
			      double *ln_rho_eos, double *ln_t_eos, double *y_e_eos,
			      double *p_eos, double *e_eos,
			      double *entropy_eos, double *mu_nu_eos,
			      double *dp_drho_eos, double *dp_deps_eos,
			      double *c_sound_squared_newtonian_eos,
				int * eflag, int *error, int* kflag)
{
  hid_t file;
  herr_t status;
  hid_t dataset;

  int m_ln_rho_check, n_ln_t_check, o_y_e_check;
  
  
  // open file
  if(*kflag ==1) {
    file = H5Fopen ("EOS/eosls180.h5", H5F_ACC_RDONLY, H5P_DEFAULT);}
  if(*kflag ==2) {
    file = H5Fopen ("EOS/eosls220.h5", H5F_ACC_RDONLY, H5P_DEFAULT);}
  if(*kflag ==3) {
    file = H5Fopen ("EOS/eosls375.h5", H5F_ACC_RDONLY, H5P_DEFAULT);}

   
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

  if (*m_ln_rho == m_ln_rho_check && *n_ln_t == n_ln_t_check && *o_y_e == o_y_e_check)
    {
      *error = 0;
      
      dataset = H5Dopen (file, "logrho");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ln_rho_eos);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "logt");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ln_t_eos);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "ye");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, y_e_eos);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "pressure");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, p_eos);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "energy");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, e_eos);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "entropy");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, entropy_eos);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "munu");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mu_nu_eos);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "dpdrhoe");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dp_drho_eos);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "dpderho");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dp_deps_eos);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "cs2");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, c_sound_squared_newtonian_eos);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "eosflag");
      status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, eflag);
      H5Dclose (dataset);
    }
  else
    {
      *error = 1;
    }
  
  
  
  // close file
  status = H5Fclose (file);
}


nuclear_array nuclear_table;

int main(int argc, const char *argv[])
{
    // read EOS table
    // read abundances
    // browse EOS table and compute rates

    read_nuclear_data("data/mass.mas12", nuclear_table);
    nucleus_scattering_cross_section(1, 1, 10, 0);

    return 0;
}
