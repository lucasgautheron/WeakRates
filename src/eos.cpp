#include "common.h"

void short_EOS_table::bind_full_table(full_EOS_table &t)
{
    m_ln_rho = t.m_ln_rho;
    n_ln_t = t.n_ln_t;
    o_y_e = t.o_y_e;
    p_mu = 11;
    
    allocate();

    for(int m = 0; m < m_ln_rho; ++m) ln_rho_eos[m] = t.ln_rho_eos[m];
    for(int n = 0; n < n_ln_t; ++n) ln_t_eos[n] = t.ln_t_eos[n];
    for(int o = 0; o < o_y_e; ++o) y_e_eos[o] = t.y_e_eos[o];
    for(int p = 0; p < p_mu; ++p)
    {
        mu_nu_eos[p] = 0.1*double(p);
    }
}

int short_EOS_table::read(const char *path, int type)
{
  int error;
  hid_t file;
  herr_t status;
  hid_t dataset;

  // open file

  file = H5Fopen (path, H5F_ACC_RDONLY, H5P_DEFAULT);


  // read data sets
  dataset = H5Dopen (file, "pointsrho");
  status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->m_ln_rho);
  H5Dclose (dataset);  
  
  dataset = H5Dopen (file, "pointst");
  status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->n_ln_t);
  H5Dclose (dataset);  

  dataset = H5Dopen (file, "pointsye");
  status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->o_y_e);
  H5Dclose (dataset);  

  dataset = H5Dopen (file, "pointsmu");
  status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->p_mu);
  H5Dclose (dataset);

  this->allocate();

  dataset = H5Dopen (file, "logrho");
  status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->ln_rho_eos);
  H5Dclose (dataset);

  dataset = H5Dopen (file, "logt");
  status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->ln_t_eos);
  H5Dclose (dataset);

  dataset = H5Dopen (file, "ye");
  status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->y_e_eos);
  H5Dclose (dataset);

  dataset = H5Dopen (file, "mu_nu");
  status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->mu_nu_eos);
  H5Dclose (dataset);

  dataset = H5Dopen (file, "elec_capt_rate");
  status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->elec_rate_tab_eos);
  H5Dclose (dataset);

  // close file
  status = H5Fclose (file);
}

int full_EOS_table::read(const char *path, int type)
{
  int error;
  hid_t file;
  herr_t status;
  hid_t dataset;

  // open file

  file = H5Fopen (path, H5F_ACC_RDONLY, H5P_DEFAULT);

  // read data sets
  switch(type)
  {
      case EOS_TYPE_LOW:
      dataset = H5Dopen (file, "pointsrho");
      status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_ln_rho);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "pointst");
      status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &n_ln_t);
      H5Dclose (dataset);  

      dataset = H5Dopen (file, "pointsye");
      status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &o_y_e);
      H5Dclose (dataset);  
      break;


      case EOS_TYPE_COMPOSE:
      dataset = H5Dopen (file, "pointsnb");
      status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_ln_rho);
      H5Dclose (dataset);  
  
      dataset = H5Dopen (file, "pointst");
      status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &n_ln_t);
      H5Dclose (dataset);  

      dataset = H5Dopen (file, "pointsyq");
      status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &o_y_e);
      H5Dclose (dataset);
      break;


      default: return 1;
  }

 
  this->allocate();

  error = 0;
  switch(type)
  {
      case EOS_TYPE_LOW:
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

      dataset = H5Dopen (file, "mub");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mub_eos);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "muq");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, muq_eos);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "mul");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mul_eos);
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

      dataset = H5Dopen (file, "xp");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xp_eos);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "xn");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xn_eos);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "xalpha");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xa_eos);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "xheavy");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xheavy_eos);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "aheavy");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, aheavy_eos);
      H5Dclose (dataset);

      dataset = H5Dopen (file, "zheavy");
      status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, zheavy_eos);
      H5Dclose (dataset);
      dataset = H5Dopen (file, "eosflag");
      status = H5Dread (dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, eflg_eos);
  break;
  case EOS_TYPE_COMPOSE:
dataset = H5Dopen (file, "nb");
status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->ln_rho_eos);
H5Dclose (dataset);

dataset = H5Dopen (file, "t");
status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->ln_t_eos);
H5Dclose (dataset);

dataset = H5Dopen (file, "yq");
status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->y_e_eos);
H5Dclose (dataset);

dataset = H5Dopen (file, "thermo");
status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->thermo_eos);
H5Dclose (dataset);

dataset = H5Dopen (file, "yi");
status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->yi_eos);
H5Dclose (dataset);

dataset = H5Dopen (file, "aav");
status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->aheavy_eos);
H5Dclose (dataset);

dataset = H5Dopen (file, "yav");
status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->xheavy_eos);
H5Dclose (dataset);

dataset = H5Dopen (file, "zav");
status = H5Dread (dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->zheavy_eos);
H5Dclose (dataset);

// compose : rho and t stored instead of ln_rho and ln_t
  for (int i = 0; i < m_ln_rho; i++) {
    ln_rho_eos[i] = log(ln_rho_eos[i]) ;
  }

  for (int j = 0; j < n_ln_t; j++) {
    ln_t_eos[j] = log(ln_t_eos[j]) ;
  }

  // compose split of large arrays
  //   split up of thermo_eos and y_i_eos
  for (int i = 0; i < m_ln_rho; i++) {
    for (int j = 0; j < n_ln_t; j++) {
      for (int k = 0; k < o_y_e; k++) {
	int index = i + j * m_ln_rho + k * m_ln_rho * n_ln_t;
	xheavy_eos [index] = xheavy_eos [index] * aheavy_eos [index] ;
	p_eos [index] = thermo_eos [index] ;

	// just a shift
	int extended_index = m_ln_rho * n_ln_t * o_y_e + index ;
	entropy_eos [index] = thermo_eos [extended_index] ;

	extended_index = /*2.*/1 * m_ln_rho * n_ln_t * o_y_e + index ;
	mub_eos [index] = thermo_eos [extended_index] ;

	extended_index = /*3.*/2 * m_ln_rho * n_ln_t * o_y_e + index ;
	muq_eos [index] = thermo_eos [extended_index] ;

	extended_index = /*4.*/3 * m_ln_rho * n_ln_t * o_y_e + index ;
	mul_eos [index] = thermo_eos [extended_index] ;

	extended_index = /*6.*/4 * m_ln_rho * n_ln_t * o_y_e + index ;
	e_eos [index] = thermo_eos [extended_index] ;

	extended_index = /*7.*/5 * m_ln_rho * n_ln_t * o_y_e + index ;
	c_sound_squared_newtonian_eos [index] = thermo_eos [extended_index] ;

	xn_eos [index] = yi_eos [index] ;

	extended_index = m_ln_rho * n_ln_t * o_y_e + index ;
	xp_eos [index] = yi_eos [extended_index] ;
      }
    }
  }
  break;
  default: return 1;
  }
  // close file
  status = H5Fclose (file);
}

int short_EOS_table::write(const char *path)
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
  H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->m_ln_rho);
  H5Dclose (dataset);
  H5Sclose(dataspace);

  dataspace = H5Screate_simple(1, &one, NULL);
  dataset = H5Dcreate2(file, "pointst", H5T_NATIVE_INT, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->n_ln_t);
  H5Dclose (dataset);
  H5Sclose(dataspace);


  dataspace = H5Screate_simple(1, &one, NULL);
  dataset = H5Dcreate2(file, "pointsye", H5T_NATIVE_INT, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->o_y_e);
  H5Dclose (dataset);
  H5Sclose(dataspace);

  dataspace = H5Screate_simple(1, &one, NULL);
  dataset = H5Dcreate2(file, "pointsmu", H5T_NATIVE_INT, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->p_mu);
  H5Dclose (dataset);
  H5Sclose(dataspace);


  hsize_t sz[1] = { this->m_ln_rho };
  dataspace = H5Screate_simple(1, sz, NULL);
  dataset = H5Dcreate2(file, "logrho", H5T_NATIVE_DOUBLE, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->ln_rho_eos);
  H5Dclose (dataset);
  H5Sclose(dataspace);


  sz[0] = this->n_ln_t;
  dataspace = H5Screate_simple(1, sz, NULL);
  dataset = H5Dcreate2(file, "logt", H5T_NATIVE_DOUBLE, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->ln_t_eos);
  H5Dclose (dataset);
  H5Sclose(dataspace);

  sz[0] = this->o_y_e;
  dataspace = H5Screate_simple(1, sz, NULL);
  dataset = H5Dcreate2(file, "ye", H5T_NATIVE_DOUBLE, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->y_e_eos);
  H5Dclose (dataset);
  H5Sclose(dataspace);

  sz[0] = this->p_mu;
  dataspace = H5Screate_simple(1, sz, NULL);
  dataset = H5Dcreate2(file, "mu_nu", H5T_NATIVE_DOUBLE, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->mu_nu_eos);
  H5Dclose (dataset);
  H5Sclose(dataspace);

  hsize_t dims[4] = { this->m_ln_rho, this->n_ln_t, this->o_y_e, this->p_mu };

  sz[0] = this->size();
  dataspace = H5Screate_simple(1, sz, NULL);
  dataset = H5Dcreate2(file, "elec_capt_rate", H5T_NATIVE_DOUBLE, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->elec_rate_tab_eos);
  H5Dclose (dataset);

  dataset = H5Dcreate2(file, "elec_capt_single_rate", H5T_NATIVE_DOUBLE, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->elec_rate_single_eos);
  H5Dclose (dataset);

  dataset = H5Dcreate2(file, "elec_capt_fast_rate", H5T_NATIVE_DOUBLE, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->elec_rate_fast_eos);
  H5Dclose (dataset);


  dataset = H5Dcreate2(file, "sigma_scattering_nuclei", H5T_NATIVE_DOUBLE, dataspace, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->scattering_xs_eos);
  H5Dclose (dataset);
  H5Sclose(dataspace);

  // close file
  status = H5Fclose (file);
  return 0;
}
