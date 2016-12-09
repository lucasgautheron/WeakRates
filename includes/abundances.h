int read_abundance_data(const char *path_abundances, const char *path_temperatures, const char *path_densities, const char *path_fractions);
double element_abundance(int A, int Z, double params[3], double *vlow = NULL, double *vhigh = NULL);
