struct element
{
    int nucleus;
    int A,Z;
    double abundance;
    double v[4];

    element() : nucleus(0), A(0), Z(0) { for(int i = 0; i < 4; ++i) v[i] = 0; }
};

int read_abundance_data(const char *path_abundances, const char *path_temperatures, const char *path_densities, const char *path_fractions);
void get_abundances(double params[3], std::vector<element> &elements);
