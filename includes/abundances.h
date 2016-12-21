struct element
{
    int nucleus;
    int A,Z;
    double abundance;
    double v[4];

    element() : nucleus(0), A(0), Z(0) { for(int i = 0; i < 4; ++i) v[i] = 0; }
};

int read_abundance_data(const char *path);
void get_abundances(double params[3], std::vector<element> &elements);
