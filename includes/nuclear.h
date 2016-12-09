struct nuclear_data
{
    int A, Z;
    double m; // MeV
    double beta_q; // beta decay Q

    nuclear_data() {}
    nuclear_data(int _A, int _Z, double _m, double _bq)
    {
        A = _A; Z = _Z; m = _m; beta_q = _bq;
    }
};
typedef std::map< std::array<int, 2>, nuclear_data *> nuclear_array;

extern nuclear_array nuclear_table;

inline double nucleus_mass(int A, int Z)
{
    std::array<int, 2> elem = { A, Z };
    return (nuclear_table.count(elem)) ? nuclear_table[elem]->m : double(A)*MASS_UNIT;
}

inline double beta_decay_Q(int A, int Z)
{
    std::array<int, 2> mother = { A, Z }, daughter = { A, Z-1 };
    if(!nuclear_table.count(mother)) return -1e10;
    if(abs(nuclear_table[mother]->beta_q) >= 0.0001) return nuclear_table[mother]->beta_q;
    if(!nuclear_table.count(daughter)) return -1e10;
    return nuclear_table[mother]->m - nuclear_table[daughter]->m;
}

int read_nuclear_data(const char *path, nuclear_array &table)
{
    std::ifstream infile(path);
    int count = 0;
    std::string line;
    std::istringstream iss;
    while (std::getline(infile, line)) {
        const char *str = line.c_str();
        int Z = atoi(str+11);
        int A = atoi(str+15);
        int base = atoi(str+96);
        double m = double(base) + 1e-6 * strtod(str+100, NULL);
        m *= MASS_UNIT;
        double beta_q = strtod(str+79, NULL)/1e3;
        //printf("%d %d %e %.3f\n", Z, A, m, beta_q);
        std::array<int, 2> elem = { A, Z };
        table[elem] = new nuclear_data(A, Z, m, beta_q);
        ++count;
    }
    return count;
}
