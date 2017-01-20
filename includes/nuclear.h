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

inline double SEMF(int A, int Z)
{
    const double a_V = 15.8, a_S = 18.3, a_C = 0.714, a_A = 23.2, a_P = 12;
    const double delta_0 = a_P/sqrt(A);
    const double delta = A&1 ? 0 : (Z&1 ? -delta_0 : delta_0);
    const double B = a_V * A - a_S * pow(A, 2./3.) - a_C * Z*Z/pow(A, 1./3.) - a_A*(A-2*Z)*(A-2*Z)/A+delta;
    return Z*M_PROTON+(A-Z)*M_NEUTRON-B;
}

inline double nucleus_mass(int A, int Z)
{
    std::array<int, 2> elem = { A, Z };
    return (nuclear_table.count(elem)) ? nuclear_table[elem]->m : SEMF(A, Z);
}

inline double beta_decay_Q(int A, int Z)
{
    std::array<int, 2> mother = { A, Z }, daughter = { A, Z-1 };
    if(!nuclear_table.count(daughter))
    {
        return nucleus_mass(A, Z)-nucleus_mass(A, Z-1);
    }
    if(std::abs(nuclear_table[daughter]->beta_q) >= 0.0001)
    {
        return -nuclear_table[daughter]->beta_q;
    }
    if(!nuclear_table.count(mother))
    {
        return nucleus_mass(A, Z)-nucleus_mass(A, Z-1);
    }
    return nuclear_table[mother]->m - nuclear_table[daughter]->m;
}

int read_nuclear_data(const char *path);
