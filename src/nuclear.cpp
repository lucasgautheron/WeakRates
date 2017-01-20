#include "common.h"

nuclear_array nuclear_table;

int read_nuclear_data(const char *path)
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
        nuclear_table[elem] = new nuclear_data(A, Z, m, beta_q);
#ifdef DEBUG
        std::cout << "Nucleus (" << A << "," << Z << "), m = " << m << " MeV, betaQ = " << beta_q << "\n";
#endif
        ++count;
    }
    return count;
}
