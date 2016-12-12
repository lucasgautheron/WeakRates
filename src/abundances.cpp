#include "common.h"

struct element
{
    int nucleus;
    int A,Z;
    double abundance;

    element() : nucleus(0), A(0), Z(0) { }
};

struct abundance_data
{
    int idx[3];
    double param[3];
    double np, nn;
    
    std::vector<element> elements;
};


typedef std::map< std::array<int, 3>, abundance_data *> abundance_array;

struct abundance_table
{
    std::vector<double> parameters[3];
    abundance_array abundances;
};

abundance_table table;

struct index_entry
{
    int idx;
    double val;
};

inline void nucleus_to_AZ(int nucleus, int &A, int &Z)
{
    A = nucleus/1000;
    Z = nucleus%1000;
}

int read_index(const char *path, std::vector<double> &index)
{
    std::ifstream infile(path);
    double val;
    int count = 0;
    std::string line;
    std::istringstream iss;
    while (std::getline(infile, line)) {
        iss.clear();
        iss.str(line);
        iss >> val;
        ++count;
        if(count <= 2) continue;
        index.push_back(val);
    }
    return count-2;
}

int read_abundance_data(const char *path_abundances, const char *path_temperatures, const char *path_densities, const char *path_fractions)
{
    read_index(path_temperatures, table.parameters[0]);
    read_index(path_densities, table.parameters[1]);
    read_index(path_fractions, table.parameters[2]);

    std::ifstream infile(path_abundances);
    
    int int_ph;
    int count = 0;
    std::string line;
    std::istringstream iss;
    while (std::getline(infile, line)) {
        abundance_data *ab = new abundance_data();
        iss.clear();
        iss.str(line);
        for(int k = 0; k < 3; ++k) iss >> ab->idx[k];
        iss >> int_ph;
        iss >> int_ph;
        iss >> int_ph;
        iss >> ab->nn;
        iss >> int_ph;
        iss >> ab->np;

        for(int k = 0; k < 3; ++k) ab->param[k] = table.parameters[k][ab->idx[k]];

        while(true)
        {
            element e;
            iss >> e.nucleus;
            iss >> e.abundance;
            if(!e.nucleus) break;

            nucleus_to_AZ(e.nucleus, e.A, e.Z);
            ab->elements.push_back(e);
            //printf("%d %d %f %e %e %e\n", e.A, e.Z, ab->param[0], bruenn_electron_capture_rate(e.A, e.Z, ab->param[0]), bruenn_electron_capture_rate(e.A, e.Z, ab->param[0], 100), bruenn_electron_capture_rate(e.A, e.Z, ab->param[0], 10000));
        }

        element neutron, proton;
        neutron.A = 1;
        neutron.Z = 0;
        neutron.abundance = ab->nn;
        proton.A = 1;
        proton.Z = 1;
        proton.abundance = ab->np;
        ab->elements.push_back(neutron);
        ab->elements.push_back(proton);

        std::array<int, 3> conditions = {ab->idx[0], ab->idx[1], ab->idx[2]};
        table.abundances[conditions] = ab;
        ++count;
    }
    return count;
}

int get_lower_key(std::vector<double> &arr, double value)
{
    auto it = std::lower_bound(arr.begin(), arr.end(), value);
    if(it == arr.end()) return -1;
    return it-arr.begin()-1;
}

// (table, A, Z, {T [MeV], nb [fm^{-3}], Ye})
double element_abundance(int A, int Z, double params[3], double *vlow, double *vhigh)
{
    int keys[3];
    for(int k = 0; k < 3; ++k) keys[k] = get_lower_key(table.parameters[k], params[k]);
 
    std::array<int, 3> lower = { keys[0], keys[1], keys[2] }, upper = { keys[0]+1, keys[1]+1, keys[2]+1 };
    std::array<int, 3> next[3] = { { keys[0]+1, keys[1], keys[2] }, { keys[0], keys[1]+1, keys[2] }, { keys[0], keys[1], keys[2]+1 } };
    double dv[3], dx[3], x[3];

    abundance_data *ptr = table.abundances[lower];
    double v = 0;
    if(ptr)
    {
        abundance_data ab = *ptr;
        for(int i = 0; i < ab.elements.size(); ++i) if(ab.elements[i].A == A && ab.elements[i].Z == Z) { v = ab.elements[i].abundance; break; }
    }

    for(int k = 0; k < 3; ++k)
    {
        dx[k] = keys[k]+1 < table.parameters[k].size() ? table.parameters[k][keys[k]+1]-table.parameters[k][keys[k]] : 1;
        dv[k] = 0;
        x[k] = table.parameters[k][keys[k]];

        if(!table.abundances.count(next[k]) || !table.abundances[next[k]]) continue;

        abundance_data ab_next = *table.abundances[next[k]];
        for(int i = 0; i < ab_next.elements.size(); ++i) if(ab_next.elements[i].A == A && ab_next.elements[i].Z == Z) { dv[k] = ab_next.elements[i].abundance - v; break; }
    }

    double v_next = 0;
    {
        abundance_data *ptr = table.abundances.count(upper) ? table.abundances[upper] : NULL;
        if(!ptr) v_next = v;
        else
        {
            abundance_data ab = *ptr;
            for(int i = 0; i < ab.elements.size(); ++i) if(ab.elements[i].A == A && ab.elements[i].Z == Z) { v_next = ab.elements[i].abundance; break; }
        }
    }

    double vals[] = { v+dv[0], v+dv[1], v+dv[2], v_next };

    if(vlow)
    {
        *vlow = *std::min_element(std::begin(vals), std::end(vals));
    }
    if(vhigh)
    {
        *vhigh = *std::max_element(std::begin(vals), std::end(vals));
    }

    for(int k = 0; k < 3; ++k) v += (params[k]-x[k]) * dv[k] / dx[k];

    return v;
}
