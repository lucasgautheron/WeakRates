#include "common.h"

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

int read_abundance_data(const char *path)
{
    read_index((std::string(path)+"/eos.t").c_str(), table.parameters[0]);
    read_index((std::string(path)+"/eos.nb").c_str(), table.parameters[1]);
    read_index((std::string(path)+"/eos.yq").c_str(), table.parameters[2]);

    std::ifstream infile((std::string(path)+"/eos.compo").c_str());
    
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
void get_abundances(double params[3], std::vector<element> &elements)
{
    int keys[4];
    for(int k = 1; k < 4; ++k) keys[k] = get_lower_key(table.parameters[k-1], params[k-1]);

    std::array<int, 3> cells[4] = { {keys[1], keys[2], keys[3]}, { keys[1]+1, keys[2], keys[3] }, { keys[1], keys[2]+1, keys[3] }, { keys[1], keys[2], keys[3]+1 } };
    double dx[4], x[4];

    for(int k = 1; k < 4; ++k)
    {
        dx[k] = keys[k]+1 < table.parameters[k-1].size() ? table.parameters[k-1][keys[k]+1]-table.parameters[k-1][keys[k]] : 1;
        x[k] = table.parameters[k-1][keys[k]];
    }

    for(int k = 0; k < 4; ++k)
    {
        if(!table.abundances.count(cells[k]) || !table.abundances[cells[k]])  continue;
        abundance_data ab = *table.abundances[cells[k]];

        for(int i = 0; i < ab.elements.size(); ++i)
        {
            bool found = false;
            for(int j = 0; j < elements.size(); ++j) if(elements[j].A == ab.elements[i].A && elements[j].Z == ab.elements[i].Z)
            {
                found = true;
                elements[j].v[k] = ab.elements[i].abundance;
                break;
            }
            if(!found)
            {
                element e;
                e.A = ab.elements[i].A;
                e.Z = ab.elements[i].Z;
                e.v[k] = ab.elements[i].abundance;
                elements.push_back(e);
            }
        }
    }

    for(int i = 0; i < elements.size(); ++i)
    {
        elements[i].abundance = elements[i].v[0];
        for(int k = 1; k < 4; ++k)
        {
            elements[i].abundance += (params[k-1]-x[k]) * (elements[i].v[k]-elements[i].v[0]) / dx[k];
        }
    }
}

