#include <iostream>
#include <fstream>
#include <vector>
#include "Isingsystems.hpp"

using namespace std;

int main() {
    const vector<int> system_sizes = {3,2,2};
    const vector<vector<vector<int>>> linklist={{{0,1},{0,2},{-1,2,2},{2,1}},{{0,0},{0,2},{-2,0},{-1,2}},{{0,0},{0,1},{1,1},{1,-2,0}}};
    vector<double> temperatures;
    int n_sp = system_sizes[0] * system_sizes[1] * system_sizes[2];
    for (double temp = 0.05; temp <= 4.0; temp += 0.05) {
        temperatures.push_back(temp);
    }

    Isingsystems ising_system(system_sizes,linklist);

    ofstream outfile("Ising_results.txt");
    outfile<<"System size: ";
    for(size_t i=0;i<system_sizes.size()-1;i++){
        outfile << system_sizes[i+1] << "x";
    }
    outfile<<endl;
    outfile << "Temperature\tMagnetization^2\tEnergy\tEnergySquare\tSpecific Heat\tUM" << endl;

    vector<vector<double>> results = ising_system.average_eval_combined(temperatures);

    for (size_t j = 0; j < temperatures.size(); j++) {
        double temp = temperatures[j];

        double specific_heat = (results[3][j] * n_sp - results[2][j] * results[2][j] * n_sp) / temp / temp;
        double UM = results[1][j] / results[0][j] / results[0][j];

        outfile << temperatures[j] << "\t" << results[0][j] << "\t" << results[2][j] << "\t" << results[3][j] << "\t" << specific_heat << "\t" << UM << endl;
    }

    outfile.close();

    return 0;
}
