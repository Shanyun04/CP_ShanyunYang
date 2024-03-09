#include <iostream>
#include <vector>
#include "IsingSystem.hpp"

int main(int argc, const char * argv[]) {
	int n_spins_spec;
	long long rep_state;

	std::cout << "Please input the number of the spin." << std::endl;
	std::cin >> n_spins_spec;
	std::cout << "Please input the rep_state. (Don't more than " 
		<< std::pow(2, n_spins_spec) -1 << " .)" << std::endl;
	IsingSystem spin(n_spins_spec);
	std::cin >> rep_state;
	spin.set_state_by_code(rep_state);
	double mz = spin.eval_mz();
    std::cout << "Magnetization of the whole IsingSystem is " << mz << std::endl;
	double Energy = spin.eval_energy_1D();
	std::cout << "Energy of the whole IsingSystem is " << Energy << std::endl;
	return 0;
}