#ifndef IsingSystem_hpp
#define IsingSystem_hpp

#include <iostream>
#include <cassert>
#include <vector>

class IsingSpin
{
private:
	int sz;

public:
	IsingSpin() { sz = 1; };
	~IsingSpin() {};

	int _sz() const { return sz; };
	void set_up() { sz = 1; };
	void set_dw() { sz = -1; };
	void set_sz(int sz_spec) {
		assert(sz_spec == 1 || sz_spec == -1);
		sz = sz_spec;
	};
	void flip() { sz *= -1; };
};

class IsingSystem
{
private:
	const double J;
	const int n_spins;
	const long long maxrep_state;
	std::vector<IsingSpin> spin;

public:
	IsingSystem(const int n_spins_spec) : J(-1.0), n_spins(n_spins_spec),
		maxrep_state(static_cast<long long>(std::pow(2, n_spins)) - 1)
	{
		spin.resize(n_spins);
	};
	virtual ~IsingSystem() {};

	double _J() const { return J; };
	int _n_spins() const { return n_spins; };
	long long _maxrep_state() const { return maxrep_state; };

	int _sz(const int site_idx) const { return spin[site_idx]._sz(); };
	void set_up_spin(const int site_idx) { spin[site_idx].set_up(); };
	void set_dw_spin(const int site_idx) { spin[site_idx].set_dw(); };
	void set_spin(const int site_idx, int s_spec) { spin[site_idx].set_sz(s_spec); };
	void flip_spin(const int site_idx) { spin[site_idx].flip(); };

	void set_state_by_code(long long rep_state) {
		int i = n_spins;
		while (rep_state != 0) {
			if (rep_state % 2 == 1) { 
				set_up_spin(i - 1);
			}
			else { 
				set_dw_spin(i - 1);
			}
			-- i;
			rep_state = rep_state / 2;
		}
		for (int i = 0; i <= n_spins - 1; i++) {
			std::cout << spin[i]._sz() << " ";
		}
		std::cout << std::endl;
	};
	double eval_mz() const {
		double mz = 0;
		for (int i = 0; i <= n_spins - 1; i++) {
			mz += spin[i]._sz();
		}
		return mz;
	};
	double eval_energy_1D() const {
		double eval_energy_1D = 0;
		for (int i = 1; i < n_spins; i++) {
			eval_energy_1D += spin[i - 1]._sz() * spin[i]._sz();
		}
		eval_energy_1D += spin[n_spins - 1]._sz() * spin[0]._sz();
		eval_energy_1D *= J;
		return eval_energy_1D;
	};
};   
#endif