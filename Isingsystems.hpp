#ifndef Isingsystems_hpp
#define Isingsystems_hpp

#include <omp.h>
#include <vector>
#include <cassert>
#include "Isingspins.hpp"
#include <cmath>
#include <bitset>
using namespace std;
//dim的第0个是在晶格中的位置,linklist三层分别是本体位置，与之连接的晶格，维度差的序号(负数代表反方向，可以有多个但必须补零对齐)与对方在晶格中的位置
//自身晶格内是第0个维度
class Isingsystems {
protected:
	const double J;
	int n_spins;
    const vector<int> dim;
	long long maxrep_state;
	vector <Isingspins> spin;

public:
	Isingsystems(vector<int> dim_spec,vector<vector<vector<int>>> linklist) : J(-1.0),dim(dim_spec) {
        int n_spins_spec=1;
        for(size_t i = 0; i < dim_spec.size(); i++) n_spins_spec *= dim_spec[i];
        n_spins = n_spins_spec;
		maxrep_state=static_cast<long long>(pow(2, n_spins)) - 1;
		spin.resize(n_spins);
        for(int i=0;i<n_spins;i++){
            spin[i].set_idx(i);
            vector<int> pos;
            int temp=i;
            for(size_t j=0;j<dim_spec.size();j++){
                pos.push_back(temp%dim_spec[j]);
                temp=temp/dim_spec[j];
            }
            spin[i].set_pos(pos);
            vector<vector<int>> linklist1=linklist[i%(dim_spec[0])];
            vector<int> nn;
            for(size_t j=0;j<linklist1.size();j++){
                vector<int> pos1=pos;
				for(size_t k=0;k<linklist1[j].size()-1;k++){
					int sign=linklist1[j][k]>0?1:-1;
                	if(linklist1[j][k]==0)sign=0;
					pos1[abs(linklist1[j][k])]= (pos1[abs(linklist1[j][k])]+sign+dim_spec[abs(linklist1[j][k])])%dim_spec[abs(linklist1[j][k])];

				}
                pos1[0]=linklist1[j][linklist1[j].size()-1];

                int idx1=0,n=1;
                for(size_t k=0;k<dim_spec.size();k++){
                    if(k!=0)
                        n=n*dim_spec[k-1];
                    idx1=idx1+n*pos1[k];
                }
                nn.push_back(idx1);
            }
            spin[i].set_NN(nn);
        }
	}
	virtual ~Isingsystems() {};
	double _J() const { return J; };
	int _n_spins() const { return n_spins; };
	long long _maxrep_state() const { return maxrep_state; };
	vector<int> _position(const int site_idx) { return spin[site_idx]._position(); };
	void set_up_spin(const int site_idx) { spin[site_idx].set_up(); }
	void set_down_spin(const int site_idx) { spin[site_idx].set_down(); }
	void set_pos(const int site_idx, vector<int> s_spec){spin[site_idx].set_pos(s_spec);}
    void set_NN(const int site_idx, vector<int> nn_spec) {spin[site_idx].set_NN(nn_spec);}
    void set_NN(const int site_idx, const int bond_idx, const int site_idx2) {spin[site_idx].set_NN(bond_idx, site_idx2);}
    void set_idx(const int site_idx, int idx_spec) {spin[site_idx].set_idx(idx_spec);}
    int _spin(const int site_idx) const { return spin[site_idx]._spin(); };
	void flip_spin(const int site_idx) { spin[site_idx].flip(); }
	int _NN(const int site_idx, const int bond_idx) { return spin[site_idx]._NN(bond_idx); };
	void set_state_by_code(long long rep_state) {
		assert(rep_state <= maxrep_state);
		for (int i = 0; i < n_spins; i++) {
			if (rep_state & 1) {  
				spin[i].set_up();
			}
			else {
				spin[i].set_down();
			}
			rep_state >>= 1;  
		}
	};
    double eval_mz() {
        double mz_total = 0;
        for (int i = 0; i < n_spins; i++) { mz_total += spin[i]._spin(); }
        return mz_total;
    };
    double eval_energy() {
        double energy_total = 0;
        for (int i = 0; i < n_spins; i++) {
            for(size_t j=0;j<spin[i]._NN().size();j++){
                energy_total += J * spin[i]._spin() * spin[spin[i]._NN(j)]._spin();
            }
        }
        return energy_total/2;
    };
    vector<vector<double>> average_eval_combined(const vector<double> temperatures) {
		vector<vector<double>> results(4, vector<double>(temperatures.size(), 0));
		int mz,energy, previous_energy, previous_mz,state;
		vector <double> Z;
		Z.resize(temperatures.size(),0);
		for (int i = 0; i <= maxrep_state; i++) {
			if (i == 0)
			{
				set_state_by_code(i);
				energy = eval_energy();
				mz = eval_mz();
			}
			else if (i != 0)
			{
				state = which_state_to_change(i);
				mz = previous_mz - 2  * spin[state]._spin();
				int energy_changed = 0;
				for (size_t j = 0; j < spin[state]._NN().size(); j++)
					energy_changed = energy_changed - spin[spin[state]._NN(j)]._spin() * spin[state]._spin() * J;
				energy = previous_energy + energy_changed * 2;
				flip_spin(state);
			}
			previous_energy = energy;
			previous_mz = mz;


			for (size_t j = 0; j < temperatures.size(); j++) {
				double weight = exp(-energy / temperatures[j]);

				results[0][j] += weight * mz * mz;
				results[1][j] += weight * mz * mz * mz * mz;
				results[2][j] += weight * energy;
				results[3][j] += weight * energy * energy;
				Z[j] += weight;
			}
		}

		for (size_t k = 0; k < temperatures.size(); k++) {

			results[0][k] /= n_spins * n_spins * Z[k];
			results[1][k] /= n_spins * n_spins * n_spins * n_spins * Z[k];
			results[2][k] /= n_spins * Z[k];
			results[3][k] /= n_spins * n_spins * Z[k];
		}

		return results;
	}

	int which_state_to_change(int state_idx)
	{
		int i = 0;
		int divisor = 2;
		while (true)
		{
			if (state_idx % divisor == (divisor / 2))
			{
				return i;
			}
			else
			{
				i++;
				divisor *= 2;
			}
		}
	};
    vector<vector<double>> average_eval_paralelled(const vector<double> temperatures) {
	vector<vector<double>> results(4, vector<double>(temperatures.size()*16, 0));
	omp_set_num_threads(16);
	vector <double> Z;
    int mz,energy;
	Z.resize(temperatures.size()*16,0);
	#pragma omp parallel
	{

	#pragma omp for		
	for (int j1 = 0; j1 < 16; j1++) {

		for (int i1 = 0; i1 < (maxrep_state + 1) / 2 / 16; i1++) {
			int i = j1 * (maxrep_state + 1) / 2 / 16 + i1;
			set_state_by_code(i);
			mz = eval_mz();
			energy = eval_energy();
			for (size_t j = 0; j < temperatures.size(); j++) {

				double weight = exp(-energy / temperatures[j]);
				results[0][j1 * temperatures.size() + j] += weight * mz * mz;
				results[1][j1 * temperatures.size() + j] += weight * mz * mz * mz * mz;
				results[2][j1 * temperatures.size() + j] += weight * energy;
				results[3][j1 * temperatures.size() + j] += weight * energy * energy;
				Z[j1 * temperatures.size() + j] += weight;
			}
		}

	}


	}
	for (size_t k = 0; k < temperatures.size(); k++) {
		for (size_t k1 = 1; k1 < 16; k1++) {
			for (size_t k2 = 0; k2 < 4; k2++) {
				results[k2][k] += results[k2][k + k1 * temperatures.size()];
			}
			Z[k] += Z[k + k1 * temperatures.size()];

		}


	}
	for (size_t k = 0; k < temperatures.size(); k++) {

		results[0][k] /= n_spins * n_spins * Z[k];
		results[1][k] /= n_spins * n_spins * n_spins * n_spins * Z[k];
		results[2][k] /= n_spins * Z[k];
		results[3][k] /= n_spins * n_spins * Z[k];
	}
	for (size_t k = 0; k < 4; k++) {
		results[k].resize(temperatures.size());
	}


	return results;
}
};



#endif
