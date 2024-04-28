ifndef Isingspins_hpp
#define Isingspins_hpp

#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
using namespace std;
class Isingspins {
private:
    int spin;
	vector<int> position;
    vector<int> NN;
    int idx;
public:
	Isingspins() { spin = 1; };
	~Isingspins() {};

	void set_up() {spin = 1;};
	void set_down() {spin = -1;};
	void set_spin(int sz_spec) {
		assert(sz_spec == 1 || sz_spec == -1);
		spin = sz_spec;
	}
    void set_pos(vector<int> pos) {position = pos;};
    void set_NN(vector<int> nn) {NN = nn;};
    void set_idx(int idx_spec) {idx = idx_spec;};
	int _spin() const { return spin; };
    vector<int> _position() { return position; };
    vector<int> _NN()  { return NN; };
    int _NN(const int bond_idx)  { return NN[bond_idx]; };
    void set_NN(const int bond_idx, const int site_idx) { NN[bond_idx] = site_idx; };
	void flip() { spin *= -1; };
};

#endif
