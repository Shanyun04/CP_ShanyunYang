#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "Isingsystems.hpp"
using namespace std;

TEST_CASE("IsingSystemCubic","[examples of 2 x 2 x 2 spins]") {
  const vector<int> system_size = {1,2,2,2};
  const vector<double> temperature = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0};
  const vector<vector<vector<int>>> linklist={{{1,0},{-1,0},{2,0},{-2,0},{3,0},{-3,0}}};
  Isingsystems model(system_size,linklist);
  SECTION("basics") {
	REQUIRE(model._n_spins() ==8);
	REQUIRE_THAT(model._J(),Catch::Matchers::WithinULP(-1.0,4));
  };

  SECTION("connectivity") {
	constexpr int i =0;
	REQUIRE(model._NN(i,0) ==1);
	REQUIRE(model._NN(i,1) ==1);
	REQUIRE(model._NN(i,2) ==2);
	REQUIRE(model._NN(i,3) ==2);
  REQUIRE(model._NN(i,4) ==4);
  REQUIRE(model._NN(i,5) ==4);

  };


  SECTION("magetization and energy"){
	model.set_state_by_code(0);
	REQUIRE(model.eval_mz() ==-8);
	REQUIRE(model.eval_energy() ==-24);
  };

  SECTION("results"){
        vector <vector<double>> results = model.average_eval_combined(temperature);
        double tolerance = 0.000000001;
        int J=-1;
//Honeycomb (2 x 2)
        for(size_t i=0;i<temperature.size();i++){
          long double beta=1/temperature[i];
          // Cubic (2 2 2 2 2)
          long double E3 = -24.0/8.0 * (-1 + 3 * exp(4 * beta * abs(J)) - 5 * exp(8 * beta * abs(J)) + 3 * exp(12 * beta * abs(J)) - 3 * exp(16 * beta * abs(J)) + 5 * exp(20 * beta * abs(J)) - 3 * exp(24 * beta * abs(J)) + exp(28 * beta * abs(J))) / ((1 + exp(4 * beta * abs(J))) * (1 + exp(8 * beta * abs(J))) * (1 - 4 * exp(4 * beta * abs(J)) + 8 * exp(8 * beta * abs(J)) - 4 * exp(12 * beta * abs(J)) + exp(16 * beta * abs(J))));
          REQUIRE_THAT(results[2][i], Catch::Matchers::WithinAbs(E3, tolerance));
        }
  }

};

/*
TEST_CASE("IsingSystemCubic", "[calculate functions]")
        const vector<int> system_size = {2,2,2};
        const vector<double> temperature = {0.1};
	IsingSystemCubic model(system_size,temperature);
    SECTION("Partition Function Calculation") {
	int index = 0;
        double partition_func = model.partition_functionsreturn(index);
        REQUIRE_THAT(partition_func,Catch::Matchers::WithinULP(,4));
    }
    SECTION("Average Magnetization Square Calculation") {
        double avg_mag_square = model.average_eval_Msquare();
        REQUIRE(avg_mag_square == Approx(1.0).epsilon(0.01));
    }
    SECTION("Average Energy Calculation") {
        double avg_energy = model.average_eval_energy();
        REQUIRE(avg_energy == Approx(-4.0).epsilon(0.01));
    }
    SECTION("Specific Heat Calculation") {
        double specific_heat = model.specific_heat();
        REQUIRE(specific_heat == Approx(0.0).epsilon(0.01));
    }
}
*/
