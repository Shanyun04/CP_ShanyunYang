#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "Isingsystems.hpp"
using namespace std;

TEST_CASE("IsingSystemCubic","[examples of 2 x 2 x 2 spins]") {
  const vector<int> system_size = {3,2,2};
  const vector<double> temperature = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0};
  const vector<vector<vector<int>>> linklist={{{0,1},{0,2},{-1,2,2},{2,1}},{{0,0},{0,2},{-2,0},{-1,2}},{{0,0},{0,1},{1,1},{1,-2,0}}};
  Isingsystems model(system_size,linklist);

  SECTION("basics") {
	REQUIRE(model._n_spins() ==12);
	REQUIRE_THAT(model._J(),Catch::Matchers::WithinULP(-1.0,4));
  };



  SECTION("connectivity") {
	constexpr int i =0;
	REQUIRE(model._NN(i,0) ==1);
	REQUIRE(model._NN(i,1) ==2);
	REQUIRE(model._NN(i,2) ==11);
  REQUIRE(model._NN(i,3) ==7);
	};



  SECTION("magetization and energy"){
	  model.set_state_by_code(0);
	  REQUIRE(model.eval_mz() ==-12);
	  REQUIRE(model.eval_energy() ==-24);
  };

  SECTION("results"){
        vector <vector<double>> results = model.average_eval_combined(temperature);
        double tolerance = 0.000000001;
        int J=-1;
//kagome (2 x 2)
        for(size_t i=0;i<temperature.size();i++){
          long double beta=1/temperature[i];
          // Kagome (2 x 2)
          long double E2 = -24.0/12.0 * (-15 - 10 * exp(4 * beta * abs(J)) + 7 * exp(8 * beta * abs(J)) + 12 * exp(12 * beta * abs(J)) + 7 * exp(16 * beta * abs(J)) - 2 * exp(20 * beta * abs(J)) + exp(24 * beta * abs(J))) / ((5 + 2 * exp(4 * beta * abs(J)) + exp(8 * beta * abs(J))) * (9 + 12 * exp(4 * beta * abs(J)) + 14 * exp(8 * beta * abs(J)) - 4 * exp(12 * beta * abs(J)) + exp(16 * beta * abs(J))));
          REQUIRE_THAT(results[2][i], Catch::Matchers::WithinAbs(E2, tolerance));
        }




  };


};
