#include <catch2/catch_test_macros.hpp>
#include <catch2/mathchers/catch_matchers_floating_point.hpp>
#include <iostream>
#include "IsingSystem.hpp"

TEST_CASE("IsingSpin", "[single spin]") {
	IsingSpin spin;

	SECTION("spin state (initial)") {
		REQUIRE(spin._sz() == 1);
	}
	SECTION("set spin state as up (1)") {
		spin.set_up();
		REQUIRE(sin._sz() == 1);
	}
	SECTION("set spin state as down (1)") {
		spin.set_dw();
		REQUIRE(sin._sz() == -1);
	}
	SECTION("set spin state as up (2)") {
		spin.set_sz(1);
		REQUIRE(sin._sz() == 1);
	}
	SECTION("set spin state as down (1)") {
		spin.set_sz(-1);
		REQUIRE(sin._sz() == -1);
	}
	SECTION("spin flip once") {
		spin.flip();
		REQUIRE(sin._sz() == -1);
	}
	SECTION("spin flip twice") {
		spin.flip();
		spin.flip();
		REQUIRE(sin._sz() == 1);
	}
};

TEST_CASE("IsingSystem", "spin") {
	IsingSystem spin;

	SECTION("set spin state as up (1)") {
		spin.set_up_spin(int i);
		REQUIRE(spin._sz(i) == 1);
	}
	SECTION("set spin state as down (1)") {
		spin.set_dw_spin(int i);
		REQUIRE(spin._sz(i) == -1);
	}
	SECTION("set spin state as up (2)") {
		spin.set_spin(int i, 1);
		REQUIRE(spin._sz(i) == 1);
	}
	SECTION("set spin state as down (2)") {
		spin.set_spin(int i, -1);
		REQUIRE(spin._sz(i) == -1);
	}
	SECTION("spin flip once") {
		spin.flip_spin(int i);
		REQUIRE(spin._sz(i) == -1);
	}
	SECTION("spin flip twice") {
		spin.flip_spin(int i);
		spin.flip_spin(int i);
		REQUIRE(spin._sz(i) == -1);
	}
};
