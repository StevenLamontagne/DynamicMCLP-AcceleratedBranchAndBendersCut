//NOTE: This file is mostly used for common functions and definitions, though most of the 
//enums here are deprecated

#pragma once
#include <ilcplex/ilocplex.h>
#include <cmath>
#include <algorithm>
#include <random>
#include  <iterator>

#include "Data.h"


ILOSTLBEGIN
typedef IloArray<IloNumArray>    Num2D;
typedef IloArray<IloBoolVarArray> BoolVar2D;
typedef IloArray<BoolVar2D> BoolVar3D;
typedef IloArray<IloArray<IloArray<IloNum>>> Num3D;
typedef IloArray<IloArray<IloArray<IloInt>>> Int3D;

#define EPS 1.0e-6

enum class CALLBACK_STATUS {
	Unknown,
	Feasible,
	Optimal,
	Infeasible,
	Unbounded,
	InfeasibleOrUnbounded,
	Error,
	Bounded
};

enum class SINGLE_CUTS :char {
	SingleB1,
	SingleB2
};

enum class BUDGET_TYPE
{
	Knapsack,
	OutletCount
};


template <typename RandomGenerator = std::default_random_engine>
struct random_selector
{
	//On most platforms, you probably want to use std::random_device("/dev/urandom")()
	random_selector(RandomGenerator g = RandomGenerator(std::random_device()()))
		: gen(g) {}

	template <typename Iter>
	Iter select(Iter start, Iter end) {
		std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
		std::advance(start, dis(gen));
		return start;
	}

	//convenience function
	template <typename Iter>
	Iter operator()(Iter start, Iter end) {
		return select(start, end);
	}

	//convenience function that works on anything with a sensible begin() and end(), and returns with a ref to the value type
	template <typename Container>
	auto operator()(const Container& c) -> decltype(*begin(c))& {
		return *select(begin(c), end(c));
	}

private:
	RandomGenerator gen;
};