#include "steepest_descent.h"
#include "functions.h"
#include <gtest/gtest.h>

TEST(SD, ROTATED_HYPER_ELLIPSOID)
{
	ROTATED_HYPER_ELLIPSOID functor;

	Options options;

	Vector x(4);
	x << RandomNumber<FLOAT>(-1000, 1000), RandomNumber<FLOAT>(-1000, 1000), RandomNumber<FLOAT>(-1000, 1000),
			RandomNumber<FLOAT>(-1000, 1000);
	options.init_x = x;

	options.optimized_performance = true;

	SteepestDescent sd(&functor, &options);
	Vector res = sd.Solve();
	ASSERT_LE(res.norm(), 10e-5);
}

TEST(SD, THREE_HUMP_CAMEL)
{

	THREE_HUMP_CAMEL functor;

	Options options;
	options.termination_type = TerminationCriterionType::GK_NORM;
	options.termination_value = 10e-6;

	Vector x(2);
	x << RandomNumber<FLOAT>(-1, 1), RandomNumber<FLOAT>(-1, 1);
	options.init_x = x;

	//std::cout<<functor.Hesse(x)<<std::endl;
	options.optimized_performance = true;

	SteepestDescent sd(&functor, &options);
	Vector res = sd.Solve();

	ASSERT_LE(res.norm(), 10e-5);

}