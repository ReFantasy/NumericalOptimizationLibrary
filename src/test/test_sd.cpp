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