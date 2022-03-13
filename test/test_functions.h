//
// Created by ReFantasy on 2022/3/13.
//

#ifndef OPTIMIZATION_TEST_FUNCTIONS_H
#define OPTIMIZATION_TEST_FUNCTIONS_H

#include "global.h"

using namespace NOL;

class ROTATED_HYPER_ELLIPSOID : public TargetFunctor
{
  public:
    FLOAT operator()(const Vector &xk) const override;

    [[nodiscard]] Vector Gradient(const Vector &xk) const override;

    [[nodiscard]] Matrix Hesse(const Vector &xk) const override;
};

#endif // OPTIMIZATION_TEST_FUNCTIONS_H
