//
// Created by ReFantasy on 2022/3/13.
//

#ifndef OPTIMIZATION_TEST_FUNCTIONS_H
#define OPTIMIZATION_TEST_FUNCTIONS_H

#include "global.h"

using namespace NOL;
/**
 * The Rotated Hyper-Ellipsoid function is continuous, convex and unimodal.
 * It is an extension of the Axis Parallel Hyper-Ellipsoid function, also referred to as the Sum Squares function.
 * The plot shows its two-dimensional form.
 *
 * Dimensions: d
 *
 * \f$ f(x) = \sum_{i=1}^{d} \sum_{j=1}^i x_j^2 \f$
 *
 * http://www.sfu.ca/~ssurjano/rothyp.html
 */
class ROTATED_HYPER_ELLIPSOID : public TargetFunctor
{
  public:
    FLOAT operator()(const Vector &xk) const override;

    [[nodiscard]] Vector Gradient(const Vector &xk) const override;

    [[nodiscard]] Matrix Hesse(const Vector &xk) const override;
};

#endif // OPTIMIZATION_TEST_FUNCTIONS_H
