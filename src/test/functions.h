//
// Created by ReFantasy on 2022/3/13.
//

#ifndef OPTIMIZATION_FUNCTIONS_H
#define OPTIMIZATION_FUNCTIONS_H

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


/**
 *The plot on the left shows the three-hump Camel function on its recommended input domain,
 * and the plot on the right shows only a portion of this domain,
 * to allow for easier viewing of the function's key characteristics. The function has three local minima.
 *
 * Dimensions: 2
 *
 * \f$ f(x) = 2 x_1^2-1.05 x_1^4+\frac{x_1^6}{6}+x_1 x_2+x_2^2 \f$
 *
 * http://www.sfu.ca/~ssurjano/camel3.html
 */
class THREE_HUMP_CAMEL : public TargetFunctor
{
public:
	FLOAT operator()(const Vector &xk) const override;

	Vector Gradient(const Vector &xk) const override;

	Matrix Hesse(const Vector &xk) const override;
};
#endif // OPTIMIZATION_FUNCTIONS_H
