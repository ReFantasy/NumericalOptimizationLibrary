/**
 *
 * Solving positive definite quadratic function by steepest descent method
 * f(x) = 1/2 * x * G * x + b * x + c;
 *
 * Author  : ReFantasy
 * Date    : 2022-02-19
 */

#ifndef __SD_H__
#define __SD_H__
#include "Eigen/Dense"
#include "global.h"

namespace NOL
{
class SteepestDescent : public UnconstrainedOptimizationLineSearchBase
{
  public:
    SteepestDescent(TargetFunctor &function, Options &options)
        : UnconstrainedOptimizationLineSearchBase(function, options)
    {
    }

    bool IsTermination(const Vector &xk, int k) const override;

    Vector DescentDirection(const Vector &xk) const override;

    FLOAT StepSize(const Vector &xk, const Vector &dk) const override;
};
} // namespace NOL

#endif //__SD_H__
