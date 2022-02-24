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
class SteepestDescent : public OptimizationBase
{
  public:
    Vector Solve(TargetFunctor &fucntor, Options &options) override;
};
} // namespace NOL

#endif //__SD_H__
