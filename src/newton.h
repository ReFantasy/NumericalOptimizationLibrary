#ifndef __NEWTON_H__
#define __NEWTON_H__
#include "Eigen/Dense"
#include "global.h"
#include "line_search.h"

namespace NOL
{
class NewtonBase : public UnconstrainedOptimizationBase
{
  public:
    Vector Solve(TargetFunctor &fucntor, Options &options) override;
};

class DampedNewton : public NewtonBase
{
  public:
    virtual Vector Solve(TargetFunctor &fucntor, Options &options);
};
} // namespace NOL
#endif //__NEWTON_H__
