#ifndef __NEWTON_H__
#define __NEWTON_H__
#include "Eigen/Dense"
#include "global.h"
#include "line_search.h"

class NewtonBase :public OptimizationBase
{
public:
    TVector Solve(TargetFunctor& fucntor, Options& options) override;
};

class DampedNewton : public NewtonBase
{
public:
    virtual TVector Solve(TargetFunctor& fucntor, Options& options);
};

#endif //__NEWTON_H__
