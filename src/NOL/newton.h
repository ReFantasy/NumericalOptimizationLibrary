#ifndef __NEWTON_H__
#define __NEWTON_H__
#include "Eigen/Dense"
#include "global.h"
#include "line_search.h"

namespace NOL
{
class NewtonBase : public UnconstrainedOptimizationLineSearchBase
{
  public:
	using UnconstrainedOptimizationLineSearchBase::UnconstrainedOptimizationLineSearchBase;

    bool IsTerminated(const Vector &xk, int k) const override;

    Vector SearchDirection(const Vector &xk) const override;

    FLOAT Step(const Vector &xk, const Vector &dk) const override
    {
        return 1.0;
    };
};

class DampedNewton : public NewtonBase
{
  public:
	using NewtonBase::NewtonBase;

    FLOAT Step(const Vector &xk, const Vector &dk) const override;
};

class QuasiNewton : public DampedNewton
{
  public:
	using DampedNewton::DampedNewton;

    Vector Solve() override;

    Vector SearchDirection(const Vector &xk) const override;

  protected:
    Matrix CorrectHk(Matrix Hk, Vector sk, Vector yk);

  private:
    Matrix _Hk;
};

} // namespace NOL
#endif //__NEWTON_H__
