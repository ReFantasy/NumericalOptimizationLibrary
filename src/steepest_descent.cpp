#include "steepest_descent.h"
#include "line_search.h"
#include <iostream>

namespace NOL
{
// page 19
// search g(k+1)^T*dk=0 -> alpha
class LineSearchForSD : public LineSearch
{
  public:
    LineSearchForSD(TargetFunctor *functor = nullptr) : LineSearch(functor)
    {
    }

  protected:
    FLOAT phi(FLOAT a) override
    {
        return _functor->FirstOrderDerivatives(xk + a * dk).transpose() * dk;
    }
    FLOAT dphi_da(FLOAT a) override
    {
        return (_functor->SecondOrderDerivatives(xk + a * dk) * dk).transpose() * dk;
    }
};

SteepestDescent::SteepestDescent()
{
    _line_search = new LineSearchForSD;
}

bool SteepestDescent::IsTerminated(const Vector &xk, int k) const
{
    FLOAT gk_norm = _functor->FirstOrderDerivatives(xk).norm();
    if (gk_norm < _options->gk_norm)
        return true;

    // <--------
    *_options << "k:" << k << " "
              << "  xk:(" << xk.transpose() << ") "
              << "  ||gk||: " << gk_norm << "\n";
    // -------->

    return false;
}

NOL::Vector SteepestDescent::SearchDirection(const Vector &xk) const
{
    Vector _gk = _functor->FirstOrderDerivatives(xk);
    Vector dk = -_gk;
    return dk;
}

FLOAT SteepestDescent::Step(const Vector &xk, const Vector &dk) const
{
    _line_search->_functor = _functor;
    _line_search->xk = xk;
    _line_search->dk = dk;

    FLOAT alpha = _line_search->Search(1.0, *_options);
    return alpha;
}

} // namespace NOL
