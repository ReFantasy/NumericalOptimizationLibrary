#include "steepest_descent.h"
#include "line_search.h"
#include <iostream>

namespace NOL
{
// page 19
// search g(k+1)^T*dk=0 -> alpha
class LineSearchForSD : public LinearSearch
{
  public:
    LineSearchForSD(TargetFunctor *functor = nullptr) : LinearSearch(functor)
    {
    }

  protected:
    FLOAT phi(FLOAT a) override
    {
        return TargetFunctorPointer()->Gradient(Xk() + a * Dk()).transpose() * Dk();
    }
    FLOAT dphi_da(FLOAT a) override
    {
        return (TargetFunctorPointer()->Hesse(Xk() + a * Dk()) * Dk()).transpose() * Dk();
    }
};

SteepestDescent::SteepestDescent()
{
    _line_search = new LineSearchForSD;
}

bool SteepestDescent::IsTerminated(const Vector &xk, int k) const
{
    FLOAT gk_norm = _functor->Gradient(xk).norm();
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
    Vector _gk = _functor->Gradient(xk);
    Vector dk = -_gk;
    return dk;
}

FLOAT SteepestDescent::Step(const Vector &xk, const Vector &dk) const
{
    _line_search->SetTargetFunctor(_functor);
    _line_search->SetXk(xk);
    _line_search->SetDk(dk);

    FLOAT alpha = _line_search->Search(1.0, *_options);
    return alpha;
}

} // namespace NOL
