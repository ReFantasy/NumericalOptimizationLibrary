#include "newton.h"
#include "line_search.h"
#include <iostream>

namespace NOL
{

bool NewtonBase::IsTermination(const Vector &xk, int k) const
{
    FLOAT xk_max_norm = _functor.FirstOrderDerivatives(xk).cwiseAbs().maxCoeff();
    if (xk_max_norm < _options.gk_norm)
        return true;
    // <--------
    _options << "k:" << k << " "
             << "  xk:(" << xk.transpose() << ") "
             << "  ||gk||: " << xk_max_norm << "\n";
    // -------->
    return false;
}

NOL::Vector NewtonBase::DescentDirection(const Vector &xk) const
{
    // compute dk
    Matrix Gk = _functor.SecondOrderDerivatives(xk);
    Vector gk = _functor.FirstOrderDerivatives(xk);

    // solve Gk*dk = -gk
    Vector dk;
    dk = Gk.colPivHouseholderQr().solve(-gk);

    return dk;
}

FLOAT DampedNewton::StepSize(const Vector &xk, const Vector &dk) const
{
    static FLOAT alpha = 10;
    LineSearch line_search{};
    line_search._functor = &_functor;
    line_search.xk = xk;
    line_search.dk = dk;
    alpha = line_search.QuadraticPolynomialInterpolation(alpha);
    return alpha;
}

} // namespace NOL
