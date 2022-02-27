#include "newton.h"
#include "line_search.h"
#include <iostream>

namespace NOL
{

bool NewtonBase::IsTermination(const Vector &xk, int k) const
{
    FLOAT xk_max_norm = _functor->FirstOrderDerivatives(xk).cwiseAbs().maxCoeff();
    if (xk_max_norm < _options->gk_norm)
        return true;
    // <--------
    *_options << "k:" << k << " "
              << "  xk:(" << xk.transpose() << ") "
              << "  ||gk||: " << xk_max_norm << "\n";
    // -------->
    return false;
}

NOL::Vector NewtonBase::DescentDirection(const Vector &xk) const
{
    // compute dk
    Matrix Gk = _functor->SecondOrderDerivatives(xk);
    Vector gk = _functor->FirstOrderDerivatives(xk);

    // solve Gk*dk = -gk
    Vector dk;
    dk = Gk.colPivHouseholderQr().solve(-gk);

    return dk;
}

FLOAT DampedNewton::StepSize(const Vector &xk, const Vector &dk) const
{
    static FLOAT alpha = 10;
    
    _line_search->_functor = _functor;
    _line_search->xk = xk;
    _line_search->dk = dk;
    alpha = _line_search->QuadraticPolynomialInterpolation(alpha);
    //alpha = _line_search->CubicPolynomialInterpolation(alpha);
    return alpha;
}

} // namespace NOL
