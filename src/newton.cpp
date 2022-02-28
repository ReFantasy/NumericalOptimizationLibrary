#include "newton.h"
#include "line_search.h"
#include <iostream>

namespace NOL
{

bool NewtonBase::IsTermination(const Vector &xk, int k) const
{
    // FLOAT xk_max_norm = _functor->FirstOrderDerivatives(xk).cwiseAbs().maxCoeff();
    FLOAT xk_max_norm = _functor->FirstOrderDerivatives(xk).norm();
    if (xk_max_norm <= _options->gk_norm)
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
    FLOAT alpha = 1.0;

    _line_search->_functor = _functor;
    _line_search->xk = xk;
    _line_search->dk = dk;
    alpha = _line_search->Search(alpha, _options->_line_search_type);
    return alpha;
}

NOL::Vector QuasiNewton::Solve()
{
    int k = 0;
    Vector xk = _options->init_x;
    _Hk = Matrix::Identity(xk.size(), xk.size());

    // <--------
    *_options << typeid(*this).name() << " initial x: ";
    *_options << xk.transpose() << "\n\n";
    // -------->

    while (true)
    {
        if (IsTermination(xk, k))
            break;

        Vector dk = DescentDirection(xk);

        double alpha = StepSize(xk, dk);

        Vector yk1 = _functor->FirstOrderDerivatives(xk);
        Vector yk2 = _functor->FirstOrderDerivatives(xk + alpha * dk);
        xk = xk + alpha * dk;

        // Hk->Hk+1
        Vector sk = alpha * dk;
        Vector yk = yk2 - yk1;
        _Hk = CorrectHk(_Hk, sk, yk);

        k++;
    }

    return xk;
}

NOL::Vector QuasiNewton::DescentDirection(const Vector &xk) const
{
    Vector gk = _functor->FirstOrderDerivatives(xk);
    Vector dk = -_Hk * gk;
    return dk;
}

Matrix QuasiNewton::CorrectHk(Matrix Hk, Vector sk, Vector yk)
{
    if (_options->_quasi_newton_type == QuasiNewtonType::SR1)
    {
        // SR1
        Vector tmp = sk - Hk * yk;
        return Hk + (tmp * tmp.transpose()) / (tmp.transpose() * yk);
    }
    else if (_options->_quasi_newton_type == QuasiNewtonType::DFP)
    {
        // DFP
        return Hk + (sk * sk.transpose()) / (sk.transpose() * yk) -
               (Hk * yk * yk.transpose() * Hk) / (yk.transpose() * Hk * yk);
    }
    else if (_options->_quasi_newton_type == QuasiNewtonType::BFGS)
    {
        // BFGS
        FLOAT a = 1.0 + (FLOAT)(yk.transpose() * Hk * yk) / (yk.transpose() * sk);
        Matrix b = (sk * sk.transpose()) / (yk.transpose() * sk);
        Matrix c = sk * yk.transpose() * Hk + Hk * yk * sk.transpose();
        FLOAT d = yk.transpose() * sk;

        return Hk + a * b - c / d;
    }

    throw std::invalid_argument("invalid quasi-newton type");
}

} // namespace NOL
