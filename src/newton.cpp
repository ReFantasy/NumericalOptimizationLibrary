#include "newton.h"
#include "line_search.h"
#include <iostream>

namespace NOL
{

bool NewtonBase::IsTerminated(const Vector &xk, int k) const
{
    // FLOAT xk_max_norm = _functor->FirstOrderDerivatives(xk).cwiseAbs().maxCoeff();
    FLOAT xk_max_norm = _functor->Gradient(xk).norm();
    if (xk_max_norm <= _options->gk_norm)
        return true;
    // <--------
    *_options << "k:" << k << " "
              << "  xk:(" << xk.transpose() << ") "
              << "  ||gk||: " << xk_max_norm << "\n";
    // -------->
    return false;
}

NOL::Vector NewtonBase::SearchDirection(const Vector &xk) const
{
    // compute dk
    Matrix Gk = _functor->Hesse(xk);
    Vector gk = _functor->Gradient(xk);

    // solve Gk*dk = -gk
    Vector dk;
    dk = Gk.colPivHouseholderQr().solve(-gk);

    return dk;
}

FLOAT DampedNewton::Step(const Vector &xk, const Vector &dk) const
{
    _line_search->SetTargetFunctor(_functor);
    _line_search->SetXk(xk);
    _line_search->SetDk(dk);
    return  _line_search->Search(1.0, *_options);
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
        if (IsTerminated(xk, k))
            break;

        Vector dk = SearchDirection(xk);

        double alpha = Step(xk, dk);

        Vector yk1 = _functor->Gradient(xk);
        Vector yk2 = _functor->Gradient(xk + alpha * dk);
        xk = xk + alpha * dk;

        // Hk->Hk+1
        Vector sk = alpha * dk;
        Vector yk = yk2 - yk1;
        _Hk = CorrectHk(_Hk, sk, yk);

        k++;
    }

    return xk;
}

NOL::Vector QuasiNewton::SearchDirection(const Vector &xk) const
{
    Vector gk = _functor->Gradient(xk);
    Vector dk = -_Hk * gk;
    return dk;
}

Matrix QuasiNewton::CorrectHk(Matrix Hk, Vector sk, Vector yk)
{
    if (_options->quasi_newton_type == QuasiNewtonSearchType::SR1)
    {
        // SR1
        Vector tmp = sk - Hk * yk;
        return Hk + (tmp * tmp.transpose()) / (tmp.transpose() * yk);
    }
    else if (_options->quasi_newton_type == QuasiNewtonSearchType::DFP)
    {
        // DFP
        return Hk + (sk * sk.transpose()) / (sk.transpose() * yk) -
               (Hk * yk * yk.transpose() * Hk) / (yk.transpose() * Hk * yk);
    }
    else if (_options->quasi_newton_type == QuasiNewtonSearchType::BFGS)
    {
        // BFGS
        FLOAT d = yk.transpose() * sk;
        FLOAT a = 1.0 + (FLOAT)(yk.transpose() * Hk * yk) / d;
        Matrix b = (sk * sk.transpose()) / d;
        Matrix c = sk * yk.transpose() * Hk + Hk * yk * sk.transpose();
        
        return Hk + a * b - c / d;
    }

    throw std::invalid_argument("invalid quasi-newton type");
}

} // namespace NOL
