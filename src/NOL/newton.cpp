#include "newton.h"
#include "helper.h"
#include "line_search.h"
#include <iostream>

namespace NOL
{

NOL::Vector NewtonBase::SearchDirection(const Vector &xk) const
{
    // compute dk
    Matrix Gk = _functor_ptr->Hesse(xk);
    Vector gk = _functor_ptr->Gradient(xk);

    // solve Gk*dk = -gk
    Vector dk;
    dk = Gk.colPivHouseholderQr().solve(-gk);

    return dk;
}

Vector LM::SearchDirection(const Vector &xk) const
{
    // compute dk
    Matrix Gk = _functor_ptr->Hesse(xk);
    Vector gk = _functor_ptr->Gradient(xk);

    //
    while (Gk.determinant() < MinStepSize<FLOAT>::value)
    {
        Gk += _vk * Matrix::Identity(Gk.rows(), Gk.cols());
        UpdateVk();
    }

    // solve Gk*dk = -gk
    Vector dk;
    dk = Gk.colPivHouseholderQr().solve(-gk);

    return dk;
}

FLOAT DampedNewton::Step(const Vector &xk, const Vector &dk) const
{
    _line_search_ptr->SetXk(xk);
    _line_search_ptr->SetDk(dk);
    return _line_search_ptr->Search(1.0, *_options_ptr);
}

NOL::Vector QuasiNewton::Solve()
{
    if (_functor_ptr == nullptr)
    {
        throw std::invalid_argument("functor pointer is null.");
    }

    int k = 0;
    Vector xk = _options_ptr->init_x;
    _Hk = Matrix::Identity(xk.size(), xk.size());

    // <--------
    *_options_ptr << typeid(*this).name() << " initial x: ";
    *_options_ptr << xk.transpose() << "\n\n";
    // -------->

    while (true)
    {
        if (IsTerminated(xk, k))
            break;

        Vector dk = SearchDirection(xk);

        double alpha = Step(xk, dk);

        _Hk = UpdateHk(_Hk, xk, dk, alpha);

        xk = xk + alpha * dk;
        k++;
    }

	_K = k;
    return xk;
}

NOL::Vector QuasiNewton::SearchDirection(const Vector &xk) const
{
    Vector gk = _functor_ptr->Gradient(xk);
    Vector dk = -_Hk * gk;
    return dk;
}

Matrix QuasiNewton::UpdateHk(const Matrix& Hk, const Vector& xk, const Vector& dk, FLOAT alpha)
{
    Vector yk1 = _functor_ptr->Gradient(xk);
    Vector yk2 = _functor_ptr->Gradient(xk + alpha * dk);
    Vector sk = alpha * dk;
    Vector yk = yk2 - yk1;

    if (_options_ptr->quasi_newton_type == QuasiNewtonSearchType::SR1)
    {
        // SR1
        Vector tmp = sk - Hk * yk;
        return Hk + (tmp * tmp.transpose()) / (tmp.transpose() * yk);
    }
    else if (_options_ptr->quasi_newton_type == QuasiNewtonSearchType::DFP)
    {
        // DFP
        return Hk + (sk * sk.transpose()) / (sk.transpose() * yk) -
               (Hk * yk * yk.transpose() * Hk) / (yk.transpose() * Hk * yk);
    }
    else if (_options_ptr->quasi_newton_type == QuasiNewtonSearchType::BFGS)
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
