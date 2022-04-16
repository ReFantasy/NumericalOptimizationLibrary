#include "conjugate_gradient.h"
#include "line_search.h"

namespace NOL
{
Vector ConjugateGradient::Solve()
{
    if (_functor_ptr == nullptr)
    {
        throw std::invalid_argument("functor pointer is null.");
    }

    // <--------
    *_options_ptr << typeid(*this).name() << " initial x: ";
    *_options_ptr << std::fixed << _options_ptr->init_x.transpose() << "\n\n";
    // -------->

    int k = 0;
    Vector xk = _options_ptr->init_x;        // x0
    Vector dk = -_functor_ptr->Gradient(xk); // d0

    while (true)
    {
        if (IsTerminated(xk, k))
            break;

        double alpha = Step(xk, dk);

        _gk = _functor_ptr->Gradient(xk);
        _dk = dk;
        xk += alpha * dk;

        if ((k % _options_ptr->conjugate_gradient_restart_num) == 0)
            dk = -_functor_ptr->Gradient(xk);
        else
            dk = SearchDirection(xk);

        k++;
    }

    _K = k;
    return xk;
}

Vector ConjugateGradient::SearchDirection(const Vector &xk) const
{
    Vector gk_add_1 = _functor_ptr->Gradient(xk);
    FLOAT beta_k = 1;

    switch (_options_ptr->conjugate_gradient_type)
    {
    case ConjugateGradientType::FR:
        beta_k = (FLOAT)(gk_add_1.transpose() * gk_add_1) / (_gk.transpose() * _gk);
        break;
    case ConjugateGradientType::PRP:
        beta_k = (FLOAT)(gk_add_1.transpose() * (gk_add_1 - _gk)) / (_gk.transpose() * _gk);
        break;
    case ConjugateGradientType::PRP_PLUS:
        beta_k = std::max((FLOAT)(gk_add_1.transpose() * (gk_add_1 - _gk)) / (_gk.transpose() * _gk), 0.0);
        break;
    case ConjugateGradientType::CD:
        beta_k = -(FLOAT)(gk_add_1.transpose() * gk_add_1) / (_dk.transpose() * _gk);
        break;
    case ConjugateGradientType::DY:
        beta_k = (FLOAT)(gk_add_1.transpose() * gk_add_1) / (_dk.transpose() * (gk_add_1 - _gk));
        break;
    }
    return -gk_add_1 + beta_k * _dk;
}
FLOAT ConjugateGradient::Step(const Vector &xk, const Vector &dk) const
{
    _line_search_ptr->SetXk(xk);
    _line_search_ptr->SetDk(dk);
    return _line_search_ptr->Search(1.0, *_options_ptr);
}

} // namespace NOL
