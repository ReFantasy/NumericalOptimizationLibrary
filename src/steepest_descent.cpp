#include "steepest_descent.h"
#include "line_search.h"
#include <iostream>

namespace NOL
{

SteepestDescent::SteepestDescent(TargetFunctor* functor)
{
    _functor = functor;
    _line_search = new LinearSearch{};
    _line_search->SetTargetFunctor(_functor);
}

SteepestDescent::SteepestDescent(TargetFunctor* functor, Options* options):SteepestDescent(functor)
{
    _options = options;
}

SteepestDescent::~SteepestDescent()
{
    if (_line_search)
        delete _line_search;
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
    Vector gk = _functor->Gradient(xk);
    Vector dk = -gk;
    return dk;
}

FLOAT SteepestDescent::Step(const Vector &xk, const Vector &dk) const
{
    _line_search->SetXk(xk);
    _line_search->SetDk(dk);

    _options->line_search_type = LineSearchType::GOLDENSECTION;

    FLOAT alpha = _line_search->Search(1.0, *_options);
    return alpha;
}

} // namespace NOL
