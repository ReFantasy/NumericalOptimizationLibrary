#include "steepest_descent.h"
#include "line_search.h"
#include <iostream>

namespace NOL
{

SteepestDescent::SteepestDescent(TargetFunctor *functor)
{
    _functor = functor;
    _line_search = new LinearSearch{};
    _line_search->SetTargetFunctor(_functor);
}

SteepestDescent::SteepestDescent(TargetFunctor *functor, Options *options) : SteepestDescent(functor)
{
    _options = options;
}

SteepestDescent::~SteepestDescent()
{
    if (_line_search)
        delete _line_search;
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

    FLOAT alpha = _line_search->Search(1.0, *_options);
    return alpha;
}

} // namespace NOL
