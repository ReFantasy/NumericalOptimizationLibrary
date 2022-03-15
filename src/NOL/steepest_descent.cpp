#include "steepest_descent.h"
#include "line_search.h"
#include <iostream>

namespace NOL
{

NOL::Vector SteepestDescent::SearchDirection(const Vector &xk) const
{
    Vector gk = _functor_ptr->Gradient(xk);
    Vector dk = -gk;
    return dk;
}

FLOAT SteepestDescent::Step(const Vector &xk, const Vector &dk) const
{
    _line_search_ptr->SetXk(xk);
    _line_search_ptr->SetDk(dk);

    if (_options_ptr->line_search_type != LineSearchType::GOLDENSECTION)
        _options_ptr->line_search_type = LineSearchType::GOLDENSECTION;
    FLOAT alpha = _line_search_ptr->Search(1.0, *_options_ptr);
    return alpha;
}

} // namespace NOL
