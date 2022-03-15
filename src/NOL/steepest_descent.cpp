#include "steepest_descent.h"
#include "line_search.h"
#include <iostream>

namespace NOL
{

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

	if(_options->line_search_type!= LineSearchType::GOLDENSECTION)
		_options->line_search_type = LineSearchType::GOLDENSECTION;
    FLOAT alpha = _line_search->Search(1.0, *_options);
    return alpha;
}

} // namespace NOL
