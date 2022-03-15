#include "global.h"
#include "line_search.h"

namespace NOL
{
UnconstrainedOptimizationLineSearchBase::UnconstrainedOptimizationLineSearchBase()
{
	_line_search_ptr = std::make_shared<LinearSearch>();
	_options_ptr = std::make_shared<Options>();
}
UnconstrainedOptimizationLineSearchBase::UnconstrainedOptimizationLineSearchBase(std::shared_ptr<TargetFunctor>  functor_ptr)
{
    _functor_ptr = functor_ptr;
    _line_search_ptr = std::make_shared<LinearSearch>();
    _line_search_ptr->SetTargetFunctor(functor_ptr);
	_options_ptr = std::make_shared<Options>();
}

Vector UnconstrainedOptimizationLineSearchBase::Solve()
{
    if (_functor_ptr == nullptr)
    {
        throw std::invalid_argument("functor pointer is null.");
    }

    int k = 0;
    Vector xk = _options_ptr->init_x;

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

        xk = xk + alpha * dk;
        k++;
    }

    return xk;
}

bool UnconstrainedOptimizationLineSearchBase::IsTerminated(const Vector &xk, int k) const
{
    static Vector last_xk(xk.size());
    if (k == 0)
    {
        _timer.ReSet();
        last_xk(0) = std::numeric_limits<FLOAT>::max();
    }
    bool is_termination = false;
    FLOAT epsilon = 0;
    switch (_options_ptr->termination_type)
    {
    case TerminationCriterionType::GK_NORM:
        if ((epsilon = _functor_ptr->Gradient(xk).norm()) < _options_ptr->termination_value)
            is_termination = true;
        break;
    case TerminationCriterionType::DELTA_XK:
        if ((epsilon = (xk - last_xk).norm()) < _options_ptr->termination_value)
            is_termination = true;
        break;
    case TerminationCriterionType::DELTA_F:
        if ((epsilon = std::abs((*_functor_ptr)(xk) - (*_functor_ptr)(last_xk))) < _options_ptr->termination_value)
            is_termination = true;
        break;
    default:
        break;
    }
    if (_timer.Elapse() / 1000.0 > _options_ptr->max_solver_time_in_seconds)
        is_termination = true;
    if (is_termination)
        return true;

    last_xk = xk;
    // <--------
    *_options_ptr << "k:" << k << " "
              << "  xk:(" << xk.transpose() << ") "
              << "  estimation: " << epsilon << "\n";
    // -------->

    return false;
}

	//UnconstrainedOptimizationLineSearchBase::~UnconstrainedOptimizationLineSearchBase()


} // namespace NOL
