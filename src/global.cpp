#include "global.h"
#include "line_search.h"

namespace NOL
{

Vector UnconstrainedOptimizationLineSearchBase::Solve()
{
    int k = 0;
    Vector xk = _options->init_x;

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
    switch (_options->termination_type)
    {
    case TerminationCriterionType::GK_NORM:
        if ((epsilon = _functor->Gradient(xk).norm()) < _options->termination_value)
            is_termination = true;
        break;
    case TerminationCriterionType::DELTA_XK:
        if ((epsilon = (xk - last_xk).norm()) < _options->termination_value)
            is_termination = true;
        break;
    case TerminationCriterionType::DELTA_F:
        if ((epsilon = std::abs((*_functor)(xk) - (*_functor)(last_xk))) < _options->termination_value)
            is_termination = true;
        break;
    default:
        break;
    }
    if (_timer.Elapse() / 1000.0 > _options->max_solver_time_in_seconds)
        is_termination = true;
    if (is_termination)
        return true;

    last_xk = xk;
    // <--------
    *_options << "k:" << k << " "
              << "  xk:(" << xk.transpose() << ") "
              << "  estimation: " << epsilon << "\n";
    // -------->

    return false;
}

} // namespace NOL
