#include "global.h"
#include "line_search.h"

namespace NOL
{

Vector UnconstrainedOptimizationLineSearchBase::Solve()
{
    int k = 0;
    Vector xk = _options->init_x;

    LineSearch line_search{};
    line_search._functor = _functor;

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

        xk = xk + alpha * dk;
        k++;
    }

    return xk;
}

} // namespace NOL