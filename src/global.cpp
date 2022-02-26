#include "global.h"
#include "line_search.h"

namespace NOL
{

Vector UnconstrainedOptimizationLineSearchBase::Solve(TargetFunctor &fucntor, Options &options)
{
    int k = 0;
    Vector xk = options.init_x;

    LineSearch line_search{};
    line_search._functor = &fucntor;

    // <--------
    options << typeid(*this).name() << " initial x: ";
    options << xk.transpose() << "\n\n";
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
