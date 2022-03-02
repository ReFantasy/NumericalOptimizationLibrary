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

} // namespace NOL
