#include "steepest_descent.h"
#include "line_search.h"
#include <iostream>

namespace NOL
{
// page 19
// search g(k+1)^T*dk=0 -> alpha
class LineSearchForSD : public LineSearch
{
  public:
    LineSearchForSD(TargetFunctor *functor = nullptr) : LineSearch(functor)
    {
    }

  protected:
    FLOAT Phi(FLOAT a) override
    {
        return _functor->FirstOrderDerivatives(xk + a * dk).transpose() * dk;
    }
    FLOAT dPhi_da(FLOAT a) override
    {
        return (_functor->SecondOrderDerivatives(xk + a * dk) * dk).transpose() * dk;
    }
};

Vector SteepestDescent::Solve(TargetFunctor &functor, Options &option)
{
    int k = 0;
    Vector xk = option.init_x0;

    // <--------
    option << "Steepest Descent with initial x: ";
    option << xk.transpose() << "\n\n";
    // -------->

    double gk_norm = 0;
    while (true)
    {
        gk_norm = functor.FirstOrderDerivatives(xk).norm();
        if (gk_norm < option.gk_norm)
            break;
        // <--------
        option << "k:" << k << " "
               << "  xk:(" << xk.transpose() << ") "
               << "  ||gk||: " << gk_norm << "\n";
        // -------->

        Vector _gk = functor.FirstOrderDerivatives(xk);
        Vector dk = -_gk;

        static LineSearchForSD line_search(&functor);
        line_search.xk = xk;
        line_search.dk = dk;
        FLOAT alpha = line_search.Zerosixeight(1.0);

        xk = xk + alpha * dk;

        k++;
    }

    // <--------
    option << "k:" << k << " "
           << "  xk:(" << xk.transpose() << ") "
           << "  ||gk||: " << gk_norm << "\n";
    // -------->
    return xk;
}

} // namespace NOL
