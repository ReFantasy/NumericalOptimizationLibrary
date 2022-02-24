#include "sd.h"
#include "line_search.h"
#include <iostream>

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
    FLOAT dPhi_dx(FLOAT a) override
    {
        return (_functor->SecondOrderDerivatives(xk + a * dk) * dk).transpose() * dk;
    }
};

Vector SteepestDescent(TargetFunctor &functor, Vector x0, FLOAT _gk_norm)
{
    int k = 0;
    Vector xk = x0;

    double gk_norm;
    while ((gk_norm = functor.FirstOrderDerivatives(xk).norm()) >= _gk_norm)
    {
        Vector _gk = functor.FirstOrderDerivatives(xk);
        Vector dk = -_gk;

        static LineSearchForSD line_search(&functor);
        line_search.xk = xk;
        line_search.dk = dk;
        FLOAT alpha = line_search.Zerosixeight(1.0);

        xk = xk + alpha * dk;

        k++;
    }

    return xk;
}
