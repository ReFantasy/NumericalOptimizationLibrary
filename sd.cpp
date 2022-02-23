#include "pqf_by_sd.h"
#include <iostream>
#include "line_search.h"


Vector SteepestDescent(TargetFunctor& functor, Vector x0, FLOAT _gk_norm)
{
    int k = 0;
    Vector xk = x0;

    double gk_norm;
    while ((gk_norm = functor.FirstOrderDerivatives(xk).norm()) >= _gk_norm)
    {
#ifdef _DEBUG
        std::cout << k << ": (" << xk.transpose() << ")"
                  << "    " << gk_norm << std::endl;
#endif

        Vector _gk = functor.FirstOrderDerivatives(xk);
        Vector dk = -_gk;

        // page 19
        // search g(k+1)^T*dk=0 -> alpha
        class Line :public LineSearch
        {
        public:
            Line(TargetFunctor& _functor) :functor(_functor) {}
            virtual FLOAT Phi(FLOAT a)
            {
                return functor.FirstOrderDerivatives(xk + a * dk).transpose() * dk;
            }

            virtual FLOAT dPhi_dx(FLOAT a)
            {
                return (functor.SecondOrderDerivatives(xk + a * dk) * dk).transpose() * dk;
            }

            Vector xk;
            Vector dk;
        protected:

            bool Criterion(FLOAT a){return true;}

            TargetFunctor& functor;
            
        };
        Line line_search(functor);
        line_search.xk = xk;
        line_search.dk = dk;
        FLOAT alpha = line_search.Zerosixeight(1.0);
       

        xk = xk + alpha * dk;

        k++;
    }

#ifdef _DEBUG
    std::cout << k << ": (" << xk.transpose() << ")"
              << "    " << gk_norm << std::endl;
#endif

    return xk;
}

