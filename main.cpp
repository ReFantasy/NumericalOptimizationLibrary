#include <iostream>

#include "example.h"
#include <sstream>
#include "global.h"

using namespace NOL;



int main(int argc, char **argv)
{
 
    //std::cout << QuadraticInterpolationMinimum(-4, 2) << std::endl;
    //example_3_1();
    //example_3_2();
    //example_3_1_by_dampednewton();


    struct Functor : public TargetFunctor
    {
    public:
        virtual FLOAT operator()(const Vector& x) const override
        {
            //return 100 * std::pow(2 * x(0) * x(0) + 2, 2) + x(1) * x(1);
            return 100 *  std::pow((x(0) * x(0) - x(1)),2) + std::pow(x(0) - 1,2);
        }
        virtual Vector FirstOrderDerivatives(const Vector& x) const override
        {
            Vector v(2);
            v(0) = 400 * x(0) * (x(0)*x(0) - x(1)) + 2 * (x(0) - 1);
            v(1) = -200 * (x(0) *x(0) - x(1));
            return v;
        }
        virtual Matrix SecondOrderDerivatives(const Vector& x) const override
        {
            Matrix m = Matrix::Identity(2, 2);
            m(0, 0) = 16 * (3 * x(0) * x(0) + 1);
            m(1, 1) = 2;
            return m;
        }

    };

    Functor functor;
    Vector x0(2);
    x0 << -1.2, 1;

    Options option;
    option.gk_norm = 1e-5;
    option.init_x = x0;
    //option.optimized_performance = true;
    option.quasi_newton_type = QuasiNewtonType::BFGS;
    option.line_search_type = LineSearchType::STRONGWOLFE;

    NewtonBase *newton = new QuasiNewton;

    newton->_functor = &functor;
    newton->_options = &option;
    LineSearch line_search;
    newton->_line_search = &line_search;


    std::cout << functor(newton->Solve()) << std::endl;

    return 0;
}
