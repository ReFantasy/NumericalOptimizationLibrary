#include "example.h"
#include "global.h"
#include "line_search.h"

using namespace NOL;

void example_3_1()
{
    struct Functor : TargetFunctor
    {
        virtual FLOAT operator()(const Vector& x) const override
        {
            return true;
        }
        virtual Vector FirstOrderDerivatives(const Vector& x) const override
        {
            return (G + G.transpose()) * x / 2 + b;
        }
        virtual Matrix SecondOrderDerivatives(const Vector& x) const override
        {
            return G;
        }

        Matrix G{ 2, 2 };
        Vector b{ 2 };
        float c = 10;
    };

    Functor functor;

    Vector x0(2);
    x0 << -30, 100;

    functor.G << 21, 4, 4, 15;
    //functor.G << 21, 4, 4, 1;
    functor.b << 2, 3;


    Options options;
    options.init_x = x0;
    options.line_search_type = LineSearchType::GOLDENSECTION;
    options.parameter_line_search_advance_and_retreat_h = 1.5;
    options.parameter_line_search_advance_and_retreat_t = 1.5;
    options.optimized_performance = true;
    SteepestDescent sd;
    sd._functor = &functor;
    sd._options = &options;

    std::cout << sd.Solve() << std::endl;
}

void example_3_2()
{
    struct Functor : TargetFunctor
    {
        virtual FLOAT operator()(const Vector& x) const override
        {
            return 3 * x(0) * x(0) + 3 * x(1) * x(1) - x(0) * x(0) * x(1);
        }
        virtual Vector FirstOrderDerivatives(const Vector& x) const override
        {
            Vector dx(2);
            dx(0) = 6 * x(0) - 2 * x(0) * x(1);
            dx(1) = 6 * x(1) - x(0) * x(0);
            return dx;
        }
        virtual Matrix SecondOrderDerivatives(const Vector& x) const override
        {
            Matrix J(2, 2);
            J(0, 0) = 6 - 2 * x(1);
            J(0, 1) = -2 * x(0);
            J(1, 0) = -2 * x(0);
            J(1, 1) = 6;
            return J;
        }
    };

    Functor functor;


    Vector x0(2);
    // x0(0) = x0(1) = 1.5;
    x0(0) = -2;
    x0(1) = 4;

    Options option;
    option.gk_norm = 10e-6;
    option.init_x = x0;

    NewtonBase newton;

    newton._functor = &functor;
    newton._options = &option;
    LineSearch line_search;
    newton._line_search = &line_search;
    std::cout << "res:   " << newton.Solve().transpose() << std::endl;
}

void example()
{
    struct Functor : TargetFunctor
    {
        virtual FLOAT operator()(const Vector& x) const override
        {
            int n = x.size();
            FLOAT sum = 0;
            for (int i = 1; i < n; i += 2)
            {
                sum = sum + 100 * std::pow(x(i - 1) * x(i - 1) - x(i), 2) + std::pow(x(i - 1) - 1, 2);
            }
            return sum;
        }
        virtual Vector FirstOrderDerivatives(const Vector& x) const override
        {
            Vector v(x.size());

            int n = x.size();
            for (int i = 1; i <= n; i++)
            {
                if (i % 2 == 0)
                {
                    v(i - 1) = 200 * (x(i - 1) - x(i - 2) * x(i - 2));
                }
                else
                {
                    v(i - 1) = 400 * (std::pow(x(i - 1), 2) - x(i)) * x(i - 1) + 2 * x(i - 1) - 2;
                }
            }


            return v;
        }
        virtual Matrix SecondOrderDerivatives(const Vector& x) const override
        {
            return Matrix{};
        }

    };

    Functor functor;

    int n = 4;
    Vector x(n);
    for (int i = 0; i < x.size(); i++)
    {
        if (i % 2 == 0)
        {
            x(i) = -2.2;
        }
        else
        {
            x(i) = 3.1;
        }
    }

    /*Vector x(2);
    x(0) =  1;
    x(1) = 2;
    std::cout << functor.FirstOrderDerivatives(x) << std::endl;*/
    Options options;
    options.line_search_type = LineSearchType::STRONGWOLFE;
    options.quasi_newton_type = QuasiNewtonType::DFP;
    options.parameter_line_search_armijo_rho = 10e-4;
    options.parameter_line_search_strong_wolfe_sigma = 0.1;
    //options.parameter_line_search_strong_wolfe_alpha_max = 1;
    options.init_x = x;

    LineSearch line_search;


    class TQuasiNewton :public QuasiNewton
    {
    public:
        bool IsTerminated(const Vector& xk, int k) const override
        {
            
            static Vector last_xk = xk;
            
            if (k != 0)
            {
                if ((*_functor)(last_xk) - (*_functor)(xk) <= 10e-8)
                    return true;
                last_xk = xk;
            }

            FLOAT xk_max_norm = _functor->FirstOrderDerivatives(xk).cwiseAbs().maxCoeff();

            // <--------
            *_options << "k:" << k << " "
                << "  xk:(" << xk.transpose() << ") "
                << "  ||gk_max_norm||: " << xk_max_norm << "\n";
            // -------->
            return false;
        }
    };

    NewtonBase* newton = new TQuasiNewton;
    newton->_functor = &functor;
    newton->_line_search = &line_search;
    newton->_options = &options;




    std::cout << functor(newton->Solve()) << std::endl;
}
