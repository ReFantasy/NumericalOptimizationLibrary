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
        virtual Vector Gradient(const Vector& x) const override
        {
            return (G + G.transpose()) * x / 2 + b;
        }
        virtual Matrix Hesse(const Vector& x) const override
        {
            return G;
        }

        Matrix G{ 2, 2 };
        Vector b{ 2 };
        float c = 10;
    };

    Functor functor;

    Vector init_x(2);
    init_x << -30, 100;

    functor.G << 21, 4, 4, 15;
    //functor.G << 21, 4, 4, 1;
    functor.b << 2, 3;


    Options options;
    options.init_x = init_x;
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
        virtual Vector Gradient(const Vector& x) const override
        {
            Vector dx(2);
            dx(0) = 6 * x(0) - 2 * x(0) * x(1);
            dx(1) = 6 * x(1) - x(0) * x(0);
            return dx;
        }
        virtual Matrix Hesse(const Vector& x) const override
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
    LinearSearch line_search;
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

            //return 100 * std::pow((x(1) - x(0)* x(0)),2)  + std::pow((1 - x(0)),2);
        }
        virtual Vector Gradient(const Vector& x) const override
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
            
            //v(0) = -400 * x(0) * x(1) + 400 * std::pow(x(0), 3) + 2 * x(0) - 2;
            //v(1) = 200 * x(1) - 200 * x(0) * x(0);

            return v;
        }
        virtual Matrix Hesse(const Vector& x) const override
        {
            return Matrix{};
        }

    };

    Functor functor;

    int n = 8;
    Vector x(n);
    //x << -10, 10;
    for (int i = 0; i < x.size(); i++)
    {
        if (i % 2 == 0)
        {
            x(i) = rand()%100;
        }
        else
        {
            x(i) = -rand()%100;
        }
        std::cout << x(i) << ", " ;
    }
    std::cout << std::endl;
    

    /*Vector x(2);
    x(0) =  1;
    x(1) = 2;
    std::cout << functor.FirstOrderDerivatives(x) << std::endl;*/
    Options options;
    options.line_search_type = LineSearchType::STRONGWOLFE;
    options.quasi_newton_type = QuasiNewtonSearchType::BFGS;
    options.parameter_line_search_armijo_rho = 1e-3;
    options.parameter_line_search_armijo_t = 3.0;
    options.parameter_line_search_wolfe_rho = 1e-3;
    options.parameter_line_search_wolfe_sigma = 0.01;
    //options.parameter_line_search_wolfe_alpha_max = 2000;
    //options.optimized_performance = true;
    options.init_x = x;

    LinearSearch line_search;


    class TQuasiNewton :public QuasiNewton
    {
    public:
        bool IsTerminated(const Vector& xk, int k) const override
        {
            
            static Vector last_xk = xk;
            
            if (k != 0)
            {
                /*if ((*_functor)(last_xk) - (*_functor)(xk) <= 10e-8)
                    return true;*/
                /*if ((*_functor)(xk) < 1e-5)
                    return true;*/
                if (_functor->Gradient(xk).norm() < 1e-7)
                    return true;
                last_xk = xk;
            }

            FLOAT xk_max_norm = _functor->Gradient(xk).cwiseAbs().maxCoeff();

            // <--------
            *_options << "k:" << k << " "
                << "  xk:(" << xk.transpose() << ") "
                << "  ||gk_max_norm||: " << xk_max_norm <<"   fvalue:"<< (*_functor)(xk) << "\n";
            // -------->
            return false;
        }
    };

    NewtonBase* newton = new TQuasiNewton;
    newton->_functor = &functor;
    newton->_line_search = &line_search;
    newton->_options = &options;

	Timer timer;
	Vector res = newton->Solve();
	std::cout<<timer.Elapse()<<"ms"<<std::endl;
    //std::cout << functor(newton->Solve()) << std::endl;
}
