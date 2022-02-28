#include "example.h"
#include "global.h"
#include "line_search.h"

using namespace NOL;

void example_3_1()
{
    struct Functor : TargetFunctor
    {
        virtual FLOAT operator()(const Vector &x) const override
        {
            return true;
        }
        virtual Vector FirstOrderDerivatives(const Vector &x) const override
        {
            return (G + G.transpose()) * x / 2 + b;
        }
        virtual Matrix SecondOrderDerivatives(const Vector &x) const override
        {
            return G;
        }

        Matrix G{2, 2};
        Vector b{2};
        float c = 10;
    };

    Functor functor;

    Vector x0(2);
    x0 << -30, 100;

    functor.G << 21, 4, 4, 15;
    functor.b << 2, 3;

    // SteepestDescent(functor, x0);

    // printf("\n\n");

    // functor.G << 21, 4, 4, 1;

    Options option;
    option.init_x = x0;
    SteepestDescent sd;
    sd._functor = &functor;
    sd._options = &option;
    LineSearch line_search;
    sd._line_search = &line_search;
    std::cout << sd.Solve() << std::endl;
}

void example_3_2()
{
    struct Functor : TargetFunctor
    {
        virtual FLOAT operator()(const Vector &x) const override
        {
            return 3 * x(0) * x(0) + 3 * x(1) * x(1) - x(0) * x(0) * x(1);
        }
        virtual Vector FirstOrderDerivatives(const Vector &x) const override
        {
            Vector dx(2);
            dx(0) = 6 * x(0) - 2 * x(0) * x(1);
            dx(1) = 6 * x(1) - x(0) * x(0);
            return dx;
        }
        virtual Matrix SecondOrderDerivatives(const Vector &x) const override
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

void example_3_1_by_dampednewton()
{
    struct Functor : TargetFunctor
    {
        virtual FLOAT operator()(const Vector &x) const override
        {
            return (FLOAT)(x.transpose() * G * x) * 0.5 + (FLOAT)(b.transpose() * x) + c;
        }
        virtual Vector FirstOrderDerivatives(const Vector &x) const override
        {
            return (G + G.transpose()) * x / 2 + b;
        }
        virtual Matrix SecondOrderDerivatives(const Vector &x) const override
        {
            return G;
        }

        Matrix G{2, 2};
        Vector b{2};
        float c = 10;
    };
    Functor functor;

    Vector b(2);
    b << 2, 3;
    float c = 10;
    Matrix G(2, 2);
    //G << 21, 4, 4, 1;
    G << 21, 4, 4, 15;

    Vector x0(2);
    x0 << -30, 100;

    functor.G = G;
    functor.b = b;
    functor.c = c;

    // 阻尼牛顿法求解
    
    Options option;
    option.gk_norm = 10e-6;
    option.init_x = x0;
    DampedNewton newton;
    newton._functor = &functor;
    newton._options = &option;
    LineSearch line_search;
    newton._line_search = &line_search;
    Vector res = newton.Solve();
    std::cout << "res:   " << res.transpose() << std::endl;
}

void example_test_1()
{
    struct Functor : TargetFunctor
    {
        virtual FLOAT operator()(const Vector &x) const override
        {
            return (x(0) - 4) * (x(0) - 4) - 3;
        }
        virtual Vector FirstOrderDerivatives(const Vector &x) const override
        {
            Vector dx(1);
            dx(0) = 2 * (x(0) - 4);
            return dx;
        }
        virtual Matrix SecondOrderDerivatives(const Vector &x) const override
        {
            Matrix d2(1, 1);
            d2(0, 0) = 2;
            return d2;
        }
    };

    Functor functor;

    Options option;
    Vector x0(1);
    x0(0) = 10;
    option.init_x = x0;
    option.gk_norm = 10e-8;
    //option.optimized_performance = true;
    SteepestDescent sd;
    sd._functor = &functor;
    sd._options = &option;
    LineSearch line_search;
    sd._line_search = &line_search;
    sd.Solve();
}
