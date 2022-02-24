#include "example.h"
#include "global.h"
#include "line_search.h"

void example_3_1()
{
    struct Functor : TargetFunctor
    {
        virtual TFLOAT operator()(const TVector &x) const override
        {
            return true;
        }
        virtual TVector FirstOrderDerivatives(const TVector &x) const override
        {
            return (G + G.transpose()) * x / 2 + b;
        }
        virtual TMatrix SecondOrderDerivatives(const TVector &x) const override
        {
            return G;
        }

        TMatrix G{2, 2};
        TVector b{2};
        float c = 10;
    };

    Functor functor;

    TVector x0(2);
    x0 << -30, 100;

    functor.G << 21, 4, 4, 15;
    functor.b << 2, 3;

    // SteepestDescent(functor, x0);

    // printf("\n\n");

    // functor.G << 21, 4, 4, 1;

    Options option;
    option.init_x0 = x0;
    SteepestDescent sd;
    sd.Solve(functor, option);
}

void example_3_2()
{
    struct Functor : TargetFunctor
    {
        virtual TFLOAT operator()(const TVector &x) const override
        {
            return 3 * x(0) * x(0) + 3 * x(1) * x(1) - x(0) * x(0) * x(1);
        }
        virtual TVector FirstOrderDerivatives(const TVector &x) const override
        {
            TVector dx(2);
            dx(0) = 6 * x(0) - 2 * x(0) * x(1);
            dx(1) = 6 * x(1) - x(0) * x(0);
            return dx;
        }
        virtual TMatrix SecondOrderDerivatives(const TVector &x) const override
        {
            TMatrix J(2, 2);
            J(0, 0) = 6 - 2 * x(1);
            J(0, 1) = -2 * x(0);
            J(1, 0) = -2 * x(0);
            J(1, 1) = 6;
            return J;
        }
    };

    Functor functor;

    NewtonBase newton;
    TVector x0(2);
    // x0(0) = x0(1) = 1.5;
    x0(0) = -2;
    x0(1) = 4;

    Options option;
    option.gk_norm = 10e-6;
    option.init_x0 = x0;

    std::cout << "res:   " << newton.Solve(functor, option).transpose() << std::endl;
}

void example_3_1_by_dampednewton()
{
    struct Functor : TargetFunctor
    {
        virtual TFLOAT operator()(const TVector &x) const override
        {
            return (TFLOAT)(x.transpose() * G * x) * 0.5 + (TFLOAT)(b.transpose() * x) + c;
        }
        virtual TVector FirstOrderDerivatives(const TVector &x) const override
        {
            return (G + G.transpose()) * x / 2 + b;
        }
        virtual TMatrix SecondOrderDerivatives(const TVector &x) const override
        {
            return G;
        }

        TMatrix G{2, 2};
        TVector b{2};
        float c = 10;
    };
    Functor functor;

    TVector b(2);
    b << 2, 3;
    float c = 10;
    TMatrix G(2, 2);
    // G << 21, 4, 4, 1;
    G << 21, 4, 4, 15;

    TVector x0(2);
    x0 << -30, 100;

    functor.G = G;
    functor.b = b;
    functor.c = c;

    // 阻尼牛顿法求解
    DampedNewton newton;
    Options option;
    option.gk_norm = 10e-6;
    option.init_x0 = x0;
    std::cout << "res:   " << newton.Solve(functor, option).transpose() << std::endl;
}



void example_test_1()
{
    struct Functor : TargetFunctor
    {
        virtual TFLOAT operator()(const TVector& x) const override
        {
            return (x(0) - 4) * (x(0) - 4) - 3;
        }
        virtual TVector FirstOrderDerivatives(const TVector& x) const override
        {
            TVector dx(1);
            dx(0) = 2 * (x(0) - 4);
            return dx;
        }
        virtual TMatrix SecondOrderDerivatives(const TVector& x) const override
        {
            TMatrix d2(1, 1);
            d2(0, 0) = 2;
            return d2;
        }
    };

    Functor functor;

    Options option;
    TVector x0(1);
    x0(0) = 10;
    option.init_x0 = x0;
    option.gk_norm = 10e-8;
    SteepestDescent sd;
    sd.Solve(functor, option);
}
