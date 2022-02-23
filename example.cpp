#include "example.h"
#include "global.h"
#include "line_search.h"

void example_3_1()
{

    struct Functor :TargetFunctor
    {
        virtual FLOAT operator()(const Vector& x) const override
        {
            return true;
        }
        virtual Vector FirstOrderDerivatives(const Vector& x) const override
        {
            return (G + G.transpose())* x / 2 + b;
        }
        virtual Matrix SecondOrderDerivatives(const Vector& x) const override
        {
            return G;
        }

        Matrix G{ 2, 2 };
        Vector b{2};
        float c = 10;
    };

    Functor functor;

    Vector x0(2);
    x0 << -30, 100;

    functor.G << 21, 4, 4, 15;
    functor.b << 2, 3;
   
 

    //SteepestDescent(functor, x0);


    printf("\n\n");

    functor.G << 21, 4, 4, 1;

    std::cout << SteepestDescent(functor, x0) << std::endl;

    
}

void example_3_2()
{
    class NewtonMethod : public NewtonBase
    {
      public:
        Vector gk(Vector xk) override
        {
            Vector _gk(2);
            _gk(0) = 6 * xk(0) - 2 * xk(0) * xk(1);
            _gk(1) = 6 * xk(1) - xk(0) * xk(0);
            return _gk;
        }
        Matrix Gk(Vector xk) override
        {
            Matrix _GK(2, 2);
            _GK << 6 - 2 * xk(1), -2 * xk(0), -2 * xk(0), 6;
            return _GK;
        }
    };

    NewtonMethod newton;
    Vector x0(2);
    x0(0) = x0(1) = 1.5;
    // x0(0) = -2;
    // x0(1) = 4;

    Vector x(2);
    x = newton.Solve(x0, 10e-6);
    std::cout << "res:   " << x.transpose() << std::endl;
}

void Newton_example_3_1()
{
    Vector b(2);
    b << 2, 3;
    float c = 10;
    Matrix G(2, 2);

    Vector x0(2);
    x0 << -30, 100;

    // G << 21, 4, 4, 1;
    G << 21, 4, 4, 15;

    // 阻尼牛顿法求解

    class LS : public LineSearch
    {
      public:
        FLOAT Phi(FLOAT a) override
        {
            return f(xk + a * dk);
        }

        FLOAT dPhi_dx(FLOAT a) override
        {
            return (G * (xk + a * dk) + b).transpose() * dk;
        }

        FLOAT f(Vector x)
        {
            return (FLOAT)(x.transpose() * G * x) * 0.5 + (FLOAT)(b.transpose() * x) + c;
        }

        Vector xk;
        Vector dk;

        Matrix G{2, 2};
        Vector b{2};
        float c = 10;

      protected:
        bool Criterion(FLOAT a)
        {
            // Armijo
            /*FLOAT p = (10e-3)/2.0;
            FLOAT left = f(xk + a * dk);
            FLOAT right = f(xk) + p * (FLOAT)(  (G * xk+b).transpose() * dk  ) * a;
            return left <= right;*/

            // Goldstein
            FLOAT p = 0.3;
            FLOAT left = f(xk + a * dk);
            FLOAT right = f(xk) + p * (FLOAT)((G * xk).transpose() * dk) * a;
            bool b1 = left <= right;

            FLOAT right2 = f(xk) + (1 - p) * (FLOAT)((G * xk).transpose() * dk) * a;
            bool b2 = left >= right2;

            return b1 && b2;
        }
    };

    class NewtonMethod : public NewtonBase
    {
      public:
        Vector gk(Vector xk) override
        {
            return _G * xk + _b;
        }
        Matrix Gk(Vector xk) override
        {

            return _G;
        }
        Vector Solve(Vector x0, FLOAT _gk_norm) override
        {
            int k = 0;
            Vector xk = x0;

            double xk_max_norm;

            while ((xk_max_norm = gk(xk).cwiseAbs().maxCoeff()) >= _gk_norm)
            {

                // compute dk
                Matrix _GK = Gk(xk);
                Vector _n_gk = -gk(xk);

#ifdef _DEBUG
                std::cout << "||gk||2: " << _n_gk.norm() << std::endl;
#endif

                // solve _GK*dk = _n_gk
                Vector dk;
                dk = _GK.colPivHouseholderQr().solve(_n_gk);

                // 阻尼牛顿法 添加线搜索
                static FLOAT a = 10;
                ls.xk = xk;
                ls.dk = dk;
                a = ls.QuadraticPolynomialInterpolation(a);
                xk = xk + a * dk;
                k++;
            }

            return xk;
        }

      public:
        LS ls;
        Matrix _G{2, 2};
        Vector _b{2};
    };

    NewtonMethod newton;
    newton._G = G;
    newton._b = b;
    newton.ls.G = G;
    newton.ls.b = b;
    newton.ls.c = c;

    std::cout << newton.Solve(x0, 10e-3).transpose() << std::endl;
}
