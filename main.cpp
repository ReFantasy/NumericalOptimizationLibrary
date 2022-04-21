#include <iostream>
#include "example.h"
#include "least_square_method.h"
using namespace NOL;

class F1:public TargetFunctor
{
public:
	FLOAT operator()(const Vector &xk) const override
	{
		return xk(0)+1;
	}

	Vector Gradient(const Vector &xk) const override
	{
		Vector g(1);
		g(0) = 1;
		return g;
	}
};
class F2:public TargetFunctor
{
public:
	FLOAT operator()(const Vector &xk) const override
	{

		return lambda*xk(0)*xk(0)+xk(0)-1;
	}

	Vector Gradient(const Vector &xk) const override
	{
		Vector g(1);
		g(0) = 2*lambda*xk(0)+1;
		return g;
	}
	FLOAT lambda = 0.9;
};

int main(int argc, char *argv[])
{
    try
    {
        // example_3_1();
        // example_3_2();
        // example_cg();
        // example_3_3();
        // example();

//		Vector  x(1);
//		x(0) = 300.0;
//
//		auto functor = std::make_shared<LeastSquareFunctorRx>();
//		(*functor).push_back_rx(std::make_shared<F1>());
//		(*functor).push_back_rx(std::make_shared<F2>());
//
//		LSGaussNewton solver;
//		solver.SetFunctor(functor);
//
//		auto options = solver.GetOptions();
//
//
//		options->optimized_performance = false;
//		options->init_x = x;
//
//		LinearEquationSolver::solver_type = LinearEquationSolver::SOLVER_TYPE::GAUSS_SEIDE;
//
//		auto res = solver.Solve();
//
//		std::cout << std::fixed << "Optimal solution : " << res.transpose() << std::endl
//				  << "min f(x): " << (*functor)(res) << std::endl
//				  << "Iteration Number: " << solver.NumOfIteration() << std::endl
//				  << "time : " << time << "ms" << std::endl
//				  << std::endl;

    }
    catch (std::exception &e)
    {
        std::cout << std::endl << e.what() << std::endl;
    }
    catch (...)
    {
        std::cout << "error" << std::endl;
    }
}
