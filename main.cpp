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

        // 求解线性方程组
        Eigen::Matrix3f A;
//        Eigen::Vector3f b;
        A<<8,-3,2,4,11,-1,6,3,12;
//        b <<20,33,36;
//		A.row(0) = b;
        // Eigen::Vector3f  x = SOR(A, b);
        //std::cout<<A.row(0)<<std::endl;

// QR分解
		Eigen::HouseholderQR<Eigen::MatrixXf> qr;
		qr.compute(A);
		Eigen::MatrixXf R = qr.matrixQR().triangularView<Eigen::Upper>();
		Eigen::MatrixXf Q = qr.householderQ();
// 显示
		std::cout << "A = " << std::endl << A << std::endl << std::endl;
		std::cout << "Q = " << std::endl << Q << std::endl << std::endl;
		std::cout << "R = " << std::endl << R << std::endl << std::endl;


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