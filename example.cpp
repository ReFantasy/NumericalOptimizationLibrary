#include "example.h"
#include "global.h"
#include "line_search.h"
using namespace NOL;

void example_3_1()
{
    // 定义目标函数
    struct Functor : TargetFunctor
    {
        virtual FLOAT operator()(const Vector &x) const override
        {
            return x.transpose() * G * x + (FLOAT)(b.transpose() * x) + c;
        }
        virtual Vector Gradient(const Vector &x) const override
        {
            return (G + G.transpose()) * x / 2 + b;
        }
        virtual Matrix Hesse(const Vector &x) const override
        {
            return G;
        }

        Matrix G{2, 2};
        Vector b{2};
        float c = 10;
    };

    // 构造目标函数对象
    // Functor functor;
    std::shared_ptr<Functor> functor = std::make_shared<Functor>();
    (*functor).b << 2, 3;
    (*functor).G << 21, 4, 4, 15;
    // functor.G << 21, 4, 4, 1;

    // 构造求解器对象
    auto sd = OptimizationFactory::CreateSolver(OptimizationMethodType::SD, functor);

    // 获取求解器参数指针
    auto options = sd->GetOptions();

    // 设置求解参数
    Vector x(2);
    x << RandomNumber<FLOAT>(-10, 10), RandomNumber<FLOAT>(-10, 10);
    options->init_x = x;
    options->termination_value = 10e-5;
    options->optimized_performance = true;

    // 求解
    Vector res = sd->Solve();

    // 输出结果
    std::cout << "example 3.1" << std::endl;
    std::cout << "initial x: " << x.transpose() << std::endl;
    std::cout << "Optimal solution : " << res.transpose() << "  time : " << sd->GetTimer()->Elapse() << "ms"
              << std::endl
              << std::endl;
}

void example_3_2()
{
    struct Functor : TargetFunctor
    {
        virtual FLOAT operator()(const Vector &x) const override
        {
            return 3 * x(0) * x(0) + 3 * x(1) * x(1) - x(0) * x(0) * x(1);
        }
        virtual Vector Gradient(const Vector &x) const override
        {
            Vector dx(2);
            dx(0) = 6 * x(0) - 2 * x(0) * x(1);
            dx(1) = 6 * x(1) - x(0) * x(0);
            return dx;
        }
        virtual Matrix Hesse(const Vector &x) const override
        {
            Matrix J(2, 2);
            J(0, 0) = 6 - 2 * x(1);
            J(0, 1) = -2 * x(0);
            J(1, 0) = -2 * x(0);
            J(1, 1) = 6;
            return J;
        }
    };

    auto functor = std::make_shared<Functor>();

    auto solver = OptimizationFactory::CreateSolver(OptimizationMethodType::NEWTON, functor);

    auto option = solver->GetOptions();

    option->termination_value = 10e-6;
    option->optimized_performance = true;
    Vector x0(2);
    x0(0) = x0(1) = 1.5;
    // x0(0) = -2;
    // x0(1) = 4;
    // x0(0) = 0;
    // x0(1) = 3;
    option->init_x = x0;

    std::cout << "example 3.2" << std::endl;
    std::cout << "initial x: " << x0.transpose() << std::endl;
    std::cout << "Optimal solution : " << solver->Solve().transpose() << "  time : " << solver->GetTimer()->Elapse()
              << "ms" << std::endl
              << std::endl;
}
void example_cg()
{
	struct Functor : TargetFunctor
	{
		virtual FLOAT operator()(const Vector &x) const override
		{
			int n = x.size();
			FLOAT sum = 0;
			for (int i = 1; i < n; i += 2)
			{
				sum = sum + 100 * std::pow(x(i - 1) * x(i - 1) - x(i), 2) + std::pow(x(i - 1) - 1, 2);
			}
			return sum;

			// return 100 * std::pow((x(1) - x(0)* x(0)),2)  + std::pow((1 - x(0)),2);
		}
		virtual Vector Gradient(const Vector &x) const override
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
		virtual Matrix Hesse(const Vector &x) const override
		{
			return Matrix{};
		}
	};

	auto functor = std::make_shared<Functor>();

	auto solver = OptimizationFactory::CreateSolver(OptimizationMethodType::ConjugateGradient, functor);

	auto options = solver->GetOptions();

	options->line_search_type = LineSearchType::STRONGWOLFE;
	options->conjugate_gradient_type = ConjugateGradientType::FR;
	options->conjugate_gradient_restart_num = 10;
	options->optimized_performance = false;

	int n = 8;
	Vector x(n);
	for (int i = 0; i < x.size(); i++)
	{
		if (i % 2 == 0)
		{
			// x(i) = RandomNumber<FLOAT>(-100.0, 100.0);
			x(i) = -1.2;
		}
		else
		{
			// x(i) = RandomNumber<FLOAT>(-100.0, 100.0);
			x(i) = 1.0;
		}
		// std::cout << x(i) << ", " ;
	}
	options->init_x = x;

	std::cout << "TEST ROSENBROCK FUNCTION By Conjugate Gradient Method" << std::endl;
	std::cout << std::fixed << "dim: " << x.size() << "  initial x: " << x.transpose() << std::endl;
	auto res = solver->Solve();
	auto time = solver->GetTimer()->Elapse();
	std::cout << std::fixed << "Optimal solution : " << res.transpose() << std::endl
			  << "min f(x): " << (*functor)(res) << std::endl
			  << "Iteration Number: " << solver->NumOfIteration() << std::endl
			  << "time : " << time << "ms" << std::endl
			  << std::endl;



}
void example_3_3()
{
    struct Functor : TargetFunctor
    {
        virtual FLOAT operator()(const Vector &x) const override
        {
            return 3 * x(0) * x(0) + 3 * x(1) * x(1) - x(0) * x(0) * x(1);
        }
        virtual Vector Gradient(const Vector &x) const override
        {
            Vector dx(2);
            dx(0) = 6 * x(0) - 2 * x(0) * x(1);
            dx(1) = 6 * x(1) - x(0) * x(0);
            return dx;
        }
        virtual Matrix Hesse(const Vector &x) const override
        {
            Matrix J(2, 2);
            J(0, 0) = 6 - 2 * x(1);
            J(0, 1) = -2 * x(0);
            J(1, 0) = -2 * x(0);
            J(1, 1) = 6;
            return J;
        }
    };

    auto functor = std::make_shared<Functor>();

    auto solver = OptimizationFactory::CreateSolver(OptimizationMethodType::LM, functor);

    auto option = solver->GetOptions();
    option->line_search_type = LineSearchType::WOLFE;

    option->termination_value = 10e-6;
    option->optimized_performance = true;
    Vector x0(2);
    // x0(0) = x0(1) = 1.5;
    // x0(0) = -2;
    // x0(1) = 4;
    x0(0) = 2;
    x0(1) = 3;
    option->init_x = x0;

    //	std::cout.setf(std::ios::right); // 设置对齐方式
    //	std::cout.width(8); //设置输出宽度
    //	std::cout.fill('0'); //将多余的空格用0填充
    //	std::cout.precision(4); //设置输出精度，保留有效数字

    std::cout << "example 3.3" << std::endl;
    std::cout << std::fixed << "initial x: " << x0.transpose() << std::endl;
    std::cout << std::fixed << "Optimal solution : " << solver->Solve().transpose()
              << "  time : " << solver->GetTimer()->Elapse() << "ms" << std::endl
              << std::endl;
}

void example()
{
    struct Functor : TargetFunctor
    {
        virtual FLOAT operator()(const Vector &x) const override
        {
            int n = x.size();
            FLOAT sum = 0;
            for (int i = 1; i < n; i += 2)
            {
                sum = sum + 100 * std::pow(x(i - 1) * x(i - 1) - x(i), 2) + std::pow(x(i - 1) - 1, 2);
            }
            return sum;

            // return 100 * std::pow((x(1) - x(0)* x(0)),2)  + std::pow((1 - x(0)),2);
        }
        virtual Vector Gradient(const Vector &x) const override
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
        virtual Matrix Hesse(const Vector &x) const override
        {
            return Matrix{};
        }
    };

    auto functor = std::make_shared<Functor>();

    auto solver = OptimizationFactory::CreateSolver(OptimizationMethodType::QUASI_NEWTON, functor);

    auto option = solver->GetOptions();
    option->line_search_type = LineSearchType::STRONGWOLFE;
    option->quasi_newton_type = QuasiNewtonSearchType::BFGS;
    option->termination_type = TerminationCriterionType::DELTA_F;
    option->termination_value = 10e-8;
    option->max_solver_time_in_seconds = 10;
    option->parameter_line_search_wolfe_rho = 10e-4;
    option->parameter_line_search_wolfe_sigma = 0.9;
    option->optimized_performance = true;

    int n = 40;
    Vector x(n);
    for (int i = 0; i < x.size(); i++)
    {
        if (i % 2 == 0)
        {
            // x(i) = RandomNumber<FLOAT>(-100.0, 100.0);
            x(i) = -1.2;
        }
        else
        {
            // x(i) = RandomNumber<FLOAT>(-100.0, 100.0);
            x(i) = 1.0;
        }
        // std::cout << x(i) << ", " ;
    }
    option->init_x = x;

    std::cout << "example 3.3" << std::endl;
    std::cout << std::fixed << "dim: " << x.size() << "  initial x: " << x.transpose() << std::endl;
    auto res = solver->Solve();
    auto time = solver->GetTimer()->Elapse();
    std::cout << std::fixed << "Optimal solution : " << res.transpose() << std::endl
              << "min f(x): " << (*functor)(res) << std::endl
              << "Iteration Number: " << solver->NumOfIteration() << std::endl
              << "time : " << time << "ms" << std::endl
              << std::endl;
}
