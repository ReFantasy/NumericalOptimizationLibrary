#include "example.h"

void example_3_1()
{
	Timer timer;

	Vector b(2);
	b << 2, 3;
	float c = 10;
	Matrix G(2, 2);

	Vector x0(2);
	x0 << -30, 100;

	G << 21, 4, 4, 15;
	timer.ReSet();
	SteepestDescent(G, b, c, x0);
	int t1 = timer.Elapse();
	//std::cout << t1 << "ms" << std::endl;

	printf("\n\n");

	G << 21, 4, 4, 1;
	timer.ReSet();
	SteepestDescent(G, b, c, x0);
	int t2 = timer.Elapse();
	//std::cout << t2 << "ms" << std::endl;
}

void example_3_2()
{
	class NewtonMethod :public NewtonBase
	{
	public:
		Vector gk(Vector xk)override
		{
			Vector _gk(2);
			_gk(0) = 6 * xk(0) - 2 * xk(0) * xk(1);
			_gk(1) = 6 * xk(1) - xk(0) * xk(0);
			return _gk;
		}
		Matrix Gk(Vector xk)override
		{
			Matrix _GK(2, 2);
			_GK << 6 - 2 * xk(1), -2 * xk(0), -2 * xk(0), 6;
			return _GK;
		}
	};

	NewtonMethod newton;
	Vector x0(2);
	x0(0) = x0(1) = 1.5;
	//x0(0) = -2;
	//x0(1) = 4;

	Vector x(2);
	x = newton.Solve(x0, 10e-6);
	std::cout <<"res:   " << x.transpose() << std::endl;
}
