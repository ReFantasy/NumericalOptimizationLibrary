#include "pqf_by_sd.h"
#include <iostream>

/**
 * @brief Positive definite quadratic function
 * @return f(x) = 1/2 * x * G * x + b * x + c;
*/
FLOAT pqf(Vector x, Matrix G, Vector b, FLOAT c)
{
	return (FLOAT)(0.5 * x.transpose() * G * x) + b.transpose() * x + c;
}

/**
 * @brief First derivative of positive definite quadratic function
 * @return f'(x)
*/
Vector gk(Vector x, Matrix G, Vector b)
{
	return (G + G.transpose()) * x / 2 + b;
}

Vector SD(Vector x0, FLOAT _gk_norm, Matrix G, Vector b, FLOAT c)
{
	int k = 0;
	Vector xk = x0;

	double gk_norm;
	while ((gk_norm = gk(xk, G, b).norm()) >= _gk_norm)
	{
#ifdef _DEBUG
		std::cout << k << ": (" << xk.transpose() << ")" << "    " << gk_norm << std::endl;
#endif

		/*page 19 (2.7)
		Eigen::Vector2d dk = -gk(xk, G, b);
		Eigen::Vector2d G_dot_xk = (G * xk);
		float alpha = (float)((-b.transpose() - (G * xk).transpose()) * dk) / ((G * dk).transpose() * dk);
		xk = xk + alpha * dk;*/

		// page 43 (3.3)
		Vector _gk = gk(xk, G, b);
		FLOAT alpha = (FLOAT)(_gk.transpose() * _gk) / (_gk.transpose() * G * _gk);
		xk = xk - alpha * _gk;

		k++;
	}
	 
#ifdef _DEBUG
	std::cout << k << ": (" << xk.transpose() << ")" << "    " << gk_norm << std::endl;
#endif

	return xk;
}


Vector SteepestDescent(Matrix G, Vector b, FLOAT c, Vector x0, FLOAT gk_norm)
{
	//Eigen::EigenSolver<Matrix> eigen_solver(G);
	//std::cout << "matrix values = \n" << eigen_solver.eigenvalues() << std::endl;//形式为二维向量(4,0)和(-1,0)。真实值为4,-1。

	return SD(x0, gk_norm, G, b, c);
}


