#include "newton.h"
#include <iostream>

Vector NewtonBase::gk(Vector xk)
{
	return {};
}

Matrix NewtonBase::Gk(Vector xk)
{
	return {};
}

Vector NewtonBase::Solve(Vector x0, FLOAT _gk_norm)
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
		std::cout <<"||gk||2: " << _n_gk.norm() << std::endl;
#endif

		// solve _GK*dk = _n_gk
		Vector dk;
		dk = _GK.colPivHouseholderQr().solve(_n_gk);

		// iteration
		xk = xk + dk;
		k++;

	}


	return xk;
}
