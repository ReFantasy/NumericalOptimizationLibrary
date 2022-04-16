#pragma once
#include "Eigen/Dense"
#include "global.h"
#include "line_search.h"
#include <vector>
#include <memory>
namespace NOL
{
	class LeastSquareFunctorRx:public TargetFunctor
	{
	public:
		FLOAT operator()(const Vector &xk) const override
		{
			FLOAT res = 0;
			for(int i=0; i < rix.size(); i++)
			{
				res += std::pow((*rix[i])(xk),2) * 0.5;
			}
			return res;
		}
		Vector Gradient(const Vector &xk) const override
		{
			return J(xk).transpose() * r(xk);
		}

	public:
		Matrix J(const Vector& x)const override
		{
			Matrix jacobi(rix.size(), x.size());
			for (int i = 0; i < rix.size(); i++)
			{
				jacobi.row(i) = (*rix[i]).Gradient(x);
			}
			return jacobi;
		}

		Vector r(const Vector& x)const override
		{
			Vector res(rix.size());
			for (int i = 0; i < rix.size(); i++)
			{
				res(i) = (*rix[i])(x);
			}
			return res;
		}

		void push_back_rx(std::shared_ptr<TargetFunctor> f)
		{
			rix.push_back(f);
		}


	private:
		std::vector<std::shared_ptr<TargetFunctor>> rix;
	};

	class LSGaussNewton : public UnconstrainedOptimizationLineSearchBase
	{
	public:
		using UnconstrainedOptimizationLineSearchBase::UnconstrainedOptimizationLineSearchBase;

		Vector SearchDirection(const Vector& xk) const override
		{
			Matrix Jk = _functor_ptr->J(xk);
			Vector rk = _functor_ptr->r(xk);
			Matrix A = Jk.transpose()*Jk;
			Vector b = -Jk.transpose()*rk;
			return LinearEquationSolver::Solver(A,b);
		}

		FLOAT Step(const Vector& xk, const Vector& dk) const override
		{
			return 1.0;
		};
	};
}
