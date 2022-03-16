#ifndef __DECORATOR_H__
#define __DECORATOR_H__
#include "global.h"
#include <memory>
#include <utility>
namespace NOL
{
	class Decorator:public UnconstrainedOptimizationLineSearchBase
	{
	public:
		void SetSolver(std::shared_ptr<UnconstrainedOptimizationLineSearchBase> solver_ptr)
		{
			_solver_ptr = std::move(solver_ptr);
		}
		Vector Solve()override
		{
			return _solver_ptr->Solve();
		}

	private:
		std::shared_ptr<UnconstrainedOptimizationLineSearchBase> _solver_ptr = nullptr;
	};
}


#endif//__DECORATOR_H__