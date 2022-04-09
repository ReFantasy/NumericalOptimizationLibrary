#pragma once

#include "global.h"

namespace NOL
{
	class ConjugateGradient : public UnconstrainedOptimizationLineSearchBase
	{
	public:
		using UnconstrainedOptimizationLineSearchBase::UnconstrainedOptimizationLineSearchBase;

		Vector Solve() override;

		Vector SearchDirection(const Vector& xk) const override;


		FLOAT Step(const Vector& xk, const Vector& dk) const override;

	private:
		Vector _gk;
		mutable Vector _dk;
	};
}