#ifndef __NEWTON_H__
#define __NEWTON_H__
#include "global.h"
#include "Eigen/Dense"

class NewtonBase
{
public:
	virtual Vector gk(Vector xk);
	virtual Matrix Gk(Vector xk);
	
	virtual Vector Solve(Vector x0, FLOAT _gk_norm);

private:
};

#endif//__NEWTON_H__