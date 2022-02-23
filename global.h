#pragma once

#include "Eigen/Dense"

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;
using FLOAT = double;


struct TargetFunctor 
{
    virtual FLOAT operator()(const Vector &x) const = 0;
    
    virtual Vector FirstOrderDerivatives(const Vector& x) const = 0;
    virtual Matrix SecondOrderDerivatives(const Vector& x) const = 0;

};

