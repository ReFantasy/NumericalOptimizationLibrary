#pragma once

#include "Eigen/Dense"

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;
using FLOAT = double;

class TargetFunctor
{
  public:
    virtual FLOAT operator()(const Vector &x) const
    {
        return FLOAT{};
    };

    virtual Vector FirstOrderDerivatives(const Vector &x) const
    {
        return Vector{0, 0};
    };
    virtual Matrix SecondOrderDerivatives(const Vector &x) const
    {
        return Matrix{};
    };
};
