//
// Created by ReFantasy on 2022/3/13.
//

#include "test_functions.h"

FLOAT ROTATED_HYPER_ELLIPSOID::operator()(const Vector &xk) const
{
    int d = xk.size();
    FLOAT sum = 0;
    for (int i = 1; i <= d; i++)
    {
        sum += (d + 1 - i) * xk(i - 1) * xk(i - 1);
    }
    return sum;
}

Vector ROTATED_HYPER_ELLIPSOID::Gradient(const Vector &xk) const
{
    int d = xk.size();
    Vector df(d);
    for (int i = 1; i <= d; i++)
    {
        df(i - 1) = 2.0 * (d + 1 - i) * xk(i - 1);
    }
    return df;
}

Matrix ROTATED_HYPER_ELLIPSOID::Hesse(const Vector &xk) const
{
    int d = xk.size();
    Matrix hesse = Matrix::Identity(d, d);
    for (int i = 1; i <= d; i++)
    {
        hesse(i - 1, i - 1) = 2 * (d + 1 - i);
    }
    return hesse;
}
