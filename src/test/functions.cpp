//
// Created by ReFantasy on 2022/3/13.
//

#include "functions.h"

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

FLOAT THREE_HUMP_CAMEL::operator()(const Vector &x) const
{
    return 2 * x(0) * x(0) - 1.05 * std::pow(x(0), 4) + std::pow(x(0), 6) / 6.0 + x(0) * x(1) + x(1) * x(1);
}

Vector THREE_HUMP_CAMEL::Gradient(const Vector &x) const
{
    Vector v(2);
    v(0) = std::pow(x(0), 5) - 4.2 * std::pow(x(0), 3) + 4 * x(0) + x(1);
    v(1) = x(0) + 2 * x(1);
    return v;
}

Matrix THREE_HUMP_CAMEL::Hesse(const Vector &x) const
{
    Matrix m(2, 2);
    m(1, 0) = m(0, 1) = 1;
    m(1, 1) = 2;
    m(0, 0) = 5 * std::pow(x(0), 4) - 12.6 * x(0) * x(0) + 4;
    return m;
}