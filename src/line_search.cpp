#include "line_search.h"
#include <iostream>

namespace NOL
{
FLOAT LineSearch::Zerosixeight(FLOAT a0, FLOAT r0, FLOAT epsilon, FLOAT t)
{
    FLOAT secton_a, secton_b;
    AdvanceandRetreat(a0, r0, t, secton_a, secton_b);
    return GoldenSection(secton_a, secton_b, epsilon);
}

FLOAT LineSearch::QuadraticPolynomialInterpolation(FLOAT a0)
{
    // 满足准则直接返回
    if (Criterion(a0))
        return a0;

    // 不满足准则
    FLOAT a1 = -dphi_da(0) * a0 * a0 / (phi(a0) - phi(0) - dphi_da(0) * a0) / 2;

    return QuadraticPolynomialInterpolation(a1);
}

FLOAT LineSearch::CubicPolynomialInterpolation(FLOAT a0)
{
	if (Criterion(a0))
		return a0;

	// 不满足准则
	FLOAT a1 = -dphi_da(0) * a0 * a0 / (phi(a0) - phi(0) - dphi_da(0) * a0) / 2;

    if (Criterion(a1))
        return a1;

	return CubicPolynomial(a0, a1);
}

FLOAT LineSearch::phi(FLOAT a)
{
    return (*_functor)(xk + a * dk);
}

FLOAT LineSearch::dphi_da(FLOAT a)
{
    return _functor->FirstOrderDerivatives(xk + a * dk).transpose() * dk;
}

FLOAT LineSearch::CubicPolynomial(FLOAT a0, FLOAT a1)
{
    FLOAT c = dphi_da(0);
    FLOAT d = phi(0);

    Matrix m(2, 2);
    m << a0 * a0, -a1 * a1, -a0 * a0 * a0, a1* a1* a1;
    Vector v(2);
    v << (phi(a1) - phi(0) - dphi_da(0) * a1), (phi(a0) - phi(0) - dphi_da(0) * a0);
    std::cout << m << std::endl;
    std::cout << v << std::endl;
    Vector res = (m * v) / (a0 * a0 * a1 * a1 * (a1 - a0));
    FLOAT a = res(0);
    FLOAT b = res(1);

    FLOAT zero1;
    FLOAT zero2;
    if (std::abs(a) < 10e-10)//a==0
    {
        zero1 = zero2 = -c / 2.0 / b;
    }
    else
    {
        FLOAT zero1 = (-b - std::sqrt(b * b - 3 * a * c)) / 3.0 / a;
        FLOAT zero2 = (-b + std::sqrt(b * b - 3 * a * c)) / 3.0 / a;
    }

    

    FLOAT z = std::min(zero1, zero2);
    if (z < 0)
        z = std::max(zero1, zero2);

	if (Criterion(z))
		return z;
	else
		return CubicPolynomial(a1, z);
}

bool LineSearch::Criterion(FLOAT a)
{
    bool res = false;

    FLOAT left;
    FLOAT right;
    FLOAT right2;

    switch (_criterion_type)
    {
    case CriterionType::Armijo:
        if ((_armijo_p <= 0) || (_armijo_p >= 1))
        {
            throw std::out_of_range("_armijo_p is out of range (0,1)");
        }
        left = (*_functor)(xk + a * dk);
        right = (*_functor)(xk) + _armijo_p * (FLOAT)(_functor->FirstOrderDerivatives(xk).transpose() * dk) * a;
        res = (left <= right);
        break;

    case CriterionType::Goldstein:
        if ((_goldstein_p <= 0) || (_armijo_p >= 0.5))
        {
            throw std::out_of_range("_goldstein_p is out of range (0,1)");
        }
        left = (*_functor)(xk + a * dk);
        right = (*_functor)(xk) + _goldstein_p * (FLOAT)(_functor->FirstOrderDerivatives(xk).transpose() * dk) * a;
        right2 = (*_functor)(xk) + (1 - _goldstein_p) * (FLOAT)(_functor->FirstOrderDerivatives(xk).transpose() * dk) * a;
        res = (left <= right) && (left >= right2);
        break;

    case CriterionType::Wolfe:
        if ((_wolfe_p > 0) && (_wofe_sigma > _wolfe_p) && (_wofe_sigma < 1))
        {
            return (((*_functor)(xk + a * dk)) <= ((*_functor)(xk) + _wolfe_p * (FLOAT)(_functor->FirstOrderDerivatives(xk).transpose() * dk) * a))
                &&
                ((_functor->FirstOrderDerivatives(xk + a * dk).transpose() * dk) >= (_wofe_sigma * _functor->FirstOrderDerivatives(xk).transpose() * dk));
        }
        else
        {
            throw std::out_of_range("wofe criterion is fail");
        }
        break;

    case CriterionType::StrongWolfe:
        if ((_wolfe_p > 0) && (_wofe_sigma > _wolfe_p) && (_wofe_sigma < 1))
        {
            return (((*_functor)(xk + a * dk)) <= ((*_functor)(xk) + _wolfe_p * (FLOAT)(_functor->FirstOrderDerivatives(xk).transpose() * dk) * a))
                &&
                (std::abs((_functor->FirstOrderDerivatives(xk + a * dk).transpose() * dk)) <= (-(FLOAT)(_wofe_sigma * _functor->FirstOrderDerivatives(xk).transpose() * dk)));
        }
        else
        {
            throw std::out_of_range("wofe criterion is fail");
        }
        break;
    default:
        break;
    }

    return res;
}

void LineSearch::AdvanceandRetreat(FLOAT a0, FLOAT r0, FLOAT t, FLOAT &secton_a, FLOAT &secton_b)
{
    assert(a0 >= 0);
    assert(r0 > 0);
    assert(t > 1);

    int i = 0;
    FLOAT a, ai, ai1, ri;
    a = ai = a0;
    ri = r0;

    FLOAT _a, _b;
    while (true)
    {
        ai1 = ai + ri;
        if (ai1 <= 0)
        {
            ai1 = 0;
        }
        else if (phi(ai1) < phi(ai))
        {
            // step 3
            ri = t * ri;
            a = ai;
            ai = ai1;
            i++;
            continue;
        }

        // step 4
        if (i == 0)
        {
            ri = -ri;
            a = ai1;
            i++;
        }
        else
        {
            _a = std::min(a, ai1);
            _b = std::max(a, ai1);
            break;
        }
    }

    secton_a = _a;
    secton_b = _b;
}

FLOAT LineSearch::GoldenSection(FLOAT secton_a, FLOAT secton_b, FLOAT epsilon /*= 10e-3*/)
{
    FLOAT a = secton_a;
    FLOAT b = secton_b;
    static const FLOAT r = (std::sqrt(5) - 1) / 2; // 0.618

    while ((b - a) > epsilon)
    {
        FLOAT al = a + (1.0 - r) * (b - a);
        FLOAT ar = a + r * (b - a);

        if (phi(al) < phi(ar))
        {
            b = ar;
        }
        else
        {
            a = al;
        }
    }

    return (a + b) / 2.0;
}
} // namespace NOL
