#include "line_search.h"
#include <iostream>

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
    FLOAT a1 = -dPhi_da(0) * a0 * a0 / (Phi(a0) - Phi(0) - dPhi_da(0) * a0) / 2;

    return QuadraticPolynomialInterpolation(a1);
}

FLOAT LineSearch::Phi(FLOAT a)
{
    return (*_functor)(xk + a * dk);
}

FLOAT LineSearch::dPhi_da(FLOAT a)
{
    return _functor->FirstOrderDerivatives(xk + a * dk).transpose() * dk;
}

bool LineSearch::Criterion(FLOAT a)
{
    bool res = false;
    FLOAT p;
    FLOAT left;
    FLOAT right;
    FLOAT right2;

    switch (_criterion_type)
    {
    case CriterionType::Armijo:
        p = (10e-3) / 2.0;
        left = (*_functor)(xk + a * dk);
        right = (*_functor)(xk) + p * (FLOAT)(_functor->FirstOrderDerivatives(xk).transpose() * dk) * a;
        res = (left <= right);
        break;

    case CriterionType::Goldstein:
        p = 0.3;
        left = (*_functor)(xk + a * dk);
        right = (*_functor)(xk) + p * (FLOAT)(_functor->FirstOrderDerivatives(xk).transpose() * dk) * a;
        right2 = (*_functor)(xk) + (1 - p) * (FLOAT)(_functor->FirstOrderDerivatives(xk).transpose() * dk) * a;
        res = (left <= right) && (left >= right2);
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
        else if (Phi(ai1) < Phi(ai))
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

        if (Phi(al) < Phi(ar))
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
