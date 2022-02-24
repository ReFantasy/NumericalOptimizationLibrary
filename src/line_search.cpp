#include "line_search.h"
#include <iostream>

TFLOAT LineSearch::Zerosixeight(TFLOAT a0, TFLOAT r0, TFLOAT epsilon, TFLOAT t)
{
    TFLOAT secton_a, secton_b;
    AdvanceandRetreat(a0, r0, t, secton_a, secton_b);
    return GoldenSection(secton_a, secton_b, epsilon);
}

TFLOAT LineSearch::QuadraticPolynomialInterpolation(TFLOAT a0)
{
    // 满足准则直接返回
    if (Criterion(a0))
        return a0;

    // 不满足准则
    TFLOAT a1 = -dPhi_da(0) * a0 * a0 / (Phi(a0) - Phi(0) - dPhi_da(0) * a0) / 2;

    return QuadraticPolynomialInterpolation(a1);
}

TFLOAT LineSearch::Phi(TFLOAT a)
{
    return (*_functor)(xk + a * dk);
}

TFLOAT LineSearch::dPhi_da(TFLOAT a)
{
    return _functor->FirstOrderDerivatives(xk + a * dk).transpose() * dk;
}

bool LineSearch::Criterion(TFLOAT a)
{
    bool res = false;
    TFLOAT p;
    TFLOAT left;
    TFLOAT right;
    TFLOAT right2;

    switch (_criterion_type)
    {
    case CriterionType::Armijo:
        p = (10e-3) / 2.0;
        left = (*_functor)(xk + a * dk);
        right = (*_functor)(xk) + p * (TFLOAT)(_functor->FirstOrderDerivatives(xk).transpose() * dk) * a;
        res = (left <= right);
        break;

    case CriterionType::Goldstein:
        p = 0.3;
        left = (*_functor)(xk + a * dk);
        right = (*_functor)(xk) + p * (TFLOAT)(_functor->FirstOrderDerivatives(xk).transpose() * dk) * a;
        right2 = (*_functor)(xk) + (1 - p) * (TFLOAT)(_functor->FirstOrderDerivatives(xk).transpose() * dk) * a;
        res = (left <= right) && (left >= right2);
        break;

    default:
        break;
    }

    return res;
}

void LineSearch::AdvanceandRetreat(TFLOAT a0, TFLOAT r0, TFLOAT t, TFLOAT &secton_a, TFLOAT &secton_b)
{
    assert(a0 >= 0);
    assert(r0 > 0);
    assert(t > 1);

    int i = 0;
    TFLOAT a, ai, ai1, ri;
    a = ai = a0;
    ri = r0;

    TFLOAT _a, _b;
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

TFLOAT LineSearch::GoldenSection(TFLOAT secton_a, TFLOAT secton_b, TFLOAT epsilon /*= 10e-3*/)
{
    TFLOAT a = secton_a;
    TFLOAT b = secton_b;
    static const TFLOAT r = (std::sqrt(5) - 1) / 2; // 0.618

    while ((b - a) > epsilon)
    {
        TFLOAT al = a + (1.0 - r) * (b - a);
        TFLOAT ar = a + r * (b - a);

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
