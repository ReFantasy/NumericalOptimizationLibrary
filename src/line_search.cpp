#include "line_search.h"
#include <iostream>

namespace NOL
{
FLOAT LineSearch::GoldenMethod(FLOAT a0, const Options &options)
{
    FLOAT secton_a, secton_b;
    AdvanceAndRetreat(a0, options.parameter_line_search_advance_and_retreat_h,
                      options.parameter_line_search_advance_and_retreat_t, secton_a, secton_b);
    return GoldenSection(secton_a, secton_b, options.parameter_line_search_golden_section_size);
}

FLOAT LineSearch::QuadraticInterpolation(FLOAT alpha, const Options &options)
{
    FLOAT secton_a, secton_b;
    AdvanceAndRetreat(alpha, options.parameter_line_search_advance_and_retreat_h,
                      options.parameter_line_search_advance_and_retreat_t, secton_a, secton_b);
    return QuadraticInterpolationMinimum(std::max(0.0, secton_a), (secton_b < 0) ? 1.0 : secton_b);
}

FLOAT LineSearch::Armijo(FLOAT alpha, const Options &options)
{
    FLOAT phi0 = phi(0.0);
    FLOAT dphi0 = dphi_da(0.0);

    while (true)
    {
        FLOAT phi_alpha = phi(alpha);

        if (phi_alpha > (phi0 + options.parameter_line_search_armijo_rho * dphi0 * alpha))
        {
            alpha = alpha / options.parameter_line_search_armijo_t;
        }
        break;
    }
    return alpha;
}

FLOAT LineSearch::Goldstein(FLOAT alpha, const Options &options)
{
    FLOAT a = 0.0;
    FLOAT b = std::numeric_limits<FLOAT>::max();
    FLOAT phi0 = phi(0.0);
    FLOAT dphi0 = dphi_da(0.0);
    int k = 0;

    while (true)
    {
        FLOAT phi_alpha = phi(alpha);
        if (phi_alpha > (phi0 + options.parameter_line_search_goldstein_p * dphi0 * alpha))
        {
            b = alpha;
            alpha = (a + b) / 2.0;
            k++;
            continue;
        }
        if (phi_alpha < (phi0 + (1 - options.parameter_line_search_goldstein_p) * dphi0 * alpha))
        {
            a = alpha;
            if (b == std::numeric_limits<FLOAT>::max())
                alpha = a * 2;
            else
                alpha = (a + b) / 2.0;

            k++;
            continue;
        }
        break;
    }

    return alpha;
}

FLOAT LineSearch::StrongWolfe(FLOAT alpha, const Options &options)
{
    auto choose = [](FLOAT a, FLOAT b) { return (a + b) / 2.0; };

    FLOAT a0 = 0;
    FLOAT amax = options.parameter_line_search_strong_wolfe_alpha_max;
    FLOAT a1 = alpha;
    if ((a1 < 0) || a1 >= amax)
        a1 = choose(0, amax);

    FLOAT c1 = options.parameter_line_search_strong_wolfe_rho;
    FLOAT c2 = options.parameter_line_search_strong_wolfe_sigma;

    int i = 1;
    while (true)
    {
        if ((phi(a1) > ( phi(0) + c1 * a1 * dphi_da(0) )) || ( (phi(a1) >= phi(a0)) && (i > 1)))
            return Zoom(a0, a1, options);
        if (std::abs(dphi_da(a1)) <= (-c2 * dphi_da(0)))
            return a1;
        if (dphi_da(a1) >= 0)
            return Zoom(a1, a0, options);

        a0 = a1;
        a1 = choose(a1, amax);
        i++;
    }
}

FLOAT LineSearch::phi(FLOAT a)
{
    return (*_functor)(xk + a * dk);
}

FLOAT LineSearch::dphi_da(FLOAT a)
{
    return _functor->FirstOrderDerivatives(xk + a * dk).transpose() * dk;
}

void LineSearch::AdvanceAndRetreat(FLOAT alpha0, FLOAT h0, FLOAT t, FLOAT &secton_a, FLOAT &secton_b)
{
    int i = 0;
    FLOAT a, ai, ai1, ri;
    a = ai = alpha0;
    ri = h0;

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

FLOAT LineSearch::QuadraticInterpolationMinimum(FLOAT a1, FLOAT a2)
{
    return a1 - (a1 - a2) / 2.0 / (1 - (phi(a1) - phi(a2)) / ((a1 - a2) * dphi_da(a1)));
}

FLOAT LineSearch::Zoom(FLOAT alpha_lo, FLOAT alpha_hi, const Options &options)
{
    FLOAT c1 = options.parameter_line_search_strong_wolfe_rho;
    FLOAT c2 = options.parameter_line_search_strong_wolfe_sigma;

    while (true)
    {
        //FLOAT alpha_j = QuadraticInterpolationMinimum(alpha_lo, alpha_hi);
        FLOAT alpha_j = (alpha_lo+ alpha_hi)/2.0;

        FLOAT phi_aj = phi(alpha_j);

        FLOAT tmp = phi(0) + c1 * alpha_j * dphi_da(0);
        /*if ((phi_aj > tmp) && (tmp > phi(alpha_lo)))
            alpha_hi = alpha_j;*/
        if ((phi_aj > tmp) || (phi_aj >= phi(alpha_lo)))
            alpha_hi = alpha_j;
        else
        {
            FLOAT dphi_aj = dphi_da(alpha_j);
            if (std::abs(dphi_aj) <= (-c2 * dphi_da(0)))
                return alpha_j;
            if (dphi_aj * (alpha_hi - alpha_lo) >= 0)
                alpha_hi = alpha_lo;
            alpha_lo = alpha_j;
        }
    }
}

} // namespace NOL
