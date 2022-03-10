#include "line_search.h"
#include <iostream>

namespace NOL
{

LinearSearch::LinearSearch(TargetFunctor *functor) : _functor(functor)
{
}

FLOAT LinearSearch::Search(FLOAT alpha, const Options &options)
{
    FLOAT step_length = alpha;

    switch (options.line_search_type)
    {
    case LineSearchType::GOLDENSECTION:
        step_length = GoldenMethod(alpha, options);
        break;
    case LineSearchType::ARMIJO:
        step_length = Armijo(alpha, options);
        break;
    case LineSearchType::GOLDSTEIN:
        step_length = Goldstein(alpha, options);
        break;
    case LineSearchType::WOLFE:
        step_length = Wolfe(alpha, options);
        break;
    case LineSearchType::STRONGWOLFE:
        step_length = StrongWolfe(alpha, options);
        break;
    default:
        step_length = 1.0;
    }
    return step_length;
}

FLOAT LinearSearch::GoldenMethod(FLOAT alpha, const Options &options)
{
    FLOAT secton_a, secton_b;
    AdvanceAndRetreat(alpha, options.parameter_line_search_advance_and_retreat_h,
                      options.parameter_line_search_advance_and_retreat_t, secton_a, secton_b);
    return GoldenSection(secton_a, secton_b, MinStepSize<FLOAT>::value);
}

FLOAT LinearSearch::Armijo(FLOAT alpha, const Options &options)
{

    FLOAT rho = options.parameter_line_search_armijo_rho;

    // p58
    FLOAT phi0 = phi(0.0);
    FLOAT dphi0 = dphi_da(0.0);
    FLOAT last_alpha = 0;
    int k = 0;
    while (true)
    {
        FLOAT phi_alpha = phi(alpha);

        if (phi_alpha < (phi0 + rho * dphi0 * alpha))
        {
            break;
        }
        FLOAT tmp = alpha;
        alpha = Interpolation(last_alpha, alpha, k);
        last_alpha = tmp;
        k++;
    }
    return alpha;
}

FLOAT LinearSearch::CubicInterpolationMinimum(FLOAT last_alpha, FLOAT alpha)
{
    FLOAT c = dphi_da(0.0);
    Matrix m(2, 2);
    m << last_alpha * last_alpha, -alpha * alpha, -last_alpha * last_alpha * last_alpha, alpha * alpha * alpha;

    Vector v(2);
    v << (phi(alpha) - phi(0) - dphi_da(0) * alpha), (phi(last_alpha) - phi(0) - dphi_da(0) * last_alpha);
    Vector res = (m * v) / (last_alpha * last_alpha * alpha * alpha * (alpha - last_alpha));
    FLOAT a = res(0);
    FLOAT b = res(1);
    last_alpha = alpha;
    alpha = (-b + std::sqrt(b * b - 3 * a * c)) / 3.0 / a;
    return alpha;
}

FLOAT LinearSearch::Interpolation(FLOAT last_alpha, FLOAT alpha, int k /*= 0*/)
{
    if (k == 0)
    {

        return QuadraticInterpolationMinimum(0, alpha);
    }
    else
    {
        return CubicInterpolationMinimum(last_alpha, alpha);
    }
}

FLOAT LinearSearch::Goldstein(FLOAT alpha, const Options &options)
{
    FLOAT a = 0.0;
    FLOAT b = std::numeric_limits<FLOAT>::max();
    FLOAT phi0 = phi(0.0);
    FLOAT dphi0 = dphi_da(0.0);
    int k = 0;

    while (true)
    {
        FLOAT phi_alpha = phi(alpha);
        if (phi_alpha > (phi0 + options.parameter_line_search_goldstein_rho * dphi0 * alpha))
        {
            b = alpha;
            alpha = (a + b) / 2.0;
            k++;
            continue;
        }
        if (phi_alpha < (phi0 + (1 - options.parameter_line_search_goldstein_rho) * dphi0 * alpha))
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

FLOAT LinearSearch::Wolfe(FLOAT alpha, const Options &options)
{
    FLOAT a = 0.0;
    FLOAT b = std::numeric_limits<FLOAT>::max();
    FLOAT phi0 = phi(0.0);
    FLOAT dphi0 = dphi_da(0.0);
    int k = 0;

    while (true)
    {
        if (phi(alpha) > (phi0 + options.parameter_line_search_wolfe_rho * dphi0 * alpha))
        {
            b = alpha;
            alpha = (a + b) / 2.0;
            k++;
            continue;
        }
        if (dphi_da(alpha) < options.parameter_line_search_wolfe_sigma * dphi0)
        {
            a = alpha;
            if (b == std::numeric_limits<FLOAT>::max())
                alpha = a * 2;
            else
                alpha = std::min(2.0 * alpha, (a + b) / 2.0);

            k++;
            continue;
        }

        break;
    }

    return alpha;
}

FLOAT LinearSearch::phi(FLOAT a)
{
    return (*_functor)(_xk + a * _dk);
}

FLOAT LinearSearch::dphi_da(FLOAT a)
{
    return _functor->Gradient(_xk + a * _dk).transpose() * _dk;
}

void LinearSearch::AdvanceAndRetreat(FLOAT alpha0, FLOAT h0, FLOAT t, FLOAT &secton_a, FLOAT &secton_b)
{
    int i = 0;
    FLOAT alpha, alpha_i, alpha_i_1, hi;
    alpha = alpha_i = alpha0;
    hi = h0;

    FLOAT _a, _b;
    while (true)
    {
        alpha_i_1 = alpha_i + hi;
        if (alpha_i_1 <= 0)
        {
            alpha_i_1 = 0;
        }
        else if (phi(alpha_i_1) < phi(alpha_i))
        {
            // step 3
            hi = t * hi;
            alpha = alpha_i;
            alpha_i = alpha_i_1;
            i++;
            continue;
        }

        // step 4
        if (i == 0)
        {
            hi = -hi;
            alpha = alpha_i_1;
            i++;
        }
        else
        {
            _a = std::min(alpha, alpha_i_1);
            _b = std::max(alpha, alpha_i_1);
            break;
        }
    }

    secton_a = _a;
    secton_b = _b;
}

/**
 * @brief convergence rate: k > log_r{epsilon/(b-a)}
 */
FLOAT LinearSearch::GoldenSection(FLOAT secton_a, FLOAT secton_b, FLOAT epsilon)
{
    FLOAT a = secton_a;
    FLOAT b = secton_b;
    static const FLOAT r = (std::sqrt(5) - 1) / 2; // 0.618

    // TODO Optimization
    while ((b - a) > epsilon)
    {
        FLOAT alpha_left = a + (1.0 - r) * (b - a);
        FLOAT alpha_right = a + r * (b - a);

        if (phi(alpha_left) < phi(alpha_right))
            b = alpha_right;
        else
            a = alpha_left;
    }

    return (a + b) / 2.0;
}

FLOAT LinearSearch::QuadraticInterpolationMinimum(FLOAT a1, FLOAT a2)
{
    return a1 - (a1 - a2) / 2.0 / (1 - (phi(a1) - phi(a2)) / ((a1 - a2) * dphi_da(a1)));
}

FLOAT LinearSearch::StrongWolfe(FLOAT alpha, const Options &options)
{
    FLOAT c1 = options.parameter_line_search_wolfe_rho;
    FLOAT c2 = options.parameter_line_search_wolfe_sigma;
    auto choose = [](FLOAT a, FLOAT b) {
        if (std::max(a, b) != std::numeric_limits<FLOAT>::max())
        {
            return (a + b) / 2.0;
        }
        else
        {
            return std::min(a, b) * 2.0;
        }
    };

    FLOAT a0 = 0;
    FLOAT amax = options.parameter_line_search_wolfe_alpha_max;
    FLOAT a1 = alpha;
    if ((a1 < 0) || a1 >= amax)
        a1 = choose(0, amax);

    int i = 1;
    while (true)
    {
        if ((phi(a1) > (phi(0) + c1 * a1 * dphi_da(0))) || ((phi(a1) >= phi(a0)) && (i > 1)))
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

FLOAT LinearSearch::Zoom(FLOAT alo, FLOAT ahi, const Options &options)
{
    FLOAT c1 = options.parameter_line_search_wolfe_rho;
    FLOAT c2 = options.parameter_line_search_wolfe_sigma;

    while (true)
    {
        if (std::abs(ahi - alo) < options.min_step_size)
        {
            // assert(phi(ahi) < phi(0));
            return std::max(ahi, alo);
        }
        FLOAT aj = (alo + ahi) / 2.0;
        FLOAT tmp = phi(0) + c1 * aj * dphi_da(0);
        if ((phi(aj) > tmp) || (phi(aj) >= phi(alo)))
        {
            ahi = aj;
        }
        else
        {
            FLOAT dphi_aj = dphi_da(aj);
            if (std::abs(dphi_aj) <= (-c2 * dphi_da(0)))
                return aj;
            if ((dphi_aj * (ahi - alo)) >= 0)
                ahi = alo;
            alo = aj;
        }
    }
}

} // namespace NOL
