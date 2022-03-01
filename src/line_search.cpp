#include "line_search.h"
#include <iostream>

namespace NOL
{
    FLOAT LineSearch::ZeroSixOneEight(FLOAT a0, FLOAT h0, FLOAT epsilon, FLOAT t)
    {
        FLOAT secton_a, secton_b;
        AdvanceAndRetreat(a0, h0, t, secton_a, secton_b);
        return GoldenSection(secton_a, secton_b, epsilon);
    }

    FLOAT LineSearch::QuadraticInterpolation(FLOAT alpha, FLOAT h0, FLOAT t)
    {
        FLOAT secton_a, secton_b;
        AdvanceAndRetreat(alpha, h0, t, secton_a, secton_b);
        return  QuadraticInterpolationMinimum(std::max(0.0, secton_a), (secton_b < 0) ? 1.0 : secton_b);
    }

    FLOAT LineSearch::CubicPolynomialInterpolation(FLOAT alpha, FLOAT h0, FLOAT t)
    {
        // TODO
        return 1.0;
    }

    FLOAT LineSearch::Armijo(FLOAT alpha, const Options& options)
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

    FLOAT LineSearch::Goldstein(FLOAT alpha, const Options& options)
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

    FLOAT LineSearch::phi(FLOAT a)
    {
        return (*_functor)(xk + a * dk);
    }

    FLOAT LineSearch::dphi_da(FLOAT a)
    {
        return _functor->FirstOrderDerivatives(xk + a * dk).transpose() * dk;
    }



    // Reference https://blog.csdn.net/cclethe/article/details/77621440?utm_source=blogxgwz2
    void LineSearch::AdvanceAndRetreat(FLOAT alpha0, FLOAT h0, FLOAT t, FLOAT& secton_a, FLOAT& secton_b)
    {
        //assert(h0 > 0);
        assert(t > 1);

        FLOAT h1;
        FLOAT alpha;
        FLOAT alpha1;

        // step 1
        FLOAT f0 = phi(alpha0);
        FLOAT f1;
        int k = 0;

        while (true)
        {
            // step 2
            alpha1 = alpha0 + h0;
            f1 = phi(alpha1);
            if (f1 < f0)
            {
                //step 3
                h1 = t * h0;
                alpha = alpha0;
                alpha0 = alpha1;
                f0 = f1;
                k++;
                h0 = h1;
            }
            else
            {
                // step 4
                if (k == 0)
                {
                    h1 = -h0;
                    alpha = alpha1;
                    alpha1 = alpha0;
                    f1 = f0;
                    k = 1;
                    h0 = h1;
                }
                else
                {
                    secton_a = std::min(alpha, alpha1);
                    secton_b = std::max(alpha, alpha1);
                    break;
                }
            }
        }

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
        return  a1 - (a1 - a2) / 2.0 / (1 - (phi(a1) - phi(a2)) / ((a1 - a2) * dphi_da(a1)));
    }

} // namespace NOL
