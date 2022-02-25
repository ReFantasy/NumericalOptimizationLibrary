#include "newton.h"
#include "line_search.h"
#include <iostream>

namespace NOL
{
Vector NewtonBase::Solve(TargetFunctor &fucntor, Options &options)
{
    int k = 0;
    Vector xk = options.init_x;

    // <--------
    options << "Base Newton Method with initial x: ";
    options << xk.transpose() << "\n\n";
    // -------->

    double xk_max_norm = 0;
    while (true)
    {
        xk_max_norm = fucntor.FirstOrderDerivatives(xk).cwiseAbs().maxCoeff();
        if (xk_max_norm < options.gk_norm)
            break;

        // <--------
        options << "k:" << k << " "
                << "  xk:(" << xk.transpose() << ") "
                << "  ||gk||: " << xk_max_norm << "\n";
        // -------->

        // compute dk
        Matrix Gk = fucntor.SecondOrderDerivatives(xk);
        Vector gk = fucntor.FirstOrderDerivatives(xk);

        // solve Gk*dk = -gk
        Vector dk;
        dk = Gk.colPivHouseholderQr().solve(-gk);

        // advance
        xk = xk + dk;
        k++;
    }

    // <--------
    options << "k:" << k << " "
            << "  xk:(" << xk.transpose() << ") "
            << "  ||gk||: " << xk_max_norm << "\n";
    // -------->

    return xk;
}

Vector DampedNewton::Solve(TargetFunctor &fucntor, Options &options)
{
    int k = 0;
    Vector xk = options.init_x;

    LineSearch line_search{};
    line_search._functor = &fucntor;

    // <--------
    options << "Damped Newton Method with initial x: ";
    options << xk.transpose() << "\n\n";
    // -------->

    double xk_max_norm;
    while (true)
    {
        xk_max_norm = fucntor.FirstOrderDerivatives(xk).cwiseAbs().maxCoeff();
        if (xk_max_norm < options.gk_norm)
            break;

        // <--------
        options << "k:" << k << " "
                << "  xk:(" << xk.transpose() << ") "
                << "  ||gk||: " << xk_max_norm << "\n";
        // -------->

        // compute dk
        Matrix Gk = fucntor.SecondOrderDerivatives(xk);
        Vector gk = fucntor.FirstOrderDerivatives(xk);

        // solve Gk*dk = -gk
        Vector dk;
        dk = Gk.colPivHouseholderQr().solve(-gk);

        static FLOAT alpha = 10;
        line_search.xk = xk;
        line_search.dk = dk;
        alpha = line_search.QuadraticPolynomialInterpolation(alpha);
        // advance
        xk = xk + alpha * dk;
        k++;
    }

    // <--------
    options << "k:" << k << " "
            << "  xk:(" << xk.transpose() << ") "
            << "  ||gk||: " << xk_max_norm << "\n";
    // -------->

    return xk;
}
} // namespace NOL
