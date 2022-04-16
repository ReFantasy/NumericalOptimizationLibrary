#include "helper.h"

namespace NOL
{
/**
 * @brief Based iterative cycle to solve linear equations
 * @param B
 * @param f
 * @param epsilon
 * @return
 */
Eigen::VectorXf IterativeCycle(Eigen::MatrixXf B, Eigen::MatrixXf f, double epsilon, size_t max_iter_num)
{
    Eigen::VectorXf x = f;

    size_t n = 0;
    while (true)
    {
        if (n >= max_iter_num)
            return x;

        Eigen::VectorXf new_x = B * x + f;
        n++;

        if ((new_x - x).lpNorm<Eigen::Infinity>() < epsilon)
        {
            x = new_x;
            break;
        }

        x = new_x;
    }
    return x;
}

Eigen::VectorXf NOL::Jacobi(Eigen::MatrixXf A, Eigen::VectorXf b, double epsilon, size_t max_iter_num)
{
    if ((A.rows() != A.cols()) || (A.diagonal().prod() == 0) || (A.rows() != b.size()))
        throw std::logic_error(std::string(__func__) + std::string("Matrix A is not square"));

    Eigen::MatrixXf U = A.triangularView<Eigen::StrictlyUpper>();
    Eigen::MatrixXf L = A.triangularView<Eigen::StrictlyLower>();
    Eigen::MatrixXf D = A - L - U;

    Eigen::MatrixXf B = D.inverse() * (-L - U);
    Eigen::MatrixXf f = D.inverse() * b;

    return IterativeCycle(B, f, epsilon, max_iter_num);
}

Eigen::VectorXf GaussSeide(Eigen::MatrixXf A, Eigen::VectorXf b, double epsilon, size_t max_iter_num)
{
    if ((A.rows() != A.cols()) || (A.diagonal().prod() == 0) || (A.rows() != b.size()))
        throw std::logic_error(std::string(__func__) + std::string("Matrix A is not square"));

    Eigen::MatrixXf U = A.triangularView<Eigen::StrictlyUpper>();
    Eigen::MatrixXf L = A.triangularView<Eigen::StrictlyLower>();
    Eigen::MatrixXf D = A - L - U;

    Eigen::MatrixXf B = (D + L).inverse() * (-U);
    Eigen::MatrixXf f = (D + L).inverse() * b;

    return IterativeCycle(B, f, epsilon, max_iter_num);
}

Eigen::VectorXf SOR(Eigen::MatrixXf A, Eigen::VectorXf b, double epsilon, size_t max_iter_num, double w)
{
    if ((A.rows() != A.cols()) || (A.diagonal().prod() == 0) || (A.rows() != b.size()))
        throw std::logic_error(std::string(__func__) + std::string("Matrix A is not square"));

    Eigen::MatrixXf U = A.triangularView<Eigen::StrictlyUpper>();
    Eigen::MatrixXf L = A.triangularView<Eigen::StrictlyLower>();
    Eigen::MatrixXf D = A - L - U;

    Eigen::MatrixXf B = (D + w * L).inverse() * ((1 - w) * D - w * U);
    Eigen::MatrixXf f = w * (D + w * L).inverse() * b;

    return IterativeCycle(B, f, epsilon, max_iter_num);
}
} // namespace NOL
