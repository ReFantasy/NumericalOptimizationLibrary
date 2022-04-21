#include "helper.h"

using namespace NOL;
/**
 * @brief Based iterative cycle to solve linear equations
 * @param B
 * @param f
 * @param epsilon
 * @return
 */
Vector IterativeCycle(Matrix B, Matrix f, FLOAT epsilon, size_t max_iter_num)
{
    auto x = f;

    size_t n = 0;
    while (true)
    {
        if (n >= max_iter_num)
            return x;

        Vector new_x = B * x + f;
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

LinearEquationSolver::SOLVER_TYPE LinearEquationSolver::solver_type = LinearEquationSolver::SOLVER_TYPE::GAUSS_SEIDE;
FLOAT LinearEquationSolver::EPSILON = 0.00001;
size_t LinearEquationSolver::MAX_ITERATION = 100000;
FLOAT LinearEquationSolver::W = 1.4;

Vector LinearEquationSolver::Solver(Matrix A, Vector b)
{
    Vector x;
    
    // 优化 Ax=b 的求解
    // MAx = Mb
    auto diag = A.diagonal();
    Matrix M(diag.size(),diag.size());
    if(diag.dot(diag)!=0)
    {
        for(int i = 0;i<diag.size();i++)
        {
            M(i,i) = 1.0/diag(i);
        }
    }
    A = M*A;
    b = M*b;
    
    switch (solver_type)
    {
    case SOLVER_TYPE::EIGEN_SOLVER:
        x = A.colPivHouseholderQr().solve(b);
        break;
    case SOLVER_TYPE::JACOBI:
        x = Jacobi(A, b, EPSILON, MAX_ITERATION);
        break;
    case SOLVER_TYPE::GAUSS_SEIDE:
        x = GaussSeide(A, b, EPSILON, MAX_ITERATION);
        break;
    case SOLVER_TYPE::SOR:
        x = Sor(A, b, EPSILON, MAX_ITERATION, W);
        break;
    default:
        x = A.colPivHouseholderQr().solve(b);
    }
    return x;
}
Vector LinearEquationSolver::Jacobi(Matrix A, Vector b, FLOAT epsilon, size_t max_iter_num)
{
    if ((A.rows() != A.cols()) || (A.diagonal().prod() == 0) || (A.rows() != b.size()))
        throw std::logic_error(std::string(__func__) + std::string("Matrix A is not square"));

    Matrix U = A.triangularView<Eigen::StrictlyUpper>();
    Matrix L = A.triangularView<Eigen::StrictlyLower>();
    Matrix D = A - L - U;

    Matrix B = D.inverse() * (-L - U);
    Matrix f = D.inverse() * b;

    return IterativeCycle(B, f, epsilon, max_iter_num);
}

Vector LinearEquationSolver::GaussSeide(Matrix A, Vector b, FLOAT epsilon, size_t max_iter_num)
{
    if ((A.rows() != A.cols()) || (A.diagonal().prod() == 0) || (A.rows() != b.size()))
        throw std::logic_error(std::string(__func__) + std::string("Matrix A is not square"));

    Matrix U = A.triangularView<Eigen::StrictlyUpper>();
    Matrix L = A.triangularView<Eigen::StrictlyLower>();
    Matrix D = A - L - U;

    Matrix B = (D + L).inverse() * (-U);
    Matrix f = (D + L).inverse() * b;

    return IterativeCycle(B, f, epsilon, max_iter_num);
}

Vector LinearEquationSolver::Sor(Matrix A, Vector b, FLOAT epsilon, size_t max_iter_num, FLOAT w)
{
    if ((A.rows() != A.cols()) || (A.diagonal().prod() == 0) || (A.rows() != b.size()))
        throw std::logic_error(std::string(__func__) + std::string("Matrix A is not square"));

    Matrix U = A.triangularView<Eigen::StrictlyUpper>();
    Matrix L = A.triangularView<Eigen::StrictlyLower>();
    Matrix D = A - L - U;

    Matrix B = (D + w * L).inverse() * ((1 - w) * D - w * U);
    Matrix f = w * (D + w * L).inverse() * b;

    return IterativeCycle(B, f, epsilon, max_iter_num);
}
