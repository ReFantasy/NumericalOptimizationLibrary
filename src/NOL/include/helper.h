#ifndef __HELPER_H__
#define __HELPER_H__
#include "Eigen/Dense"
#include <chrono>
#include <iostream>
#include <random>
class Timer
{
  public:
    Timer()
    {
        _start_time = std::chrono::high_resolution_clock::now();
    }

    void Start()
    {
        ReSet();
    }

    std::chrono::high_resolution_clock::time_point StartTime()
    {
        return _start_time;
    }
    std::chrono::high_resolution_clock::time_point CurrentTime()
    {
        return std::chrono::high_resolution_clock::now();
    }

    template <typename T = std::chrono::milliseconds> int Elapse() const
    {
        return std::chrono::duration_cast<T>(std::chrono::high_resolution_clock::now() - _start_time).count();
    }

    void ReSet()
    {
        _start_time = std::chrono::high_resolution_clock::now();
    }

  private:
    std::chrono::high_resolution_clock::time_point _start_time;
};

template <typename T, typename Distribution = std::uniform_real_distribution<T>>
T RandomNumber(const T &lo = 0, const T &hi = 1)
{
    std::random_device dev;
    std::mt19937 rng(dev());
    Distribution dist(lo, hi);
    return dist(rng);
}

template <typename T = bool> T ALL_OF(T b)
{
    return b;
}
template <typename T, typename... args> bool ALL_OF(T b1, args... b2)
{
    return b1 && ALL_OF(b2...);
}

template <typename T = bool> T ANY_OF(T b)
{
    return b;
}
template <typename T, typename... args> bool ANY_OF(T b1, args... b2)
{
    return b1 || ANY_OF(b2...);
}

/*****************************************************************************************************
 *
 *                             线性方程组 Ax=b 的数值求解（迭代法, A 为方阵）
 *
 *****************************************************************************************************/

namespace NOL
{
Eigen::VectorXf Jacobi(Eigen::MatrixXf A, Eigen::VectorXf b, double epsilon = 0.000001, size_t max_iter_num = 100);

Eigen::VectorXf GaussSeide(Eigen::MatrixXf A, Eigen::VectorXf b, double epsilon = 0.000001, size_t max_iter_num = 100);

Eigen::VectorXf SOR(Eigen::MatrixXf A, Eigen::VectorXf b, double epsilon = 0.000001, size_t max_iter_num = 100,
                    double w = 1.4);

} // namespace NOL
#endif //__HELPER_H__