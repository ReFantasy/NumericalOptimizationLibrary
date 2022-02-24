#pragma once
#include "Eigen/Dense"
#include <chrono>
#include <iostream>
#include <sstream>

using TVector = Eigen::VectorXd;
using TMatrix = Eigen::MatrixXd;
using TFLOAT = double;

class TargetFunctor
{
  public:
    virtual TFLOAT operator()(const TVector &x) const
    {
        return TFLOAT{};
    };

    virtual TVector FirstOrderDerivatives(const TVector &x) const
    {
        return TVector{0, 0};
    };
    virtual TMatrix SecondOrderDerivatives(const TVector &x) const
    {
        return TMatrix{};
    };
};

class Options
{
  public:
    /**
     * Steepest Descent Options
     */
    TVector init_x0;
    double gk_norm = 10e-5;

    void Summary() const
    {
        std::cout << ss.str() << std::endl;
    }
    void ClearLog()
    {
        ss.clear();
        ss.str("");
    }

  public:
    bool print_log_to_stdio = true;
    bool optimized_performance = false;

    std::stringstream ss;

  public:
    template <typename T> Options &operator<<(T content)
    {
        if (optimized_performance)
            return *this;
        if (print_log_to_stdio)
            std::cout << content;
        ss << content;

        return *this;
    }
};

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
