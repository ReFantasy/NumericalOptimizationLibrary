/**
 * MIT License
 * Copyright (c) 2022 ReFantasy
 * https://github.com/ReFantasy/NumericalOptimizationLibrary
 *
 * This header file defines the basic data structure of the optimization library,
 * including basic data types, base classes of optimization function, optimization option and other auxiliary utility.
 */
#pragma once
#include "Eigen/Dense"
#include <chrono>
#include <iostream>
#include <sstream>

namespace NOL
{
using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;
using FLOAT = double;

/**
 * @brief the base class of the function to be optimized
 */
class TargetFunctor
{
  public:
    /**
     * @brief Overloading bracket operators to evaluate function values
     * @param x Input value of function
     * @return Return value of function
     */
    virtual FLOAT operator()(const Vector &x) const = 0;

    /**
     * @brief Compute the first derivative of a function
     * @param x The point of derivative
     * @return Derivative value
     */
    virtual Vector FirstOrderDerivatives(const Vector &x) const = 0;

    /**
     * @brief Compute the second derivative of a function
     * @param x The point of derivative
     * @return Derivative value (Jacobian)
     */
    virtual Matrix SecondOrderDerivatives(const Vector &x) const = 0;
};

/**
 * @brief Options for optimizing algorithms
 */
class Options
{
  public:
    Vector init_x;
    double gk_norm = 10e-5;

    /**
     * @brief Output iteration record of optimization process
     */
    void Summary() const
    {
        std::cout << ss.str() << std::endl;
    }

    /**
     * @brief Clear all records and logs
     */
    void ClearLog()
    {
        ss.clear();
        ss.str("");
    }

  public:
    /**
     * @brief Output to standard I / O device
     */
    bool print_log_to_stdio = true;

    /**
     * @brief if true, no optimization process data will be output and recorded
     */
    bool optimized_performance = false;

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

  private:
    std::stringstream ss;
};

/**
 * @brief Base class of all unconstrained optimization algorithm classes
 */
class UnconstrainedOptimizationBase
{
  public:
    /**
     * @brief Find the optimal solution of function
     * @param fucntor Function object
     * @param options Optimization parameters
     * @return the point of the optimal solution
     */
    virtual Vector Solve(TargetFunctor &fucntor, Options &options) = 0;
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

} // namespace NOL
