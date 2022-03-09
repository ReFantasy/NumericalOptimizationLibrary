/**
 * MIT License
 * Copyright (c) 2022 ReFantasy
 * https://github.com/ReFantasy/NumericalOptimizationLibrary
 *
 * This header file defines the basic data structure of the optimization library,
 * including basic data types, base classes of optimization function, optimization option and other auxiliary utility.
 */
#ifndef __GLOBAL_H__
#define __GLOBAL_H__
#include "Eigen/Dense"
#include <chrono>
#include <iostream>
#include <limits>
#include <sstream>

namespace NOL
{
using FLOAT = double;
using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

template<typename T>
class MinStepSize
{
public:
    static constexpr T value = 1e-9;
};

template<>
class MinStepSize<float>
{
public:
    static constexpr float value = 1e-5;
};

class LinearSearch;

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
    virtual FLOAT operator()(const Vector &xk) const = 0;

    /**
     * @brief Compute the first derivative of a function
     * @param x The point of derivative
     * @return Derivative value
     */
    virtual Vector Gradient(const Vector &xk) const = 0;

    /**
     * @brief Compute the second derivative of a function
     * @param x The point of derivative
     * @return Derivative value (Jacobian)
     */
    virtual Matrix Hesse(const Vector &xk) const = 0;
};

enum class LineSearchType
{
    GOLDENSECTION,
    ARMIJO,
    GOLDSTEIN,
    WOLFE,
    STRONGWOLFE
};


enum class QuasiNewtonSearchType
{
    SR1,
    DFP,
    BFGS
};

/**
 * @brief Options for optimizing algorithms
 */
class Options
{
  public:
    Vector init_x;
    FLOAT gk_norm = 1e-6;

    // If during the line search, the step_size falls below this
    // value, it is truncated to zero.
    FLOAT min_step_size = MinStepSize<FLOAT>::value;

    LineSearchType line_search_type = LineSearchType::GOLDSTEIN;
    QuasiNewtonSearchType quasi_newton_type = QuasiNewtonSearchType::DFP;


  public:
    FLOAT parameter_line_search_advance_and_retreat_alpha = 1.0;
    FLOAT parameter_line_search_advance_and_retreat_h = 1.0;
    FLOAT parameter_line_search_advance_and_retreat_t = 2;

    FLOAT parameter_line_search_golden_section_size = 1e-3;

    FLOAT parameter_line_search_armijo_rho = 0.001;
    FLOAT parameter_line_search_armijo_t = 2.0; // >1

    FLOAT parameter_line_search_goldstein_p = 0.25;

    FLOAT parameter_line_search_wolfe_rho = 0.1;
    FLOAT parameter_line_search_wolfe_sigma = 0.5;
    FLOAT parameter_line_search_wolfe_alpha_max = std::numeric_limits<FLOAT>::max();

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
class UnconstrainedOptimizationLineSearchBase
{
  public:
    /**
     * @brief Find the optimal solution of function
     * @param fucntor Function object
     * @param options Optimization parameters
     * @return the point of the optimal solution
     */
    virtual Vector Solve();

    virtual bool IsTerminated(const Vector &xk, int k) const = 0;

    virtual Vector SearchDirection(const Vector &xk) const = 0;

    virtual FLOAT Step(const Vector &xk, const Vector &dk) const = 0;

    TargetFunctor *_functor;
    Options *_options;
    LinearSearch *_line_search;
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

#endif //__GLOBAL_H__
