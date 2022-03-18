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
#include "helper.h"
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <utility>

namespace NOL
{
using FLOAT = double;
using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

/**
 * If during the line search, the step_size falls below this value, it is truncated to zero.
 */
template <typename T> class MinStepSize
{
  public:
    static constexpr T value = 1e-9;
};

template <> class MinStepSize<float>
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

enum class TerminationCriterionType
{
    GK_NORM,
    DELTA_XK,
    DELTA_F
};
/**
 * @brief Options for optimizing algorithms
 */
class Options
{
  public:
    Vector init_x;
    TerminationCriterionType termination_type = TerminationCriterionType::GK_NORM;
    FLOAT termination_value = 1e-6;

    FLOAT min_step_size = 0.001;

    FLOAT max_solver_time_in_seconds = 100000000.0;
    FLOAT max_line_search_time_in_milliseconds = 50;

    LineSearchType line_search_type = LineSearchType::GOLDSTEIN;
    QuasiNewtonSearchType quasi_newton_type = QuasiNewtonSearchType::DFP;

  public:
    FLOAT parameter_line_search_advance_and_retreat_h = 1.0;
    FLOAT parameter_line_search_advance_and_retreat_t = 1.5;

    FLOAT parameter_line_search_min_gold_section = 10e-3;

    FLOAT parameter_line_search_armijo_rho = 10e-3;

    FLOAT parameter_line_search_goldstein_rho = 0.15;

    FLOAT parameter_line_search_wolfe_rho = 10e-4;
    FLOAT parameter_line_search_wolfe_sigma = 0.6;
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
    friend class Decorator;

  public:
    UnconstrainedOptimizationLineSearchBase();
    /**
     * @param functor Function object pointer
     */
    explicit UnconstrainedOptimizationLineSearchBase(const std::shared_ptr<TargetFunctor>& functor_ptr);
    virtual ~UnconstrainedOptimizationLineSearchBase() = default;

    /**
     * Find the optimal solution of function
     * @return the point of the optimal solution
     */
    virtual Vector Solve();

    /**
     * whether the iteration stop condition is satisfied
     * @param xk current iteration point
     * @param k number of iteration
     * @return if it is satisfied, return true, otherwise false
     */
    virtual bool IsTerminated(const Vector &xk, int k) const;

    /**
     * search the direction at point xk to be sure that function value is reduced
     * @param xk the point of where search is started
     * @return descend direction
     */
    virtual Vector SearchDirection(const Vector &xk) const = 0;

    /**
     * compute the step length at point xk along the direction dk
     * @param xk point location
     * @param dk descend direction
     * @return  step length
     */
    virtual FLOAT Step(const Vector &xk, const Vector &dk) const = 0;

    virtual std::shared_ptr<Options> GetOptions() const
    {
        return _options_ptr;
    }

    virtual std::shared_ptr<Timer> GetTimer() const
    {
        return _timer_ptr;
    }

    virtual void SetFunctor(std::shared_ptr<TargetFunctor> functor_ptr)
    {
        _functor_ptr = std::move(functor_ptr);
    };

    virtual std::shared_ptr<TargetFunctor> GetFunctor() const
    {
        return _functor_ptr;
    }
	unsigned int NumOfIteration()const{return _K;}
  protected:
    std::shared_ptr<TargetFunctor> _functor_ptr;
    std::shared_ptr<Options> _options_ptr;
    std::shared_ptr<LinearSearch> _line_search_ptr;
    mutable std::shared_ptr<Timer> _timer_ptr = std::make_shared<Timer>();
	unsigned int _K = 0;
};

enum class OptimizationMethodType
{
    SD,
    NEWTON,
    DAMPED_NEWTON,
    LM,
    QUASI_NEWTON
};

class OptimizationFactory
{
  public:
    static std::shared_ptr<UnconstrainedOptimizationLineSearchBase> CreateSolver(
        OptimizationMethodType type, const std::shared_ptr<TargetFunctor> &functor = nullptr);
};
} // namespace NOL

#endif //__GLOBAL_H__
