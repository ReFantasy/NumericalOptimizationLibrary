/**
 *
 * 线搜索准则及其求步长的算法实现
 *
 * 原理：假定已得下降方向 dk，求步长 alpha 的问题为一维搜索或线搜索问题，
 *      具体来说就是求解优化的目标函数 f(x) 在 xk 位置，沿着 dk 方向，迭代到新的位置 x(k+1) = xk + alpha*dk，
 *      在 x(k+1)处时，使得 alpha 满足 f(xk(+1)) <= f(xk)
 *
 * Author  : ReFantasy
 * Date    : 2022-02-21
 */
#ifndef __LINE_SEARCH_H__
#define __LINE_SEARCH_H__
#include "global.h"
#include <limits>

namespace NOL
{

class LineSearch
{
  public:
    LineSearch(TargetFunctor *functor = nullptr) : _functor(functor)
    {
    }

    FLOAT Search(FLOAT alpha, const Options &options)
    {
        FLOAT step_length = alpha;

        switch (options.line_search_type)
        {
        case LineSearchType::GOLDENSECTION:
            step_length = GoldenMethod(alpha, options);
            break;
        case LineSearchType::QUADRATIC:
            step_length = QuadraticInterpolation(alpha, options);
            break;
        case LineSearchType::ARMIJO:
            step_length = Armijo(alpha, options);
            break;
        case LineSearchType::GOLDSTEIN:
            step_length = Goldstein(alpha, options);
            break;
        case LineSearchType::STRONGWOLFE:
            step_length = StrongWolfe(alpha, options);
            break;
        default:
            step_length = 1.0;
        }
        return step_length;
    }

  public:
    Vector xk;
    Vector dk;
    TargetFunctor *_functor = nullptr;

  protected:
    /**
     * @brief 【0.618方法】[精确线搜索]求线搜索步长，
     * 该方法用于求近似满足精确线搜索准则的步长
     * @param a0 线搜索初始步长
     * @param h0 0.618方法搜索步长，h0>0
     * @param epsilon 搜索停止条件，当搜索区间小于 epsilon 时，停止搜索
     * @param t 搜索步长增长系数 assert(t>1)
     * @return 线搜索函数 FLOAT Phi(FLOAT a) 到达极值点时的搜索步长
     */
    FLOAT GoldenMethod(FLOAT a0, const Options &options);
    FLOAT QuadraticInterpolation(FLOAT a0, const Options &options);

    FLOAT Armijo(FLOAT alpha, const Options &options);
    FLOAT Goldstein(FLOAT alpha, const Options &options);
    FLOAT StrongWolfe(FLOAT alpha, const Options &options);

  protected:
    /**
     * @brief 线搜索函数 phi(a) = f( xk + a*dk ), a>0,
     * 其中，xk 为优化函数 f(x) 当前迭代点，dk 为迭代下降方向，xk,dk可以通过类继承的方式定义为成员变量，
     * 然后针对具体问题，自定义并重载 float Phi(float) 函数
     * @param a 线搜索步长
     * @return 搜索函数在步长 a 时的函数值，即待优化函数的下一个迭代点处的函数值
     */
    virtual FLOAT phi(FLOAT a);

    virtual FLOAT dphi_da(FLOAT a);

  protected:
    void AdvanceAndRetreat(FLOAT a0, FLOAT h0, FLOAT t, FLOAT &secton_a, FLOAT &secton_b);
    FLOAT GoldenSection(FLOAT secton_a, FLOAT secton_b, FLOAT epsilon = 10e-3);

    // 二次插值极小值点
    FLOAT QuadraticInterpolationMinimum(FLOAT a1, FLOAT a2);

    FLOAT Zoom(FLOAT alpha_lo, FLOAT alpha_hi, const Options &options);
};
} // namespace NOL
#endif //__LINE_SEARCH_H__
