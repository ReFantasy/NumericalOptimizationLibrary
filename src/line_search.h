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

enum class CriterionType
{
    Armijo,
    Goldstein
};

class LineSearch
{
  public:
    LineSearch(TargetFunctor *functor = nullptr) : _functor(functor)
    {
    }
    /**
     * @brief 【0.618方法】求线搜索步长，
     * 该方法用于求近似满足精确线搜索准则的步长
     * @param a0 线搜索初始步长
     * @param r0 0.618方法搜索步长，r0>0
     * @param epsilon 搜索停止条件，当搜索区间小于 epsilon 时，停止搜索
     * @param t 搜索步长增长系数 assert(t>1)
     * @return 线搜索函数 FLOAT Phi(FLOAT a) 到达极值点时的搜索步长
     */
    TFLOAT Zerosixeight(TFLOAT a0, TFLOAT r0 = 1.0, TFLOAT epsilon = 10e-3, TFLOAT t = 1.2);

    TFLOAT QuadraticPolynomialInterpolation(TFLOAT a0);

  public:
    TVector xk;
    TVector dk;
    TargetFunctor *_functor = nullptr;
    CriterionType _criterion_type = CriterionType::Goldstein;

  protected:
    /**
     * @brief 线搜索函数 phi(a) = f( xk + a*dk ), a>0,
     * 其中，xk 为优化函数 f(x) 当前迭代点，dk 为迭代下降方向，xk,dk可以通过类继承的方式定义为成员变量，
     * 然后针对具体问题，自定义并重载 float Phi(float) 函数
     * @param a 线搜索步长
     * @return 搜索函数在步长 a 时的函数值，即待优化函数的下一个迭代点处的函数值
     */
    virtual TFLOAT Phi(TFLOAT a);

    virtual TFLOAT dPhi_da(TFLOAT a);

    /**
     * @brief 非精确线搜索准则
     * @param a 线搜索步长
     * @return 满足准则返回 true
     */
    virtual bool Criterion(TFLOAT a);

  private:
    void AdvanceandRetreat(TFLOAT a0, TFLOAT r0, TFLOAT t, TFLOAT &secton_a, TFLOAT &secton_b);
    TFLOAT GoldenSection(TFLOAT secton_a, TFLOAT secton_b, TFLOAT epsilon = 10e-3);
};

#endif //__LINE_SEARCH_H__