/**
 *
 * 线搜索准则及其求步长的算法实现
 *
 * Author  : ReFantasy
 * Date    : 2022-02-21
 */
#ifndef __LINE_SEARCH_H__
#define __LINE_SEARCH_H__
#include "global.h"

class LineSearch
{
  public:
    /**
     * @brief 【0.618方法】求线搜索步长，
     * 该方法用于求近似满足精确线搜索准则的步长
     * @param a0 线搜索初始步长
     * @param r0 0.618方法搜索步长，r0>0
     * @param epsilon 搜索停止条件，当搜索区间小于 epsilon 时，停止搜索
     * @param t 搜索步长增长系数 assert(t>1)
     * @return 线搜索函数 FLOAT Phi(FLOAT a) 到达极值点时的搜索步长
     */
    FLOAT Zerosixeight(FLOAT a0, FLOAT r0 = 1.0, FLOAT epsilon = 10e-3, FLOAT t = 1.2);

    FLOAT QuadraticPolynomialInterpolation(FLOAT a0);

  public:
    /**
     * @brief 线搜索函数 phi(a) = f( xk + a*dk ), a>0,
     * 其中，xk 为优化函数 f(x) 当前迭代点，dk 为迭代下降方向，xk,dk可以通过类继承的方式定义为成员变量，
     * 然后针对具体问题，自定义并重载 float Phi(float) 函数
     * @param a 线搜索步长
     * @return 搜索函数在步长 a 时的函数值，即待优化函数的下一个迭代点处的函数值
     */
    virtual FLOAT Phi(FLOAT a) = 0;

    virtual FLOAT dPhi_dx(FLOAT a) = 0;

  protected:
    /**
     * @brief 非精确线搜索准则
     * @param a 线搜索步长
     * @return 满足准则返回 true
     */
    virtual bool Criterion(FLOAT a) = 0;

  private:
    void AdvanceandRetreat(FLOAT a0, FLOAT r0, FLOAT t, FLOAT &secton_a, FLOAT &secton_b);
    FLOAT GoldenSection(FLOAT secton_a, FLOAT secton_b, FLOAT epsilon = 10e-3);
};

#endif //__LINE_SEARCH_H__
