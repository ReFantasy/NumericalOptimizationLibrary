/**
 *
 * 一维线搜索
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
/**
 * 假定已得下降方向 \f$ d_k \f$，求步长 \f$ \alpha \f$ ,  使得\f$ \varphi (\alpha) = f(x_k + \alpha d_k) < f(x_k)
 * \f$的问题为一维搜索或线搜索问题。
 *
 */
class LinearSearch
{
  public:
    /**
     * @brief 构造函数
     * @param functor 目标函数\f$ f \f$ 的指针
     */
    explicit LinearSearch(std::shared_ptr<TargetFunctor> functor = nullptr);

    /**
     * @brief 线搜索
     * @param alpha 初始步长
     * @param options 搜索选项
     * @return 线搜索的步长结果
     */
    FLOAT Search(FLOAT alpha, const Options &options);

    /**
     * @brief 设置当前迭代点
     * @param \f$ x_k \f$ 迭代点位置
     */
    void SetXk(Vector xk)
    {
        _xk = xk;
    }

    /**
     * @brief 设置线目标函数下降方向
     * @param \f$ d_k \f$ 下降方向
     */
    void SetDk(Vector dk)
    {
        _dk = std::move(dk);
    }

    /**
     * @brief 获取当前迭代位置
     * @return
     */
    Vector Xk() const
    {
        return _xk;
    }

    /**
     * @brief 获取当前目标函数下降方向
     * @return
     */
    Vector Dk() const
    {
        return _dk;
    }

    /**
     * @brief 设置目标函数
     * @param \f$ functor \f$ 目标函数指针
     */
    void SetTargetFunctor(std::shared_ptr<TargetFunctor> functor)
    {
        _functor_ptr = std::move(functor);
    }

    /**
     * @brief 获取目标函数
     * @return 目标函数指针\f$ f \f$
     */
    std::shared_ptr<TargetFunctor> TargetFunctorPointer() const
    {
        return _functor_ptr;
    }

  protected:
    /**
     * @brief 黄金分割精确线搜索，
     * @param alpha 线搜索初始步长
     * @param options 参数选项
     * @return 线搜索函数 FLOAT Phi(FLOAT a) 到达极值点时的步长
     */
    FLOAT GoldenMethod(FLOAT alpha, const Options &options);

    /**
     * @brief 非精确线搜索 Armijo 准则。通过不断缩小参数范围，获取满足准则的可选参数区间
     *        在区间的迭代过程中，首先使用二次插值，若不满足条件，则进行三次插值
     * @param alpha 步长
     * @param options 参数选项
     * @return 一维搜索的步长
     */
    FLOAT Armijo(FLOAT alpha, const Options &options);

    /**
     * @brief 非精确线搜索 Goldstein 准则
     * @param alpha 步长
     * @param options 参数选项
     * @return 一维搜索的步长
     */
    FLOAT Goldstein(FLOAT alpha, const Options &options);

    /**
     * @brief 非精确线搜索 Wolfe 准则
     * @param alpha 步长
     * @param options 参数选项
     * @return 一维搜索的步长
     */
    FLOAT Wolfe(FLOAT alpha, const Options &options);

    /**
     * @brief Strong Wolfe condition line search.  This implementation
     *        is based on the pseudo-code algorithm presented in Nocedal & Wright [1] (p60-61)
     *        [1] Nocedal J., Wright S., Numerical Optimization, 2nd Ed., Springer, 1999.
     *
     * @param alpha 步长
     * @param options 参数选项
     * @return 一维搜索的步长
     */
    FLOAT StrongWolfe(FLOAT alpha, const Options &options);

  protected:
    void AdvanceAndRetreat(FLOAT a0, FLOAT h0, FLOAT t, FLOAT &secton_a, FLOAT &secton_b);
    FLOAT GoldenSection(FLOAT secton_a, FLOAT secton_b, FLOAT epsilon = 10e-4);
    FLOAT QuadraticInterpolationMinimum(FLOAT a1, FLOAT a2);
    FLOAT CubicInterpolationMinimum(FLOAT last_alpha, FLOAT alpha);
    FLOAT Interpolation(FLOAT last_alpha, FLOAT alpha, int k = 0);

  protected:
    /**
     * @brief 线搜索函数 \f$ \varphi (\alpha) = f( x_k + \alpha * d_k ), \alpha >0 \f$
     * 其中，\f$ x_k \f$ 为优化函数 \f$ f(x)\f$ 当前迭代点，\f$ d_k \f$ 为迭代下降方向
     * @param \f$ \alpha \f$ 线搜索步长
     * @return 搜索函数在步长 \f$ \alpha \f$ 时的函数值
     */
    virtual FLOAT phi(FLOAT a);

    virtual FLOAT dphi_da(FLOAT a);

    FLOAT Zoom(FLOAT alpha_lo, FLOAT alpha_hi, const Options &options);

  private:
    Vector _xk;
    Vector _dk;
    std::shared_ptr<TargetFunctor> _functor_ptr;
};
} // namespace NOL
#endif //__LINE_SEARCH_H__
