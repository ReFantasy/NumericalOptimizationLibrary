/**
 *
 * 线搜索求步长
 * 针对具体问题的求解，需要自定义并重载 float Phi(float) 函数
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
    virtual FLOAT Phi(FLOAT);

    FLOAT Zerosixeight(FLOAT a0, FLOAT r0, FLOAT epsilon = 10e-3, FLOAT t=1.2);

private:
    void AdvanceandRetreat(FLOAT a0, FLOAT r0, FLOAT t, FLOAT &secton_a, FLOAT& secton_b);
    FLOAT GoldenSection(FLOAT secton_a, FLOAT secton_b, FLOAT epsilon = 10e-3);
};

#endif//__LINE_SEARCH_H__
