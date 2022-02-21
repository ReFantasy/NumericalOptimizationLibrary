#include "line_search.h"

FLOAT LineSearch::Phi(FLOAT)
{
    return 0;
}

FLOAT LineSearch::Zerosixeight(FLOAT a0, FLOAT r0, FLOAT epsilon, FLOAT t )
{
    FLOAT secton_a, secton_b;
    AdvanceandRetreat(a0, r0, t, secton_a, secton_b);
    return GoldenSection(secton_a, secton_b, epsilon);
}

void LineSearch::AdvanceandRetreat(FLOAT a0, FLOAT r0, FLOAT t, FLOAT& secton_a, FLOAT& secton_b)
{
    assert(a0 >= 0);
    assert(r0 > 0);
    assert(t > 1);

    int i = 0;
    FLOAT a, ai, ai1, ri;
    a = ai = a0;
    ri = r0;

    FLOAT _a, _b;
    while (true)
    {
        ai1 = ai + ri;
        if (ai1 <= 0)
        {
            ai1 = 0;
        }
        else if(Phi(ai1)<Phi(ai))
        {
            //step 3
            ri = t * ri;
            a = ai;
            ai = ai1;
            i++;
            continue;
        }

        // step 4
        if (i == 0)
        {
            ri = -ri;
            a = ai1;
            i++;
        }
        else
        {
            _a = std::min(a, ai1);
            _b = std::max(a, ai1);
            break;
        }
    }
    
    secton_a = _a;
    secton_b = _b;
}

FLOAT LineSearch::GoldenSection(FLOAT secton_a, FLOAT secton_b, FLOAT epsilon /*= 10e-3*/)
{
    FLOAT a = secton_a;
    FLOAT b = secton_b;
    static const FLOAT r = (std::sqrt(5)-1)/2; //0.618

    while ((b - a) > epsilon)
    {
        FLOAT al = a + (1.0 - r) * (b - a);
        FLOAT ar = a + r * (b - a);

        if (Phi(al) < Phi(ar))
        {
            b = ar;
        }
        else
        {
            a = al;
        }
    }
    
    return (a + b) / 2.0;
}
