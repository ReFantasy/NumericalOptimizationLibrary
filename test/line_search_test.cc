#include <gtest/gtest.h>
#include <iostream>
#include "global.h"
#include "line_search.h"

using namespace NOL;


TEST(TestLineSearch, AdvanceAndRetreat)
{
    class TLineSearch :public LineSearch
    {
    public:
        void advance_and_retreat(FLOAT a0, FLOAT h0, FLOAT t, FLOAT& secton_a, FLOAT& secton_b)
        {
            AdvanceAndRetreat(a0, h0, t, secton_a, secton_b);
        }
    protected:
        // phi(x) = (x-3)^2-4
        FLOAT phi(FLOAT x)override
        {
            return (x - 3) * (x - 3) - 4;
        }

        FLOAT dphi_da(FLOAT x)override
        {
            return 2 * (x - 3);
        }
    };


    TLineSearch line_search;


    FLOAT a0 = 0;
    FLOAT h0 = 1;
    FLOAT t = 1.5;
    FLOAT a, b;
    line_search.advance_and_retreat(a0, h0, t, a, b);
    ASSERT_LE(a, 3.0);
    ASSERT_GE(b, 3.0);

    a0 = 1;
    h0 = 4;
    line_search.advance_and_retreat(a0, h0, t, a, b);
    ASSERT_LE(a, 3.0);
    ASSERT_GE(b, 3.0);

    a0 = 6;
    line_search.advance_and_retreat(a0, h0, t, a, b);
    ASSERT_LE(a, 3.0);
    ASSERT_GE(b, 3.0);

    h0 = -1;
    line_search.advance_and_retreat(a0, h0, t, a, b);
    ASSERT_LE(a, 3.0);
    ASSERT_GE(b, 3.0);
}


