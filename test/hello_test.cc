#include <gtest/gtest.h>
#include <iostream>
#include "global.h"
#include "steepest_descent.h"

using namespace NOL;
// Demonstrate some basic assertions.
//TEST(suit1, test1)
//{
    //struct Functor : TargetFunctor
    //{
    //    virtual FLOAT operator()(const Vector& x) const override
    //    {
    //        return true;
    //    }
    //    virtual Vector FirstOrderDerivatives(const Vector& x) const override
    //    {
    //        return (G + G.transpose()) * x / 2 + b;
    //    }
    //    virtual Matrix SecondOrderDerivatives(const Vector& x) const override
    //    {
    //        return G;
    //    }

    //    Matrix G{ 2, 2 };
    //    Vector b{ 2 };
    //    float c = 10;
    //};

    //Functor functor;

    //Vector x0(2);
    //x0 << -30, 100;

    //functor.G << 21, 4, 4, 15;
    //functor.b << 2, 3;

    //Options option;
    //option.init_x = x0;
    //option.optimized_performance = true;
    //SteepestDescent sd(functor, option);
    //Vector xk =  sd.Solve(functor, option);
    ////std::cout << xk << std::endl;
    //Vector valid(2);
    //valid << -0.06, -0.18;

    //EXPECT_LE(xk(0), valid(0));
    //EXPECT_LE(xk(1), valid(1));
//}


class A
{
public:
    void f1() {}
};

class ATest :public ::testing::Test
{
protected:
    void SetUp()override
    {
        
    }

    A a1;
    A a2;
    A a3;
};

TEST_F(ATest, test1)
{
    a1.f1();
}
