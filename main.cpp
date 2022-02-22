#include "example.h"
#include "line_search.h"
#include <Eigen/Core>
#include <cmath>
#include <iostream>

int main(int argc, char *argv[])
{
    /*class Search :public LineSearch
    {
    public:
        FLOAT Phi(FLOAT a)override
        {
            return -sin(a);
        }
        FLOAT dPhi_dx(FLOAT a)override { return -cos(a); };

    protected:
        bool Armijo(FLOAT a)override { return false; };
    };*/

    // example_3_1();
    Newton_example_3_1();
}
