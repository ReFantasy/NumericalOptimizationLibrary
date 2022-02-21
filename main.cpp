#include "example.h"
#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <format>
#include "line_search.h"

int main(int argc, char *argv[])
{
    std::string message = std::format("The answer is {}.", 42);
    assert(message == "The answer is 42.");

    /*class Search :public LineSearch
    {
    public:
        FLOAT Phi(FLOAT a)override
        {
            return std::pow(2 * a - 3, 2) - 4;
        }
    };

    Search search;
    std::cout << search.Zerosixeight(0.5, 0.2, 10e-6) << std::endl;*/
}
