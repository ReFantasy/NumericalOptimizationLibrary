#include <iostream>

#include "example.h"
#include <sstream>
#include "global.h"

int main(int argc, char **argv)
{
    //example_3_1();
    //example_3_2();
    example_3_1_by_dampednewton();

    /*NOL::Matrix m = NOL::Matrix::Identity(3,3);
    std::cout << m.cols() << std::endl;*/

    //example_test_1();

    /*NOL::Vector v(4);
    NOL::Matrix m(v.size(), v.size());
    std::cout << v.size() << std::endl;*/


    //NOL::Matrix m(3, 3);
    //std::cout << NOL::Matrix::Identity(4, 4) << std::endl;

    return 0;
}
