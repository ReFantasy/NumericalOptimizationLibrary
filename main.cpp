#include "example.h"

using namespace NOL;

int main(int argc, char *argv[])
{
    try
    {
        //example_3_1();
        //example_3_2();
		//example_cg();
        //example_3_3();
        // example();

		// 求解线性方程组
		Eigen::Matrix3f A;
		Eigen::Vector3f b;
		A<<8,-3,2,4,11,-1,6,3,12;
		b <<20,33,36;
		Eigen::Vector3f  x = SOR(A, b);
		std::cout<<x<<std::endl;

    }
    catch (std::exception &e)
    {
        std::cout << std::endl << e.what() << std::endl;
    }
    catch (...)
    {
        std::cout << "error" << std::endl;
    }
}