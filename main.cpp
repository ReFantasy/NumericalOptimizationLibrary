#include "example.h"

int main(int argc, char *argv[])
{
	try
	{
		//example_3_1();
		//example_3_2();
		example_3_3();
	}
	catch (std::exception &e)
	{
		std::cout<<std::endl<<e.what()<<std::endl;
	}


}