#include "example.h"


class Base
{
public:
	virtual void run()
	{
		std::cout<<"base run"<<std::endl;
	}
};

class Dec:public Base
{
public:
	void SetBase(Base base)
	{
		this->base = base;
	}
	void run()override
	{
		std::cout<<"do some decorator"<<std::endl;
		base.run();
	}
private:
	Base base;
};
int main(int argc, char *argv[])
{
	try
	{
		//example_3_1();
		example_3_2();
	}
	catch (std::exception &e)
	{
		std::cout<<std::endl<<e.what()<<std::endl;
	}


}