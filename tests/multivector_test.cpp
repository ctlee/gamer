#include "multivector.h"
#include <iostream>

int main()
{
	multivector<std::string, 2> v({5,5});
	v[{1,2}] = "Hello, sir!";
	std::cout << v[{1,2}] << std::endl;
	std::cout << "dude" << std::endl;

	for(std::string curr : v)
	{
		std::cout << curr << std::endl;
	}

	v.resize({6,6});

	int m = 0;
	for(auto curr = v.index_begin(); curr != v.index_end() && m < 100; ++curr, ++m)
	{
		for(auto x : (*curr))
		{
			std::cout << x << " ";
		}
		std::cout << v[*curr] << std::endl;
	}

	std::cout << "Test Complete" << std::endl;

	return 0;
}
