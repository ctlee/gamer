#include <iostream>
#include <memory>
#include <map>
#include <array>
#include <vector>
#include <queue>
#include <set>
#include <tuple>
#include <initializer_list>
#include <math.h>
#include <algorithm>
#include <tuple>
#include <future>
#include <experimental/optional>
#include "SurfaceMesh.h"
#include "SimplicialComplex.h"
#include <cstring>


int main(int argc, char* argv[])
{
	return 0;
	if(argc != 3)
	{
		std::cerr << "Wrong arguments passed" << std::endl;
		return -1;
	}
	auto pF = OFF_to_SimplicialComplex(argv[1]);

	init_orientation(*pF);
	clear_orientation(*pF);
	auto rval = compute_orientation(*pF);

	std::cout << std::get<0>(rval) << " " << std::get<1>(rval) << " " << std::get<2>(rval) << std::endl;

	if(strcmp(argv[2], "true") == 0)
	{
		if(std::get<1>(rval))
		{
			return 0;
		}
		else
		{
			std::cerr << "False negative. Algorithm decided that \"" << argv[1] << "\" was not orientable, but should have decided that it was." << std::endl;
			return -1;
		}
	}
	else
	{
		if(std::get<1>(rval))
		{
			std::cerr << "False positive. Algorithm decided that \"" << argv[1] << "\" was orientable, but should have decided that it was not." << std::endl;
			return -1;
		}
		else
		{
			return 0;
		}
	}
}

