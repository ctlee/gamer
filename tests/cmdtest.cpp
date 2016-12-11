#include "index_tracker.h"
#include <iostream>
#include <assert.h>
#include <array>
#include <vector>
#include <cstdlib>
#include <limits>
#include <sys/time.h>


using BT = BTreeNode<int,12>;

int main(int argc, char *argv[])
{
	auto A = Interval<int>(0,2);
	auto B = Interval<int>(3,5);
	auto C = Interval<int>(2,3);

	std::cout << (A < B) << std::endl;
	std::cout << (A < C) << std::endl;
	std::cout << (C < B) << std::endl;
	std::cout << (B > A) << std::endl;

	std::cout << (A > B) << std::endl;
	std::cout << (B < A) << std::endl;

	for(auto i = A.lower(); i < A.upper(); ++i)
	{
		std::cout << i << " ";
	}
	std::cout << std::endl;

    timespec start, stop;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);

	index_tracker<int,4> indices;
	std::vector<int> active;
	constexpr size_t N = 50000;
	for(int i = 0; i < N; ++i)
	{
		if(std::rand() % 200 < 10)
		{
			for(int k = 0; k < 500; ++k)
			{
				active.push_back(indices.pop());
			}
		}
		else
		{
			for(int k = 0; k < 10; ++k)
			{
				if(!active.empty())
				{
					int index = std::rand() % active.size();
					indices.insert(active[index]);
					active.erase(active.begin() + index);
				}
			}
		}
//		int add = std::rand() % N;
//		remove_scalar<BT>(head, add);
//		head = insert_scalar<BT>(head, add);
//		check_order<BT>(head,-1);
//		added_list.push_back(add);
//		std::cout << indices << std::endl;
	}

	clock_gettime(CLOCK_THREAD_CPUTIME_ID, &stop);

	std::cout << ((stop.tv_sec - start.tv_sec) * 10e8 + stop.tv_nsec - start.tv_nsec)/10e8	 << std::endl;

//	std::cout << head << std::endl;
/*
	for(int x : added_list)
	{
		remove_scalar<BT>(head, x);
//		check_order<BT>(head,-1);
		std::cout << x << " - " << head << std::endl;
	}
*/
	/*
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
	while(head != nullptr)
	{
		pop<BT>(head);
		std::cout << head << std::endl;
//		std::cout << pop<BT>(head) << std::endl;
	}
	clock_gettime(CLOCK_THREAD_CPUTIME_ID, &stop);

	std::cout << ((stop.tv_sec - start.tv_sec) * 10e8 + stop.tv_nsec - start.tv_nsec)/10e8	 << std::endl;

	std::cout << "Done." << std::endl;
	*/
/*
	Data<BT> x;
	std::cout << " 1: ";
	print<BT>(split<BT>(head, x));
	std::cout << x << std::endl;
	std::cout << " 2: ";
	print<BT>(head);
	*/
}
