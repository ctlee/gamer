#include "index_tracker.h"
#include <iostream>
#include <assert.h>
#include <array>
#include <vector>
#include <cstdlib>
#include <limits>
#include <cmath>
#include <sys/time.h>
#include "SurfaceMesh.h"


template <typename TypeGenerator>
struct lazy_tuple
{
	using T = typename TypeGenerator::head;

	T data;
	lazy_tuple<typename TypeGenerator::tail>* next;
};

/*
double compute_risk(double p, int n, int N)
{
	return p*pow(1-p,n);
}

double compute_expected(double V, int n, int N)
{
	return V*(N-n)/N;
}
*/


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

	index_tracker<long> indices;
	std::vector<int> active;
	constexpr size_t N = 5;
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
			for(int k = 0; k < 100; ++k)
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
		std::cout << indices << std::endl;
	}

	clock_gettime(CLOCK_THREAD_CPUTIME_ID, &stop);

	std::cout << ((stop.tv_sec - start.tv_sec) * 10e8 + stop.tv_nsec - start.tv_nsec)/10e8	 << std::endl;

//	auto pF = readOBJ("../data/neuron/neuron.obj");
//	compute_orientation(*pF);
//	writeOFF("neuron.off", *pF);
/*
	{
		int N = 730;
		double p = 1./(365*17);
		double V = 3000.;
		double sum = 0.;
		for(int i = 0; i < N; ++i)
		{
			double curr = compute_risk(p,i,N) * compute_expected(V,i,N);
			sum += curr;
			std::cout << " : " << curr << std::endl;
		}
		std::cout << sum << std::endl;
	}
*/
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
