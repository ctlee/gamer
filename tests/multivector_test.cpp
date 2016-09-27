#include "multivector.h"
#include "util.h"
#include <iostream>
#include <array>
#include <random>


template <typename Fn, typename... Ts>
struct flattenH {};

template <typename Fn, typename T, typename... Ts>
struct flattenH<Fn,T,Ts...> {
	template <std::size_t N>
	static void apply(Fn f, T head, Ts... tail)
	{
		f(N, head);
		flattenH<Fn,Ts...>::template apply<N+1>(f,tail...);
	}
};

template <typename Fn>
struct flattenH<Fn> {
	template <std::size_t N>
	static void apply(Fn f) {}
};

template <typename Fn, std::size_t K, typename T, typename... Ts>
struct flattenH<Fn, std::array<T,K>, Ts...> {
	template <std::size_t N>
	static void apply(Fn f, const std::array<T,K>& head, Ts... tail)
	{
		for(std::size_t k = 0; k < K; ++k)
		{
			f(N+k,head[k]);
		}
		flattenH<Fn,Ts...>::template apply<N+K>(f,tail...);
	}
};

template <typename Fn, typename... Ts>
void flatten(Fn f, Ts... args)
{
	flattenH<Fn,Ts...>::template apply<0>(f,args...);
}

template <typename... Ts>
void testme(Ts... args)
{
	std::cout << sizeof...(Ts) << std::endl;
}

int main()
{
	multivector<double, 1> v0(3);
	multivector<double, 1> v1(3);
	multivector<double, 1> v2(3);
/*
	std::default_random_engine generator;
	std::normal_distribution<double> distribution;

	for(auto& curr : v0)
	{
		curr = distribution(generator);
	}
	for(auto& curr : v1)
	{
		curr = distribution(generator);
	}
*/
	v0[0] = 1; v0[1] = 0; v0[2] = 0;
	v1[0] = 0; v1[1] = 1; v1[2] = 0;
	v2[0] = 0; v2[1] = 0; v2[2] = 1;

	std::cout << "Test Complete" << std::endl;

	std::array<std::size_t,5> x;
	std::array<std::size_t,3> a{{2,3,4}};
	std::array<std::size_t,2> b{{7,11}};
	flatten([&x](std::size_t i, std::size_t y){ x[i] = y; std::cout << i << " " << y << std::endl; },a,b);

	for(auto a : x)
	{
		std::cout << a << " ";
	}
	std::cout << std::endl;
	return -1;
}
