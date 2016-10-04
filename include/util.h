#pragma once

#include <utility>
#include <array>

namespace util
{
    template<typename T> struct range
    {
        template <class C>
        range(C &&c) : _begin(c.begin()), _end(c.end()) {}

        range(T b, T e) : _begin(b), _end(e) {}

        T begin() { return _begin; }
        T end() { return _end; }

    private:
        T _begin;
        T _end;
    };

    template<typename T> range<T> make_range(T b, T e)
    {
        return range<T>(std::move(b), std::move(e));
    }

    template<typename T> range<T> make_range(std::pair<T,T> p)
    {
        return range<T>(std::move(p.first), std::move(p.second));
    }

	template <typename... Ts>
	struct type_holder
	{
		static const size_t size = sizeof...(Ts);
	};

	template <typename T, typename... Ts>
	struct type_holder<T, Ts...>
	{
		using head = T;
		using tail = type_holder<Ts...>;
		static const size_t size = 1 + type_holder<Ts...>::size;
	};

	template <size_t k, typename T>
	struct type_get {};

	template <typename... Ts>
	struct type_get<0, type_holder<Ts...>>
	{
		using type = typename type_holder<Ts...>::head;
	};

	template <size_t k, typename... Ts>
	struct type_get<k, type_holder<Ts...>>
	{
		using type = typename type_get<k-1, typename type_holder<Ts...>::tail>::type;
	};



	template <class C, template <typename> class V, typename... Rs>
	struct type_map_helper {};

	template <template <class...> class G, template <typename> class V, typename... Rs>
	struct type_map_helper<G<>, V, Rs...>
	{
		using type = G<Rs...>;
	};

	template <typename... Rs, template < class...> class G, typename T, typename... Ts, template <typename> class V>
	struct type_map_helper<G<T, Ts...>, V, Rs...>
	{
		using type = typename type_map_helper<G<Ts...>, V, Rs..., V<T>>::type;
	};

	template <class C, template <typename> class V>
	struct type_map
	{
		using type = typename type_map_helper<C, V>::type;
	};

	template <class IntegerType, template <class...> class OutHolder, class IntegerSequence, template <IntegerType> class F, typename... Accumulators>
	struct int_type_map_helper {};

	template <class Integer, template <class...> class OutHolder, template <class, Integer...> class InHolder, template <Integer> class F, class... Accumulator>
	struct int_type_map_helper<Integer, OutHolder, InHolder<Integer>, F, Accumulator...>
	{
		using type = OutHolder<Accumulator...>;
	};

	template <class Integer, template <class...> class OutHolder, template <class, Integer...> class InHolder, Integer I, Integer... Is, template <Integer> class F, class... Accumulator>
	struct int_type_map_helper<Integer, OutHolder, InHolder<Integer, I, Is...>, F, Accumulator...>
	{
		using type = typename int_type_map_helper<Integer, OutHolder, InHolder<Integer,Is...>, F, Accumulator..., F<I>>::type;
	};

	template <class IntegerType, template <class...> class OutHolder, class IntegerSequence, template <IntegerType> class F>
	struct int_type_map
	{
		using type = typename int_type_map_helper<IntegerType, OutHolder, IntegerSequence, F>::type;
	};

	template <template <class...> class TUPLE, typename HOLDER_FULL>
	struct type_swap
	{};

	template <template <class...> class TUPLE, template <class...> class HOLDER, typename... Ts>
	struct type_swap<TUPLE, HOLDER<Ts...>>
	{
		using type = TUPLE<Ts...>;
	};



	template <typename S, std::size_t depth, std::size_t N, typename T>
	void fill_arrayH(std::array<S,N>& arr, S arg)
	{
		static_assert(depth + 1 == N, "Size of array must match number of input arguments");
		arr[depth] = arg;
	}

	template <typename S, std::size_t depth, std::size_t N, typename T, typename... Ts>
	void fill_arrayH(std::array<S,N>& arr, S head, Ts... tail)
	{
		arr[depth] = head;
		fill_arrayH<S,depth+1,N,Ts...>(arr, tail...);
	}

	template <typename S, std::size_t N, typename... Ts>
	void fill_array(std::array<S,N>& arr, Ts... args)
	{
		static_assert(sizeof...(args) == N, "Size of array must match number of input arguments");
		fill_arrayH<S,0,N,Ts...>(arr, args...);
	}


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
}
