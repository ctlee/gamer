/*
 * ***************************************************************************
 * This file is part of the GAMer software.
 * Copyright (C) 2016-2018
 * by Christopher Lee, John Moody, Rommie Amaro, J. Andrew McCammon,
 *    and Michael Holst
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * ***************************************************************************
 */

/**
 * @file 	tensor.h
 * @brief Basic tensor library
 */

#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>
#include <sstream>

/// Namespace for all things gamer
namespace gamer
{
/// Namespace for tensor array mangement utilities
namespace array_util
{
/// @cond detail
/// Namespace for tensor array management details
namespace detail
{
template <typename S, std::size_t depth, std::size_t N, typename T>
void fill_arrayH(std::array<S, N> &arr, S arg)
{
    static_assert(depth + 1 == N, "Size of array must match number of input arguments");
    arr[depth] = arg;
}

template <typename S, std::size_t depth, std::size_t N, typename T, typename ... Ts>
void fill_arrayH(std::array<S, N> &arr, S head, Ts... tail)
{
    arr[depth] = head;
    fill_arrayH<S, depth+1, N, Ts...>(arr, tail ...);
}
} // end namespace detail
/// @endcond

/**
 * @brief      Fill an array with a list of values.
 *
 * @param      arr   Array to fill
 * @param[in]  args  List of values
 *
 * @tparam     S     Typename of array elements
 * @tparam     N     Number of array elements
 * @tparam     Ts    Typename of values
 */
template <typename S, std::size_t N, typename ... Ts>
void fill_array(std::array<S, N> &arr, Ts... args)
{
    static_assert(sizeof ... (args) == N, "Size of array must match number of input arguments");
    detail::fill_arrayH<S, 0, N, Ts...>(arr, args ...);
}

/// @cond detail
namespace detail
{
template <typename Fn, typename ... Ts>
struct flattenH {};

template <typename Fn, typename T, typename ... Ts>
struct flattenH<Fn, T, Ts...> {
    template <std::size_t N>
    static void apply(Fn f, T head, Ts... tail)
    {
        f(N, head);
        flattenH<Fn, Ts...>::template apply<N+1>(f, tail ...);
    }
};

template <typename Fn>
struct flattenH<Fn> {
    template <std::size_t N>
    static void apply(Fn f) {}
};

template <typename Fn, std::size_t K, typename T, typename ... Ts>
struct flattenH<Fn, std::array<T, K>, Ts...> {
    template <std::size_t N>
    static void apply(Fn f, const std::array<T, K> &head, Ts... tail)
    {
        for (std::size_t k = 0; k < K; ++k)
        {
            f(N+k, head[k]);
        }
        flattenH<Fn, Ts...>::template apply<N+K>(f, tail ...);
    }
};
} // end namespace detail
/// @endcond

/**
 * @brief      Flatten an array
 *
 * @param[in]  f     { parameter_description }
 * @param[in]  args  The arguments
 *
 * @tparam     Fn    { description }
 * @tparam     Ts    { description }
 */
template <typename Fn, typename ... Ts>
void flatten(Fn f, Ts... args)
{
    detail::flattenH<Fn, Ts...>::template apply<0>(f, args ...);
}
}



namespace detail {
	template <std::size_t k>
	struct factorial {
		constexpr static std::size_t value = k*factorial<k-1>::value;
	};

	template <>
	struct factorial<0> {
		constexpr static std::size_t value = 1;
	};

	template <std::size_t x, std::size_t n>
	struct pow {
		constexpr static std::size_t value = x * pow<x,n-1>::value;
	};

	template <std::size_t x>
	struct pow<x,0> {
		constexpr static std::size_t value = 1;
	};
}
/// @endcond

/**
 * @brief      Class for tensor.
 *
 * @tparam     _ElemType          { description }
 * @tparam     _vector_dimension  { description }
 * @tparam     _index_dimension   { description }
 */
template <typename _ElemType, std::size_t _vector_dimension, std::size_t _index_dimension>
class tensor {
public:
	constexpr static std::size_t index_dimension = _index_dimension;
	constexpr static std::size_t vector_dimension = _vector_dimension;
	constexpr static std::size_t total_dimension = detail::pow<vector_dimension, index_dimension>::value;
	using ElemType  = _ElemType;
	using IndexType = std::array<std::size_t, index_dimension>;
	using DataType  = std::array<ElemType, total_dimension>;

	/**
	 * @brief      { struct_description }
	 */
	struct index_iterator : public std::iterator<std::bidirectional_iterator_tag, IndexType>
	{
		using super = std::iterator<std::bidirectional_iterator_tag, IndexType>;
		/**
		 * @brief      { function_description }
		 *
		 * @param[in]  iter  The iterator
		 */
		index_iterator(const index_iterator& iter)
			: i(iter.i)
		{}
		index_iterator(const index_iterator&& iter)
			: i(std::move(iter.i))
		{}
		index_iterator()
		{
			i.fill(0);
			i[0] = vector_dimension;
		}
		index_iterator(int)
		{
			i.fill(0);
		}

		/**
		 * @brief      { operator_description }
		 *
		 * @return     { description_of_the_return_value }
		 */
		index_iterator& operator++()
		{
			std::size_t k = i.size() - 1;
			++(i[k]);
			while(k > 0)
			{
				if(i[k] == vector_dimension)
				{
					i[k] = 0;
					--k;
					++(i[k]);
				}
				else
				{
					break;
				}
			}
			return *this;
		}
		/**
		 * @brief      { operator_description }
		 *
		 * @param[in]  <unnamed>  { parameter_description }
		 *
		 * @return     { description_of_the_return_value }
		 */
		index_iterator operator++(int) { auto tmp = *this; ++(*this); return tmp; }

		/**
		 * @brief      { operator_description }
		 *
		 * @return     { description_of_the_return_value }
		 */
		index_iterator& operator--()
		{
			std::size_t k = index_dimension - 1;
			--(i[k]);
			while(k > 0)
			{
				if(i[k] < 0)
				{
					i[k] = 0;
					--k;
					--(i[k]);
				}
				else
				{
					break;
				}
			}
			return *this;
		}

		/**
		 * @brief      { operator_description }
		 *
		 * @param[in]  <unnamed>  { parameter_description }
		 *
		 * @return     { description_of_the_return_value }
		 */
		index_iterator operator--(int) { auto tmp = *this; --(*this); return tmp; }

		/**
		 * @brief      { operator_description }
		 *
		 * @param[in]  j     { parameter_description }
		 *
		 * @return     { description_of_the_return_value }
		 */
		bool operator==(index_iterator j) const
		{
			for(std::size_t k = 0; k < index_dimension; ++k)
			{
				if(i[k] != j.i[k])
				{
					return false;
				}
			}
			return true;
		}

		/**
		 * @brief      { operator_description }
		 *
		 * @param[in]  j     { parameter_description }
		 *
		 * @return     { description_of_the_return_value }
		 */
		bool operator!=(index_iterator j) const { return !(*this == j); }

		typename super::reference operator*() { return i; }
		typename super::pointer operator->() { return i; }
	protected:
		IndexType i;
	};
/*
	template <typename... Ts>
	tensor(Ts... args)
	{
		array_util::fill_array(_dimensions, args...);
		resize(_dimensions);
	}
*/
	tensor() { _data.fill(0); }

	tensor(const ElemType (&s)[total_dimension])
	{
		std::copy(std::begin(s), std::end(s), std::begin(_data));
	}

	tensor(const tensor& x) : _data(x._data) {}

	/**
	 * @brief      Constructs the object.
	 *
	 * @param[in]  x     { parameter_description }
	 */
	tensor(const tensor&& x) : _data(std::move(x._data)) {}

/*
	const IndexType& size()
	{
		return _dimensions;
	}
*/
	template <typename NumType,
			typename = std::enable_if_t<std::is_arithmetic<NumType>::value>,
			typename = std::enable_if_t<std::is_arithmetic<ElemType>::value>>
	operator tensor<NumType, _vector_dimension, _index_dimension>() const {
		tensor<NumType, _vector_dimension, _index_dimension> t;
		std::copy(_data.cbegin(), _data.cend(), t.begin());
		return t;
	}

	friend std::ostream& operator<<(std::ostream& output, const tensor& t){
		output  << "Tensor(";
		bool first = true;
		for(auto curr = t.index_begin(); curr != t.index_end(); ++curr){
			if(first) {output << "{"; first = false;}
			else output << "; {";

			bool first2 = true;
			for(auto x : (*curr)){
				if(first2) {output << x; first2 = false;}
				else output << "," << x;
			}
			output << "}:" << t[*curr];
		}
		output << ")";
		return output;
	}

	/**
	 * @brief      Returns a string representation of the object.
	 *
	 * @return     String representation of the object.
	 */
	std::string to_string() const {
		std::ostringstream output;
		output << *this;
		return output.str();
	}

	/**
	 * @brief      { function_description }
	 *
	 * @param      index  The index
	 *
	 * @tparam     Ts     { description }
	 *
	 * @return     { description_of_the_return_value }
	 */
	template <typename... Ts>
	const _ElemType& get(Ts&&... index) const
	{
		return _data[get_index(std::forward<Ts>(index)...)];
	}

	/**
	 * @brief      { function_description }
	 *
	 * @param      index  The index
	 *
	 * @tparam     Ts     { description }
	 *
	 * @return     { description_of_the_return_value }
	 */
	template <typename... Ts>
	_ElemType& get(Ts&&... index)
	{
		return _data[get_index(std::forward<Ts>(index)...)];
	}

	/**
	 * @brief      { operator_description }
	 *
	 * @param[in]  index  The index
	 *
	 * @return     { description_of_the_return_value }
	 */
	const _ElemType& operator[](const IndexType& index) const
	{
		return _data[get_index(index)];
	}

	/**
	 * @brief      { operator_description }
	 *
	 * @param[in]  index  The index
	 *
	 * @return     { description_of_the_return_value }
	 */
	_ElemType& operator[](const IndexType& index)
	{
		return _data[get_index(index)];
	}

	/**
	 * @brief      { operator_description }
	 *
	 * @param[in]  index  The index
	 *
	 * @return     { description_of_the_return_value }
	 */
	const _ElemType& operator[](std::size_t index) const
	{
		static_assert(_index_dimension == 1, "operator[] with integer index only allowed on 1-tensors");
		return _data[index];
	}

	/**
	 * @brief      { operator_description }
	 *
	 * @param[in]  index  The index
	 *
	 * @return     { description_of_the_return_value }
	 */
	_ElemType& operator[](std::size_t index)
	{
		static_assert(_index_dimension == 1, "operator[] with integer index only allowed on 1-tensors");
		return _data[index];
	}

  /**
   * @brief      { operator_description }
   *
   * @param[in]  rhs   The right hand side
   *
   * @return     { description_of_the_return_value }
   */
  bool operator==(const tensor& rhs) const
  {
      auto rhs_curr = rhs.begin();
      for(auto tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr)
      {
           if(*tcurr != *rhs_curr) return false;
      }
      return true;
  }

  /**
   * @brief      { operator_description }
   *
   * @param[in]  rhs   The right hand side
   *
   * @return     { description_of_the_return_value }
   */
  bool operator!=(const tensor& rhs) const
  {
  	return !(*this == rhs);
  }

  /**
   * @brief      { operator_description }
   *
   * @param[in]  rhs   The right hand side
   */
  void operator=(const tensor& rhs)
  {
      auto rhs_curr = rhs.begin();
      for(auto tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr)
      {
          *tcurr = *rhs_curr;
      }
  }

	/**
	 * @brief      { operator_description }
	 *
	 * @param[in]  rhs   The right hand side
	 *
	 * @return     { description_of_the_return_value }
	 */
	tensor& operator+=(const tensor& rhs)
	{
		typename DataType::const_iterator rhs_curr = rhs.begin();
		for(typename DataType::iterator tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr)
		{
			*tcurr += *rhs_curr;
		}
		return *this;
	}

	/**
	 * @brief      { operator_description }
	 *
	 * @param[in]  rhs   The right hand side
	 *
	 * @return     { description_of_the_return_value }
	 */
	tensor& operator-=(const tensor& rhs)
	{
		typename DataType::const_iterator rhs_curr = rhs.begin();
		for(typename DataType::iterator tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr)
		{
			*tcurr -= *rhs_curr;
		}
		return *this;
	}

	/**
	 * @brief      { operator_description }
	 *
	 * @param[in]  x     { parameter_description }
	 *
	 * @return     { description_of_the_return_value }
	 */
	tensor& operator*=(ElemType x)
	{
		for(auto& a : *this)
		{
			a *= x;
		}
		return *this;
	}

	/**
	 * @brief      Scalar division
	 *
	 * @param[in]  x     { parameter_description }
	 *
	 * @return     { description_of_the_return_value }
	 */
	tensor& operator/=(ElemType x)
	{
		for(auto& a : *this)
		{
			a /= x;
		}
		return *this;
	}

	/**
	 * @brief      Negation operator
	 *
	 * @return     { description_of_the_return_value }
	 */
	tensor operator-() const
	{
		tensor rhs;
		typename DataType::iterator rhs_curr = rhs.begin();
		for(typename DataType::const_iterator tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr)
		{
			*rhs_curr = - (*tcurr);
		}
		return rhs;
	}

	/**
	 * @brief      { function_description }
	 *
	 * @param[in]  rhs   The right hand side
	 *
	 * @return     { description_of_the_return_value }
	 */
	tensor ElementwiseProduct(const tensor& rhs) const{
		tensor ret;
		typename DataType::iterator ret_curr = ret.begin();
		typename DataType::const_iterator rhs_curr = rhs.begin();
		for(typename DataType::const_iterator tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr, ++ret_curr)
		{
			*ret_curr = *tcurr * (*rhs_curr);
		}
		return ret;
	}

	/**
	 * @brief      { function_description }
	 *
	 * @param[in]  rhs   The right hand side
	 *
	 * @return     { description_of_the_return_value }
	 */
	tensor ElementwiseDivision(const tensor& rhs) const{
		tensor ret;
		typename DataType::iterator ret_curr = ret.begin();
		typename DataType::const_iterator rhs_curr = rhs.begin();
		for(typename DataType::const_iterator tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr, ++ret_curr)
		{
			*ret_curr = *tcurr / (*rhs_curr);
		}
		return ret;
	}

	/**
	 * @brief      { function_description }
	 *
	 * @return     { description_of_the_return_value }
	 */
	auto data()
	{
		return _data.data();
	}

	/**
	 * @brief      { function_description }
	 *
	 * @return     { description_of_the_return_value }
	 */
	index_iterator index_begin() { return index_iterator(0); }

	/**
	 * @brief      { function_description }
	 *
	 * @return     { description_of_the_return_value }
	 */
	index_iterator index_end()   { return index_iterator(); }

	/**
	 * @brief      { function_description }
	 *
	 * @return     { description_of_the_return_value }
	 */
	index_iterator index_begin() const { return index_iterator(0); }

	/**
	 * @brief      { function_description }
	 *
	 * @return     { description_of_the_return_value }
	 */
	index_iterator index_end()   const { return index_iterator(); }

	/**
	 * @brief      { function_description }
	 *
	 * @return     { description_of_the_return_value }
	 */
	typename DataType::iterator       begin()       { return _data.begin(); }

	/**
	 * @brief      { function_description }
	 *
	 * @return     { description_of_the_return_value }
	 */
	typename DataType::iterator       end()         { return _data.end();   }

	/**
	 * @brief      { function_description }
	 *
	 * @return     { description_of_the_return_value }
	 */
	typename DataType::const_iterator begin() const { return _data.begin(); }

	/**
	 * @brief      { function_description }
	 *
	 * @return     { description_of_the_return_value }
	 */
	typename DataType::const_iterator end()   const { return _data.end();   }

private:
	/**
	 * @brief      Gets the index.
	 *
	 * @param[in]  args  The arguments
	 *
	 * @tparam     Ts    { description }
	 *
	 * @return     The index.
	 */
	template <typename... Ts>
	std::size_t get_index(Ts... args) const
	{
		std::size_t rval = 0;
		array_util::flatten([&rval](std::size_t ignored, std::size_t i){
			rval *= vector_dimension;
			rval += i;
		}, args...);
		return rval;
	}

	DataType  _data;
};

/**
 * @brief      { function_description }
 *
 * @param[in]  A         { parameter_description }
 *
 * @tparam     ElemType  { description }
 * @tparam     D         { description }
 * @tparam     N         { description }
 *
 * @return     { description_of_the_return_value }
 */
template <typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType,D,N> Sym(const tensor<ElemType,D,N>& A)
{
	tensor<ElemType,D,N> rval;

	std::array<std::size_t, N> sigma;
	for(std::size_t i = 0; i < N; ++i)
	{
		//std::cout << " : " << i << std::endl;
		sigma[i] = i;
	}

	do {
		for(auto curr = rval.index_begin(); curr != rval.index_end(); ++curr)
		{
			std::array<std::size_t, N> index = *curr;
			for(std::size_t i = 0; i < N; ++i)
			{
				index[i] = (*curr)[sigma[i]];
			}
			rval[*curr] += A[index];
		}
	} while(std::next_permutation(sigma.begin(), sigma.end()));

	rval *= 1.0/detail::factorial<N>::value;
	return rval;
}

/**
 * @brief      { function_description }
 *
 * @param[in]  arr   The arr
 *
 * @tparam     N     { description }
 *
 * @return     { description_of_the_return_value }
 */
template <std::size_t N>
int sgn(const std::array<std::size_t,N>& arr)
{
	int rval = 1;
	std::array<bool,N> visited;
	visited.fill(false);

	for(std::size_t i = 1; i < N; ++i)
	{
		if(!visited[i])
		{
			std::size_t len = 0;
			std::size_t next = i;
			do
			{
				++len;
				visited[next] = true;
				next = arr[next];
			} while(!visited[next]);

			rval *= 2*(len % 2) - 1;
		}
	}
	return rval;
}

/**
 * @brief      { function_description }
 *
 * @param[in]  A         { parameter_description }
 *
 * @tparam     ElemType  { description }
 * @tparam     D         { description }
 * @tparam     N         { description }
 *
 * @return     { description_of_the_return_value }
 */
template <typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType,D,N> Alt(const tensor<ElemType,D,N>& A)
{
	tensor<ElemType,D,N> rval;

	std::array<std::size_t, N> sigma;
	for(std::size_t i = 0; i < N; ++i)
	{
		//std::cout << " : " << i << std::endl;
		sigma[i] = i;
	}

	do {
		for(auto curr = rval.index_begin(); curr != rval.index_end(); ++curr)
		{
			std::array<std::size_t, N> index = *curr;
			for(std::size_t i = 0; i < N; ++i)
			{
				index[i] = (*curr)[sigma[i]];
			}
			rval[*curr] += sgn(sigma)*(A[index]);
		}
	} while(std::next_permutation(sigma.begin(), sigma.end()));

	rval *= 1.0/detail::factorial<N>::value;
	return rval;
}

/**
 * @brief      { operator_description }
 *
 * @param[in]  A         { parameter_description }
 * @param[in]  B         { parameter_description }
 *
 * @tparam     ElemType  { description }
 * @tparam     D         { description }
 * @tparam     N         { description }
 * @tparam     M         { description }
 *
 * @return     { description_of_the_return_value }
 */
template <typename ElemType, std::size_t D, std::size_t N, std::size_t M>
tensor<ElemType,D,N+M> operator*(const tensor<ElemType,D,N>& A, const tensor<ElemType,D,M>& B)
{
	tensor<ElemType,D,N+M> rval;

	for(auto pA = A.index_begin(); pA != A.index_end(); ++pA)
	{
		for(auto pB = B.index_begin(); pB != B.index_end(); ++pB)
		{
			rval.get(*pA,*pB) = A[*pA] * B[*pB];
		}
	}

	return rval;
}

/**
 * @brief      Elementwise sum
 *
 * @param[in]  A         { parameter_description }
 * @param[in]  B         { parameter_description }
 *
 * @tparam     ElemType  { description }
 * @tparam     D         { description }
 * @tparam     N         { description }
 *
 * @return     { description_of_the_return_value }
 */
template <typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType,D,N> operator+(const tensor<ElemType,D,N>& A, const tensor<ElemType,D,N>& B)
{
	tensor<ElemType,D,N> rval(A);
	rval += B;
	return rval;
}

/**
 * @brief      { operator_description }
 *
 * @param[in]  A         { parameter_description }
 * @param[in]  B         { parameter_description }
 *
 * @tparam     ElemType  { description }
 * @tparam     D         { description }
 * @tparam     N         { description }
 *
 * @return     { description_of_the_return_value }
 */
template <typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType,D,N> operator-(const tensor<ElemType,D,N>& A, const tensor<ElemType,D,N>& B)
{
	tensor<ElemType,D,N> rval(A);
	rval -= B;
	return rval;
}

/**
 * @brief      { operator_description }
 *
 * @param[in]  A           { parameter_description }
 * @param[in]  x           { parameter_description }
 *
 * @tparam     ScalarType  { description }
 * @tparam     ElemType    { description }
 * @tparam     D           { description }
 * @tparam     N           { description }
 *
 * @return     { description_of_the_return_value }
 */
template <typename ScalarType, typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType,D,N> operator*(const tensor<ElemType,D,N>& A, ScalarType x)
{
	auto rval(A);
	rval *= x;
	return rval;
}

/**
 * @brief      { operator_description }
 *
 * @param[in]  x           { parameter_description }
 * @param[in]  A           { parameter_description }
 *
 * @tparam     ScalarType  { description }
 * @tparam     ElemType    { description }
 * @tparam     D           { description }
 * @tparam     N           { description }
 *
 * @return     { description_of_the_return_value }
 */
template <typename ScalarType, typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType,D,N> operator*(ScalarType x, const tensor<ElemType,D,N>& A)
{
	auto rval(A);
	rval *= x;
	return rval;
}

/**
 * @brief      { operator_description }
 *
 * @param[in]  A           { parameter_description }
 * @param[in]  x           { parameter_description }
 *
 * @tparam     ScalarType  { description }
 * @tparam     ElemType    { description }
 * @tparam     D           { description }
 * @tparam     N           { description }
 *
 * @return     { description_of_the_return_value }
 */
template <typename ScalarType, typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType,D,N> operator/(const tensor<ElemType,D,N>& A, ScalarType x)
{
	auto rval(A);
	rval /= x;
	return rval;
}

/**
 * @brief      { operator_description }
 *
 * @param[in]  A         { parameter_description }
 * @param[in]  B         { parameter_description }
 *
 * @tparam     ElemType  { description }
 * @tparam     D         { description }
 * @tparam     N         { description }
 * @tparam     M         { description }
 *
 * @return     { description_of_the_return_value }
 */
template <typename ElemType, std::size_t D, std::size_t N, std::size_t M>
tensor<ElemType,D,N+M> operator^(const tensor<ElemType,D,N>& A, const tensor<ElemType,D,M>& B)
{
	auto rval = Alt(A*B);
	ElemType num = detail::factorial<N+M>::value;
	ElemType den = detail::factorial<N>::value * detail::factorial<M>::value;
	rval *= num / den;
	//std::cout << " ~ " << N << " " << M << " " << num << " " << den << std::endl;
	return rval;
}

/**
 * @brief      Inner product
 *
 * @param[in]  A         { parameter_description }
 * @param[in]  B         { parameter_description }
 *
 * @tparam     ElemType  { description }
 * @tparam     D         { description }
 * @tparam     N         { description }
 *
 * @return     { description_of_the_return_value }
 */
template <typename ElemType, std::size_t D, std::size_t N>
ElemType dot(const tensor<ElemType,D,N>& A, const tensor<ElemType,D,N>& B)
{
	ElemType rval = 0;
	auto Acurr = A.begin();
	auto Bcurr = B.begin();
	for(; Acurr != A.end(); ++Acurr, ++Bcurr)
	{
		rval += (*Acurr)*(*Bcurr);
	}
	double scale = detail::factorial<N>::value;
	return rval / scale;
}

/**
 * @brief      Inner product
 *
 * @param[in]  A         { parameter_description }
 * @param[in]  B         { parameter_description }
 *
 * @tparam     ElemType  { description }
 * @tparam     D         { description }
 * @tparam     N         { description }
 *
 * @return     { description_of_the_return_value }
 */
template <typename ElemType, std::size_t D, std::size_t N>
ElemType operator|(const tensor<ElemType,D,N>& A, const tensor<ElemType,D,N>& B)
{
	ElemType rval = 0;
	auto Acurr = A.begin();
	auto Bcurr = B.begin();
	for(; Acurr != A.end(); ++Acurr, ++Bcurr)
	{
		rval += (*Acurr)*(*Bcurr);
	}
	double scale = detail::factorial<N>::value;
	return rval / scale;
}

/**
 * @brief      Normalize tensor
 *
 * @param[in]  A         { parameter_description }
 *
 * @tparam     ElemType  { description }
 * @tparam     D         { description }
 * @tparam     N         { description }
 *
 * @return     { description_of_the_return_value }
 */
template <typename ElemType, std::size_t D, std::size_t N>
double norm(const tensor<ElemType,D,N>& A){
	return std::sqrt(A|A);
}

/**
 * @brief      Squared norm
 *
 * @param[in]  A         { parameter_description }
 *
 * @tparam     ElemType  { description }
 * @tparam     D         { description }
 * @tparam     N         { description }
 *
 * @return     { description_of_the_return_value }
 */
template <typename ElemType, std::size_t D, std::size_t N>
double norm_sq(const tensor<ElemType,D,N>& A){
	return A|A;
}

/**
 * @brief      Cross product for vectors only
 *
 * @param[in]  x         { parameter_description }
 * @param[in]  y         { parameter_description }
 *
 * @tparam     ElemType  { description }
 *
 * @return     { description_of_the_return_value }
 */
template <typename ElemType>
tensor<ElemType,3,1> cross(const tensor<ElemType,3,1>& x, const tensor<ElemType,3,1>& y)
{
	tensor<ElemType,3,1> rval;
	rval[0] = x[1]*y[2] - y[1]*x[2];
	rval[1] = x[2]*y[0] - y[2]*x[0];
	rval[2] = x[0]*y[1] - y[0]*x[1];
	return rval;
}
} // end namespace gamer