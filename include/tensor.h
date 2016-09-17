#include <iostream>
#include <array>
#include <vector>
#include <algorithm>
#include <cassert>
#include "util.h"

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

template <typename _ElemType, std::size_t _vector_dimension, std::size_t _index_dimension>
class tensor {
public:
	constexpr static std::size_t index_dimension = _index_dimension;
	constexpr static std::size_t vector_dimension = _vector_dimension;
	constexpr static std::size_t total_dimension = detail::pow<vector_dimension, index_dimension>::value;
	using ElemType  = _ElemType;
	using IndexType = std::array<std::size_t, index_dimension>;
	using DataType  = std::vector<ElemType>;

	struct index_iterator : public std::iterator<std::bidirectional_iterator_tag, IndexType>
	{
		using super = std::iterator<std::bidirectional_iterator_tag, IndexType>;
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
		index_iterator operator++(int) { auto tmp = *this; ++(*this); return tmp; }
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
		index_iterator operator--(int) { auto tmp = *this; --(*this); return tmp; }
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
		util::fill_array(_dimensions, args...);
		resize(_dimensions);
	}
*/
	tensor() { resize(); }

	tensor(const tensor& x) : _data(x._data) {}

	tensor(const tensor&& x) : _data(std::move(x._data)) {}

	void resize()
	{
		_data.resize(total_dimension);
	}
/*
	const IndexType& size()
	{
		return _dimensions;
	}
*/
	const _ElemType& operator[](const IndexType& index) const
	{
		return _data[get_index(index)];
	}

	_ElemType& operator[](const IndexType& index)
	{
		return _data[get_index(index)];
	}

	tensor& operator+=(const tensor& rhs)
	{
		typename DataType::const_iterator rhs_curr = rhs.begin();
		for(typename DataType::iterator tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr)
		{
			*tcurr += *rhs_curr;
		}
		return *this;
	}

	tensor& operator-=(const tensor& rhs)
	{
		typename DataType::const_iterator rhs_curr = rhs.begin();
		for(typename DataType::iterator tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr)
		{
			*tcurr -= *rhs_curr;
		}
		return *this;
	}

	tensor& operator*=(ElemType x)
	{
		for(auto& a : *this)
		{
			a *= x;
		}
		return *this;
	}

	index_iterator index_begin() { return index_iterator(0); }
	index_iterator index_end()   { return index_iterator(); }

	index_iterator index_begin() const { return index_iterator(0); }
	index_iterator index_end()   const { return index_iterator(); }

	typename DataType::iterator       begin()       { return _data.begin(); }
	typename DataType::iterator       end()         { return _data.end();   }

	typename DataType::const_iterator begin() const { return _data.begin(); }
	typename DataType::const_iterator end()   const { return _data.end();   }

private:
	std::size_t get_index(const IndexType& index) const
	{
		std::size_t rval = 0;
		for(auto i : index)
		{
			rval *= vector_dimension;
			rval += i;
		}
		return rval;
	}
private:
	DataType  _data;
};





template <typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType,D,N> Sym(const tensor<ElemType,D,N>& A)
{
	tensor<ElemType,D,N> rval;

	std::array<std::size_t, N> sigma;
	for(std::size_t i = 0; i < N; ++i)
	{
		std::cout << " : " << i << std::endl;
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
	return std::move(rval);
}

template <std::size_t N>
int sgn(const std::array<std::size_t,N>& arr)
{
	int rval = 1;
	std::array<bool,N> visited;
	visited.fill(false);

	for(int i = 1; i < N; ++i)
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

template <typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType,D,N> Alt(const tensor<ElemType,D,N>& A)
{
	tensor<ElemType,D,N> rval;

	std::array<std::size_t, N> sigma;
	for(std::size_t i = 0; i < N; ++i)
	{
		std::cout << " : " << i << std::endl;
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
	return std::move(rval);
}


template <typename ElemType, std::size_t D, std::size_t N, std::size_t M>
tensor<ElemType,D,N+M> operator*(const tensor<ElemType,D,N>& A, const tensor<ElemType,D,M>& B)
{
	std::array<std::size_t,N+M> dimensions;
	tensor<ElemType,D,N+M> rval;

	for(auto pA = A.index_begin(); pA != A.index_end(); ++pA)
	{
		for(auto pB = B.index_begin(); pB != B.index_end(); ++pB)
		{
			std::size_t i = 0;
			for(auto k : *pA)
			{
				dimensions[i++] = k;
			}
			for(auto k : *pB)
			{
				dimensions[i++] = k;
			}
			rval[dimensions] = A[*pA] * B[*pB];
		}
	}

	return std::move(rval);
}

template <typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType,D,N> operator-(const tensor<ElemType,D,N>& A, const tensor<ElemType,D,N>& B)
{
	tensor<ElemType,D,N> rval(A);
	rval -= B;
	return std::move(rval);
}

template <typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType,D,N> operator*(const tensor<ElemType,D,N>& A, ElemType x)
{
	auto rval(A);
	for(auto& a : rval)
	{
		a *= x;
	}
	return std::move(rval);
}

template <typename ElemType, std::size_t D, std::size_t N, std::size_t M>
tensor<ElemType,D,N+M> operator^(const tensor<ElemType,D,N>& A, const tensor<ElemType,D,M>& B)
{
	auto rval = Alt(A*B);
	ElemType num = detail::factorial<N+M>::value;
	ElemType den = detail::factorial<N>::value * detail::factorial<M>::value;
	rval *= num / den;
	std::cout << " ~ " << N << " " << M << " " << num << " " << den << std::endl;
	return std::move(rval);
}

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
	return rval;
}
