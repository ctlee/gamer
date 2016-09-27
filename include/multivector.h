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
}

template <typename _ElemType, std::size_t _index_dimension>
class multivector {
public:
	constexpr static std::size_t index_dimension = _index_dimension;
	using ElemType  = _ElemType;
	using IndexType = std::array<std::size_t, index_dimension>;
	using DataType  = std::vector<ElemType>;

	struct index_iterator : public std::iterator<std::bidirectional_iterator_tag, IndexType>
	{
		using super = std::iterator<std::bidirectional_iterator_tag, IndexType>;
		index_iterator(const index_iterator& iter)
			: i(iter.i), n(iter.n)
		{}
		index_iterator(const index_iterator&& iter)
			: i(std::move(iter.i)), n(std::move(iter.n))
		{}
		index_iterator(const IndexType& m)
			: n(m)
		{
			i.fill(0);
		}
		index_iterator(int, const IndexType& m)
			: n(m)
		{
			i.fill(0);
			i[0] = m[0];
		}
		index_iterator(const IndexType& j, const IndexType& m)
			: i(j), n(m)
		{}
		index_iterator& operator++()
		{
			std::size_t k = i.size() - 1;
			++(i[k]);
			while(k > 0)
			{
				if(i[k] == n[k])
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
			std::size_t k = i.size() - 1;
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
			for(std::size_t k = 0; k < n.size(); ++k)
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
		IndexType n;
	};

	template <typename... Ts>
	multivector(Ts... args)
	{
		util::fill_array(_dimensions, args...);
		resize(_dimensions);
	}

	multivector(const multivector& x) : _dimensions(x._dimensions), _data(x._data) {}

	multivector(const multivector&& x) : _dimensions(std::move(x._dimensions)), _data(std::move(x._data)) {}

	multivector(std::initializer_list<std::size_t> init)
	{
		assert(init.size() == index_dimension);
		int i = 0;
		for(auto curr : init)
		{
			_dimensions[i++] = curr;
		}
		std::size_t n = 1;
		std::for_each(_dimensions.begin(), _dimensions.end(), [&n](auto x){n *= x;});
		_data.resize(n);
	}

	multivector(const IndexType& dims)
		: _dimensions(dims)
	{
		assert(dims.size() == index_dimension);
		std::size_t n = 1;
		std::for_each(_dimensions.begin(), _dimensions.end(), [&n](auto x){n *= x;});
		_data.resize(n);
	}

	multivector(const IndexType&& dims)
		: _dimensions(std::move(dims))
	{
		assert(dims.size() == index_dimension);
		std::size_t n = 1;
		std::for_each(_dimensions.begin(), _dimensions.end(), [&n](auto x){n *= x;});
		_data.resize(n);
	}

	const IndexType& dims() const
	{
		return _dimensions;
	}

	void resize(const IndexType& size)
	{
		_data.resize(compute_size(size));
		_dimensions = size;
	}

	const IndexType& size()
	{
		return _dimensions;
	}

	const _ElemType& operator[](const IndexType& index) const
	{
		return _data[get_index(index)];
	}

	_ElemType& operator[](const IndexType& index)
	{
		return _data[get_index(index)];
	}

	const _ElemType& operator[](std::size_t index) const
	{
		static_assert(_index_dimension == 1, "operator[] with integer index only allowed on 1-tensors");
		return _data[index];
	}

	_ElemType& operator[](std::size_t index)
	{
		static_assert(_index_dimension == 1, "operator[] with integer index only allowed on 1-tensors");
		return _data[index];
	}

	multivector& operator+=(const multivector& rhs)
	{
		assert(_dimensions == rhs._dimensions);

		typename DataType::const_iterator rhs_curr = rhs.begin();
		for(typename DataType::iterator tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr)
		{
			*tcurr += *rhs_curr;
		}
		return *this;
	}

	multivector& operator-=(const multivector& rhs)
	{
		assert(_dimensions == rhs._dimensions);

		typename DataType::const_iterator rhs_curr = rhs.begin();
		for(typename DataType::iterator tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr)
		{
			*tcurr -= *rhs_curr;
		}
		return *this;
	}

	multivector& operator*=(ElemType x)
	{
		for(auto& a : *this)
		{
			a *= x;
		}
		return *this;
	}

	multivector& Sym()
	{
		
	}

	index_iterator index_begin() { return index_iterator(_dimensions); }
	index_iterator index_end()   { return index_iterator(0,_dimensions); }

	index_iterator index_begin() const { return index_iterator(_dimensions); }
	index_iterator index_end()   const { return index_iterator(0,_dimensions); }

	typename DataType::iterator       begin()       { return _data.begin(); }
	typename DataType::iterator       end()         { return _data.end();   }

	typename DataType::const_iterator begin() const { return _data.begin(); }
	typename DataType::const_iterator end()   const { return _data.end();   }

private:
	std::size_t compute_size(const IndexType& size)
	{
		std::size_t n = 1;
		std::for_each(size.begin(), size.end(), [&n](auto x){n *= x;});
		return n;
	}

/*
	std::size_t get_index(const IndexType& index) const
	{
		std::size_t rval = 0;
		for(typename IndexType::const_iterator i = index.begin(), j = _dimensions.begin();
			(i != index.end()) || (j != _dimensions.end());
			++i, ++j)
		{
			rval *= *j;
			rval += *i;
		}

		return rval;
	}
*/
	template <typename... Ts>
	std::size_t get_index(Ts... args) const
	{
		std::size_t rval = 0;
		util::flatten([&rval,this](std::size_t i, std::size_t j){
			rval *= this->_dimensions[i];
			rval += j;
		}, args...);
		return rval;
	}


private:
	IndexType _dimensions;
	DataType  _data;
};


template <typename ElemType, std::size_t N>
multivector<ElemType,N> operator-(const multivector<ElemType,N>& A, const multivector<ElemType,N>& B)
{
	multivector<ElemType,N> rval(A);
	rval -= B;
	return std::move(rval);
}

template <typename ElemType, std::size_t N>
multivector<ElemType,N> operator*(const multivector<ElemType,N>& A, ElemType x)
{
	auto rval(A);
	for(auto& a : rval)
	{
		a *= x;
	}
	return std::move(rval);
}
