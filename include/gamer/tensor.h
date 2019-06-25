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
 * @file    tensor.h
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
#include <casc/casc>

/// Namespace for all things gamer
namespace gamer
{
/// Namespace for tensor array management utilities
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
 * @param[in]  args  Sequence of values
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
/**
 * @brief      Base case
 *
 * @tparam     Fn    Functor
 * @tparam     Ts    Arguments list
 */
template <typename Fn, typename ... Ts>
struct flattenH {};

/**
 * @brief      Pass each element of sequence into functor(N, head)
 *
 * @tparam     Fn    Functor
 * @tparam     T     Typename of head which is passed to the functor
 * @tparam     Ts    Arguments list
 */
template <typename Fn, typename T, typename ... Ts>
struct flattenH<Fn, T, Ts...> {
    template <std::size_t N>
    static void apply(Fn f, T head, Ts... tail)
    {
        f(N, head);
        flattenH<Fn, Ts...>::template apply<N+1>(f, tail ...);
    }
};

/**
 * @brief      Overload for call with no arguments
 *
 * @tparam     Fn    Functor
 */
template <typename Fn>
struct flattenH<Fn> {
    template <std::size_t N>
    static void apply(Fn f) {}
};

/**
 * @brief      Overload for processing array plus Ts...
 */
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
 * @brief      Apply a functor to a sequence
 *
 * @param[in]  f     Functor to apply with prototype void f(type, type)
 * @param[in]  args  The arguments, array, or sequence
 *
 * @tparam     Fn    Typename of the functor
 * @tparam     Ts    Typenames of arguments
 */
template <typename Fn, typename ... Ts>
void flatten(Fn f, Ts... args)
{
    detail::flattenH<Fn, Ts...>::template apply<0>(f, args ...);
}
} // end namespace array_util

/// @cond detail
namespace detail
{
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
    constexpr static std::size_t value = x * pow<x, n-1>::value;
};

template <std::size_t x>
struct pow<x, 0> {
    constexpr static std::size_t value = 1;
};
}
/// @endcond

/**
 * @brief      General multidimensional tensor class
 *
 * The implementation is based on a flattened array. Utilities are provided
 * to assist with indexing and basic tensor mathematics.
 *
 * @tparam     _ElemType          Typename of the tensor elements
 * @tparam     _vector_dimension  Dimension of the vector space
 * @tparam     _tensor_rank       Rank of the tensor
 */
template <typename _ElemType, std::size_t _vector_dimension, std::size_t _tensor_rank>
class tensor
{
    public:
        /// Tensor rank of the tensor
        constexpr static std::size_t tensor_rank = _tensor_rank;
        /// Dimension of the vector space
        constexpr static std::size_t vector_dimension = _vector_dimension;
        /// Total number of tensor components
        constexpr static std::size_t total_dimension = detail::pow<vector_dimension, tensor_rank>::value;
        /// Alias for typename of tensor elements
        using ElemType = _ElemType;
        /// Typename of an index for the tensor
        using IndexType = std::array<std::size_t, tensor_rank>;
        /// Typename of the underlying flat data representation
        using DataType = std::array<ElemType, total_dimension>;

        /**
         * @brief      Iterator over tensor indices
         */
        struct index_iterator : public std::iterator<std::bidirectional_iterator_tag, IndexType>
        {
            using super = std::iterator<std::bidirectional_iterator_tag, IndexType>;
            /**
             * @brief      Copy constructor
             *
             * @param[in]  iter  The iterator to copy from
             */
            index_iterator(const index_iterator &iter)
                : i(iter.i)
            {}

            /**
             * @brief      Move constructor
             *
             * @param[in]  iter  The iterator to move from
             */
            index_iterator(const index_iterator &&iter)
                : i(std::move(iter.i))
            {}

            /**
             * @brief      Constructor initialized at end of indices
             */
            index_iterator()
            {
                i.fill(0);
                i[0] = vector_dimension;
            }

            /**
             * @brief      Constructor initialized at start of indices
             */
            index_iterator(int)
            {
                i.fill(0);
            }

            /**
             * @brief      Prefix incrementation of the index
             *
             * @return     Iterator to the next index
             */
            index_iterator &operator++()
            {
                std::size_t k = tensor_rank - 1;
                ++(i[k]);
                while (k > 0)
                {
                    // incrementing overflows, advance the next index
                    if (i[k] == vector_dimension)
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
             * @brief      Postfix incrementation of the index
             *
             * @return     Iterator to the next index
             */
            index_iterator operator++(int) { auto tmp = *this; ++(*this); return tmp; }

            /**
             * @brief      Prefix decrementation of the index
             *
             * @return     Iterator to the previous index
             */
            index_iterator &operator--()
            {
                std::size_t k = tensor_rank - 1;

                if (i[k] > 0)
                {
                    --(i[k]);
                }
                else
                {
                    std::size_t p = 1;
                    while (k > 0)
                    {
                        --k;
                        ++p;
                        if (i[k] > 0)
                        {
                            --(i[k]);
                            for (std::size_t j = 1; j < p; ++j)
                            {
                                i[k+j] = vector_dimension - 1;
                            }
                            break;
                        }
                    }
                }

                return *this;
            }

            /**
             * @brief      Postfix decrementation of the index
             *
             * @return     Iterator to the previous index
             */
            index_iterator operator--(int) { auto tmp = *this; --(*this); return tmp; }

            /**
             * @brief      Equality operator for indices
             *
             * @param[in]  j     Other index iterator to compare to
             *
             * @return     True if indices are equal
             */
            bool operator==(index_iterator j) const
            {
                for (std::size_t k = 0; k < tensor_rank; ++k)
                {
                    if (i[k] != j.i[k])
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

            /**
             * @brief      Dereference the iterator
             *
             * @return     Index array of indices
             */
            typename super::reference operator*() { return i; }

            /**
             * @brief      Dereference the iterator
             *
             * @return     Index array of indices
             */
            typename super::pointer operator->() { return i; }

            protected:
                /// Array of indices
                IndexType i;
        };

        /**
         * @brief     Default constructor initializes to zero
         */
        tensor() { _data.fill(0); }

        /**
         * @brief      Constructor fill with same value
         *
         * @param[in]  s     Value to fill
         */
        tensor(const ElemType &s) { _data.fill(s); }

        /**
         * @brief      Constructor from flat array
         *
         * @param[in]  s     Row major array to copy from
         */
        tensor(const ElemType (&s)[total_dimension])
        {
            std::copy(std::begin(s), std::end(s), std::begin(_data));
        }

        /**
         * @brief      Copy constructor
         *
         * @param[in]  x     Other tensor to copy
         */
        tensor(const tensor &x) : _data(x._data) {}

        /**
         * @brief      Move constructor
         *
         * @param[in]  x     Other tensor to copy from
         */
        tensor(const tensor &&x) : _data(std::move(x._data)) {}


        /**
         * @brief      Type casting for numerical typed tensors
         *
         * @tparam     NumType    Typename of new tensor data
         * @tparam     <unnamed>  Check that resulting type is numeric
         * @tparam     <unnamed>  Check that current type is numeric
         */
        template <typename NumType,
                  typename = std::enable_if_t<std::is_arithmetic<NumType>::value>,
                  typename = std::enable_if_t<std::is_arithmetic<ElemType>::value> >
        operator tensor<NumType, _vector_dimension, _tensor_rank>() const
        {
            tensor<NumType, _vector_dimension, _tensor_rank> t;
            std::copy(_data.cbegin(), _data.cend(), t.begin());
            return t;
        }

        /**
         * @brief      Print operator overload
         *
         * @param      output  The output stream
         * @param[in]  t       Tensor to print
         *
         * @return     Output stream
         */
        friend std::ostream &operator<<(std::ostream &output, const tensor &t)
        {
            output << "Tensor(";
            bool first = true;
            for (auto curr = t.index_begin(); curr != t.index_end(); ++curr)
            {
                if (first) {output << "{"; first = false;}
                else output << "; {";

                bool first2 = true;
                for (auto x : (*curr))
                {
                    if (first2) {output << x; first2 = false;}
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
        std::string to_string() const
        {
            std::ostringstream output;
            output << *this;
            return output.str();
        }

        /**
         * @brief      Get element by index
         *
         * @param      index  Sequence of indices
         *
         * @tparam     Ts     Typenames of indices
         *
         * @return     Reference to value stored at the index
         */
        template <typename ... Ts>
        const _ElemType &get(Ts && ... index) const
        {
            return _data[get_index(std::forward<Ts>(index)...)];
        }

        /**
         * @brief      Get element by index
         *
         * @param      index  Sequence of indices
         *
         * @tparam     Ts     Typenames of indices
         *
         * @return     Reference to value stored at the index
         */
        template <typename ... Ts>
        _ElemType &get(Ts && ... index)
        {
            return _data[get_index(std::forward<Ts>(index)...)];
        }

        /**
         * @brief      Get element by index
         *
         * @param      index  Sequence of indices
         *
         * @return     Reference to value stored at the index
         */
        const _ElemType &operator[](const IndexType &index) const
        {
            return _data[get_index(index)];
        }

        /**
         * @brief      Get element by index
         *
         * @param      index  Sequence of indices
         *
         * @return     Reference to value stored at the index
         */
        _ElemType &operator[](const IndexType &index)
        {
            return _data[get_index(index)];
        }

        /**
         * @brief      Special overload for data accession of 1-tensors
         *
         * @param      index  Sequence of indices
         *
         * @return     Reference to value stored at the index
         */
        const _ElemType &operator[](std::size_t index) const
        {
            static_assert(_tensor_rank == 1, "operator[] with integer index only allowed on 1-tensors");
            return _data[index];
        }

        /**
         * @brief      Special overload for data accession of 1-tensors
         *
         * @param      index  Sequence of indices
         *
         * @return     Reference to value stored at the index
         */
        _ElemType &operator[](std::size_t index)
        {
            static_assert(_tensor_rank == 1, "operator[] with integer index only allowed on 1-tensors");
            return _data[index];
        }

        /**
         * @brief      Equality comparison of tensors
         *
         * @param[in]  rhs   The right hand side
         *
         * @return     True if all elements equal or False otherwise
         */
        bool operator==(const tensor &rhs) const
        {
            auto rhs_curr = rhs.begin();
            for (auto tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr)
            {
                if (*tcurr != *rhs_curr) return false;
            }
            return true;
        }

        /**
         * @brief      Inquality comparison of tensors
         *
         * @param[in]  rhs   The right hand side
         *
         * @return     False if all elements equal or True otherwise
         */
        bool operator!=(const tensor &rhs) const
        {
            return !(*this == rhs);
        }

        /**
         * @brief      Assignment operator
         *
         * @param[in]  rhs   The right hand side
         */
        void operator=(const tensor &rhs)
        {
            auto rhs_curr = rhs.begin();
            for (auto tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr)
            {
                *tcurr = *rhs_curr;
            }
        }

        /**
         * @brief      Elementwise tensor sum operator
         *
         * @param[in]  rhs   The right hand side
         *
         * @return     Reference to this
         */
        tensor &operator+=(const tensor &rhs)
        {
            typename DataType::const_iterator rhs_curr = rhs.begin();
            for (typename DataType::iterator tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr)
            {
                *tcurr += *rhs_curr;
            }
            return *this;
        }

        /**
         * @brief      Elementwise tensor difference operator
         *
         * @param[in]  rhs   The right hand side
         *
         * @return     Reference to this
         */
        tensor &operator-=(const tensor &rhs)
        {
            typename DataType::const_iterator rhs_curr = rhs.begin();
            for (typename DataType::iterator tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr)
            {
                *tcurr -= *rhs_curr;
            }
            return *this;
        }

        /**
         * @brief      Scalar multiplication
         *
         * @param[in]  x     Scalar to multiply by
         *
         * @return     Reference to this
         */
        tensor &operator*=(ElemType x)
        {
            for (auto &a : *this)
            {
                a *= x;
            }
            return *this;
        }

        /**
         * @brief      Scalar division
         *
         * @param[in]  x     Scalar to divide by
         *
         * @return     Reference to this
         */
        tensor &operator/=(ElemType x)
        {
            for (auto &a : *this)
            {
                a /= x;
            }
            return *this;
        }

        /**
         * @brief      Unary negation operator
         *
         * @return     Tensor with negated values
         */
        tensor operator-() const
        {
            tensor rhs;
            typename DataType::iterator rhs_curr = rhs.begin();
            for (typename DataType::const_iterator tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr)
            {
                *rhs_curr = -(*tcurr);
            }
            return rhs;
        }

        /**
         * @brief      Elementwise product with another tensor
         *
         * @param[in]  rhs   The right hand side
         *
         * @return     Tensor with elementwise product
         */
        tensor ElementwiseProduct(const tensor &rhs) const
        {
            tensor ret;
            typename DataType::iterator       ret_curr = ret.begin();
            typename DataType::const_iterator rhs_curr = rhs.begin();
            for (typename DataType::const_iterator tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr, ++ret_curr)
            {
                *ret_curr = *tcurr * (*rhs_curr);
            }
            return ret;
        }

        /**
         * @brief      Elementwise division with another tensor
         *
         * @param[in]  rhs   The right hand side
         *
         * @return     Tensor with elementwise division
         */
        tensor ElementwiseDivision(const tensor &rhs) const
        {
            tensor ret;
            typename DataType::iterator       ret_curr = ret.begin();
            typename DataType::const_iterator rhs_curr = rhs.begin();
            for (typename DataType::const_iterator tcurr = _data.begin(); tcurr != _data.end(); ++tcurr, ++rhs_curr, ++ret_curr)
            {
                *ret_curr = *tcurr / (*rhs_curr);
            }
            return ret;
        }

        /**
         * @brief      Direct access to the underlying array data
         *
         * @return     Direct access to the underlying array
         */
        auto data()
        {
            return _data.data();
        }

        /**
         * @brief      Iterator to first index
         *
         * @return     Iterator to first index
         */
        index_iterator index_begin() { return index_iterator(0); }

        /**
         * @brief      Iterator to past the end
         *
         * @return     Iterator to past the end
         */
        index_iterator index_end()   { return index_iterator(); }

        /**
         * @brief      Iterator to first index
         *
         * @return     Iterator to first index
         */
        index_iterator index_begin() const { return index_iterator(0); }

        /**
         * @brief      Iterator to past the end
         *
         * @return     Iterator to past the end
         */
        index_iterator index_end()   const { return index_iterator(); }

        /**
         * @brief      Direct iterator to beginning of data
         *
         * @return     Direct iterator to beginning of data
         */
        typename DataType::iterator       begin()       { return _data.begin(); }

        /**
         * @brief      Direct iterator to end of data
         *
         * @return     Direct iterator to end of data
         */
        typename DataType::iterator       end()         { return _data.end();   }

        /**
         * @brief      Direct iterator to beginning of data
         *
         * @return     Direct iterator to beginning of data
         */
        typename DataType::const_iterator begin() const { return _data.begin(); }

        /**
         * @brief      Direct iterator to end of data
         *
         * @return     Direct iterator to end of data
         */
        typename DataType::const_iterator end()   const { return _data.end();   }

    private:
        /**
         * @brief      Gets the index.
         *
         * @param[in]  args  List of indices
         *
         * @tparam     Ts    Typename of indices
         *
         * @return     The flat index to check
         */
        template <typename ... Ts>
        std::size_t get_index(Ts... args) const
        {
            std::size_t rval = 0;
            array_util::flatten([&rval](std::size_t ignored, std::size_t i){
                rval *= vector_dimension;
                rval += i;
            }, args ...);
            return rval;
        }

        DataType _data;
};

/**
 * @brief      Symmetric part of a tensor
 *
 * @param[in]  A         Tensor to symmetrize
 *
 * @tparam     ElemType  Typename of tensor elements
 * @tparam     D         Vector dimension of the tensor
 * @tparam     N         Rank of the tensor
 *
 * @return     Symmetric part of the tensor
 */
template <typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType, D, N> Sym(const tensor<ElemType, D, N> &A)
{
    tensor<ElemType, D, N>     rval;

    std::array<std::size_t, N> sigma;
    for (std::size_t i = 0; i < N; ++i)
    {
        sigma[i] = i;
    }

    do
    {
        for (auto curr = rval.index_begin(); curr != rval.index_end(); ++curr)
        {
            std::array<std::size_t, N> index = *curr;
            for (std::size_t i = 0; i < N; ++i)
            {
                index[i] = (*curr)[sigma[i]];
            }
            rval[*curr] += A[index];
        }
    }
    while (std::next_permutation(sigma.begin(), sigma.end()));

    rval *= 1.0/detail::factorial<N>::value;
    return rval;
}

/**
 * @brief      Parity of a permutation
 *
 * @param[in]  arr   Sequence
 *
 * @tparam     N     Length of sequence
 *
 * @return     1 if even or -1 if odd number of inversions
 */
template <std::size_t N>
int sgn(const std::array<std::size_t, N> &arr)
{
    int                 rval = 1;
    std::array<bool, N> visited;
    visited.fill(false);

    for (std::size_t i = 1; i < N; ++i)
    {
        if (!visited[i])
        {
            std::size_t len  = 0;
            std::size_t next = i;
            do
            {
                ++len;
                visited[next] = true;
                next = arr[next];
            }
            while (!visited[next]);

            rval *= 2*(len % 2) - 1;
        }
    }
    return rval;
}

/**
 * @brief      Alternize a tensor
 *
 * @param[in]  A         Tensor to alternatize
 *
 * @tparam     ElemType  Typename of tensor elements
 * @tparam     D         Vector dimension of the tensor
 * @tparam     N         Rank of the tensor
 *
 * @return    Alternating tensor
 */
template <typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType, D, N> Alt(const tensor<ElemType, D, N> &A)
{
    tensor<ElemType, D, N>     rval;

    std::array<std::size_t, N> sigma;
    for (std::size_t i = 0; i < N; ++i)
    {
        sigma[i] = i;
    }

    do
    {
        for (auto curr = rval.index_begin(); curr != rval.index_end(); ++curr)
        {
            std::array<std::size_t, N> index = *curr;
            for (std::size_t i = 0; i < N; ++i)
            {
                index[i] = (*curr)[sigma[i]];
            }
            rval[*curr] += sgn(sigma)*(A[index]);
        }
    }
    while (std::next_permutation(sigma.begin(), sigma.end()));

    rval *= 1.0/detail::factorial<N>::value;
    return rval;
}

/**
 * @brief      Tensor product
 *
 * @param[in]  A         First tensor
 * @param[in]  B         Second tensor
 *
 * @tparam     ElemType  Typename of elements
 * @tparam     D         Vector dimension
 * @tparam     N         Rank of first tensor
 * @tparam     M         Rank of second tensor
 *
 * @return     Tensor with dimension D and rank N+M
 */
template <typename ElemType, std::size_t D, std::size_t N, std::size_t M>
tensor<ElemType, D, N+M> operator*(const tensor<ElemType, D, N> &A, const tensor<ElemType, D, M> &B)
{
    tensor<ElemType, D, N+M> rval;

    for (auto pA = A.index_begin(); pA != A.index_end(); ++pA)
    {
        for (auto pB = B.index_begin(); pB != B.index_end(); ++pB)
        {
            rval.get(*pA, *pB) = A[*pA] * B[*pB];
        }
    }

    return rval;
}

/**
 * @brief      Elementwise sum
 *
 * @param[in]  A         First tensor
 * @param[in]  B         Second tensor
 *
 * @tparam     ElemType  Typename of elements
 * @tparam     D         Vector dimension
 * @tparam     N         Rank of tensors
 *
 * @return     Elementwise sum of tensors
 */
template <typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType, D, N> operator+(const tensor<ElemType, D, N> &A, const tensor<ElemType, D, N> &B)
{
    tensor<ElemType, D, N> rval(A);
    rval += B;
    return rval;
}

/**
 * @brief      Elementwise difference
 *
 * @param[in]  A         First tensor
 * @param[in]  B         Second tensor
 *
 * @tparam     ElemType  Typename of elements
 * @tparam     D         Vector dimension
 * @tparam     N         Rank of tensors
 *
 * @return     Elementwise difference of tensors
 */
template <typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType, D, N> operator-(const tensor<ElemType, D, N> &A, const tensor<ElemType, D, N> &B)
{
    tensor<ElemType, D, N> rval(A);
    rval -= B;
    return rval;
}

/**
 * @brief      Scalar tensor product
 *
 * @param[in]  A           Tensor to multiply
 * @param[in]  x           Scalar to multiply by
 *
 * @tparam     ScalarType  Typename of the scalar
 * @tparam     ElemType    Typename of tensor elements
 * @tparam     D           Vector dimension
 * @tparam     N           Rank
 *
 * @return     Scalar tensor product
 */
template <typename ScalarType, typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType, D, N> operator*(const tensor<ElemType, D, N> &A, ScalarType x)
{
    auto rval(A);
    rval *= x;
    return rval;
}

/**
 * @brief      Scalar tensor product
 *
 * @param[in]  A           Tensor to multiply
 * @param[in]  x           Scalar to multiply by
 *
 * @tparam     ScalarType  Typename of the scalar
 * @tparam     ElemType    Typename of tensor elements
 * @tparam     D           Vector dimension
 * @tparam     N           Rank
 *
 * @return     Scalar tensor product
 */
template <typename ScalarType, typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType, D, N> operator*(ScalarType x, const tensor<ElemType, D, N> &A)
{
    auto rval(A);
    rval *= x;
    return rval;
}

/**
 * @brief      Scalar division of tensor
 *
 * @param[in]  A           Tensor to divide
 * @param[in]  x           Scalar denominator
 *
 * @tparam     ScalarType  Typename of the scalar
 * @tparam     ElemType    Typename of tensor elements
 * @tparam     D           Vector dimension
 * @tparam     N           Rank
 *
 * @return     Divided tensor
 */
template <typename ScalarType, typename ElemType, std::size_t D, std::size_t N>
tensor<ElemType, D, N> operator/(const tensor<ElemType, D, N> &A, ScalarType x)
{
    auto rval(A);
    rval /= x;
    return rval;
}

/**
 * @brief      Wedge product of two tensors
 *
 * @param[in]  A         First tensor
 * @param[in]  B         Second tensor
 *
 * @tparam     ElemType  Typename of elements
 * @tparam     D         Vector dimension
 * @tparam     N         Rank of first tensor
 * @tparam     M         Rank of second tensor
 *
 * @return     Tensor with dimension D and rank N+M
 */
template <typename ElemType, std::size_t D, std::size_t N, std::size_t M>
tensor<ElemType, D, N+M> operator^(const tensor<ElemType, D, N> &A, const tensor<ElemType, D, M> &B)
{
    auto     rval = Alt(A*B);
    ElemType num  = detail::factorial<N+M>::value;
    ElemType den  = detail::factorial<N>::value * detail::factorial<M>::value;
    rval *= num / den;
    return rval;
}

/**
 * @brief      Inner product
 *
 * @param[in]  A         First tensor
 * @param[in]  B         Second tensor
 *
 * @tparam     ElemType  Typename of elements
 * @tparam     D         Vector dimension
 * @tparam     N         Rank of first tensor
 *
 * @return     Tensor inner product
 */
template <typename ElemType, std::size_t D, std::size_t N>
ElemType dot(const tensor<ElemType, D, N> &A, const tensor<ElemType, D, N> &B)
{
    return A|B;
}

/**
 * @brief      Inner product
 *
 * @param[in]  A         First tensor
 * @param[in]  B         Second tensor
 *
 * @tparam     ElemType  Typename of elements
 * @tparam     D         Vector dimension
 * @tparam     N         Rank of first tensor
 *
 * @return     Tensor inner product
 */
template <typename ElemType, std::size_t D, std::size_t N>
ElemType operator|(const tensor<ElemType, D, N> &A, const tensor<ElemType, D, N> &B)
{
    ElemType rval  = 0;
    auto     Acurr = A.begin();
    auto     Bcurr = B.begin();
    for (; Acurr != A.end(); ++Acurr, ++Bcurr)
    {
        rval += (*Acurr)*(*Bcurr);
    }
    double scale = detail::factorial<N>::value;
    return rval / scale;
}

/**
 * @brief      Cross product for vectors only
 *
 * @param[in]  x         First vector
 * @param[in]  y         Second vector
 *
 * @tparam     ElemType  Typename of data
 *
 * @return     Cross product
 */
template <typename ElemType>
tensor<ElemType, 3, 1> cross(const tensor<ElemType, 3, 1> &x, const tensor<ElemType, 3, 1> &y)
{
    tensor<ElemType, 3, 1> rval;
    rval[0] = x[1]*y[2] - y[1]*x[2];
    rval[1] = x[2]*y[0] - y[2]*x[0];
    rval[2] = x[0]*y[1] - y[0]*x[1];
    return rval;
}
} // end namespace gamer
