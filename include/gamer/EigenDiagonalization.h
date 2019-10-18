
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

#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "gamer.h"


namespace gamer
{

template <typename T = REAL, std::size_t dim = 3>
class EigenDiagonalizeTraits
{
public:
using CovarianceMatrix = std::array<T, (dim * (dim + 1) / 2)>;

private:
/// Eigen 3x3 Matrix of type REAL
using _EigenMatrix = Eigen::Matrix<T, dim, dim>;
/// Eigen 3 Vector of type REAL
using _EigenVector = Eigen::Matrix<T, dim, 1>;

/**
 * @brief     Construct a covariance matrix
 */
static _EigenMatrix construct_covariance_matrix(const CovarianceMatrix &cov)
{
    _EigenMatrix m;

    for (std::size_t i = 0; i < dim; ++i)
    {
        for (std::size_t j = i; j < dim; ++j)
        {
            m(i, j) = static_cast<REAL>(cov[(dim * i) + j - ( (i * (i + 1) ) / 2)]);

            if (i != j)
                m(j, i) = m(i, j);
        }
    }
    return m;
}

public:
/**
 * @brief      Diagonalize a Self Adjoint Matrix
 *
 * @param[in]  mat           Self adjoint matrix to diagonalize
 * @param[out] eigenvalues   Resulting eigenvalues
 * @param[out] eigenvectors  Resulting eigenvectors
 *
 * @return     True on success
 */
static bool diagonalizeSelfAdjointMatrix(const _EigenMatrix &mat,
                                         _EigenVector       &eigenvalues,
                                         _EigenMatrix       &eigenvectors                                                 )
{
    Eigen::SelfAdjointEigenSolver<_EigenMatrix> eigensolver;

    if (dim == 2 || dim == 3)
        eigensolver.computeDirect(mat);
    else
        eigensolver.compute(mat);         // More accurate but slower

    if (eigensolver.info() != Eigen::Success)
    {
        return false;
        // std::stringstream ss;
        // ss << "getEigenvalues has encountered Eigen error " <<
        // eigensolver.info();
        // throw std::runtime_error(ss.str());
    }

    eigenvalues  = eigensolver.eigenvalues();
    eigenvectors = eigensolver.eigenvectors();
    return true;
}

/**
 * @brief      Diagonalize an upper triangular covariance matrix,
 *
 * @param[in]  cov           Upper triangular covariance matrix
 * @param[out] eigenvalues   Resulting eigenvalues
 * @param[out] eigenvectors  Resulting eigenvectors
 *
 * @return     True on success
 */
static bool diagonalizeSelfAdjointCovMatrix(const CovarianceMatrix &cov,
                                            _EigenVector           &eigenvalues,
                                            _EigenMatrix           &eigenvectors)
{
    _EigenMatrix m   = construct_covariance_matrix(cov);
    bool res = diagonalizeSelfAdjointMatrix(m, eigenvalues, eigenvectors);
    return res;
}
}; // end class EigenDiagonalizeTraits
}  // end namespace gamer
