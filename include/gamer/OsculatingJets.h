
/*****************************************************************************
 * This file is part of the GAMer software.
 * Copyright (C) 2019
 * by Christopher Lee, John Moody, Rommie Amaro, J. Andrew McCammon,
 *    and Michael Holst
 *
 * Copyright (c) 2007  INRIA Sophia-Antipolis (France), INRIA Lorraine LORIA.
 * all rights reserved.
 *
 * You can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free SoREALware Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the
 * software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author(s)          : Marc Pouget and Frédéric Cazals
 * Christopher T. Lee : Adapted from CGAL for GAMER
 * ***************************************************************************
 */

#pragma once

#include <Eigen/Eigenvalues>
#include <Eigen/SVD>

#include "gamer/gamer.h"
#include "gamer/EigenDiagonalization.h"

namespace gamer
{

/**
 * @brief      Iterative computation of factorial
 *
 * @param[in]  n     Integer to get factorial of
 *
 * @return     Value of the factorial
 */
inline unsigned int fact(unsigned int n)
{
    unsigned int i, p = 1;
    for (i = 2; i <= n; i++)
        p *= i;
    return p;
}

class EigenSVD
{
public:
static REAL solve(const EigenMatrixN &M, EigenVectorN &B)
{
    Eigen::JacobiSVD<EigenMatrixN> jacobiSvd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
    B = jacobiSvd.solve(EigenVectorN(B));
    return jacobiSvd.singularValues().array().abs().maxCoeff() /
           jacobiSvd.singularValues().array().abs().minCoeff();
}
};

class Monge_via_jet_fitting
{
public:
class MongeForm
{
protected:
Vector m_origin_pt;
Vector m_d1;                 // Maximal principal direction
Vector m_d2;                 // Minimal principal direction
Vector m_n;                  // Normal direction

//coeff = (k1, k2, //ppal curv
//         b0, b1, b2, b3, //third order
//         c0, c1, c2, c3, c4) //fourth order
//     if (degree==1) no coeff needed
std::vector<REAL> m_coefficients;

public:
/**
 * @brief      Default constructor
 */
MongeForm()
{
    m_origin_pt = Vector({0., 0., 0.});
    m_d1 = Vector({0., 0., 0.});
    m_d2 = Vector({0., 0., 0.});
    m_n  = Vector({0., 0., 0.});
    m_coefficients = std::vector<REAL>();
}

/**
 * @brief      Destroys the object.
 */
~MongeForm() {}

//access
const Vector origin() const {
    return m_origin_pt;
}
Vector &origin() {
    return m_origin_pt;
}
const Vector maximal_principal_direction() const {
    return m_d1;
}
Vector &maximal_principal_direction() {
    return m_d1;
}
const Vector minimal_principal_direction() const {
    return m_d2;
}
Vector &minimal_principal_direction() {
    return m_d2;
}
const Vector normal_direction() const {
    return m_n;
}
Vector &normal_direction() {
    return m_n;
}

const std::vector<REAL> coefficients() const {
    return m_coefficients;
}

std::vector<REAL> &coefficients() {
    return m_coefficients;
}

const REAL principal_curvatures(size_t i) const
{
    if ((i == 0 || i == 1) && coefficients().size() >= 2)
    {
        return coefficients()[i];
    }
    else
    {
        throw std::runtime_error("Index out of bounds for principal curvatures.");
    }
}
const REAL third_order_coefficients(size_t i) const
{
    if (i <= 3 && coefficients().size() >= 6)
    {
        return coefficients()[i+2];
    }
    else
    {
        throw std::runtime_error("Index out of bounds for third order coefficients.");
    }
}
const REAL fourth_order_coefficients(size_t i) const
{
    if (i <= 4 && coefficients().size() >= 11)
    {
        return coefficients()[i+6];
    }
    else
    {
        throw std::runtime_error("Index out of bounds for fourth order coefficients.");
    }
}

//if d>=2, number of coeffs = (d+1)(d+2)/2-4.
//we remove cst, linear and the xy coeff which vanish
void set_up(std::size_t degree);

//switch min-max ppal curv/dir wrt a given normal orientation.
// if given_normal.monge_normal < 0 then change the orientation
// if z=g(x,y) in the basis (d1,d2,n) then in the basis
// (d2,d1,-n)
// z=h(x,y)=-g(y,x)
void comply_wrt_given_normal(const Vector &given_normal);

void dump_verbose(std::ostream &out_stream) const;
void dump_4ogl(std::ostream &out_stream, const REAL scale);
};         // end nested class MongeForm

public:
using Aff_transformation = Eigen::Affine3d;
using LAVector = EigenVectorN;
using LAMatrix = EigenMatrixN;


public:
/// Default constructor
Monge_via_jet_fitting();

/**
 * @brief      Function call
 *
 * @param[in]  begin          Iterator to first vertex and origin of fitting
 * @param[in]  end            Iterator just past the end
 * @param[in]  dJet           Order of jet fitting
 * @param      dPrime         Order of differentials to compute
 *
 * @tparam     InputIterator  Typename of iterator
 *
 * @return     Monge form at vertex
 */
template <class InputIterator>
MongeForm operator()(InputIterator begin, InputIterator end,
                     std::size_t dJet, std::size_t dPrime){
    // Precondition verifying that the differential
    if (!((dJet >= 1) && (dPrime >= 1) && (dPrime <= 4) && (dPrime <= dJet))) {
        std::stringstream ss;
        ss << "Cannot compute " << dPrime
           << "-order differential property using "
           << dJet << "-jet.";
        throw std::runtime_error(ss.str());
    }

    this->deg = static_cast<int>(dJet);
    this->deg_monge = static_cast<int>(dPrime);
    // Number of points required to fit dJet-jet
    this->nb_d_jet_coeff = static_cast<int>((dJet+1)*(dJet+2)/2);
    this->nb_input_pts   = static_cast<int>(end - begin);
    // Check that enough points are given to fit d-jet
    if (nb_input_pts < nb_d_jet_coeff)
        throw std::runtime_error("Insufficient points provided to perform jet fitting.");

    // Initialize MongeForm
    MongeForm monge_form;
    monge_form.set_up(dPrime);

    // Assemble the linear system
    LAMatrix M(nb_input_pts, nb_d_jet_coeff);
    LAVector Z(nb_input_pts);

    // Compute
    compute_PCA(begin, end);
    fill_matrix(begin, end, dJet, M, Z);         //with precond

    // std::cout << "M:" << std::endl << M << std::endl;
    // std::cout << "Z:" << std::endl << Z << std::endl;

    solve_linear_system(M, Z);                //correct with precond
    // std::cout << "Solved Z: " << Z << std::endl;
    compute_Monge_basis(Z.data(), monge_form);
    if (dPrime >= 3)
        compute_Monge_coefficients(Z.data(), dPrime, monge_form);
    return monge_form;
}

const REAL condition_number() const {
    return condition_nb;
}

const std::pair<REAL, Vector> pca_basis(size_t i) const
{
    if (i >= 3) throw std::runtime_error("Out of bounds for PCA basis...");
    return m_pca_basis[i];
}

protected:
int deg;
int deg_monge;
int nb_d_jet_coeff;
int nb_input_pts;
REAL preconditionning;
REAL condition_nb;

std::vector< std::pair<REAL, Vector> > m_pca_basis;

//translate_p0 changes the origin of the world to p0 the first point
//  of the input data points
//change_world2fitting (coord of a vector in world) = coord of this
//  vector in fitting. The matrix tranform has as lines the coord of
//  the basis vectors of fitting in the world coord.
//idem for change_fitting2monge
Aff_transformation translate_p0, change_world2fitting,
                   change_fitting2monge;

//eigen val and vect stored in m_pca_basis
// change_world2fitting is computed
template <class InputIterator>
void compute_PCA(InputIterator begin, InputIterator end){
    int n = this->nb_input_pts;
    REAL x, y, z,
         sumX = 0., sumY = 0., sumZ = 0.,
         sumX2 = 0., sumY2 = 0., sumZ2 = 0.,
         sumXY = 0., sumXZ = 0., sumYZ = 0.,
         xx, yy, zz, xy, xz, yz;

    for (; begin != end; begin++)
    {
        // TODO: clean up this messiness...
        Vector lp = (**begin).position;
        // std::cout << lp << std::endl;
        x = lp[0];
        y = lp[1];
        z = lp[2];
        // std::cout << x << " " << y << " " << z << std::endl;
        sumX  += x / n;
        sumY  += y / n;
        sumZ  += z / n;
        sumX2 += x * x / n;
        sumY2 += y * y / n;
        sumZ2 += z * z / n;
        sumXY += x * y / n;
        sumXZ += x * z / n;
        sumYZ += y * z / n;
    }
    xx = sumX2 - sumX * sumX;
    yy = sumY2 - sumY * sumY;
    zz = sumZ2 - sumZ * sumZ;
    xy = sumXY - sumX * sumY;
    xz = sumXZ - sumX * sumZ;
    yz = sumYZ - sumY * sumZ;

    // assemble covariance matrix as a
    // semi-definite matrix.
    // Matrix numbering:
    // 0 1 2
    //   3 4
    //     5
    std::array<REAL, 6> covariance = {{ xx, xy, xz, yy, yz, zz }};
    EigenVector eigenvalues;
    EigenMatrix eigenvectors;

    // solve for eigenvalues and eigenvectors.
    // eigen values are sorted in ascending order,
    // eigen vectors are sorted in accordance.
    EigenDiagonalizeTraits<REAL, 3>::diagonalizeSelfAdjointCovMatrix
        (covariance, eigenvalues, eigenvectors);
    // std::cout << "Eigenvalues: " << eigenvalues << std::endl;
    // std::cout << "Eigenvectors: " << std::endl << eigenvectors << std::endl;

    //store in m_pca_basis
    for (int i = 0; i < 3; i++)
    {
        // std::cout << "Eigenvalue[" << i << "]: " << eigenvalues[2-i] << std::endl;;
        m_pca_basis[i].first = eigenvalues[2-i];
    }

    Vector v1({eigenvectors(6), eigenvectors(7), eigenvectors(8)});
    m_pca_basis[0].second = v1;
    std::cout << v1 << std::endl;
    Vector v2({eigenvectors(3), eigenvectors(4), eigenvectors(5)});
    m_pca_basis[1].second = v2;
    std::cout << v2 << std::endl;
    Vector v3({eigenvectors(0), eigenvectors(1), eigenvectors(2)});
    m_pca_basis[2].second = v3;
    std::cout << v3 << std::endl;

    switch_to_direct_orientation(m_pca_basis[0].second,
                                 m_pca_basis[1].second,
                                 m_pca_basis[2].second);

    std::cout << m_pca_basis[0].second << std::endl
              << m_pca_basis[1].second << std::endl
              << m_pca_basis[2].second << std::endl;


    EigenMatrix tmp;
    tmp << m_pca_basis[0].second[0], m_pca_basis[0].second[1], m_pca_basis[0].second[2],
        m_pca_basis[1].second[0], m_pca_basis[1].second[1], m_pca_basis[1].second[2],
        m_pca_basis[2].second[0], m_pca_basis[2].second[1], m_pca_basis[2].second[2];

    // std::cout << "tmp: " << std::endl << tmp << std::endl;

    //Store the change of basis W->F
    Aff_transformation change_basis(tmp);
    // std::cout << change_basis.matrix() << std::endl;
    this->change_world2fitting = change_basis;
}

//Coordinates of input points are computed in the fitting basis with
//  p0 as origin.
//Preconditionning is computed, M and Z are filled
template <class InputIterator>
void fill_matrix(InputIterator begin, InputIterator end,
                 std::size_t d, LAMatrix &M, LAVector &Z)
{
    //origin of fitting coord system = first input data point
    Vector point0 = (**begin).position;
    // std::cout << "point0: " << point0 << std::endl;
    //transform coordinates of sample points with a
    //translation ($-p$) and multiplication by $ P_{W\rightarrow F}$.
    Vector orig({0., 0., 0.});
    Vector v_point0_orig(orig - point0);
    // std::cout << "v_point0_orig: " << v_point0_orig << std::endl;

    this->translate_p0 = Eigen::Translation3d(v_point0_orig);
    Aff_transformation transf_points = this->change_world2fitting *
                                       this->translate_p0;
    // std::cout << "transf_points: " << std::endl << transf_points.matrix() << std::endl;

    //compute and store transformed points
    std::vector<Vector> pts_in_fitting_basis;
    pts_in_fitting_basis.reserve(this->nb_input_pts);

    for(auto it = begin; it != end; ++it) {
        Vector cur_pt = (**it).position;
        cur_pt = transf_points*EigenMap(cur_pt);
        pts_in_fitting_basis.push_back(cur_pt);
    }
    // CGAL_For_all(begin, end){
    //     Vector cur_pt = transf_points(D2L_converter(*begin));
    //     pts_in_fitting_basis.push_back(cur_pt);
    // }

    //Compute preconditionning
    REAL precond = 0.;
    typename std::vector<Vector>::iterator itb = pts_in_fitting_basis.begin(),
                                           ite = pts_in_fitting_basis.end();

    for(auto it = itb; it != ite; ++it) {
        precond += std::abs( (*it)[0] ) + std::abs( (*it)[1] );
    }
    // std::cout << "Precond: " << precond << std::endl;
    // CGAL_For_all(itb, ite) precond += CGAL::abs(itb->x()) + CGAL::abs(itb->y());
    precond /= 2 * this->nb_input_pts;
    this->preconditionning = precond;
    //fill matrices M and Z
    int line_count = 0;
    REAL x, y;

    for(auto it = itb; it != ite; ++it) {
        // CGAL_For_all(itb, ite) {
        x = (*it)[0];         // x
        y = (*it)[1];         // y
        Z.coeffRef(line_count) = (*it)[2];         //itb->z());
        for (std::size_t k = 0; k <= d; k++)
        {
            for (std::size_t i = 0; i <= k; i++)
            {
                M.coeffRef(line_count, k * (k + 1) / 2 + i) =
                    std::pow( x, static_cast<int>(k - i) )
                    * std::pow( y, static_cast<int>(i) )
                    / (fact( static_cast<unsigned int>(i) ) *
                       fact( static_cast<unsigned int>(k - i) )
                       * std::pow( this->preconditionning, static_cast<int>(k) ) );
            }
        }
        line_count++;
    }
}

//A is computed, solving MA=Z in the ls sense, the solution A is stored
// in Z
//Preconditionning is needed
void solve_linear_system(LAMatrix &M, LAVector &Z);

//Classical differential geometric calculus
//change_fitting2monge is computed
//if deg_monge =1 only 1st order info
//if deg_monge >= 2 2nd order info are computed
void compute_Monge_basis(const REAL* A, MongeForm &monge_form);

//if deg_monge >=3 then 3rd (and 4th) order info are computed
void compute_Monge_coefficients(REAL* A, std::size_t dprime,
                                MongeForm &monge_form);

//for a trihedron (v1,v2,v3) switches v1 to -v1 if det(v1,v2,v3) < 0
void switch_to_direct_orientation(Vector &v1, const Vector &v2,
                                  const Vector &v3);

friend
std::ostream &
operator<<(std::ostream                                     &out_stream,
           const typename Monge_via_jet_fitting::MongeForm &monge)
{
    monge.dump_verbose(out_stream);
    return out_stream;
}
}; // end class Monge_via_jet_fitting
} // end namespace gamer
