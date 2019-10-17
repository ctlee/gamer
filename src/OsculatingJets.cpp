// Code adapted from CGAL
// License reproduced as follows

// Copyright (c) 2007  INRIA Sophia-Antipolis (France), INRIA Lorraine LORIA.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free SoREALware Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the soREALware.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Marc Pouget and Frédéric Cazals

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>

#include "gamer/gamer.h"
#include "gamer/OsculatingJets.h"

namespace gamer {

//-------------------------------------------------------------------
// Implementation
//------------------------------------------------------------------
void Monge_via_jet_fitting::MongeForm::set_up(std::size_t degree)
{
    if (degree >= 2)
        std::fill_n(back_inserter(m_coefficients),
                    (degree + 1) * (degree + 2) / 2 - 4, 0.);
}

void Monge_via_jet_fitting::MongeForm::comply_wrt_given_normal(const Vector &given_normal)
{
    if (dot( given_normal, this->normal_direction() ) < 0)
    {
        normal_direction() = -normal_direction();
        std::swap( maximal_principal_direction(), minimal_principal_direction() );
        if (coefficients().size() >= 2)
            std::swap(coefficients()[0], coefficients()[1]);
        if (coefficients().size() >= 6)
        {
            std::swap(coefficients()[2], coefficients()[5]);
            std::swap(coefficients()[3], coefficients()[4]);
        }
        if (coefficients().size() >= 11)
        {
            std::swap(coefficients()[6], coefficients()[10]);
            std::swap(coefficients()[7], coefficients()[9]);
        }
        typename std::vector<REAL>::iterator itb = coefficients().begin(),
                                             ite = coefficients().end();
        for (; itb != ite; itb++)
        {
            *itb = -(*itb);
        }
    }
}

void Monge_via_jet_fitting::MongeForm::dump_verbose(std::ostream &out_stream) const
{
    out_stream << "origin : " << origin() << std::endl
               << "n : " << normal_direction() << std::endl;
    if (coefficients().size() >= 2)
        out_stream << "d1 : " << maximal_principal_direction() << std::endl
                   << "d2 : " << minimal_principal_direction() << std::endl
                   << "k1 : " << coefficients()[0] << std::endl
                   << "k2 : " << coefficients()[1] << std::endl;
    if (coefficients().size() >= 6)
        out_stream << "b0 : " << coefficients()[2] << std::endl
                   << "b1 : " << coefficients()[3] << std::endl
                   << "b2 : " << coefficients()[4] << std::endl
                   << "b3 : " << coefficients()[5] << std::endl;
    if (coefficients().size() >= 11)
        out_stream << "c0 : " << coefficients()[6] << std::endl
                   << "c1 : " << coefficients()[7] << std::endl
                   << "c2 : " << coefficients()[8] << std::endl
                   << "c3 : " << coefficients()[9] << std::endl
                   << "c4 : " << coefficients()[10] << std::endl
                   << std::endl;
}

void Monge_via_jet_fitting::MongeForm::
dump_4ogl(std::ostream &out_stream, const REAL scale)
{
    // CGAL_precondition( coefficients().size() >= 2 );
    if (coefficients().size() < 2)
        throw std::runtime_error("Insufficient coefficients");
    out_stream << origin()  << " "
               << maximal_principal_direction() * scale << " "
               << minimal_principal_direction() * scale << " "
               << coefficients()[0] << " "
               << coefficients()[1] << " "
               << std::endl;
}




Monge_via_jet_fitting::Monge_via_jet_fitting()
{
    m_pca_basis = std::vector< std::pair<REAL, Vector> >(3);
}

void Monge_via_jet_fitting::solve_linear_system(LAMatrix &M, LAVector &Z)
{
    condition_nb = EigenSVD::solve(M, Z);
    for (int k = 0; k <= this->deg; k++)
        for (int i = 0; i <= k; i++)
            // Z[k*(k+1)/2+i] /= std::pow(this->preconditionning,k);
            Z( k * (k + 1) / 2 + i) = Z(k * (k + 1) / 2 + i) / std::pow(this->preconditionning, k);
}

void Monge_via_jet_fitting::compute_Monge_basis(const REAL* A, MongeForm &monge_form)
{
    // std::cout << A[0] << std::endl;
    // std::cout << A[1] << std::endl;
    // std::cout << A[2] << std::endl;
    // std::cout << A[3] << std::endl;
    // std::cout << A[4] << std::endl;
    // std::cout << A[5] << std::endl;

    // only 1st order info.
    if (this->deg_monge == 1)
    {
        Vector orig_monge({0., 0., A[0]});
        Vector normal({-A[1], -A[2], 1.});
        REAL norm2 = normal | normal;
        normal /= std::sqrt(norm2);
        monge_form.origin().mapEigen() = ( this->translate_p0.inverse() * this->change_world2fitting.inverse() ) * orig_monge.mapEigen();
        monge_form.normal_direction().mapEigen() = this->change_world2fitting.inverse() * normal.mapEigen();
    }
    // else (deg_monge >= 2) then 2nd order info are computed
    else
    {
        //bi-index to uni-index conversion : A(i,j)=A[(i+j)(i+j+1)/2+j]
        Vector orig_monge({0., 0., A[0]});
        //normal = Xu crossprod Xv
        Vector Xu({1., 0., A[1]}), Xv({0., 1., A[2]}), normal({-A[1], -A[2], 1.});
        REAL norm2 = normal | normal;
        normal /= std::sqrt(norm2);

        //Surface in fitting_basis : X(u,v)=(u,v,J_A(u,v))
        //in the basis Xu=(1,0,A[1]), Xv=(0,1,A[2]), Weingarten=-I^{-1}II
        //first fond form I=(e,f,f,g)
        //                 =(Xu.Xu, Xu.Xv, Xu.Xv, Xv.Xv)
        //second fond form II=(l,m,m,n)/norm2^(1/2)
        //                   =(n.Xuu, n.Xuv, n.Xuv, n.Xvv)
        //ppal curv are the opposite of the eigenvalues of Weingarten or the
        //  eigenvalues of weingarten = -Weingarten = I^{-1}II

        using Matrix = Eigen::Matrix<REAL,2,2>;
        // typedef typename CGAL::Linear_algebraCd<REAL>::Matrix Matrix;

        REAL e = 1 + A[1] * A[1], f = A[1] * A[2], g = 1 + A[2] * A[2],
             l = A[3], m = A[4], n = A[5];
        Matrix weingarten;
        weingarten(0, 0) = (g * l - f * m) / (std::sqrt(norm2) * norm2);
        weingarten(0, 1) = (g * m - f * n) / (std::sqrt(norm2) * norm2);
        weingarten(1, 0) = (e * m - f * l) / (std::sqrt(norm2) * norm2);
        weingarten(1, 1) = (e * n - f * m) / (std::sqrt(norm2) * norm2);
        // Y, Z are normalized GramSchmidt of Xu, Xv
        // Xu->Y=Xu/||Xu||;
        // Xv->Z=Xv-(Xu.Xv)Xu/||Xu||^2;
        // Z-> Z/||Z||
        Vector Y, Z;
        REAL normXu = std::sqrt( Xu | Xu );
        Y = Xu / normXu;
        REAL XudotXv = Xu | Xv;
        Z = Xv - XudotXv * Xu / (normXu * normXu);
        REAL normZ = std::sqrt( Z | Z );
        Z /= normZ;
        Matrix change_XuXv2YZ;
        change_XuXv2YZ(0, 0) = 1 / normXu;
        change_XuXv2YZ(0, 1) = -XudotXv / (normXu * normXu * normZ);
        change_XuXv2YZ(1, 0) = 0;
        change_XuXv2YZ(1, 1) = 1 / normZ;

        // REAL     det = change_XuXv2YZ.determinant();
        //in the new orthonormal basis (Y,Z) of the tangent plane :
        // weingarten = inv *(1/det) * weingarten * change_XuXv2YZ;
        weingarten = change_XuXv2YZ.inverse() * weingarten * change_XuXv2YZ;

        // diagonalization of weingarten
        std::array<REAL, 3> W = {{ weingarten(0, 0), weingarten(1, 0), weingarten(1, 1) }};

        // std::array<REAL, 2> eval = {{ 0., 0. }};
        // std::array<REAL, 4> evec = {{ 0., 0., 0., 0. }};
        Eigen::Matrix<REAL,2,1> eval;
        Eigen::Matrix<REAL,2,2> evec;

        //eval in increasing order
        EigenDiagonalizeTraits<REAL, 2>::diagonalizeSelfAdjointCovMatrix(W, eval, evec);
        // CGAL::Default_diagonalize_traits<REAL, 2>::diagonalize_selfadjoint_covariance_matrix
        // (W, eval, evec);

        // std::cout << "Weingarten: " << casc::to_string(W) << std::endl;
        // std::cout << "eigenvalues: " << std::endl << eval << std::endl;
        // std::cout << "eigenvectors: " << std::endl << evec << std::endl;

        // std::cout << Y << "; " << Z << std::endl;

        Vector d_max = evec(2) * Y + evec(3) * Z,
               d_min = evec(0) * Y + evec(1) * Z;

        switch_to_direct_orientation(d_max, d_min, normal);
        // Aff_transformation change_basis (d_max[0], d_max[1], d_max[2],
        //                                  d_min[0], d_min[1], d_min[2],
        //                                  normal[0], normal[1], normal[2]);
        EigenMatrix tmp;
        tmp << d_max[0], d_max[1], d_max[2],
               d_min[0], d_min[1], d_min[2],
               normal[0], normal[1], normal[2];
        Aff_transformation change_basis (tmp);

        this->change_fitting2monge = change_basis;

        //store the monge basis origin and vectors with their world coord
        //store ppal curv
        monge_form.origin() = ( this->translate_p0.inverse() *
                                          this->change_world2fitting.inverse() ) * orig_monge.mapEigen();
        monge_form.maximal_principal_direction().mapEigen() = this->change_world2fitting.inverse() * d_max.mapEigen();
        monge_form.minimal_principal_direction().mapEigen() = this->change_world2fitting.inverse() * d_min.mapEigen();
        monge_form.normal_direction().mapEigen() = this->change_world2fitting.inverse() * normal.mapEigen();
        monge_form.coefficients()[0]  = eval[1];
        monge_form.coefficients()[1]  = eval[0];
    }
    //end else
}

void Monge_via_jet_fitting::compute_Monge_coefficients(REAL* A, std::size_t dprime, MongeForm &monge_form)
{
    //One has the equation w=J_A(u,v) of the fitted surface S
    // in the fitting_basis
    //Substituing (u,v,w)=change_fitting2monge^{-1}(x,y,z)
    //One has the equation f(x,y,z)=0 on this surface S in the monge
    //  basis
    //The monge form of the surface at the origin is the bivariate fct
    //   g(x,y) s.t. f(x,y,g(x,y))=0
    //voir les calculs Maple dans monge.mws
    //Notations are f123= d^3f/dxdydz
    //              g(x,y)=sum (gij x^i y^j/ i!j!) with
    //              g00=g10=g01=g11=0, g20=kmax, g02=kmin
    //
    //g(x,y)= 1/2*(k1x^2 +k2y^2)
    //       +1/6*(b0x^3 +3b1x^2y +3b2xy^2 +b3y^3)
    //       +1/24*(c0x^4 +4c1x^3y +6c2x^2y^2 +4c3xy^3 +c4y^4)
    //       +...
    // p stores change_fitting2monge^{-1}=change_fitting2monge^{T}
    REAL p[3][3];
    p[0][0] = this->change_fitting2monge(0, 0);
    p[1][0] = this->change_fitting2monge(0, 1);
    p[2][0] = this->change_fitting2monge(0, 2);
    p[0][1] = this->change_fitting2monge(1, 0);
    p[1][1] = this->change_fitting2monge(1, 1);
    p[2][1] = this->change_fitting2monge(1, 2);
    p[0][2] = this->change_fitting2monge(2, 0);
    p[1][2] = this->change_fitting2monge(2, 1);
    p[2][2] = this->change_fitting2monge(2, 2);

    // formula are designed for w=sum( Aij ui vj), but we have J_A = sum(
    // Aij/i!j! ui vj)
    for (int k = 0; k <= this->deg; k++)
        for (int i = 0; i <= k; i++)
            A[k * (k + 1) / 2 + i] /= fact(k - i) * fact(i);
    //this is A(k-i;i)

/*   //debug */
/*   std::cout << "coeff of A" << std::endl */
/*      << A[0] << " "<< A[1] << " "<< A[2] << std::endl */
/*      << A[3] << " "<< A[4] << " "<< A[5] << std::endl */
/*      << A[6] << " "<< A[7] << " "<< A[8] << " "<< A[9]<< std::endl */
/*      << A[10] << " "<< A[11] << " "<< A[12] << " "<< A[13]<< " " << A[14] <<
   std::endl; */



    //     note f1 = f2 = f12 = 0
    //     REAL f1 = A[1] * p[0][0] + A[2] * p[1][0] - p[2][0];
    //     REAL f2 = A[2] * p[1][1] + A[1] * p[0][1] - p[2][1];
    //     REAL f12 =
    //     2 * A[3] * p[0][0] * p[0][1]
    //     + 2 * A[5] * p[1][0] * p[1][1]
    //     + A[4] * p[0][1] * p[1][0]
    //     + A[4] * p[0][0] * p[1][1];
    //         -f11 / f3 = kmax
    //         -f22 / f3 = kmin

    REAL f3  = A[1] * p[0][2] + A[2] * p[1][2] - p[2][2];
    REAL f11 =
        2 * A[4] * p[0][0] * p[1][0]
        + 2 * A[5] * p[1][0] * p[1][0]
        + 2 * A[3] * p[0][0] * p[0][0];
    REAL f13 =
        A[4] * p[0][0] * p[1][2]
        + A[4] * p[0][2] * p[1][0]
        + 2 * A[5] * p[1][0] * p[1][2]
        + 2 * A[3] * p[0][0] * p[0][2];
    REAL f22 =
        2 * A[4] * p[0][1] * p[1][1]
        + 2 * A[5] * p[1][1] * p[1][1]
        + 2 * A[3] * p[0][1] * p[0][1];
    REAL f23 =
        A[4] * p[0][1] * p[1][2]
        + 2 * A[5] * p[1][1] * p[1][2]
        + A[4] * p[0][2] * p[1][1]
        + 2 * A[3] * p[0][1] * p[0][2];
    REAL f33 =
        2 * A[5] * p[1][2] * p[1][2]
        + 2 * A[3] * p[0][2] * p[0][2]
        + 2 * A[4] * p[0][2] * p[1][2];
    REAL f111 =
        6 * A[8] * p[0][0] * p[1][0] * p[1][0]
        + 6 * A[7] * p[0][0] * p[0][0] * p[1][0]
        + 6 * A[6] * p[0][0] * p[0][0] * p[0][0]
        + 6 * A[9] * p[1][0] * p[1][0] * p[1][0];
    REAL f222 =
        6 * A[7] * p[0][1] * p[0][1] * p[1][1]
        + 6 * A[8] * p[0][1] * p[1][1] * p[1][1]
        + 6 * A[9] * p[1][1] * p[1][1] * p[1][1]
        + 6 * A[6] * p[0][1] * p[0][1] * p[0][1];
    REAL f112 =
        2 * A[7] * p[0][0] * p[0][0] * p[1][1]
        + 6 * A[6] * p[0][0] * p[0][0] * p[0][1]
        + 2 * A[8] * p[0][1] * p[1][0] * p[1][0]
        + 4 * A[8] * p[0][0] * p[1][0] * p[1][1]
        + 6 * A[9] * p[1][0] * p[1][0] * p[1][1]
        + 4 * A[7] * p[0][0] * p[0][1] * p[1][0];
    REAL f122 =
        4 * A[8] * p[0][1] * p[1][0] * p[1][1]
        + 2 * A[8] * p[0][0] * p[1][1] * p[1][1]
        + 6 * A[6] * p[0][0] * p[0][1] * p[0][1]
        + 2 * A[7] * p[0][1] * p[0][1] * p[1][0]
        + 4 * A[7] * p[0][0] * p[0][1] * p[1][1]
        + 6 * A[9] * p[1][0] * p[1][1] * p[1][1];
    REAL f113 =
        6 * A[6] * p[0][0] * p[0][0] * p[0][2]
        + 6 * A[9] * p[1][0] * p[1][0] * p[1][2]
        + 2 * A[7] * p[0][0] * p[0][0] * p[1][2]
        + 2 * A[8] * p[0][2] * p[1][0] * p[1][0]
        + 4 * A[7] * p[0][0] * p[0][2] * p[1][0]
        + 4 * A[8] * p[0][0] * p[1][0] * p[1][2];
    REAL f223 =
        2 * A[8] * p[0][2] * p[1][1] * p[1][1]
        + 6 * A[6] * p[0][1] * p[0][1] * p[0][2]
        + 6 * A[9] * p[1][1] * p[1][1] * p[1][2]
        + 2 * A[7] * p[0][1] * p[0][1] * p[1][2]
        + 4 * A[7] * p[0][1] * p[0][2] * p[1][1]
        + 4 * A[8] * p[0][1] * p[1][1] * p[1][2];
    REAL f123 =
        2 * A[8] * p[0][2] * p[1][0] * p[1][1]
        + 2 * A[7] * p[0][0] * p[0][1] * p[1][2]
        + 2 * A[7] * p[0][0] * p[0][2] * p[1][1]
        + 6 * A[9] * p[1][0] * p[1][1] * p[1][2]
        + 2 * A[7] * p[0][1] * p[0][2] * p[1][0]
        + 6 * A[6] * p[0][0] * p[0][1] * p[0][2]
        + 2 * A[8] * p[0][0] * p[1][1] * p[1][2]
        + 2 * A[8] * p[0][1] * p[1][0] * p[1][2];

    REAL b0 = 1 / (f3 * f3) * (-f111 * f3 + 3 * f13 * f11);
    REAL b1 = 1 / (f3 * f3) * (-f112 * f3 + f23 * f11);
    REAL b2 = 1 / (f3 * f3) * (-f122 * f3 + f13 * f22);
    REAL b3 = -1 / (f3 * f3) * (f222 * f3 - 3 * f23 * f22);

    monge_form.coefficients()[2] = b0;
    monge_form.coefficients()[3] = b1;
    monge_form.coefficients()[4] = b2;
    monge_form.coefficients()[5] = b3;

    if (dprime == 4)
    {
        REAL f1111 =
            24 * A[13] * p[0][0] * p[1][0] * p[1][0] * p[1][0]
            + 24 * A[12] * p[0][0] * p[0][0] * p[1][0] * p[1][0]
            + 24 * A[11] * p[0][0] * p[0][0] * p[0][0] * p[1][0]
            + 24 * A[14] * p[1][0] * p[1][0] * p[1][0] * p[1][0]
            + 24 * A[10] * p[0][0] * p[0][0] * p[0][0] * p[0][0];
        REAL f1112 =
            6 * A[13] * p[0][1] * p[1][0] * p[1][0] * p[1][0]
            + 18 * A[13] * p[0][0] * p[1][0] * p[1][0] * p[1][1]
            + 24 * A[10] * p[0][0] * p[0][0] * p[0][0] * p[0][1]
            + 12 * A[12] * p[0][0] * p[0][1] * p[1][0] * p[1][0]
            + 18 * A[11] * p[0][0] * p[0][0] * p[0][1] * p[1][0]
            + 24 * A[14] * p[1][0] * p[1][0] * p[1][0] * p[1][1]
            + 6 * A[11] * p[0][0] * p[0][0] * p[0][0] * p[1][1]
            + 12 * A[12] * p[0][0] * p[0][0] * p[1][0] * p[1][1];
        REAL f1122 =
            12 * A[11] * p[0][0] * p[0][0] * p[0][1] * p[1][1]
            + 12 * A[13] * p[0][0] * p[1][0] * p[1][1] * p[1][1]
            + 12 * A[13] * p[0][1] * p[1][0] * p[1][0] * p[1][1]
            + 16 * A[12] * p[0][0] * p[0][1] * p[1][0] * p[1][1]
            + 12 * A[11] * p[0][0] * p[0][1] * p[0][1] * p[1][0]
            + 24 * A[10] * p[0][0] * p[0][0] * p[0][1] * p[0][1]
            + 4 * A[12] * p[0][1] * p[0][1] * p[1][0] * p[1][0]
            + 4 * A[12] * p[0][0] * p[0][0] * p[1][1] * p[1][1]
            + 24 * A[14] * p[1][0] * p[1][0] * p[1][1] * p[1][1];
        REAL f1222 =
            6 * A[13] * p[0][0] * p[1][1] * p[1][1] * p[1][1]
            + 24 * A[10] * p[0][0] * p[0][1] * p[0][1] * p[0][1]
            + 24 * A[14] * p[1][0] * p[1][1] * p[1][1] * p[1][1]
            + 6 * A[11] * p[0][1] * p[0][1] * p[0][1] * p[1][0]
            + 18 * A[11] * p[0][0] * p[0][1] * p[0][1] * p[1][1]
            + 12 * A[12] * p[0][0] * p[0][1] * p[1][1] * p[1][1]
            + 12 * A[12] * p[0][1] * p[0][1] * p[1][0] * p[1][1]
            + 18 * A[13] * p[0][1] * p[1][0] * p[1][1] * p[1][1];
        REAL f2222 =
            24 * A[13] * p[0][1] * p[1][1] * p[1][1] * p[1][1]
            + 24 * A[11] * p[0][1] * p[0][1] * p[0][1] * p[1][1]
            + 24 * A[12] * p[0][1] * p[0][1] * p[1][1] * p[1][1]
            + 24 * A[10] * p[0][1] * p[0][1] * p[0][1] * p[0][1]
            + 24 * A[14] * p[1][1] * p[1][1] * p[1][1] * p[1][1];

        REAL c0 =
            -1 / (f3 * f3 * f3) * (f1111 * (f3 * f3) - 4 * f13 * f3 * f111 + 12 * f13 * f13 * f11 - 6 * f113 * f3 * f11 + 3 * f33 * f11 * f11);
        REAL c1 =
            1 / (f3 * f3 * f3) * (f23 * f3 * f111 + 3 * f3 * f123 * f11 + 3 * f13 * f3 * f112 - f1112 * (f3 * f3) - 6 * f13 * f23 * f11);
        REAL c2 =
            1 / (f3 * f3 * f3) * ( -f33 * f22 * f11 + f113 * f3 * f22 + 2 * f13 * f3 * f122 - 2 * f13 * f13 * f22 + f223 * f3 * f11 + 2 * f23 * f3 * f112 - 2 * f23 * f23 * f11 - f1122 * (f3 * f3) );
        REAL c3 =
            1 / (f3 * f3 * f3) * (-f1222 * (f3 * f3) - 6 * f13 * f23 * f22 + 3 * f123 * f3 * f22 + f13 * f3 * f222 + 3 * f23 * f3 * f122);
        REAL c4 =
            -1 / (f3 * f3 * f3) * (f2222 * (f3 * f3) + 3 * f33 * f22 * f22 - 6 * f223 * f3 * f22 - 4 * f23 * f3 * f222 + 12 * f23 * f23 * f22);

        monge_form.coefficients()[6]  = c0;
        monge_form.coefficients()[7]  = c1;
        monge_form.coefficients()[8]  = c2;
        monge_form.coefficients()[9]  = c3;
        monge_form.coefficients()[10] = c4;
    }
}

void Monge_via_jet_fitting::switch_to_direct_orientation(Vector &v1, const Vector &v2, const Vector &v3)
{
    if (dot( v1, cross(v2,v3) ) < 0)
        v1 = -v1;
}

} // end namespace gamer
