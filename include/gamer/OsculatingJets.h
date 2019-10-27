
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
 * @param[in]  n     Positive integer for factorial
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

/**
 * @brief      This class describes a monge via jet fitting.
 */
class Monge_via_jet_fitting
{
  public:
    /**
     * @brief      Representation of Monge parameterization
     */
    class MongeForm
    {
      protected:
        Vector m_origin_pt;  /// Origin of parameterization
        Vector m_d1;         /// Maximal principal direction
        Vector m_d2;         /// Minimal principal direction
        Vector m_n;          /// Normal directions
        /// Vector of differential values
        std::vector<REAL> m_coefficients;
        // coeff = (k1, k2, //ppal curv
        // b0, b1, b2, b3, //third order
        // c0, c1, c2, c3, c4) //fourth order
      public:
        /**
         * @brief      Default constructor
         */
        MongeForm(std::size_t degree)
        {
            m_origin_pt = Vector({0., 0., 0.});
            m_d1 = Vector({0., 0., 0.});
            m_d2 = Vector({0., 0., 0.});
            m_n  = Vector({0., 0., 0.});
            m_coefficients = std::vector<REAL>();

            //if d>=2, number of coeffs = (d+1)(d+2)/2-4.
            //we remove cst, linear and the xy coeff which vanish
            if (degree >= 2)
                std::fill_n(back_inserter(m_coefficients),
                            (degree + 1) * (degree + 2) / 2 - 4, 0.);
        }

        /**
         * @brief      Destroys the object.
         */
        ~MongeForm() {
        }

        /**
         * @brief      Get the origin
         *
         * @return     Vector origin
         */
        const Vector origin() const {
            return m_origin_pt;
        }

        /**
         * @brief      Get the origin
         *
         * @return     Vector origin
         */
        Vector &origin() {
            return m_origin_pt;
        }

        /**
         * @brief      Get the max principal direction
         *
         * @return     Vector max principal direction
         */
        const Vector maximal_principal_direction() const {
            return m_d1;
        }

        /**
         * @brief      Get the max principal direction
         *
         * @return     Vector max principal direction
         */
        Vector &maximal_principal_direction() {
            return m_d1;
        }

        /**
         * @brief      Get the min principal direction
         *
         * @return     Vector min principal direction
         */
        const Vector minimal_principal_direction() const {
            return m_d2;
        }

        /**
         * @brief      Get the min principal direction
         *
         * @return     Vector min principal direction
         */
        Vector &minimal_principal_direction() {
            return m_d2;
        }

        /**
         * @brief      Get the normal direction
         *
         * @return     Vector normal direction
         */
        const Vector normal_direction() const {
            return m_n;
        }

        /**
         * @brief      Get the normal direction
         *
         * @return     Vector normal direction
         */
        Vector &normal_direction() {
            return m_n;
        }

        /**
         * @brief      Access coefficients
         *
         * @return     Coefficients
         */
        const std::vector<REAL> coefficients() const {
            return m_coefficients;
        }

        /**
         * @brief      Access coefficients
         *
         * @return     Coefficients
         */
        std::vector<REAL> &coefficients() {
            return m_coefficients;
        }

        /**
         * @brief      Access principal curvature values
         *
         * @param[in]  i     Index
         *
         * @return     Value of the curvature
         */
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

        /**
         * @brief      Access third order differential properties
         *
         * @param[in]  i     Index
         *
         * @return     Values
         */
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

        /**
         * @brief      Access fourth order differential properties
         *
         * @param[in]  i     Index
         *
         * @return     Values
         */
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

        //switch min-max ppal curv/dir wrt a given normal orientation.
        // if given_normal.monge_normal < 0 then change the orientation
        // if z=g(x,y) in the basis (d1,d2,n) then in the basis
        // (d2,d1,-n)
        // z=h(x,y)=-g(y,x)
        void comply_wrt_given_normal(const Vector &given_normal){
            if (dot(given_normal, normal_direction()) < 0)
            {
                //normal_direction() = -normal_direction();
                m_n = -m_n;
                // std::swap(maximal_principal_direction(),
                // minimal_principal_direction());
                std::swap(m_d1, m_d2);
                if (m_coefficients.size() >= 2)
                    std::swap(m_coefficients[0], m_coefficients[1]);
                if (m_coefficients.size() >= 6)
                {
                    std::swap(m_coefficients[2], m_coefficients[5]);
                    std::swap(m_coefficients[3], m_coefficients[4]);
                }
                if (m_coefficients.size() >= 11)
                {
                    std::swap(m_coefficients[6], m_coefficients[10]);
                    std::swap(m_coefficients[7], m_coefficients[9]);
                }
                typename std::vector<REAL>::iterator itb = m_coefficients.begin(),
                                                     ite = m_coefficients.end();
                for (; itb != ite; itb++)
                {
                    *itb = -(*itb);
                }
            }
        }

        void dump_verbose(std::ostream &out_stream) const {
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
        void dump_4ogl(std::ostream &out_stream, const REAL scale){
            if (coefficients().size() < 2)
                throw std::runtime_error("Insufficient coefficients");
            out_stream << origin()  << " "
                       << maximal_principal_direction() * scale << " "
                       << minimal_principal_direction() * scale << " "
                       << coefficients()[0] << " "
                       << coefficients()[1] << " "
                       << std::endl;
        }
    }; // end nested class MongeForm

    /// Default constructor
    Monge_via_jet_fitting(){
        m_pca_basis = std::vector< std::pair<REAL, Vector> >(3);
    }

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
        // Degree of jet fitting
        deg = static_cast<int>(dJet);
        // Degree of differential to compute
        deg_monge = static_cast<int>(dPrime);
        // Number of points required to fit d-jet
        nb_d_jet_coeff = static_cast<int>((dJet+1)*(dJet+2)/2);
        // Number of input points
        nb_input_pts   = static_cast<int>(end - begin);
        // Check that enough points are given to fit d-jet
        if (nb_input_pts < nb_d_jet_coeff)
            throw std::runtime_error("Insufficient points provided to perform jet fitting.");

        // Initialize MongeForm
        MongeForm monge_form(dPrime);

        // Assemble the linear system
        EigenMatrixN M(nb_input_pts, nb_d_jet_coeff);
        EigenVectorN Z(nb_input_pts);

        // Compute
        compute_PCA(begin, end);
        fill_matrix(begin, end, dJet, M, Z);     //with precond

        // std::cout << "M:" << std::endl << M << std::endl;
        // std::cout << "Z:" << std::endl << Z << std::endl;
        // Solve MA=Z in the ls sense. The solution A is stored in Z.
        Eigen::JacobiSVD<Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic, Eigen::DontAlign>> jacobiSvd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Z = jacobiSvd.solve(EigenVectorN(Z));
        condition_nb = jacobiSvd.singularValues().array().abs().maxCoeff() /
                       jacobiSvd.singularValues().array().abs().minCoeff();

        for (int k = 0; k <= deg; k++)
            for (int i = 0; i <= k; i++)
                Z(k*(k+1)/2+i) /= std::pow(preconditionning, k);

        compute_Monge_basis(Z.data(), monge_form);
        if (dPrime >= 3)
            compute_Monge_coefficients(Z.data(), dPrime, monge_form);
        return monge_form;
    }

    /**
     * @brief      Access to condition number
     *
     * @return     Condition number
     */
    const REAL condition_number() const {
        return condition_nb;
    }

    /**
     * @brief      Get the PCA basis
     *
     * @param[in]  i     Index
     *
     * @return     PCA basis vector
     */
    const std::pair<REAL, Vector> pca_basis(std::size_t i) const
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

    /// Translate the points such that p0 (the first point) is at the origin
    Eigen::Affine3d p02origin;
    /// Rotate local orientation to PCA
    Eigen::Affine3d world2pca;
    /// Transform from PCA fitting to Monge
    Eigen::Affine3d pca2monge;

    /**
     * @brief      Compute PCA of points
     *
     * @param[in]  begin          Iterator to origin point
     * @param[in]  end            Iterator to just past the end.
     *
     * @tparam     InputIterator  Typename of iterator
     */
    template <class InputIterator>
    void compute_PCA(InputIterator begin, InputIterator end){
        int n = nb_input_pts;
        REAL x, y, z,
             sumX = 0., sumY = 0., sumZ = 0.,
             sumX2 = 0., sumY2 = 0., sumZ2 = 0.,
             sumXY = 0., sumXZ = 0., sumYZ = 0.,
             xx, yy, zz, xy, xz, yz;

        for (; begin != end; begin++)
        {
            Vector lp = (**begin).position;
            x = lp[0];
            y = lp[1];
            z = lp[2];
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
        // std::cout << "PCA Eigenvalues: "  << std::endl
        //           << eigenvalues << std::endl;
        // std::cout << "PCA Eigenvectors: "  << std::endl
        //           << eigenvectors << std::endl;

        // Store eigenvalues in m_pca_basis
        for (int i = 0; i < 3; i++)
        {
            m_pca_basis[i].first = eigenvalues[2-i];
        }

        Vector v1({eigenvectors(6), eigenvectors(7), eigenvectors(8)});
        m_pca_basis[0].second = v1;
        Vector v2({eigenvectors(3), eigenvectors(4), eigenvectors(5)});
        m_pca_basis[1].second = v2;
        Vector v3({eigenvectors(0), eigenvectors(1), eigenvectors(2)});
        m_pca_basis[2].second = v3;

        switch_to_direct_orientation(m_pca_basis[0].second,
                                     m_pca_basis[1].second,
                                     m_pca_basis[2].second);

        EigenMatrix tmp;
        tmp << m_pca_basis[0].second[0], m_pca_basis[0].second[1], m_pca_basis[0].second[2],
            m_pca_basis[1].second[0], m_pca_basis[1].second[1], m_pca_basis[1].second[2],
            m_pca_basis[2].second[0], m_pca_basis[2].second[1], m_pca_basis[2].second[2];

        //Store the change of basis W->F
        Eigen::Affine3d change_basis(tmp);
        world2pca = change_basis;
    }

    //Coordinates of input points are computed in the fitting basis with
    //  p0 as origin.
    //Preconditionning is computed, M and Z are filled
    template <class InputIterator>
    void fill_matrix(InputIterator begin, InputIterator end,
                     std::size_t d, EigenMatrixN &M, EigenVectorN &Z)
    {
        //origin of fitting coord system = first input data point
        Vector point0 = (**begin).position;
        // std::cout << "point0: " << point0 << std::endl;
        //transform coordinates of sample points with a
        //translation ($-p$) and multiplication by $ P_{W\rightarrow F}$.
        Vector orig({0., 0., 0.});
        Vector v_point0_orig(orig - point0);
        // std::cout << "v_point0_orig: " << v_point0_orig << std::endl;

        p02origin = Eigen::Translation3d(v_point0_orig);
        Eigen::Affine3d transf_points = world2pca *
                                        p02origin;

        //compute and store transformed points
        std::vector<Vector> pts_in_fitting_basis;
        pts_in_fitting_basis.reserve(nb_input_pts);

        for(auto it = begin; it != end; ++it) {
            Vector cur_pt = (**it).position;
            cur_pt = transf_points*EigenMap(cur_pt);
            pts_in_fitting_basis.push_back(cur_pt);
        }

        //Compute preconditionning
        REAL precond = 0.;
        for(auto it = pts_in_fitting_basis.begin(); it != pts_in_fitting_basis.end(); ++it) {
            precond += std::abs((*it)[0]) + std::abs((*it)[1]);
        }
        precond /= 2 * nb_input_pts;
        preconditionning = precond;

        //fill matrices M and Z
        int line_count = 0;
        REAL x, y;
        for(auto it = pts_in_fitting_basis.begin(); it != pts_in_fitting_basis.end(); ++it) {
            // CGAL_For_all(itb, ite) {
            x = (*it)[0];     // x
            y = (*it)[1];     // y
            Z.coeffRef(line_count) = (*it)[2];     //itb->z());
            for (std::size_t k = 0; k <= d; k++)
            {
                for (std::size_t i = 0; i <= k; i++)
                {
                    M.coeffRef(line_count, k * (k + 1) / 2 + i) =
                        std::pow( x, static_cast<int>(k - i) )
                        * std::pow( y, static_cast<int>(i) )
                        / (fact( static_cast<unsigned int>(i) ) *
                           fact( static_cast<unsigned int>(k - i) )
                           * std::pow( preconditionning, static_cast<int>(k) ) );
                }
            }
            line_count++;
        }
    }

    /**
     * @brief      Compute the Monge basis
     *
     * @param[in]  A           Coefficients of the d-jet
     * @param      monge_form  The monge form
     */
    void compute_Monge_basis(const REAL* A, MongeForm &monge_form){
        //bi-index to uni-index conversion : A(i,j)=A[(i+j)(i+j+1)/2+j]
        Vector orig_monge({0., 0., A[0]});
        Vector normal({-A[1], -A[2], 1.});
        REAL norm2 = normal | normal;
        normal /= std::sqrt(norm2);

        monge_form.origin() = p02origin.inverse() * world2pca.inverse() * EigenMap(orig_monge);
        monge_form.normal_direction() = world2pca.inverse() * EigenMap(normal);

        if (deg_monge >= 2)
        {
            // normal = cross(Xu, Xv)
            Vector Xu({1., 0., A[1]});
            Vector Xv({0., 1., A[2]});

            //Surface in fitting_basis : X(u,v)=(u,v,J_A(u,v))
            //in the basis Xu=(1,0,A[1]), Xv=(0,1,A[2]), Weingarten=-I^{-1}II
            //first fond form I=(e,f,f,g)
            //                 =(Xu.Xu, Xu.Xv, Xu.Xv, Xv.Xv)
            //second fond form II=(l,m,m,n)/norm2^(1/2)
            //                   =(n.Xuu, n.Xuv, n.Xuv, n.Xvv)
            //ppal curv are the opposite of the eigenvalues of Weingarten or the
            //  eigenvalues of weingarten = -Weingarten = I^{-1}II

            REAL e = 1 + A[1] * A[1], f = A[1] * A[2], g = 1 + A[2] * A[2],
                 l = A[3], m = A[4], n = A[5];
            Eigen::Matrix<REAL,2,2> weingarten;
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
            Eigen::Matrix<REAL,2,2> change_XuXv2YZ;
            change_XuXv2YZ(0, 0) = 1 / normXu;
            change_XuXv2YZ(0, 1) = -XudotXv / (normXu * normXu * normZ);
            change_XuXv2YZ(1, 0) = 0;
            change_XuXv2YZ(1, 1) = 1 / normZ;

            //in the new orthonormal basis (Y,Z) of the tangent plane :
            weingarten = change_XuXv2YZ.inverse() * weingarten * change_XuXv2YZ;

            // diagonalization of weingarten
            std::array<REAL, 3> W = {{weingarten(0, 0), weingarten(1, 0), weingarten(1, 1)}};
            Eigen::Matrix<REAL,2,1> eval;
            Eigen::Matrix<REAL,2,2> evec;
            EigenDiagonalizeTraits<REAL, 2>::diagonalizeSelfAdjointCovMatrix(W, eval, evec);

            // Principal directions
            Vector d_max = evec(2) * Y + evec(3) * Z,
                   d_min = evec(0) * Y + evec(1) * Z;
            switch_to_direct_orientation(d_max, d_min, normal);

            // Construct transformation to Monge basis
            EigenMatrix tmp;
            tmp << d_max[0], d_max[1], d_max[2],
                d_min[0], d_min[1], d_min[2],
                normal[0], normal[1], normal[2];
            Eigen::Affine3d change_basis (tmp);
            pca2monge = change_basis;

            // Store the monge basis origin and vectors with their world coord
            monge_form.maximal_principal_direction() = world2pca.inverse() * EigenMap(d_max);
            monge_form.minimal_principal_direction() = world2pca.inverse() * EigenMap(d_min);

            // Store principal curvatures
            monge_form.coefficients()[0]  = eval[1];
            monge_form.coefficients()[1]  = eval[0];
        }
    }

    /**
     * @brief      Calculates the 3rd and 4th order differential values
     *
     * @param      A           Vector of coefficients
     * @param[in]  dprime      Degree differential to compute
     * @param      monge_form  MongeForm object
     */
    void compute_Monge_coefficients(REAL* A, std::size_t dprime,
                                    MongeForm &monge_form){
        //One has the equation w=J_A(u,v) of the fitted surface S
        // in the fitting_basis
        //Substituing (u,v,w)=pca2monge^{-1}(x,y,z)
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
        // p stores pca2monge^{-1}=pca2monge^{T}
        REAL p[3][3];
        p[0][0] = pca2monge(0, 0);
        p[1][0] = pca2monge(0, 1);
        p[2][0] = pca2monge(0, 2);
        p[0][1] = pca2monge(1, 0);
        p[1][1] = pca2monge(1, 1);
        p[2][1] = pca2monge(1, 2);
        p[0][2] = pca2monge(2, 0);
        p[1][2] = pca2monge(2, 1);
        p[2][2] = pca2monge(2, 2);

        // formula are designed for w=sum( Aij ui vj), but we have J_A = sum(
        // Aij/i!j! ui vj)
        for (int k = 0; k <= deg; k++)
            for (int i = 0; i <= k; i++)
                A[k * (k + 1) / 2 + i] /= fact(k - i) * fact(i);
        //this is A(k-i;i)

        /*   //debug */
        /*   std::cout << "coeff of A" << std::endl */
        /*      << A[0] << " "<< A[1] << " "<< A[2] << std::endl */
        /*      << A[3] << " "<< A[4] << " "<< A[5] << std::endl */
        /*      << A[6] << " "<< A[7] << " "<< A[8] << " "<< A[9]<< std::endl */
        /*      << A[10] << " "<< A[11] << " "<< A[12] << " "<< A[13]<< " " <<
           A[14] <<
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

    /**
     * @brief      Flip the orientation if det(v1,v2,v3) < 0
     *
     * @param      v1    Vector 1
     * @param[in]  v2    Vector 2
     * @param[in]  v3    Vector 3
     */
    void switch_to_direct_orientation(Vector &v1, const Vector &v2,
                                      const Vector &v3){
        if (dot(v1, cross(v2,v3)) < 0.) {
            v1 = -v1;
        }
    }

    /**
     * @brief      Print operator overload
     *
     * @param      out_stream  The stream to write data to
     * @param[in]  monge       MongeForm object to print
     *
     * @return     The stream
     */
    friend std::ostream& operator<<(std::ostream  &out_stream,
                                    const typename Monge_via_jet_fitting::MongeForm &monge)
    {
        monge.dump_verbose(out_stream);
        return out_stream;
    }
}; // end class Monge_via_jet_fitting
} // end namespace gamer
