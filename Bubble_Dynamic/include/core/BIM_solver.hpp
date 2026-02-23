/*  __       __          _
 * |__)||\/||__) /\ |\/||__)/  \|\/|
 * |__)||  ||__)/--\|  ||__)\__/|  |
 *
 * -----------------------------------------------------------------------------
 * Copyright (C) 2023 Armand Sieber
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 * -----------------------------------------------------------------------------
 *
 */

/*! \file BIM_solver.hpp
    \brief Handles the computation of the normal and tangential velocities on the boundaries of the domain.
 */

#ifndef BIM_SOLVER_HPP
#define BIM_SOLVER_HPP

#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#ifdef _OPENMP
    #include <omp.h>
#endif

#include <armadillo>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "cubic_spline.hpp"
#include "integrands_BIM.hpp"
#include "Case_SphericalBoundary.hpp"

/*! \class BIM_solver
    \brief Class handling the computation of the normal and tangential velocities on the boundaries of the domain.

    The nornmal velocities are computed by solving the boundary integral formulation of Lapalce's equation and the
    tangential velocities are computed from the values of the potentials on the boundaries of the domain.
*/

class BIM_solver {
public:
    BIM_solver(int Nb, int Ns, std::string dumper_name) {
        un_b.resize(Nb + 1);
        un_1.resize(Ns + 1);
        un_2.resize(Ns + 1);
        phi1.resize(Ns + 1);
        phi2.resize(Ns + 1);

        nodes_position.open(dumper_name);

    }

    // Vectors for bubble and fluid-fluid interface dynamics
    std::vector<double> un_b; /*!< \brief Normal velocity on bubble surface */
    std::vector<double> un_1; /*!< \brief Normal velocity on fluid-fluid interface (fluid 1) */
    std::vector<double> un_2; /*!< \brief Normal velocity on fluid-fluid interface (fluid 2) */
    std::vector<double> phi1; /*!< \brief Potential on fluid-fluid interface (fluid 1) */
    std::vector<double> phi2; /*!<  brief Potential on fluid-fluid interface (fluid 2) */

    // Simulation parameters
    double dt; /*!< \brief Time step value */
    double time; /*!< \brief Simulation time */
    int time_step; /*!< \brief Time step index */

    std::ofstream nodes_position; /*!<\brief  Output stream to be written in dumper file */

    // Functions
    /*! \brief Computes the normal velocities on the boundaries using the boundary integral formulation of Laplace equation */
    template<typename BubbleData, typename BoundaryData, typename Inputs>
    void compute_un(BubbleData &bubble, BoundaryData &boundary, Inputs &data); //compute nodes normal velocities

    /*! \brief Computes the tangential velocities from the potentials on the boundaries and their position*/
    template<typename BubbleData, typename BoundaryData, typename Inputs>
    void compute_ut(BubbleData &bubble, BoundaryData &boundary, Inputs &data); //compute nodes tangential velocities

    /*! \brief Write solution to dumper file*/
    template<typename BubbleData, typename BoundaryData, typename Inputs>
    void write_solution(BubbleData &bubble, BoundaryData &boundary, Inputs &data); // write solution

};

using namespace std;
using namespace arma;


template<typename BubbleData, typename BoundaryData, typename Inputs>
void BIM_solver::compute_un(BubbleData &bubble, BoundaryData &boundary, Inputs &data) {

    gsl_set_error_handler_off();

    // ----- Set up quantities -----
    // -----------------------------


    int Nb = data.Nb; // elements on bubble
    int Ns = data.Ns; // elements on fluid-fluid interface
    const int N = Nb + 1 + Ns + 1;

    double alpha = data.alpha;

    // G: coefficient matrix containing the integral of Green's function
    mat G_nodes(N, N, fill::zeros);
    // H: coefficient matrix containing the integral of Green's function normal derivative
    mat H_nodes(N, N, fill::zeros);

    vec phi = conv_to<vec>::from(bubble->phi_nodes);
    vec F = conv_to<vec>::from(boundary->F_nodes);

    vec r0b = conv_to<vec>::from(bubble->r_nodes);
    vec r0s = conv_to<vec>::from(boundary->r_nodes);
    vec z0b = conv_to<vec>::from(bubble->z_nodes);
    vec z0s = conv_to<vec>::from(boundary->z_nodes);
    vec r0_both = join_cols(r0b, r0s);
    vec z0_both = join_cols(z0b, z0s);
    double r0 = 0.0;
    double z0 = 0.0;

    // ----- Set up cubic splines -----
    // --------------------------------

    vec drb = diff(r0b);
    vec dzb = diff(z0b);
    vec drs = diff(r0s);
    vec dzs = diff(z0s);

    vec pointsb = sqrt(drb % drb + dzb % dzb);
    vec pointss = sqrt(drs % drs + dzs % dzs);

    vec distance_b = cumsum(pointsb);
    vec distance_s = cumsum(pointss);
    vec begin(1, fill::zeros);

    distance_b = join_cols(begin, distance_b);
    distance_s = join_cols(begin, distance_s);
    std::vector<double> dist_b = conv_to < std::vector < double >> ::from(distance_b);
    std::vector<double> dist_s = conv_to < std::vector < double >> ::from(distance_s);

    // Cubic splines on bubble
    cubic_spline spline_rb;
    spline_rb.set_spline(dist_b, bubble->r_nodes, 1.0, -1.0);

    cubic_spline spline_zb;
    spline_zb.set_spline(dist_b, bubble->z_nodes, 0.0, 0.0);

    // Cubic splines on fluid-fluid interface
    cubic_spline spline_rs;
    spline_rs.set_spline(dist_s, boundary->r_nodes, boundary->drds1, boundary->drds2);

    cubic_spline spline_zs;
    spline_zs.set_spline(dist_s, boundary->z_nodes, boundary->dzds1, boundary->dzds2);

    // ----- Loop over all node points -----
    // -------------------------------------

#pragma omp parallel num_threads(omp_get_max_threads())
    {
#pragma omp for private(r0, z0) schedule(dynamic)
        for (int i = 0; i < N; ++i) {
            r0 = r0_both(i);
            z0 = z0_both(i);
            for (int j = 0; j < N; ++j) {

                // Extract integration intervals and cubic splines coefficients for integration
                // ----------------------------------------------------------------------------

                double a, a2, b, b2; // two integration interval with subscript 1 and 2
                double m_ar, m_br, m_cr, m_dr, m_az, m_bz, m_cz, m_dz; // spline coefficients interval 1
                double m_ar2, m_br2, m_cr2, m_dr2, m_az2, m_bz2, m_cz2, m_dz2; // spline coefficients interval 2

                // The integration interval lies on the fluid-fluid interface
                if (j > Nb) {
                    a = 0.0;
                    a2 = 0.0;
                    b = 0.0;
                    b2 = 0.0;

                    if (j < (N - 1)) {
                        b = distance_s[1 + j - (Nb + 1)] - distance_s[j - (Nb + 1)];
                        m_ar = spline_rs.m_a[j - (Nb + 1)];
                        m_br = spline_rs.m_b[j - (Nb + 1)];
                        m_cr = spline_rs.m_c[j - (Nb + 1)];
                        m_dr = r0_both(j);
                        m_az = spline_zs.m_a[j - (Nb + 1)];
                        m_bz = spline_zs.m_b[j - (Nb + 1)];
                        m_cz = spline_zs.m_c[j - (Nb + 1)];
                        m_dz = z0_both(j);

                    }
                    if (j > Nb + 1) {
                        b2 = distance_s[j - (Nb + 1)] - distance_s[j - 1 - (Nb + 1)];
                        m_ar2 = spline_rs.m_a[j - 1 - (Nb + 1)];
                        m_br2 = spline_rs.m_b[j - 1 - (Nb + 1)];
                        m_cr2 = spline_rs.m_c[j - 1 - (Nb + 1)];
                        m_dr2 = r0_both(j - 1);
                        m_az2 = spline_zs.m_a[j - 1 - (Nb + 1)];
                        m_bz2 = spline_zs.m_b[j - 1 - (Nb + 1)];
                        m_cz2 = spline_zs.m_c[j - 1 - (Nb + 1)];
                        m_dz2 = z0_both(j - 1);
                    }
                }

                    // The integration interval lies on the fluid-fluid interface
                else {
                    a = 0.0;
                    a2 = 0.0;
                    b = 0.0;
                    b2 = 0.0;

                    if (j < Nb) {
                        b = distance_b[1 + j] - distance_b[j];
                        m_ar = spline_rb.m_a[j];
                        m_br = spline_rb.m_b[j];
                        m_cr = spline_rb.m_c[j];
                        m_dr = r0_both(j);
                        m_az = spline_zb.m_a[j];
                        m_bz = spline_zb.m_b[j];
                        m_cz = spline_zb.m_c[j];
                        m_dz = z0_both(j);
                    }
                    if (j > 0) {
                        b2 = distance_b[j] - distance_b[j - 1];
                        m_ar2 = spline_rb.m_a[j - 1];
                        m_br2 = spline_rb.m_b[j - 1];
                        m_cr2 = spline_rb.m_c[j - 1];
                        m_dr2 = r0_both(j - 1);
                        m_az2 = spline_zb.m_a[j - 1];
                        m_bz2 = spline_zb.m_b[j - 1];
                        m_cz2 = spline_zb.m_c[j - 1];
                        m_dz2 = z0_both(j - 1);
                    }

                }

                // Compute the boundary integrals
                // ------------------------------

                // Special treatment for integration intervals at the bubble and fluid-fluid interface endpoints
                if ((j == 0) || (j == (Nb + 1))) {
                    gsl_integration_workspace *w_G1 = gsl_integration_workspace_alloc(2000);
                    gsl_integration_workspace *w_H1 = gsl_integration_workspace_alloc(2000);

                    double result_G1, error_G1;
                    double result_H1, error_H1;

                    G_H_params params;

                    params.coeff_r0 = m_ar;
                    params.coeff_r1 = m_br;
                    params.coeff_r2 = m_cr;
                    params.coeff_r3 = m_dr;
                    params.coeff_z0 = m_az;
                    params.coeff_z1 = m_bz;
                    params.coeff_z2 = m_cz;
                    params.coeff_z3 = m_dz;
                    params.r0 = r0;
                    params.z0 = z0;
                    params.a = a;
                    params.b = b;
                    params.iter = j;
                    params.N_b = Nb;
                    params.N_s = Ns;

                    gsl_function F1_1;
                    F1_1.function = &G1;
                    F1_1.params = &params;
                    gsl_function F2_1;
                    F2_1.function = &H1;
                    F2_1.params = &params;

                    gsl_integration_qags(&F1_1, a, b, 1e-12, 1e-12, 2000, w_G1, &result_G1, &error_G1);
                    gsl_integration_qags(&F2_1, a, b, 1e-12, 1e-12, 2000, w_H1, &result_H1, &error_H1);

                    G_nodes(i, j) = result_G1;
                    H_nodes(i, j) = result_H1;

                    gsl_integration_workspace_free(w_G1);
                    gsl_integration_workspace_free(w_H1);

                }

                    // Special treatment for integration intervals at the bubble and fluid-fluid interface endpoints
                else if ((j == Nb) || (j == (Nb + Ns + 1))) {

                    gsl_integration_workspace *w_G2 = gsl_integration_workspace_alloc(2000);
                    gsl_integration_workspace *w_H2 = gsl_integration_workspace_alloc(2000);

                    double result_G2, error_G2;
                    double result_H2, error_H2;

                    G_H_params params;

                    params.coeff_r0 = m_ar2;
                    params.coeff_r1 = m_br2;
                    params.coeff_r2 = m_cr2;
                    params.coeff_r3 = m_dr2;
                    params.coeff_z0 = m_az2;
                    params.coeff_z1 = m_bz2;
                    params.coeff_z2 = m_cz2;
                    params.coeff_z3 = m_dz2;
                    params.r0 = r0;
                    params.z0 = z0;
                    params.a = a2;
                    params.b = b2;
                    params.iter = j;
                    params.N_b = Nb;
                    params.N_s = Ns;

                    gsl_function F1_2;
                    F1_2.function = &G2;
                    F1_2.params = &params;
                    gsl_function F2_2;
                    F2_2.function = &H2;
                    F2_2.params = &params;


                    gsl_integration_qags(&F1_2, a2, b2, 1e-12, 1e-12, 2000, w_G2, &result_G2, &error_G2);
                    gsl_integration_qags(&F2_2, a2, b2, 1e-12, 1e-12, 2000, w_H2, &result_H2, &error_H2);

                    G_nodes(i, j) = result_G2;
                    H_nodes(i, j) = result_H2;


                    gsl_integration_workspace_free(w_G2);
                    gsl_integration_workspace_free(w_H2);

                }

                    // For all other integration intervals on the bubble and fluid-fluid interface
                else {

                    gsl_integration_workspace *w_G1 = gsl_integration_workspace_alloc(2000);
                    gsl_integration_workspace *w_H1 = gsl_integration_workspace_alloc(2000);

                    double result_G1, error_G1;
                    double result_H1, error_H1;

                    G_H_params params;

                    params.coeff_r0 = m_ar;
                    params.coeff_r1 = m_br;
                    params.coeff_r2 = m_cr;
                    params.coeff_r3 = m_dr;
                    params.coeff_z0 = m_az;
                    params.coeff_z1 = m_bz;
                    params.coeff_z2 = m_cz;
                    params.coeff_z3 = m_dz;
                    params.r0 = r0;
                    params.z0 = z0;
                    params.a = a;
                    params.b = b;
                    params.iter = j;
                    params.N_b = Nb;
                    params.N_s = Ns;

                    gsl_function F1_1;
                    F1_1.function = &G1;
                    F1_1.params = &params;
                    gsl_function F2_1;
                    F2_1.function = &H1;
                    F2_1.params = &params;

                    gsl_integration_qags(&F1_1, a, b, 1e-12, 1e-12, 2000, w_G1, &result_G1, &error_G1);
                    gsl_integration_qags(&F2_1, a, b, 1e-12, 1e-12, 2000, w_H1, &result_H1, &error_H1);

                    double G1_tmp = result_G1;
                    double H1_tmp = result_H1;

                    gsl_integration_workspace_free(w_G1);
                    gsl_integration_workspace_free(w_H1);

                    gsl_integration_workspace *w_G2 = gsl_integration_workspace_alloc(2000);
                    gsl_integration_workspace *w_H2 = gsl_integration_workspace_alloc(2000);


                    double result_G2, error_G2;
                    double result_H2, error_H2;

                    G_H_params params2;

                    params2.coeff_r0 = m_ar2;
                    params2.coeff_r1 = m_br2;
                    params2.coeff_r2 = m_cr2;
                    params2.coeff_r3 = m_dr2;
                    params2.coeff_z0 = m_az2;
                    params2.coeff_z1 = m_bz2;
                    params2.coeff_z2 = m_cz2;
                    params2.coeff_z3 = m_dz2;
                    params2.r0 = r0;
                    params2.z0 = z0;
                    params2.a = a2;
                    params2.b = b2;
                    params2.iter = j;
                    params2.N_b = Nb;
                    params2.N_s = Ns;

                    gsl_function F1_2;
                    F1_2.function = &G2;
                    F1_2.params = &params2;
                    gsl_function F2_2;
                    F2_2.function = &H2;
                    F2_2.params = &params2;


                    gsl_integration_qags(&F1_2, a2, b2, 1e-12, 1e-12, 2000, w_G2, &result_G2, &error_G2);
                    gsl_integration_qags(&F2_2, a2, b2, 1e-12, 1e-12, 2000, w_H2, &result_H2, &error_H2);

                    double G2_tmp = result_G2;
                    double H2_tmp = result_H2;


                    gsl_integration_workspace_free(w_G2);
                    gsl_integration_workspace_free(w_H2);

                    G_nodes(i, j) = G1_tmp + G2_tmp;
                    H_nodes(i, j) = H1_tmp + H2_tmp;


                }
            }
        }
    }

    // ----- 4 pi method -----
    // -----------------------

    // The diagonal elements of matrix H for bublle-bubble contributions
    // can be expressed as H_ii = 4*pi - sum(H_ij) with j =! i
    // see Wang et al. (Strong interaction between a buoyancy bubble and a free surface, 1996)

    double pi = 3.14159265358979323846264338328;
    for (int i = 0; i < Nb + 1; ++i) {
        H_nodes(i, i) = 4.0 * pi;
        for (int j = 0; j < Nb + 1; ++j) {
            if (i != j) {
                H_nodes(i, i) = H_nodes(i, i) - H_nodes(i, j);
            }
        }
    }


    // ----- Solve matrix system -----
    // ------------------------------

    // Re-arange the computed matrices to extract the normal velocities (un_b, un_1, un_2)
    // and the potentials phi1 and phi2. The "subscript" b refers to quantities on the
    // bubble surface, 1 to quantities at the fluid-fluid interface on the side of the
    // first fluid and 2 to quantities at the fluid-fluid interface on the side of the
    // second fluid

    mat A1(Nb + 1, Nb + 1, fill::zeros);
    A1 = H_nodes.submat(0, 0, Nb, Nb);
    mat A2(Nb + 1, Ns + 1, fill::zeros);
    A2 = H_nodes.submat(0, Nb + 1, Nb, N - 1);
    mat A3(Ns + 1, Nb + 1, fill::zeros);
    A3 = H_nodes.submat(Nb + 1, 0, N - 1, Nb);
    mat A4(Ns + 1, Ns + 1, fill::zeros);
    A4 = H_nodes.submat(Nb + 1, Nb + 1, N - 1, N - 1);

    mat B1(Nb + 1, Nb + 1, fill::zeros);
    B1 = G_nodes.submat(0, 0, Nb, Nb);
    mat B2(Nb + 1, Ns + 1, fill::zeros);
    B2 = G_nodes.submat(0, Nb + 1, Nb, N - 1);
    mat B3(Ns + 1, Nb + 1, fill::zeros);
    B3 = G_nodes.submat(Nb + 1, 0, N - 1, Nb);
    mat B4(Ns + 1, Ns + 1, fill::zeros);
    B4 = G_nodes.submat(Nb + 1, Nb + 1, N - 1, N - 1);

    mat A11 = -B1;
    mat A12 = A2 * (1.0 - alpha);
    mat A21 = -B3;
    mat A22 = A4 + mat(Ns + 1, Ns + 1, fill::eye) * 2.0 * pi + alpha * (mat(Ns + 1, Ns + 1, fill::eye) * 2.0 * pi - A4);
    vec b_1 = -A1 * phi + A2 * F;
    vec b_2 = -A3 * phi - (mat(Ns + 1, Ns + 1, fill::eye) * 2.0 * pi - A4) * F;
    mat A_horizontal1 = join_rows(A11, A12);
    mat A_horizontal2 = join_rows(A21, A22);
    mat A = join_cols(A_horizontal1, A_horizontal2);
    
    /* Comment out regularization for testing
    // Add small regularization to improve numerical stability
    mat A_reg = A + mat(A.n_rows, A.n_cols, fill::eye) * 1e-10;
    */
    
    vec b_tot = join_cols(b_1, b_2);
    vec X_sol;
    
    /* Error handling with pseudoinverse fallback - commented out for testing
    try {
        X_sol = solve(A, b_tot, solve_opts::refine);
    } catch (const std::runtime_error& e) {
        // If the solve fails, try with pseudoinverse
        X_sol = pinv(A) * b_tot;
    }
    */
    
    // Direct solve without error handling for testing
    X_sol = solve(A, b_tot, solve_opts::refine);
    
    vec X_sol_unb(Nb + 1);
    vec X_sol_phi1(Ns + 1);
    vec uns(Ns + 1);

    for (int i = 0; i < N; ++i) {
        if (i > Nb) {
            X_sol_phi1(i - (Nb + 1)) = X_sol(i);
        } else {
            X_sol_unb(i) = X_sol(i);
        }
    }

    /* Comment out regularization for testing
    // Add small regularization to B4 for stability
    mat B4_reg = B4 + mat(B4.n_rows, B4.n_cols, fill::eye) * 1e-10;
    */
    mat eye_mat = mat(Ns + 1, Ns + 1, fill::eye) * 2.0 * pi;
    
    // Regular matrix inversion without regularization
    uns = (-inv(B4) * (eye_mat - A4)) * (alpha * X_sol_phi1 + F);

    un_b = conv_to < std::vector < double >> ::from(-X_sol_unb);
    un_1 = conv_to < std::vector < double >> ::from(-uns);
    un_2 = conv_to < std::vector < double >> ::from(uns);
    phi1 = conv_to < std::vector < double >> ::from(X_sol_phi1);
    phi2 = conv_to < std::vector < double >> ::from(F + alpha * X_sol_phi1);

}


template<typename BubbleData, typename BoundaryData, typename Inputs>
void BIM_solver::compute_ut(BubbleData &bubble, BoundaryData &boundary, Inputs &data) {


    // ----- Set up quantities -----
    // -----------------------------

    int Nb = data.Nb;
    int Ns = data.Ns;

    vec r0b = conv_to<vec>::from(bubble->r_nodes);
    vec r0s = conv_to<vec>::from(boundary->r_nodes);
    vec z0b = conv_to<vec>::from(bubble->z_nodes);
    vec z0s = conv_to<vec>::from(boundary->z_nodes);

    // ----- Set up cubic splines -----
    // --------------------------------

    vec drb = diff(r0b);
    vec dzb = diff(z0b);
    vec drs = diff(r0s);
    vec dzs = diff(z0s);
    vec pointsb = sqrt(drb % drb + dzb % dzb);
    vec pointss = sqrt(drs % drs + dzs % dzs);
    vec distance_b = cumsum(pointsb);
    vec distance_s = cumsum(pointss);
    vec begin(1, fill::zeros);
    distance_b = join_cols(begin, distance_b);
    distance_s = join_cols(begin, distance_s);
    std::vector<double> dist_b = conv_to < std::vector < double >> ::from(distance_b);
    std::vector<double> dist_s = conv_to < std::vector < double >> ::from(distance_s);


    // bubble
    cubic_spline spline_rb;
    spline_rb.set_spline(dist_b, bubble->r_nodes, 1.0, -1.0);

    cubic_spline spline_zb;
    spline_zb.set_spline(dist_b, bubble->z_nodes, 0.0, 0.0);

    cubic_spline spline_phib;
    spline_phib.set_spline(dist_b, bubble->phi_nodes, 0.0, 0.0);

    // fluid-fluid interface
    cubic_spline spline_rs;
    spline_rs.set_spline(dist_s, boundary->r_nodes, boundary->drds1, boundary->drds2);

    cubic_spline spline_zs;
    spline_zs.set_spline(dist_s, boundary->z_nodes, boundary->dzds1, boundary->dzds2);

    cubic_spline spline_phi1;
    spline_phi1.set_spline(dist_s, phi1, boundary->dphi1ds1, boundary->dphi1ds2);

    cubic_spline spline_phi2;
    spline_phi2.set_spline(dist_s, phi2, boundary->dphi2ds1, boundary->dphi2ds2);


    // ----- Compute tangential velocity -----
    // ---------------------------------------

    vector<double> drbds(Nb + 1);
    vector<double> dzbds(Nb + 1);
    vector<double> dphibds(Nb + 1);
    vector<double> drsds(Ns + 1);
    vector<double> dzsds(Ns + 1);
    vector<double> dphi1ds(Ns + 1);
    vector<double> dphi2ds(Ns + 1);

    // Bubble
    double m_ar, m_br, m_cr;
    double m_az, m_bz, m_cz;
    double m_aphi, m_bphi, m_cphi;
    double xib;

    for (int i = 0; i < Nb + 1; ++i) {

        if (i < Nb) {
            xib = 0.0;

            m_ar = spline_rb.m_a[i];
            m_br = spline_rb.m_b[i];
            m_cr = spline_rb.m_c[i];

            m_az = spline_zb.m_a[i];
            m_bz = spline_zb.m_b[i];
            m_cz = spline_zb.m_c[i];

            m_aphi = spline_phib.m_a[i];
            m_bphi = spline_phib.m_b[i];
            m_cphi = spline_phib.m_c[i];
        } else {
            xib = dist_b[i] - dist_b[i - 1];
        }
        drbds[i] = m_cr + 2.0 * m_br * xib + 3.0 * m_ar * xib * xib;
        dzbds[i] = m_cz + 2.0 * m_bz * xib + 3.0 * m_az * xib * xib;
        dphibds[i] = m_cphi + 2.0 * m_bphi * xib + 3.0 * m_aphi * xib * xib;

        bubble->ur_nodes[i] = -un_b[i] * dzbds[i] + dphibds[i] * drbds[i];
        bubble->uz_nodes[i] = un_b[i] * drbds[i] + dphibds[i] * dzbds[i];
        bubble->u_nodes[i] = sqrt(un_b[i] * un_b[i] + dphibds[i] * dphibds[i]);

    }
    bubble->ur_nodes[0] = 0.0;
    bubble->ur_nodes[Nb] = 0.0;

    // Fluid-fluid interface
    double m_ars, m_brs, m_crs;
    double m_azs, m_bzs, m_czs;
    double m_aphi1, m_bphi1, m_cphi1;
    double m_aphi2, m_bphi2, m_cphi2;
    double xis;

    /* 
    // VOLUME CONSERVATION MECHANISM COMMENTED OUT FOR TESTING
    
    // 1) Compute net flux on inner surface (bubble)
    double netFluxInner = 0.0;
    for (int i = 0; i < Nb; ++i) {
        double ri_mid = 0.5 * (bubble->r_nodes[i] + bubble->r_nodes[i+1]);
        double ds_i = dist_b[i+1] - dist_b[i];
        netFluxInner += un_b[i] * ri_mid * ds_i;
    }
    
    // 2) Compute net flux on outer surface (droplet)
    double netFluxOuter = 0.0;
    for (int i = 0; i < Ns; ++i) {
        double ro_mid = 0.5 * (boundary->r_nodes[i] + boundary->r_nodes[i+1]);
        double ds_o = dist_s[i+1] - dist_s[i];
        netFluxOuter += un_1[i] * ro_mid * ds_o;
    }
    
    // 3) Enforce ΔV_bubble - ΔV_droplet = 0 (volumes change in same direction)
    double totalFlux = netFluxOuter - netFluxInner;
    
    // Compute ∫ r ds over outer loop for normalization
    double areaCoeff = 0.0;
    for (int i = 0; i < Ns; ++i) {
        double ro_mid = 0.5 * (boundary->r_nodes[i] + boundary->r_nodes[i+1]);
        double ds_o = dist_s[i+1] - dist_s[i];
        areaCoeff += ro_mid * ds_o;
    }
    
    // Average normal-velocity shift to enforce volume conservation
    // Using ΔV_bubble - ΔV_droplet = 0 which means volumes change in same direction
    // This results in the droplet volume increasing as bubble volume increases
    double avgShift = totalFlux / areaCoeff;
    
    // Adjust outer boundary normal velocities to conserve volume
    // Negative sign on un_1 adjustment makes droplet volume change in same direction as bubble
    for (int i = 0; i < Ns + 1; ++i) {
        un_1[i] -= avgShift;
        un_2[i] += avgShift; // Adjust the opposing side by the same amount, opposite sign
    }
    */
    
    // With volume conservation commented out, we're using the directly calculated normal velocities
    // This will allow us to observe how the simulation behaves without volume enforcement
    
    // Process derivatives and compute velocities
    for (int i = 0; i < Ns + 1; ++i) {
        if (i < Ns) {
            xis = 0.0;
            
            m_ars = spline_rs.m_a[i];
            m_brs = spline_rs.m_b[i];
            m_crs = spline_rs.m_c[i];
            
            m_azs = spline_zs.m_a[i];
            m_bzs = spline_zs.m_b[i];
            m_czs = spline_zs.m_c[i];
            
            m_aphi1 = spline_phi1.m_a[i];
            m_bphi1 = spline_phi1.m_b[i];
            m_cphi1 = spline_phi1.m_c[i];
            
            m_aphi2 = spline_phi2.m_a[i];
            m_bphi2 = spline_phi2.m_b[i];
            m_cphi2 = spline_phi2.m_c[i];
        } else {
            xis = dist_s[i] - dist_s[i - 1];
        }
        
        drsds[i] = m_crs + 2.0 * m_brs * xis + 3.0 * m_ars * xis * xis;
        dzsds[i] = m_czs + 2.0 * m_bzs * xis + 3.0 * m_azs * xis * xis;
        dphi1ds[i] = m_cphi1 + 2.0 * m_bphi1 * xis + 3.0 * m_aphi1 * xis * xis;
        dphi2ds[i] = m_cphi2 + 2.0 * m_bphi2 * xis + 3.0 * m_aphi2 * xis * xis;
        
        boundary->ur_nodes[i] = -un_1[i] * dzsds[i] + dphi1ds[i] * drsds[i];
        boundary->uz_nodes[i] = un_1[i] * drsds[i] + dphi1ds[i] * dzsds[i];
        boundary->u_nodes1[i] = sqrt(un_1[i] * un_1[i] + dphi1ds[i] * dphi1ds[i]);
        boundary->u_nodes2[i] = sqrt(un_2[i] * un_2[i] + dphi2ds[i] * dphi2ds[i]);
        boundary->u_vec_prod[i] = -un_1[i] * un_2[i] + dphi1ds[i] * dphi2ds[i];

    }

    boundary->ur_nodes[Ns] = 0.0;


}

template<typename BubbleData, typename BoundaryData, typename Inputs>
void BIM_solver::write_solution(BubbleData &bubble, BoundaryData &boundary, Inputs &data) {

    std::vector<double> data_dump = bubble->r_nodes; // r-coordinate of the nodes on the bubble surface
    data_dump.insert(data_dump.end(), bubble->z_nodes.begin(), bubble->z_nodes.end()); // z-coordinate of the nodes on the bubble surface
    data_dump.insert(data_dump.end(), bubble->phi_nodes.begin(), bubble->phi_nodes.end()); // value of the potentials on the bubble surface
    data_dump.insert(data_dump.end(), un_b.begin(), un_b.end()); // value of the normal velocities on the bubble surface
    data_dump.insert(data_dump.end(), boundary->r_nodes.begin(), boundary->r_nodes.end()); // r-coordinate of the nodes on the fluid-fluid interface
    data_dump.insert(data_dump.end(), boundary->z_nodes.begin(), boundary->z_nodes.end()); // z-coordinate of the nodes on the fluid-fluid interface
    data_dump.insert(data_dump.end(), phi1.begin(), phi1.end()); // value of the potentials on the fluid-fluid interface (fluid 1)
    data_dump.insert(data_dump.end(), un_1.begin(), un_1.end()); // value of the normal velocities on the fluid-fluid interface (fluid 1)
    data_dump.insert(data_dump.end(), boundary->F_nodes.begin(), boundary->F_nodes.end()); // value of F on the fluid-fluid interface
    data_dump.insert(data_dump.end(), bubble->V); // bubble volume
    
    // Add spherical volume to output if using spherical boundary
    if (auto* sphericalBoundary = dynamic_cast<Case_SphericalBoundary*>(boundary.get())) {
        data_dump.insert(data_dump.end(), sphericalBoundary->V_spherical); // spherical droplet volume
    }

    if (time_step == 0){
        nodes_position << showpoint << "Nb" << "\t" << data.Nb << "\n";
        nodes_position << showpoint << "Ns" << "\t" << data.Ns << "\n";
        nodes_position << showpoint
                   << "gamma" << "\t" << data.gamma
                   << "\t" << "xi"    << "\t" << data.xi
                   << "\t" << "Rmax"  << "\t" << data.Rmax
                   << "\t" << "e"     << "\t" << data.e
                   << "\t" << "Omega" << "\t" << data.Omega
                   << "\t" << "A"     << "\t" << data.A << "\n";
        nodes_position << showpoint << "alpha" << "\t" << data.alpha << "\n";
        nodes_position << showpoint << "zeta" << "\t" << data.zeta << "\n";
        nodes_position << showpoint << "Specific heat ratio" << "\t" << data.k << "\n";
        nodes_position << showpoint << "Bubble dynammics" << "\t" << data.bubble_dynamics << "\n";
        nodes_position << showpoint << "epsilon (for Rayleigh-Plesset bubble)" << "\t" << data.epsilon << "\n";
        nodes_position << showpoint << "Surface elasticity (0 false, 1 true)" << "\t" << data.surface_elasticity << "\n";
        nodes_position << showpoint << "sigma_s" << "\t" << data.sigma_s << "\n";
        
        // No need to add a comment about spherical volume - SimpleMP4.py already knows to look for it
        // The header should always be exactly 10 lines
    }
    nodes_position << showpoint << time_step << "\t" << time << "\t";

    for (int i = 0; i < data_dump.size(); ++i) {
        nodes_position << showpoint << data_dump[i] << "\t";
    }
    nodes_position << "\n";

}

#endif // BIM_SOLVER_HPP
