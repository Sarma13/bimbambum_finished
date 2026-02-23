/*  __       __          __
 * |__)||\/||__) /\ |\/||__)/  \|\/|
 * |__)||  ||__)/--\|  ||__)\__/|  |
 *
 * This file is part of BIMBAMBUM.
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

/*! \file Case_FlatBoundary.cpp
  \brief Derived class: handles the dynamics of a fluid-fluid interface that has been initiated as a 'infinite' flat surface.

  This derived class inherits members from BoundaryData.cpp and BoundaryData.cpp
  */

#include "Case_FlatBoundary.hpp"

#include <armadillo>
#include <cmath>
#include <vector>

#include "cubic_spline.hpp"

using namespace std;
using namespace arma;


/*! Fluid-fluid interface initialization. The interface is initiated as flat and located at z = 0.
    Theoretically, it extends to infinity (r -> infinity). Concretely, however, it is truncated at
    a distance of r = 30R_max, where R_max is the maximum bubble radius.
 */
void Case_FlatBoundary::initialize() {

    arma::vec r1(5);
    r1 = {30.0, 20.0, 15.0, 10.0, 7.5}; // large node spacing 'far away' from the bubble

    // smaller node spacing close to the bubble (initiated as evenly spaced nodes)
    arma::vec r2 = arma::linspace(7.0, 0.0, Ns - 4);
    arma::vec r_both = arma::join_cols(r1, r2);

    for (int i = 0; i < Ns + 1; ++i) {
        r_nodes[i] = r_both(i);
        z_nodes[i] = 0.0;
        F_nodes[i] = 0.0;
        curv_nodes[i] = 0.0;
    }
}

/*! Remeshing of the of the fluid-fluid interface nodes.
    The remeshing is done by rearranging the node points along the fluid-fluid interface with a
    spacing defined by a geometric progression.
  */
void Case_FlatBoundary::remesh_boundary() {

    // generates surface vectors
    vec rs = conv_to<vec>::from(r_nodes);
    vec zs = conv_to<vec>::from(z_nodes);

    // compute the arclength
    vec drs = diff(rs);
    vec dzs = diff(zs);
    vec pointss = sqrt(drs % drs + dzs % dzs);
    vec distance_s = cumsum(pointss);
    vec begin(1, fill::zeros);
    distance_s = join_cols(begin, distance_s);
    vector<double> distance = conv_to < std::vector < double >> ::from(distance_s);

    // build the cubic splines
    cubic_spline spline_rs;
    spline_rs.set_spline(distance, r_nodes, drds1, drds2);

    cubic_spline spline_zs;
    spline_zs.set_spline(distance, z_nodes, dzds1, dzds2);

    cubic_spline spline_F;
    spline_F.set_spline(distance, F_nodes, 0.0, 0.0);

    // rearrange the nodes along the arclength
    vec d_imposed1(5);
    d_imposed1 = {0.0, 10.0, 15.0, 20.0, 22.5}; // large node spacing 'far away' from the bubble

    double n = Ns - 4; // geometric progression for nodes closer to bubble
    double r = 0.97; // this value ensures that the node spacing decreases as the nodes are located closer to the axis of symmetry
    double h = (r - 1.0) / (pow(r, (1.0 * n - 1.0)) - 1.0);

    vec e(Ns);
    for (int i = 0; i < Ns - 5; ++i) {
        e(i) = pow(r, 1.0 * i);
    }

    vec s = h * cumsum(e) * (distance[Ns] - 23.0);
    s = join_cols(begin, s);
    s(Ns) = distance[Ns];

    double si;
    for (int i = 0; i < Ns + 1; ++i) {
        if (i < 5) {
            si = d_imposed1[i];
        } else {
            si = s(i - 5) + 23.0;
        }
        r_nodes[i] = spline_rs.interpolate(si);
        z_nodes[i] = spline_zs.interpolate(si);
        F_nodes[i] = spline_F.interpolate(si);
    }

    r_nodes[Ns] = 0.0; // imposed value at the axis of symmetry

}

/*! Filtering of the of the fluid-fluid interface nodes.
    The filtering is achieved with a five-point smoothing formula.
    For details, see Longuet-Higgins and Cokelet (The deformation of steep surface waves on water, 1976)
    or Pearson (Hydrodynamics of jet impact in a collapsing bubble, 2002).
  */
void Case_FlatBoundary::filter_boundary() {

    // generates surface vectors
    vector<double> r_copy = r_nodes;
    vector<double> z_copy = z_nodes;
    vector<double> F_copy = F_nodes;

    // Additional dummy node points created to compute polynomial at endpoints of fluid-fluid interface
    r_copy.insert(r_copy.begin(), r_nodes[0] + 5.0);
    r_copy.insert(r_copy.begin(), r_nodes[0] + 10.0);
    r_copy.push_back(-r_nodes[Ns - 1]);
    r_copy.push_back(-r_nodes[Ns - 2]);

    z_copy.insert(z_copy.begin(), 0.0);
    z_copy.insert(z_copy.begin(), 0.0);
    z_copy.push_back(z_nodes[Ns - 1]);
    z_copy.push_back(z_nodes[Ns - 2]);

    F_copy.insert(F_copy.begin(), 0.0);
    F_copy.insert(F_copy.begin(), 0.0);
    F_copy.push_back(F_nodes[Ns - 1]);
    F_copy.push_back(F_nodes[Ns - 2]);

    // compute the arclength
    vec rs = conv_to<vec>::from(r_copy);
    vec zs = conv_to<vec>::from(z_copy);
    vec drs = diff(rs);
    vec dzs = diff(zs);
    vec pointss = sqrt(drs % drs + dzs % dzs);
    vec distance_s = cumsum(pointss);
    vec begin(1, fill::zeros);
    distance_s = join_cols(begin, distance_s);
    vector<double> distance = conv_to < std::vector < double >> ::from(distance_s);

    // Longuet-Higgings & Cokelet smoothing for irregular node distribution
    double s1, s2, s3, s4, s5;
    double c1, c2, c3, c4, c5;

    for (int i = 0; i < Ns + 1; ++i) {
        s1 = distance[i];
        s2 = distance[i + 1];
        s3 = distance[i + 2];
        s4 = distance[i + 3];
        s5 = distance[i + 4];
        c1 = 0.5 * ((s3 - s2) * (s3 - s4)) / ((s1 - s3) * (s1 - s5));
        c2 = 0.5 * ((s4 - s3)) / ((s4 - s2));
        c3 = 0.5 * (1.0 + ((s3 - s2) * (s3 - s4)) / ((s3 - s1) * (s3 - s5)));
        c4 = 0.5 * ((s2 - s3)) / ((s2 - s4));
        c5 = 0.5 * ((s3 - s2) * (s3 - s4)) / ((s5 - s1) * (s5 - s3));
        r_nodes[i] = c1 * r_copy[i] + c2 * r_copy[i + 1] + c3 * r_copy[i + 2] +
                     c4 * r_copy[i + 3] + c5 * r_copy[i + 4];
        z_nodes[i] = c1 * z_copy[i] + c2 * z_copy[i + 1] + c3 * z_copy[i + 2] +
                     c4 * z_copy[i + 3] + c5 * z_copy[i + 4];
        F_nodes[i] = c1 * F_copy[i] + c2 * F_copy[i + 1] + c3 * F_copy[i + 2] +
                     c4 * F_copy[i + 3] + c5 * F_copy[i + 4];
    }

    r_nodes[Ns] = 0.0;
}

/*! Computation of the of the fluid-fluid interface curvature.
    The curvature is computed using a 4th order polynomial fitted over 9 node points.
    For details, see Curtiss (Non-linear, non-spherical bubble dynamics near a two fluid interface, 2009)
  */
void Case_FlatBoundary::boundary_curvature() {

    // Solution vectors
    vector<double> r_copy = r_nodes;
    vector<double> z_copy = z_nodes;
    vector<double> F_copy = F_nodes;

    // Additional dummy node points created to compute polynomial at endpoints of fluid-fluid interface
    r_copy.insert(r_copy.begin(), r_nodes[0] + 5.0);
    r_copy.insert(r_copy.begin(), r_nodes[0] + 10.0);
    r_copy.insert(r_copy.begin(), r_nodes[0] + 15.0);
    r_copy.insert(r_copy.begin(), r_nodes[0] + 20.0);
    r_copy.push_back(-r_nodes[Ns - 1]);
    r_copy.push_back(-r_nodes[Ns - 2]);
    r_copy.push_back(-r_nodes[Ns - 3]);
    r_copy.push_back(-r_nodes[Ns - 4]);

    z_copy.insert(z_copy.begin(), 0.0);
    z_copy.insert(z_copy.begin(), 0.0);
    z_copy.insert(z_copy.begin(), 0.0);
    z_copy.insert(z_copy.begin(), 0.0);
    z_copy.push_back(z_nodes[Ns - 1]);
    z_copy.push_back(z_nodes[Ns - 2]);
    z_copy.push_back(z_nodes[Ns - 3]);
    z_copy.push_back(z_nodes[Ns - 4]);


    //compute distance
    vec rs = conv_to<vec>::from(r_copy);
    vec zs = conv_to<vec>::from(z_copy);
    vec drs = diff(rs);
    vec dzs = diff(zs);
    vec pointss = sqrt(drs % drs + dzs % dzs);
    vec distance_s = cumsum(pointss);
    vec begin(1, fill::zeros);
    distance_s = join_cols(begin, distance_s);

    int degree = 4; // polynomial degree
    int n_elem = 9; // number of nodes points over which the polynomial curve fitting is conducted
    double numerator, denominator;

    double drds_i, d2rds2_i, dzds_i, d2zds2_i;

    vec r_i(n_elem, fill::zeros);
    vec z_i(n_elem, fill::zeros);
    vec distance_i(n_elem, fill::zeros);

    for (int i = 0; i < Ns + 1; ++i) {
        for (int j = 0; j < n_elem; ++j) {
            distance_i(j) = distance_s(i + j);
            r_i(j) = r_copy[i + j];
            z_i(j) = z_copy[i + j];
        }
        double s_tmp = distance_s(i + 4);

        vec r_coeff = polyfit(distance_i, r_i, degree);
        vec z_coeff = polyfit(distance_i, z_i, degree);

        drds_i =
                4.0 * s_tmp * s_tmp * s_tmp * r_coeff(0) + 3.0 * s_tmp * s_tmp * r_coeff(1) + 2.0 * s_tmp * r_coeff(2) +
                r_coeff(3);
        d2rds2_i = 12.0 * s_tmp * s_tmp * r_coeff(0) + 6.0 * s_tmp * r_coeff(1) + 2.0 * r_coeff(2);
        dzds_i =
                4.0 * s_tmp * s_tmp * s_tmp * z_coeff(0) + 3.0 * s_tmp * s_tmp * z_coeff(1) + 2.0 * s_tmp * z_coeff(2) +
                z_coeff(3);
        d2zds2_i = 12.0 * s_tmp * s_tmp * z_coeff(0) + 6.0 * s_tmp * z_coeff(1) + 2.0 * z_coeff(2);


        numerator = drds_i * d2zds2_i - dzds_i * d2rds2_i;
        denominator = pow((drds_i * drds_i + dzds_i * dzds_i), 1.5);
        curv_nodes[i] = numerator / denominator + dzds_i / (r_nodes[i] * pow((drds_i * drds_i + dzds_i * dzds_i), 0.5));

        // at the symmetry axis
        if (i == Ns) {
            curv_nodes[i] = 2 * d2zds2_i / (drds_i * drds_i * drds_i);
        }

    }
}

/*! Computation of the derivatives at the endpoints of fluid-fluid interface.
    Needed for the evaluation of the cubic splines clamped-end conditions.

    The derivatives are known at the axis of symmetry (indicated with 2): dr/ds2 = -1,
    dz/ds2 = 0.0, dphi1/ds2 = 0.0 and dphi2/ds2 = 0.0.
    At the extremum of the domain (indicated with 1), we employ the formulation of Robinson et al.
    (Interaction of cavitation bubbles with a free surface, 2001) for dr/ds and dz/ds. We then
    assume that the derivatives of the potentials phi1 and phi2 are equal to 0 since the fluid-fluid
    interface remains essentially flat and motionless far away from the bubble.
  */

void Case_FlatBoundary::boundary_endpoints_derivatives() {
    vec rs = conv_to<vec>::from(r_nodes);
    vec zs = conv_to<vec>::from(z_nodes);

    drds1 = -1.0 / sqrt(1.0 + (3.0 * zs(0) / rs(0)) * (3.0 * zs(0) / rs(0)));
    drds2 = -1.0;
    dzds1 = 3.0 * zs(0) / rs(0) * 1.0 /
            sqrt(1.0 + (3.0 * zs(0) / rs(0)) * (3.0 * zs(0) / rs(0)));
    dzds2 = 0.0;

    dphi1ds1 = 0.0;
    dphi1ds2 = 0.0;
    dphi2ds1 = 0.0;
    dphi2ds2 = 0.0;
}
