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

/*! \file BubbleData.cpp
    \brief Base class containing function definition, properties and quantities associated with the bubble
    dynamics.
*/

#include <cassert>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include <armadillo>
#include <gsl/gsl_integration.h>

#include "BubbleData.hpp"
#include "cubic_spline.hpp"
#include "integrands_volume.hpp"

using namespace std;
using namespace arma;

/*! Filtering of the of the bubble interface nodes (to avoid the onset of saw-tooth effects).
    Filtering is achieved with a five-point smoothing formula.
    For details, see Longuet-Higgins and Cokelet (The deformation of steep surface waves on water, 1976)
    or Pearson (Hydrodynamics of jet impact in a collapsing bubble, 2002).
  */

void BubbleData::filter_bubble() {

    // generates surface vectors
    vector<double> r_copy = r_nodes;
    vector<double> z_copy = z_nodes;
    vector<double> phi_copy = phi_nodes;

    // Additional dummy node points created to compute polynomial at endpoints of bubble interface
    r_copy.insert(r_copy.begin(), -r_nodes[1]);
    r_copy.insert(r_copy.begin(), -r_nodes[2]);
    r_copy.push_back(-r_nodes[Nb - 1]);
    r_copy.push_back(-r_nodes[Nb - 2]);

    z_copy.insert(z_copy.begin(), z_nodes[1]);
    z_copy.insert(z_copy.begin(), z_nodes[2]);
    z_copy.push_back(z_nodes[Nb - 1]);
    z_copy.push_back(z_nodes[Nb - 2]);

    phi_copy.insert(phi_copy.begin(), phi_nodes[1]);
    phi_copy.insert(phi_copy.begin(), phi_nodes[2]);
    phi_copy.push_back(phi_nodes[Nb - 1]);
    phi_copy.push_back(phi_nodes[Nb - 2]);

    // compute the arclength
    vec rb = conv_to<vec>::from(r_copy);
    vec zb = conv_to<vec>::from(z_copy);
    vec drb = diff(rb);
    vec dzb = diff(zb);
    vec pointsb = sqrt(drb % drb + dzb % dzb);
    vec distance_b = cumsum(pointsb);
    vec begin(1, fill::zeros);
    distance_b = join_cols(begin, distance_b);
    vector<double> distance = conv_to < std::vector < double >> ::from(distance_b);

    double s1, s2, s3, s4, s5;
    double c1, c2, c3, c4, c5;

    // Longuet-Higgings & Cokelet smoothing for irregular node distribution
    for (int i = 0; i < Nb + 1; ++i) {
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
        r_nodes[i] = c1 * r_copy[i] + c2 * r_copy[i + 1] + c3 * r_copy[i + 2] + c4 * r_copy[i + 3] + c5 * r_copy[i + 4];
        z_nodes[i] = c1 * z_copy[i] + c2 * z_copy[i + 1] + c3 * z_copy[i + 2] + c4 * z_copy[i + 3] + c5 * z_copy[i + 4];
        phi_nodes[i] = c1 * phi_copy[i] + c2 * phi_copy[i + 1] + c3 * phi_copy[i + 2] + c4 * phi_copy[i + 3] +
                       c5 * phi_copy[i + 4];
    }

    // imposed value at the axis of symmetry
    r_nodes[0] = 0.0;
    r_nodes[Nb] = 0.0;


}

/*! Remeshing of the of the bubble interface nodes.
    The remeshing is done by rearranging the node points along the bubble interface with a
    spacing defined by a geometric progression.
  */

void BubbleData::remesh_bubble() {

    // generates surface vectors
    vec rb = conv_to<vec>::from(r_nodes);
    vec zb = conv_to<vec>::from(z_nodes);

    // compute the arclength
    vec drb = diff(rb);
    vec dzb = diff(zb);
    vec pointsb = sqrt(drb % drb + dzb % dzb);
    vec distance_b = cumsum(pointsb);
    vec begin(1, fill::zeros);
    distance_b = join_cols(begin, distance_b);
    vector<double> distance = conv_to < std::vector < double >> ::from(distance_b);

    cubic_spline spline_rb;
    spline_rb.set_spline(distance, r_nodes, 1.0, -1.0);

    cubic_spline spline_zb;
    spline_zb.set_spline(distance, z_nodes, 0.0, 0.0);

    cubic_spline spline_phib;
    spline_phib.set_spline(distance, phi_nodes, 0.0, 0.0);

    // rearrange the nodes along the arclength - geometrical progression re-griding technique
    double n = Nb + 1;
    double r = 0.9999; // this value generates practically evenly spaced nodes
    double h = (r - 1.0) / (pow(r, (1.0 * n - 1.0)) - 1.0);

    vec e(Nb);

    for (int i = 0; i < Nb; ++i) {
        e(i) = pow(r, 1.0 * i);
    }

    vec s = h * cumsum(e) * distance[Nb];
    s = join_cols(begin, s);
    s(0) = distance[0];
    s(Nb) = distance[Nb];

    for (int i = 0; i < Nb + 1; ++i) {
        r_nodes[i] = spline_rb.interpolate(s(i));
        z_nodes[i] = spline_zb.interpolate(s(i));
        phi_nodes[i] = spline_phib.interpolate(s(i));
    }

    // imposed value at the axis of symmetry
    r_nodes[0] = 0.0;
    r_nodes[Nb] = 0.0;
}

/*! Check for bubble intersection (condition for end of simulation).

    Checks if any line segment (constructed between adjacent node points on the
    bubble surface) intersects another line segment.

    The first line segment is defined by two distinct points (r1, z1) and (r2, z2)
    and the second line segment by (r3, z3) and (r4, z4) and therefore the associated
    lines by L1 = [x1; y1] + t*[(x2-x1); (y2-y1)] and L2 = [x3; y3] + u*[(x4-x3); (y4-y3)]
    The two line segments intersect if 0 < t < 1 and 0 < u < 1
 */

int BubbleData::intersect() {

    double dri;
    double drj;
    double dzi;
    double dzj;
    double t;
    double u;
    int intersection = 0;

    for (int i = 0; i < Nb; ++i) {
        for (int j = 0; j < Nb; ++j) {
            if (i != j) {
                dri = r_nodes[i + 1] - r_nodes[i];
                dzi = z_nodes[i + 1] - z_nodes[i];
                drj = r_nodes[j + 1] - r_nodes[j];
                dzj = z_nodes[j + 1] - z_nodes[j];
                t = (dzj * (r_nodes[j] - r_nodes[i]) + drj * (z_nodes[i] - z_nodes[j])) / (dzj * dri - drj * dzi);
                u = (dzi * (r_nodes[i] - r_nodes[j]) + dri * (z_nodes[j] - z_nodes[i])) / (dzi * drj - dri * dzj);
                if ((t > 0.0) && (t < 1.0) && (u > 0.0) && (u < 1.0)) {
                    intersection = 1;
                }
            }
        }
    }

    return intersection;
}

/*! Estimation of most appropriate time step based on the bubble interface dynamics.*/
double BubbleData::time_step_bubble(double epsilon, double k) {

    double u_nodes_max = u_nodes[0];
    for (int i = 0; i < Nb + 1; ++i) {
        u_nodes_max = max(u_nodes_max, u_nodes[i]);
    }

    return delta_phi / (1.0 + 0.5 * u_nodes_max * u_nodes_max + epsilon * pow(V0 / V, k));
}

/*! Computes the bubble volume.
    For details, see Pearson (Hydrodynamics of jet impact in a collapsing bubble, 2002).
*/
double BubbleData::compute_volume() {

    double volume = 0.0;

    // generates surface vectors
    vec rb = conv_to<vec>::from(r_nodes);
    vec zb = conv_to<vec>::from(z_nodes);

    // compute the arclength
    vec drb = diff(rb);
    vec dzb = diff(zb);
    vec pointsb = sqrt(drb % drb + dzb % dzb);
    vec distance_b = cumsum(pointsb);
    vec begin(1, fill::zeros);
    distance_b = join_cols(begin, distance_b);
    vector<double> distance = conv_to < std::vector < double >> ::from(distance_b);

    cubic_spline spline_rb;
    spline_rb.set_spline(distance, r_nodes, 1.0, -1.0);

    cubic_spline spline_zb;
    spline_zb.set_spline(distance, z_nodes, 0.0, 0.0);

    // perform the integration
    double a, b;
    double m_ar, m_br, m_cr, m_dr, m_az, m_bz, m_cz, m_dz;

    for (int i = 0; i < Nb; ++i) {

        m_ar = spline_rb.m_a[i];
        m_br = spline_rb.m_b[i];
        m_cr = spline_rb.m_c[i];
        m_dr = rb(i);

        m_az = spline_zb.m_a[i];
        m_bz = spline_zb.m_b[i];
        m_cz = spline_zb.m_c[i];
        m_dz = zb(i);

        gsl_integration_workspace *w_volume = gsl_integration_workspace_alloc(1000);
        double result, error;

        volume_params params;
        params.coeff_r0 = m_ar;
        params.coeff_r1 = m_br;
        params.coeff_r2 = m_cr;
        params.coeff_r3 = m_dr;
        params.coeff_z0 = m_az;
        params.coeff_z1 = m_bz;
        params.coeff_z2 = m_cz;
        params.coeff_z3 = m_dz;

        gsl_function vol;
        vol.function = &volume_b;
        vol.params = &params;

        a = 0.0;
        b = distance_b[1 + i] - distance_b[i];

        gsl_integration_qags(&vol, a, b, 0.0, 1e-8, 1000, w_volume, &result, &error);

        volume = volume + result;
        gsl_integration_workspace_free(w_volume);

    }
    return volume;
}

double BubbleData::compute_spherical_volume() {

    double volume = 0.0;

    // generates surface vectors
    vec rb = conv_to<vec>::from(r_nodes);
    vec zb = conv_to<vec>::from(z_nodes);

    // compute the arclength
    vec drb = diff(rb);
    vec dzb = diff(zb);
    vec pointsb = sqrt(drb % drb + dzb % dzb);
    vec distance_b = cumsum(pointsb);
    vec begin(1, fill::zeros);
    distance_b = join_cols(begin, distance_b);
    vector<double> distance = conv_to < std::vector < double >> ::from(distance_b);

    cubic_spline spline_rb;
    spline_rb.set_spline(distance, r_nodes, 1.0, -1.0);

    cubic_spline spline_zb;
    spline_zb.set_spline(distance, z_nodes, 0.0, 0.0);

    // perform the integration
    double a, b;
    double m_ar, m_br, m_cr, m_dr, m_az, m_bz, m_cz, m_dz;

    for (int i = 0; i < Nb; ++i) {

        m_ar = spline_rb.m_a[i];
        m_br = spline_rb.m_b[i];
        m_cr = spline_rb.m_c[i];
        m_dr = rb(i);

        m_az = spline_zb.m_a[i];
        m_bz = spline_zb.m_b[i];
        m_cz = spline_zb.m_c[i];
        m_dz = zb(i);

        gsl_integration_workspace *w_volume = gsl_integration_workspace_alloc(1000);
        double result, error;

        volume_params params;
        params.coeff_r0 = m_ar;
        params.coeff_r1 = m_br;
        params.coeff_r2 = m_cr;
        params.coeff_r3 = m_dr;
        params.coeff_z0 = m_az;
        params.coeff_z1 = m_bz;
        params.coeff_z2 = m_cz;
        params.coeff_z3 = m_dz;

        gsl_function vol;
        vol.function = &volume_spherical_droplet;  // Use the spherical droplet volume function
        vol.params = &params;

        a = 0.0;
        b = distance_b[1 + i] - distance_b[i];

        gsl_integration_qags(&vol, a, b, 0.0, 1e-8, 1000, w_volume, &result, &error);

        volume = volume + result;
        gsl_integration_workspace_free(w_volume);
    }
    
    return volume;
}
