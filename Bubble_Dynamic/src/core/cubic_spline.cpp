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

/*! \file cubic_spline.cpp */

#include <armadillo>
#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>

#include "cubic_spline.hpp"

using namespace arma;

void cubic_spline::set_spline(const std::vector<double> &x, const std::vector<double> &y, double left_derivative,
                              double right_derivative) {


    // Assign value and resize vectors
    int n = x.size();
    m_a.resize(n);
    m_b.resize(n);
    m_c.resize(n);
    m_x = x;
    m_y = y;

    // Check for vector size and monotonicity
    for (int i = 0; i < n - 1; i++) {
        assert(m_x[i] < m_x[i + 1]);
    }
    assert(x.size() == y.size());


    // ----- Setup and solve matrix equation system -----
    // --------------------------------------------------

    // Interpolation function: f(x) = m_a*(x-x_i)^3 + m_b*(x-x_i)^2 + m_c*(x-x_i) + y_i

    // Interior points
    // ---------------

    mat A_cubic(n, n, fill::zeros);
    vec b_cubic(n);
    b_cubic.zeros();
    for (int i = 1; i < n - 1; i++) {
        A_cubic(i, i - 1) = 1.0 / 3.0 * (x[i] - x[i - 1]);
        A_cubic(i, i) = 2.0 / 3.0 * (x[i + 1] - x[i - 1]);
        A_cubic(i, i + 1) = 1.0 / 3.0 * (x[i + 1] - x[i]);
        b_cubic(i) = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
    }

    // Boundary conditions
    // --------------------

    // Left derivative
    A_cubic(0, 0) = 2.0 * (x[1] - x[0]);
    A_cubic(0, 1) = 1.0 * (x[1] - x[0]);
    b_cubic(0) = 3.0 * ((y[1] - y[0]) / (x[1] - x[0]) - left_derivative);

    // Right derivative
    A_cubic(n - 1, n - 1) = 2.0 * (x[n - 1] - x[n - 2]);
    A_cubic(n - 1, n - 2) = 1.0 * (x[n - 1] - x[n - 2]);
    b_cubic(n - 1) = 3.0 * (right_derivative - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));

    // Solve the equation system A_cubic * m_b = b_cubic
    // -------------------------------------------------

    vec m_b_tmp(n);
    m_b_tmp = solve(A_cubic, b_cubic, solve_opts::refine);
    m_b = conv_to < std::vector < double >> ::from(m_b_tmp);

    // Compute m_a and m_c
    // -------------------

    for (int i = 0; i < n - 1; i++) {
        m_a[i] = 1.0 / 3.0 * (m_b[i + 1] - m_b[i]) / (x[i + 1] - x[i]);
        m_c[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - 1.0 / 3.0 * (2.0 * m_b[i] + m_b[i + 1]) * (x[i + 1] - x[i]);
    }
}


double cubic_spline::interpolate(double x) {
    // Checks if x lies within the interpolation interval
    assert(x >= m_x.front());
    assert(x <= 2.0 * m_x.back()); // Times 2.0 is a dummy correction


    // Find the closest index idx so that m_x[idx] < x
    std::vector<double>::const_iterator it;
    it = std::lower_bound(m_x.begin(), m_x.end(), x);
    int idx = std::max(int(it - m_x.begin()) - 1, 0);

    double dx = x - m_x[idx];
    double interpolation = m_a[idx] * dx * dx * dx + m_b[idx] * dx * dx + m_c[idx] * dx + m_y[idx];

    return interpolation;
}
