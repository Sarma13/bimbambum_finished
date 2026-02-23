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

/*! \file cubic_spline.hpp
    \brief Header for cubic_spline.cpp.
*/

/*! \class cubic_spline
    \brief Conduct a cubic spline interpolation through the nodes on the bubble surface and on
    the fluid-fluid interface.
*/

#ifndef CUBICSPLINE_HPP
#define CUBICSPLINE_HPP


class cubic_spline {
public:
    cubic_spline() {}

    ~cubic_spline() {}

    std::vector<double> m_a, m_b, m_c;        // spline coefficients
    void set_spline(const std::vector<double> &x, const std::vector<double> &y, double left_derivative,
                    double right_derivative);

    double interpolate(double x);

private:

    std::vector<double> m_x, m_y;            // x,y coordinates of points
    double left_derivative, right_derivative; //endpoints derivatives
};

#endif // CUBICSPLINE_HPP
