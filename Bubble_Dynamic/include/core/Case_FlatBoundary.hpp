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

/*! \file Case_FlatBoundary.hpp
    \brief Header for Case_FlatBoundary.cpp
*/

/*! \class Case_FlatBoundary
    \brief Derived class: handles the dynamics of fluid-fluid interface that has been initiated as a 'infinite' flat surface.

    The fluid-fluid interface is taken as initially flat and located at z = 0 with potentials on both side equal to zero.
    This derived class inherits members from BoundaryData.cpp and BoundaryData.hpp.


*/

#ifndef CASE_FLATBOUNDARY_HPP
#define CASE_FLATBOUNDARY_HPP

#include "BoundaryData.hpp"

class Case_FlatBoundary : public BoundaryData {
public:
    template<typename Input>
    Case_FlatBoundary(Input &data): BoundaryData(data) {}

    virtual ~Case_FlatBoundary() {}

    /*! \brief Initializes fluid-fluid interface node points.*/
    void initialize() override;

    /*! \brief Redistributes the node points position along fluid-fluid interface. */
    void remesh_boundary() override;

    /*! \brief Filters the fluid-fluid interface. */
    void filter_boundary() override;

    /*! \brief Estimates the derivatives at the endpoints on the fluid-fluid interface. */
    void boundary_curvature() override;

    /*! \brief Computes the time step value based on the fluid-fluid interface dynamics*/
    void boundary_endpoints_derivatives() override;

};


#endif //CASE_FLATBOUNDARY_HPP
